#!/usr/bin/env python
"""Super-simple converter from blasr m4 alignments to pbdagcon 'pre'
alignments. For use in the pre-assembler dagcon workflow.
"""

import sys
import heapq
import string # pylint: disable=W0402
from itertools import ifilter
from collections import namedtuple, defaultdict
import numpy as np
from pbcore.io import BaxH5Reader

# qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart
# ... tend tseqlength mapqv
#
# store only fields we need
__m4fields__ = [0, 1, 2, 5, 6, 8, 9, 10, 11]
M4RECORD = namedtuple(
    'M4RECORD', 'qname tname score qstart qend tstrand tstart tend tseqlength')

__tuplfy__ = M4RECORD._make # pylint: disable=W0212

# dna compliment
__rc__ = string.maketrans('actgACTG', 'tgacTGAC')


def parse_m4(rec):
    """Parse in the m4 file, returning a list of records"""
    return [y for (x, y) in enumerate(rec.split()) if x in __m4fields__]


def rating(rec):
    """Rates the alignment for by length and score (revisit at some point)"""
    score = -int(rec.score)
    alen = int(rec.tend) - int(rec.tstart)
    return score + alen


def schwartzian(rec):
    """Provides a schwartzian transform for the given record, used for 
    sorting
    """
    flds = rec.split()
    return (flds[1], float(flds[2]), rec)


def sort_targ_score(recs):
    """Sorts the list in place by target id (string), then score (float)"""
    recs[:] = [schwartzian(x) for x in recs]
    recs.sort()
    recs[:] = [rec for (target, score, rec) in recs] # pylint: disable=W0612


def rescore(recs):
    """Rescore alignments using coverage based statistics"""
    prev = ""
    cov = np.zeros(1)
    for idx, rec in enumerate(recs):
        fields = rec.split()
        rec = __tuplfy__(fields)
        if rec.tname != prev:
            prev = rec.tname
            cov = np.zeros(int(rec.tseqlength), dtype=np.float16)

        if rec.tstrand:
            start = int(rec.tseqlength) - int(rec.tend)
            end = int(rec.tseqlength) - int(rec.tstart)
        else:
            start = int(rec.tstart)
            end = int(rec.tend)

        cov[start:end] += 1
        score = np.sum(1/cov[start:end])
        fields[2] = str(-score)
        recs[idx] = " ".join(fields)


def bestn_true(recstr, myq):
    """Checks if the record falls inside bestn (used when blasr is chunked)"""
    rec = __tuplfy__(recstr.split())
    rate = rating(rec)
    return rate in myq[rec.qname[32:]].top


class AlnLimiter(object): # pylint: disable=R0903
    """Functor that returns alignments until some count is reached. Alignments
    should be sorted. 
    """
    def __init__(self, limit=76):
        self.count = 0
        self.target = ''
        self.limit = limit

    def __call__(self, rec):
        target = rec.split()[1]
        if target != self.target:
            self.count = 0
            self.target = target
        self.count += 1
        return self.count < self.limit


class TopAlignments(object): # pylint: disable=R0903
    """Tracks the top alignments for a given query, used for bestn calc"""
    bestn = 10

    def __init__(self):
        self.top = [0] * TopAlignments.bestn

    def __call__(self):
        return  # noop

    def add(self, aln):
        """Adds an alignment to a bounded list, kicking out another if 
        necessary
        """
        heapq.heappushpop(self.top, aln)


def main(): # pylint: disable=R0914
    """Drives the program"""
    _, mym4, allm4, basfofn, rgnfofn, TopAlignments.bestn = sys.argv
    TopAlignments.bestn = int(TopAlignments.bestn)

    # tracks bestn
    my_queries = defaultdict(TopAlignments)
    my_m4recs = []

    # load my m4 chunk
    m4h = open(mym4)
    rec_add = my_m4recs.append
    for line in m4h:
        flds = parse_m4(line)
        rec = __tuplfy__(flds)
        rate = rating(rec)
        my_queries[rec.qname[32:]].add(rate)
        rec_add(' '.join(flds))

    m4h.close()

    # if we're chunked locate relevant alignments
    if mym4 != allm4:
        # assuming fofn here
        m4files = [x.rstrip() for x in open(allm4) if x.rstrip() != mym4]
        for m4f in m4files:
            m4h = open(m4f)
            for recstr in m4h:
                rec = __tuplfy__(parse_m4(recstr))
                if rec.qname[32:] in my_queries:
                    rate = rating(rec)
                    my_queries[rec.qname[32:]].add(rate)
            m4h.close()

        # remove alignments that fall outside of bestn
        my_m4recs[:] = [x for x in my_m4recs if bestn_true(x, my_queries)]

    # sort by target name/score
    sort_targ_score(my_m4recs)

    # rescore based on coverage
    rescore(my_m4recs)

    # sort one more time be new score
    sort_targ_score(my_m4recs)

    # take a max number of alignments for each target
    limiter = AlnLimiter()
    my_m4recs[:] = [x for x in ifilter(limiter, my_m4recs)]

    # load only related sequences
    seqs = {}
    basfiles = open(basfofn).read().splitlines()
    rgnfiles = open(rgnfofn).read().splitlines()
    for bas, rgn in zip(basfiles,rgnfiles):
        reader = BaxH5Reader(bas, rgn)
        for sr in reader.subreads():
            if sr.readName[32:] in my_queries:
                seqs[sr.readName] = sr.basecalls() 

    # may or may not help
    del my_queries

    # generate pre-alignments
    for recstr in my_m4recs:
        rec = __tuplfy__(recstr.split())

        # Bug 24538, rare case missing self hit
        if rec.tname not in seqs:
            msg = "Warning: skipping query %s target %s\n"
            sys.stderr.write(msg % (rec.qname, rec.tname))
            continue

        qst = int(rec.qstart)
        qnd = int(rec.qend)
        qseq = seqs[rec.qname][qst:qnd]
        strand = '-' if rec.tstrand == '1' else '+'
        tst = int(rec.tstart)
        tnd = int(rec.tend)
        if strand == '+':
            tseq = seqs[rec.tname][tst:tnd]
        else:
            tseq = seqs[rec.tname].translate(__rc__)[::-1][tst:tnd]

        print ' '.join([rec.qname, rec.tname, strand,
                       rec.tseqlength, str(tst), str(tnd), qseq, tseq])

if __name__ == '__main__':
    sys.exit(main())
