#!/usr/bin/env python

import sys
from pbcore.io import BaxH5Reader

def main():
    """Sorts and prints a list of subread lengths from a base file."""
    basfofn = sys.argv[1]
    rgnfofn = sys.argv[2]
    cutoff = int(open(sys.argv[3]).read().strip())

    basfiles = open(basfofn).read().splitlines()
    rgnfiles = open(rgnfofn).read().splitlines()

    for basfile, rgnfile in zip(basfiles, rgnfiles):
        reader = BaxH5Reader(basfile, rgnfile)
        for read in reader.subreads():
            if len(read) >= cutoff:
                print ">%s\n%s\n" % (read.readName, read.basecalls())

if __name__ == '__main__':
    sys.exit(main())
