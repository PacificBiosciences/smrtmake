#!/usr/bin/env perl

my $chunkFile = shift;
my $idMapFile = shift;
my $lasFile = shift;
my $cutFile = shift;

sub pbidToLen {
    my $pbid = shift;
    my @v = $pbid =~ m#\d+/(\d+)_(\d+)#; 
    return $v[1] - $v[0];
}

sub selfHit {
    my $qid = shift;
    my $pbidQ = $rinfo{$qid};
    my $lenQ = pbidToLen $pbidQ;
    return "$pbidQ $pbidQ -40000 100.0000 0 0 $lenQ $lenQ 0 0 $lenQ $lenQ 0";
}

sub chimeric {
    my ($qid, $tid, $strand, $minQ, undef, $minT) = split / +/, $_[0];
    my (undef, undef, undef, undef, $maxQ, undef, $maxT) = split / +/, $_[-1]; 
    
    return $maxT < $minT;
}

sub collapse {
    # Most are simply well-behaved one line or incrementing alignment coordinates ...
    #   2     19,611 c   [ 2,829.. 3,874] x [     0.. 1,000] :   <    202 diffs  ( 10 trace pts)
    #   2     19,611 c   [ 3,867.. 7,146] x [ 1,043.. 4,182] :   <    598 diffs  ( 33 trace pts)
    #   2     19,611 c   [ 7,164.. 9,226] x [ 4,257.. 6,220] :   <    363 diffs  ( 21 trace pts)

    # but, there are more complex cases where alignments overlap or map to multiple locations.
    # These are likely chimeric and should be skipped.
    #     7     36,114 n   [ 8,555..10,341] x [ 4,718.. 6,476] :   <    298 diffs  ( 18 trace pts)
    #     7     36,114 n   [ 8,555..12,992] x [ 4,718.. 9,166] :   <    835 diffs  ( 44 trace pts)
    #     7     36,114 n   [12,930..14,674] x [     0.. 1,722] :   <    302 diffs  ( 17 trace pts)
    #     7     36,114 n   [12,930..15,948] x [     0.. 2,989] :   <    565 diffs  ( 30 trace pts)
    # --
    #     7     36,114 c   [ 9,415..12,844] x [ 6,172.. 9,584] :   <    675 diffs  ( 34 trace pts)
    #     7     36,114 c   [12,986..14,673] x [   626.. 2,360] :   <    321 diffs  ( 17 trace pts)
    #     7     36,114 c   [12,986..15,948] x [   626.. 3,607] :   <    563 diffs  ( 30 trace pts)
    # ...
    #   279     51,623 n   [     0.. 3,254] x [ 6,333.. 9,553] :   <    662 diffs  ( 32 trace pts)
    # --
    #   279     51,623 c   [     0.. 3,254] x [10,508..13,597] :   <    660 diffs  ( 32 trace pts)
    #   279     51,623 c   [   609.. 2,522] x [     0.. 1,963] :   <    421 diffs  ( 19 trace pts)
 
    my ($qid, $tid, $strand, $minQ, undef, $minT) = split / +/, $_[0];
    my (undef, undef, undef, undef, $maxQ, undef, $maxT) = split / +/, $_[-1]; 
    my $alnlenQ = 0;
    my $alnlenT = 0;
    my $diffcnt = 0;
    for my $rec (@_) {
        my @f = split / +/, $rec;
        $alnlenQ += $f[4] - $f[3];
        $alnlenT += $f[6] - $f[5];  
        $diffcnt += $f[7];
    }
    
    my $pbidQ = $rinfo{$qid};
    my $pbidT = $seeds{$tid};
    my $lenQ = pbidToLen $pbidQ;
    my $lenT = pbidToLen $pbidT;
    my $s = $strand eq "c" ? 1 : 0;
    my $alndiff = abs($alnlenQ - $alnlenT);
    my $pctid = sprintf "%.4f", (1 - ($diffcnt + $alndiff) / $alnlenQ) * 100;
    my $score = $alnlenQ;
    return "$pbidQ $pbidT -$score $pctid 0 $minQ $maxQ $lenQ $s $minT $maxT $lenT 0";
}

sub emitSet {
    my @sorted = map { $_->[0] } 
                 sort { $a->[1] <=> $b->[1] } 
                 map { [$_, (split)[2]] } @_;

    # XXX parameterize?
    my $count = 0;
    foreach my $rec (@sorted) {
        print "$rec\n";
        last if ++$count > 50;
    }
}

open CUT, $cutFile;
our $cutoff = <CUT>;
chomp $cutoff;

our %seeds = {};
open IDM, $idMapFile;
while (<IDM>) {
    chomp;
    my ($iid, $pbid) = split;
    my $len = pbidToLen $pbid;
    if ($len >= $cutoff) {
        $seeds{$iid} = $pbid;
    }
}

open CHUNK, $chunkFile; 
our %rinfo = {};
while (<CHUNK>) {
    chomp;
    my ($iid, $pbid) = split;
    $rinfo{$iid} = $pbid;
}
close CHUNK;

my $prevPair = "0:0:0";
my $prevQ = 0;
my $currQ = 0;
my $currT = 0;
my @recs = ();
my @alnset = ();
open LA, "LAshow $lasFile|";
while (<LA>) {
    # <spaces>
    # subreads.merge: 5,688,372 records
    next if /^\S*$/;
    next if /records/;

    #       2        320 c   [     0.. 2,381] x [ 1,951.. 4,455] :   <    413 diffs  ( 23 trace pts)
    chomp;
    s/[],\[<:x]//g;
    s/\.\./ /g;
    s/^ +//;

    # 2        320 c        0  2381   1951  4455        413 diffs  ( 23 trace pts)
    my @f = split;
    # currS: strand (n = forward, c = reverse)
    ($currQ, $currT, $currS) = @f[0..2];
    if (not exists $rinfo{$currQ}) {
        next;
    }

    # finished processing this query, print out the m4 records and empty the alignment set
    if ($prevQ != $currQ) {
        # first add a self-hit if the previous query is big enough to be a seed read
        push @alnset, selfHit $prevQ if exists $seeds{$prevQ};
        emitSet @alnset;
        @alnset = ();
        $prevQ = $currQ;
    }

    # handle the next pair for the current query
    # each query:target pair can appear 1 or more times in the alignment set
    my $currPair = "$currQ:$currT:$currS";
    if ($currPair ne $prevPair) {
        # only collapse and add the record set if its target is a seed read, collapse the 
        # alignment record set into an alignment super-set.
        @qtp = split /:/, $prevPair;
        if (exists $seeds{$qtp[1]} and not chimeric @recs) {
            push @alnset, collapse @recs;
        }

        @recs = ();
        $prevPair = $currPair;
    }

    # always add the next alignment record to the set.
    push @recs, $_;
}
close LA;

if (exists $seeds{$currT} and not chimeric @recs) {
    push @alnset, collapse @recs;
}
push @alnset, selfHit $currQ if exists $seeds{$currQ};
emitSet @alnset if exists $rinfo{$currQ};

