#!/usr/bin/env perl

sub pbidToLen {
    my $pbid = shift;
    my @v = $pbid =~ m#\d+/(\d+)_(\d+)#; 
    return $v[1] - $v[0];
}

sub selfHit {
    my $qid = shift;
    my $pbidQ = $rinfo{$qid};
    my $lenQ = pbidToLen $pbidQ;
    return "$qid $qid -40000 100.0000 0 0 $lenQ $lenQ 0 0 $lenQ $lenQ 0";
}

sub wellBehaved {
    my ($qid, $tid, $strand, $minQ, undef, $minT) = split / +/, $_[0];
    my (undef, undef, undef, undef, $maxQ, undef, $maxT) = split / +/, $_[-1]; 

    my $pbidQ = $rinfo{$qid};
    
    return $maxT > $minT;
}

sub collapse {
    # Most are simply well-behaved one line or incrementing alignment coordinates ...
    #   2     19,611 c   [ 2,829.. 3,874] x [     0.. 1,000] :   <    202 diffs  ( 10 trace pts)
    #   2     19,611 c   [ 3,867.. 7,146] x [ 1,043.. 4,182] :   <    598 diffs  ( 33 trace pts)
    #   2     19,611 c   [ 7,164.. 9,226] x [ 4,257.. 6,220] :   <    363 diffs  ( 21 trace pts)

    # but, there are more complex cases where alignments overlap or map to multiple locations.
    # These are likely chimeric in nature and should be skipped.
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
    my $pbidT = $rinfo{$tid};
    my $lenQ = pbidToLen $pbidQ;
    my $lenT = pbidToLen $pbidT;
    my $s = $strand eq "c" ? 1 : 0;
    my $alndiff = abs($alnlenQ - $alnlenT);
    my $pctid = sprintf "%.4f", (1 - ($diffcnt + $alndiff) / $alnlenQ) * 100;
    my $score = $alnlenQ;
    return "$qid $tid -$score $pctid 0 $minQ $maxQ $lenQ $s $minT $maxT $lenT 0";
}

sub emitSet {
    my @sorted = map { $_->[0] } 
                 sort { $a->[1] <=> $b->[1] } 
                 map { [$_, (split)[2]] } @_;

    my $count = 0;
    # take only the first 200 alignments
    foreach my $rec (grep {$_} @sorted[0..198]) {
        print "$rec\n";
    }
}

my $block = shift;
my $readFile = shift;
my $lasFile = shift;

# load the the subread mapping information
open RF, $readFile || die "Failed to open $readFile: $!";
our %seeds;
our %rinfo;
my $iidOffs = 0;
while (<RF>) {
    chomp;
    my ($blk, $iid, $pbid, $seed) = split;
    $rinfo{$iid} = $pbid;
    $seeds{$iid}++ if $seed and $blk == $block;
}
close RF;

my @recset;
my $prevPair = "0:0:0";
my $prevT = 0;
my $currQ = 0;
my $currT = 0;
my @frgset;
my @alnset;

open LA, "LAshow $lasFile | sort -sk2,2 -T/lustre/archive/tmp |";
while (<LA>) {
    # <spaces>
    # subreads.merge: 5,688,372 records
    next if /^\S*$/;
    next if /records/;

    #       2        320 c   [     0.. 2,381] x [ 1,951.. 4,455] :   <    413 diffs  ( 23 trace pts)
    chomp;
    s/[],\[<:x]//g;
    s/\.\./ /g;
    s/ diffs.*//;
    s/ +/ /g;
    s/^ +//;

    # 2 320 c 0 2381 1951 4455 413
    my @f = split;
    ($currQ, $currT, $currS) = @f[0..2];

    # screen out non-seed targets
    next unless exists $seeds{$currT};

    # finished processing this target, print out the m4 records and empty the 
    # alignment set.
    if ($prevT != $currT) {
        if ($prevT != 0) {
            # add a self-hit
            push @alnset, selfHit $prevT; 
            emitSet @alnset;
        }
        $prevT = $currT;
        @alnset = ();
        @frgset = ();
    }

    # handle the next pair for the current target.
    # each query:target pair can appear 1 or more times in the alignment set
    my $currPair = "$currQ:$currT:$currS";
    if ($currPair ne $prevPair) {
        # only collapse and add the record set if it's well behaved.
        if ($#frgset > -1 && wellBehaved @frgset) {
            push @alnset, collapse @frgset;
        }

        @frgset = ();
        $prevPair = $currPair;
    }

    # add the next alignment fragment to the current alignment fragment set
    push @frgset, $_;
}
close LA;

if ($#frgset > -1 && wellBehaved @frgset) {
    push @alnset, collapse @frgset;
}
push @alnset, selfHit $currT;
emitSet @alnset;

