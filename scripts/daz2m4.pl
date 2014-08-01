#!/usr/bin/env perl

my $dbFile = shift;
my $lasFile = shift;

sub pbidToLen {
    my $pbid = shift;
    my @v = $pbid =~ m#\d+/(\d+)_(\d+)#; 
    return $v[1] - $v[0];
}

sub prolog {
    $iid = shift;
    foreach my $partition (@dbmeta) {
        if ($iid <= $partition->[0]) {
            return $partition->[1];
        }
    } 
}

sub collapse {
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
    
    my $movie = prolog $qid;
    my $pbidQ = $movie."/".$rinfo{$qid};
    my $pbidT = $movie."/".$rinfo{$tid};
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

    # XXX limit to bestn, parameterize
    foreach my $rec (@sorted[0..24]) {
        print "$rec\n"; 
    }
}

# XXX fix this when daz partitions are better understood
open DBM, $dbFile;
our @dbmeta = ();
while (<DBM>) {
    if (/\s+(\d+) m\d+_\S+ (m\d+_\S+)/) {
        my $maxId = $1;
        my $movPfx = $2; 
        push @dbmeta, [$maxId, $movPfx];
    }
}
close DBM;

open DB, "DBshow $dbFile|";
%rinfo = {};
$iid = 1;
while (<DB>) {
    if (m#>Prolog/(\S+) RQ#) {
        my $pbid = $1;
        $rinfo{$iid++} = $pbid;
    }
}
close DB;

my $prevPair = "0:0";
my $prevQ = 0;
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
    my $currQ = $f[0];
    my $currPair = "$f[0]:$f[1]";
    if ($currPair ne $prevPair) {
        if ($prevPair ne "0:0") {
            my $m4 = collapse(@recs);
            push @alnset, $m4;
            @recs = (); 
        }
        if ($prevQ != $currQ) {
            if ($prevQ != 0) {
                emitSet @alnset;
                @alnset = ();
            }
            $prevQ = $currQ;
        }
        $prevPair = $currPair;
    }
    push @recs, $_;
}
close LA;

collapse @recs;
emitSet @alnset;

