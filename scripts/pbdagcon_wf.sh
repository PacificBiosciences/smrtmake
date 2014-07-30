#!/bin/bash

# Simple pbdagcon workflow script.  Written for the benefit of running via
# smrtpipe so I can communicate pipe errors to the task.  We're overcoming
# the limitation of smrtpipe forcing tasks to run serially, enabling a new
# level of pipelining that's extremely efficient in an imperfect world ...
# However, direct file I/O is faster by default.

tmp=${tmp-"/tmp"}

trap "rm -f $tmp/aln.$$.pre" EXIT SIGINT

echo "Generating pre-alignments"
echo "m4topre.py $mym4 $allm4 $basfofn $rgnfofn ${bestn-24} > $tmp/aln.$$.pre"

# generate pre-alignments to a tmp directory
m4topre.py $mym4 $allm4 $basfofn $rgnfofn ${bestn-24} > $tmp/aln.$$.pre || exit $?

echo "Correcting reads"
# pipe it to consensus and generate fasta
pbdagcon -c ${cov-8} -a -j ${nproc-15} $tmp/aln.$$.pre | tee ${fasta-"corrected.fa"} | \
# generate a fastq
awk '{if($0~/>/){sub(/>/,"@",$0);print;}else{l=length($0);q="";while(l--){q=q "9"}printf("%s\n+\n%s\n",$0,q)}}' > ${fastq-"corrected.fq"}


# check the status of each piped command and exit non-zero if found
for exitval in ${PIPESTATUS[*]}
do
    if [ $exitval -gt 0 ]
    then
        exit $exitval
    fi
done


exit 0;
