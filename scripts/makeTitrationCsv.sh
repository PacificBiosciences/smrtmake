
echo Condition,Coverage,Variants,ConcordanceQV

for file in Output/*/MaskedVariants/*.gff; do
    conditionName=`echo $file | awk -F/ '{print $(NF-2)}'`
    coverage=`echo $file | awk -F. '{print $(NF-1)}'`
    genomeLength=`grep sequence-region $file | awk '{print $(NF)}'`
    variantCount=`egrep -v '^#' $file | wc -l`
    errorRate=`bc -l <<< "scale=7;($variantCount+1)/($genomeLength+1)"`
    qv=`bc -l <<< "scale=1;-10*l($errorRate)/l(10)"`
    echo $conditionName,$coverage,$variantCount,$qv
done
