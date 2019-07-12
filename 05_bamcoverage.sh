for i in ./*_rmdup_picard.bam
do
base=$(basename $i "_rmdup_picard.bam")

num1=10000000
num2="$(samtools view -c  $i  2>&1 )"
res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
bamCoverage --scaleFactor  $res -b  $i   -o   ./${base}.10M.bw
done