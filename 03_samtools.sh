for sample in `ls *.sam`
do
base=$(basename $sample ".sam")
samtools view -h -b -S  ${base}.sam   |  samtools sort  -  |  samtools view -h -bq 30 -  >  ./${base}_sort_uni.bam
done
