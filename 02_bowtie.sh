for sample in `ls *_trimmed.R1.fq.gz`
do
base=$(basename $sample "_trimmed.R1.fq.gz")

bowtie2 -p 10 -X 2000 -x /database/mm10_bowtie2/mm10 -1 ./${base}_trimmed.R1.fq.gz -2 ./${base}_trimmed.R2.fq.gz -S ../02_alignment/${base}.sam  2>../02_alignment/${base}.align.log
done&