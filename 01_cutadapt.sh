for i in `ls *_1.fq.gz`
do 
base=$(basename $i "_1.fq.gz")
cutadapt -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B CTGTCTCTTATACACATCTGACGCTGCCGACGA -q 20 -O 10  --trim-n  -m 30  --max-n 0.1  -o  ../01_cutadapt/${base}_trimmed.R1.fq.gz  -p ../01_cutadapt/${base}_trimmed.R2.fq.gz   ${base}_1.fq.gz  ${base}_2.fq.gz 
done

