for i in `ls *_sort_uni.bam`
do 
base=$(basename $i "_sort_uni.bam")
java -jar ./picard-tools-2.2.4/picard.jar MarkDuplicates REMOVE_DUPLICATES=true   I=${base}_sort_uni.bam   O=./${base}_rmdup_picard.bam  M=./${base}_rmdup_picard.txt &
done
