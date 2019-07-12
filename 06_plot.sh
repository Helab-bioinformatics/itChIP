#scatter plot
multiBigwigSummary BED-file -b *bw  -out results.npz --BED  mm10_RefSeq.bed  -p  10
plotCorrelation -in  results.npz -o correlation_scatter.svg -c pearson --removeOutliers -p scatterplot  --skipZeros  --colorMap inferno

#heatmap
computeMatrix reference-point --referencePoint center -p 10  -R ./bulk_sortregions.bed\
                                         -S *.bw\
                                         -a 5000 -b 5000  -bs 50   --sortRegions no   --missingDataAsZero   -q\
                                         -o   ./all_matrix.gz   --outFileNameMatrix  ./all_matrix.tab   --outFileSortedRegions  ./all_sortregions.bed &

plotHeatmap  -m    ./all_matrix.gz   -o  ./all_50bin.svg   --sortRegions  no     --zMin 0  --zMax 30    --colorList  darkblue yellow   --heatmapHeight 20  -x  all -y  signal   --plotTitle   all &

