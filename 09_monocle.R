library(monocle)
#input count matrix
load("scRNA-seq_pseudotime.RData")
load("scitChIP-seq_pseudotime.RData")

pheno <- data.frame(colnames(df.RNA),stage)
names(pheno) = c("samplename", "stage")
rownames(pheno) <- colnames(df.RNA)

feature <- data.frame(rownames(df.RNA))
colnames(feature) <- "gene_short_name"
rownames(feature) <- rownames(df.RNA)

pd <- new("AnnotatedDataFrame", data = pheno)
fd <- new("AnnotatedDataFrame", data = feature)
heart <- newCellDataSet(as.matrix(df.RNA), phenoData = pd, featureData = fd,lowerDetectionLimit = 0.5,expressionFamily =tobit())

#Estimate size factors and dispersions
heart <- estimateSizeFactors(heart)
heart <- estimateDispersions(heart)

heart <-setOrderingFilter(heart,ordering_genes = test.gene.USE)
heart_pseudo <- reduceDimension(heart, max_components = 2,norm_method ="log",
                               method = 'DDRTree')
heart_pseudo <- orderCells(heart_pseudo)

plot_cell_trajectory(heart_pseudo, 1, 2, color_by = 'Pseudotime', branches=c(2,3),show_branch_points = F, cell_size = 3) +
  scale_colour_gradient2(name = 'Pseudotime', mid='blue', high = 'yellow')
  

heart_pseudo2 <- orderCells(heart_pseudo,root_state = 3)
#which(heart_pseudo@phenoData@data$State == 2)



###identify DEG
BEAM_res <- BEAM(heart_pseudo, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(heart_pseudo[row.names(subset(BEAM_res,
                                                  qval < 1e-13)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
#write.table(subset(BEAM_res,qval < 1e-13),"DEG.pseudotime.txt",row.names = T,col.names = T,quote = F,sep = "\t")

