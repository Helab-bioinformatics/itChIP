################
author:Qiang Shi
################

#### Receiver operating characteristic (ROC) curve 
#### delete information about call-peak in peak.xls
### sed '1,28d' 100-it_broad_q1_peaks.xls > 100-it_broad_q1_peaks.bed
### sed '1,28d' 500-it_broad_q1_peaks.xls > 500-it_broad_q1_peaks.bed
### sed '1,28d' 10E4-it_broad_q1_peaks.xls > 10E4-it_broad_q1_peaks.bed
### sed '1,28d' 1M-it_broad_q1_peaks.xls > 1M-it_broad_q1_peaks.bed
### sed '1,28d' SRR058074_broad_peaks.xls > SRR058074_broad_peaks.bed



#### reserve -log10(qvalue)
## awk -v OFS='\t' '{print $1,$2,$3,$8,$9}' 100-it_broad_q1_peaks.bed > 100-it_broad_q1_qvalue_peaks.bed
## awk -v OFS='\t' '{print $1,$2,$3,$8,$9}' 500-it_broad_q1_peaks.bed > 500-it_broad_q1_qvalue_peaks.bed
## awk -v OFS='\t' '{print $1,$2,$3,$8,$9}' 10E4-it_broad_q1_peaks.bed > 10E4-it_broad_q1_qvalue_peaks.bed
## awk -v OFS='\t' '{print $1,$2,$3,$8,$9}' 1M-it_broad_q1_peaks.bed > 1M-it_broad_q1_qvalue_peaks.bed
## awk -v OFS='\t' '{print $1,$2,$3,$8,$9}' SRR058074_broad_peaks.bed > SRR058074_broad_qvalue_peaks.bed

#### To obtain peaks overlapping with TSS±5kb regions, excute following system command
### bedtools intersect -a 100-it_broad_q1_qvalue_peaks.bed -b mm10_TSS_fl5kb_all+.bed6 -wa > 100-it_broad_q1_qvalue_TSS±5kb_peaks.bed
### bedtools intersect -a 500-it_broad_q1_qvalue_peaks.bed -b mm10_TSS_fl5kb_all+.bed6 -wa > 500-it_broad_q1_qvalue_TSS±5kb_peaks.bed
### bedtools intersect -a 10E4-it_broad_q1_qvalue_peaks.bed -b mm10_TSS_fl5kb_all+.bed6 -wa > 10E4-it_broad_q1_qvalue_TSS±5kb_peaks.bed
### bedtools intersect -a 1M-it_broad_q1_qvalue_peaks.bed -b mm10_TSS_fl5kb_all+.bed6 -wa > 1M-it_broad_q1_qvalue_TSS±5kb_peaks.bed
### bedtools intersect -a soni-1M_broad_q0.001_peaks.broadPeak -b mm10_TSS_fl5kb_all+.bed6 -wa > soni-1M_broad_q0.001_TSS±5kb_peaks.broadPeak
### bedtools intersect -a SRR058074_broad_qvalue_peaks.bed -b mm10_TSS_fl5kb_all+.bed6 -wa > SRR058074_TSS±5kb_peaks.bed

#### obtain unique peaks
### sort -k2n 100-it_broad_q1_qvalue_TSS±5kb_peaks.bed | uniq > 100-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed
### sort -k2n 500-it_broad_q1_qvalue_TSS±5kb_peaks.bed | uniq > 500-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed
### sort -k2n 10E4-it_broad_q1_qvalue_TSS±5kb_peaks.bed | uniq > 10E4-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed
### sort -k2n 1M-it_broad_q1_qvalue_TSS±5kb_peaks.bed | uniq > 1M-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed
### sort -k2n soni-1M_broad_q0.001_TSS±5kb_peaks.broadPeak | uniq > soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak
### sort -k2n SRR058074_TSS±5kb_peaks.bed | uniq > SRR058074_TSS±5kb_uniq_peaks.bed

#### define condition positive
### soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak
condition.positive <- read.table("./soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak", header=FALSE)

#### define condition negative
### bedtools subtract -a mm10_TSS_fl5kb_all+.bed6 -b soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak > soni-1M_broad_q0.001_TSS±5kb_negative_region.bed
condition.negative <- read.table("./soni-1M_broad_q0.001_TSS±5kb_negative_region.bed", header=FALSE)

#### candidate positive peak
candidate.positive.100 <- read.table("./100-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed")
candidate.positive.500 <- read.table("./500-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed")
candidate.positive.10E4 <- read.table("./10E4-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed")
candidate.positive.1M <- read.table("./1M-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed")
candidate.positive.nano <- read.table("./SRR058074_TSS±5kb_uniq_peaks.bed")

colnames(candidate.positive.100) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.positive.500) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.positive.10E4) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.positive.1M) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.positive.nano) <- c("chrom","chromStart","chromEnd","-log10(q)","name")

cutoff.100 <- sort(unique(candidate.positive.100$`-log10(q)`), decreasing=TRUE)
cutoff.500 <- sort(unique(candidate.positive.500$`-log10(q)`), decreasing=TRUE)
cutoff.10E4 <- sort(unique(candidate.positive.10E4$`-log10(q)`), decreasing=TRUE)
cutoff.1M <- sort(unique(candidate.positive.1M$`-log10(q)`), decreasing=TRUE)
cutoff.nano <- sort(unique(candidate.positive.nano$`-log10(q)`), decreasing=TRUE)


#### candidate true positive peak
### bedtools intersect -a 100-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak -wa > 100-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed
### bedtools intersect -a 500-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak -wa > 500-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed
### bedtools intersect -a 10E4-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak -wa > 10E4-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed
### bedtools intersect -a 1M-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak -wa > 1M-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed
### bedtools intersect -a SRR058074_TSS±5kb_uniq_peaks.bed -b soni-1M_broad_q0.001_TSS±5kb_uniq_peaks.broadPeak -wa > nano_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed

#### unique peak
### sort -k2n 100-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed | uniq > 100-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed
### sort -k2n 500-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed | uniq > 500-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed
### sort -k2n 10E4-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed | uniq > 10E4-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed
### sort -k2n 1M-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed | uniq > 1M-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed
### sort -k2n nano_broad_q1_qvalue_TSS±5kb_candidate.true.positive_peaks.bed | uniq > nano_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed

candidate.true.positive.100 <- read.table("./100-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed", header=FALSE)
candidate.true.positive.500 <- read.table("./500-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed", header=FALSE)
candidate.true.positive.10E4 <- read.table("./10E4-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed", header=FALSE)
candidate.true.positive.1M <- read.table("./1M-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed", header=FALSE)
candidate.true.positive.nano <- read.table("./nano_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed", header=FALSE)

colnames(candidate.true.positive.100) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.true.positive.500) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.true.positive.10E4) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.true.positive.1M) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.true.positive.nano) <- c("chrom","chromStart","chromEnd","-log10(q)","name")

#### candidate false positive peak
### bedtools subtract -a 100-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b 100-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed -wa > 100-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed
### bedtools subtract -a 500-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b 500-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed -wa > 500-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed
### bedtools subtract -a 10E4-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b 10E4-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed -wa > 10E4-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed
### bedtools subtract -a 1M-it_broad_q1_qvalue_TSS±5kb_uniq_peaks.bed -b 1M-it_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed -wa > 1M-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed
### bedtools subtract -a SRR058074_TSS±5kb_uniq_peaks.bed -b nano_broad_q1_qvalue_TSS±5kb_candidate.true.positive_uniq_peaks.bed -wa > nano_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed

candidate.false.positive.100 <- read.table("./100-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed", header=FALSE)
candidate.false.positive.500 <- read.table("./500-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed", header=FALSE)
candidate.false.positive.10E4 <- read.table("./10E4-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed", header=FALSE)
candidate.false.positive.1M <- read.table("./1M-it_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed", header=FALSE)
candidate.false.positive.nano <- read.table("./nano_broad_q1_qvalue_TSS±5kb_candidate.false.positive_peaks.bed", header=FALSE)

colnames(candidate.false.positive.100) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.false.positive.500) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.false.positive.10E4) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.false.positive.1M) <- c("chrom","chromStart","chromEnd","-log10(q)","name")
colnames(candidate.false.positive.nano) <- c("chrom","chromStart","chromEnd","-log10(q)","name")

###### 100
######
ROC.100 <- as.data.frame(matrix(NA,nrow=length(cutoff.100),ncol=3))
colnames(ROC.100) <- c("cutoff", "tp", "fp")

for(i in 1:length(cutoff.100)){
  ROC.100$cutoff[i] <- cutoff.100[i]
  
  true.positive.100 <- candidate.true.positive.100[candidate.true.positive.100$`-log10(q)`>=cutoff.100[i], ]
  ROC.100$tp[i] <- nrow(true.positive.100)
  
  false.positive.100 <- candidate.false.positive.100[candidate.false.positive.100$`-log10(q)`>=cutoff.100[i], ]
  ROC.100$fp[i] <- nrow(false.positive.100)
}

ROC.100$tpr <- ROC.100$tp/nrow(condition.positive)
ROC.100$fpr <- ROC.100$fp/nrow(condition.negative)

####### AUC
tpr.100 <- c(0,sort(ROC.100$tpr, decreasing=FALSE),1)
fpr.100 <- c(0,sort(ROC.100$fpr, decreasing=FALSE),1)
AUC.100 <- c()
for(i in 2:length(fpr.100)){
  S = ( fpr.100[i] - fpr.100[i-1] ) * ( tpr.100[i-1] + (tpr.100[i] - tpr.100[i-1])/2 )
  AUC.100 <- sum(AUC.100, S)
}

#####

###### 500
#######
ROC.500 <- as.data.frame(matrix(NA,nrow=length(cutoff.500),ncol=3))
colnames(ROC.500) <- c("cutoff", "tp", "fp")

for(i in 1:length(cutoff.500)){
  ROC.500$cutoff[i] <- cutoff.500[i]
  
  true.positive.500 <- candidate.true.positive.500[candidate.true.positive.500$`-log10(q)`>=cutoff.500[i], ]
  ROC.500$tp[i] <- nrow(true.positive.500)
  
  false.positive.500 <- candidate.false.positive.500[candidate.false.positive.500$`-log10(q)`>=cutoff.500[i], ]
  ROC.500$fp[i] <- nrow(false.positive.500)
}

ROC.500$tpr <- ROC.500$tp/nrow(condition.positive)
ROC.500$fpr <- ROC.500$fp/nrow(condition.negative)

####### AUC
tpr.500 <- c(0,sort(ROC.500$tpr, decreasing=FALSE),1)
fpr.500 <- c(0,sort(ROC.500$fpr, decreasing=FALSE),1)
AUC.500 <- c()
for(i in 2:length(fpr.500)){
  S = ( fpr.500[i] - fpr.500[i-1] ) * ( tpr.500[i-1] + (tpr.500[i] - tpr.500[i-1])/2 )
  AUC.500 <- sum(AUC.500, S)
}


#####

###### 10E4
######
ROC.10E4 <- as.data.frame(matrix(NA,nrow=length(cutoff.10E4),ncol=3))
colnames(ROC.10E4) <- c("cutoff", "tp", "fp")

for(i in 1:length(cutoff.10E4)){
  ROC.10E4$cutoff[i] <- cutoff.10E4[i]
  
  true.positive.10E4 <- candidate.true.positive.10E4[candidate.true.positive.10E4$`-log10(q)`>=cutoff.10E4[i], ]
  ROC.10E4$tp[i] <- nrow(true.positive.10E4)
  
  false.positive.10E4 <- candidate.false.positive.10E4[candidate.false.positive.10E4$`-log10(q)`>=cutoff.10E4[i], ]
  ROC.10E4$fp[i] <- nrow(false.positive.10E4)
}

ROC.10E4$tpr <- ROC.10E4$tp/nrow(condition.positive)
ROC.10E4$fpr <- ROC.10E4$fp/nrow(condition.negative)

####### AUC
tpr.10E4 <- c(0,sort(ROC.10E4$tpr, decreasing=FALSE),1)
fpr.10E4 <- c(0,sort(ROC.10E4$fpr, decreasing=FALSE),1)
AUC.10E4 <- c()
for(i in 2:length(fpr.10E4)){
  S = ( fpr.10E4[i] - fpr.10E4[i-1] ) * ( tpr.10E4[i-1] + (tpr.10E4[i] - tpr.10E4[i-1])/2 )
  AUC.10E4 <- sum(AUC.10E4, S)
}



#####

###### 1M
######
ROC.1M <- as.data.frame(matrix(NA,nrow=length(cutoff.1M),ncol=3))
colnames(ROC.1M) <- c("cutoff", "tp", "fp")

for(i in 1:length(cutoff.1M)){
  ROC.1M$cutoff[i] <- cutoff.1M[i]
  
  true.positive.1M <- candidate.true.positive.1M[candidate.true.positive.1M$`-log10(q)`>=cutoff.1M[i], ]
  ROC.1M$tp[i] <- nrow(true.positive.1M)
  
  false.positive.1M <- candidate.false.positive.1M[candidate.false.positive.1M$`-log10(q)`>=cutoff.1M[i], ]
  ROC.1M$fp[i] <- nrow(false.positive.1M)
}

ROC.1M$tpr <- ROC.1M$tp/nrow(condition.positive)
ROC.1M$fpr <- ROC.1M$fp/nrow(condition.negative)

####### AUC
tpr.1M <- c(0,sort(ROC.1M$tpr, decreasing=FALSE),1)
fpr.1M <- c(0,sort(ROC.1M$fpr, decreasing=FALSE),1)
AUC.1M <- c()
for(i in 2:length(fpr.1M)){
  S = ( fpr.1M[i] - fpr.1M[i-1] ) * ( tpr.1M[i-1] + (tpr.1M[i] - tpr.1M[i-1])/2 )
  AUC.1M <- sum(AUC.1M, S)
}

###### nano
######
ROC.nano <- as.data.frame(matrix(NA,nrow=length(cutoff.nano),ncol=3))
colnames(ROC.nano) <- c("cutoff", "tp", "fp")

for(i in 1:length(cutoff.nano)){
  ROC.nano$cutoff[i] <- cutoff.nano[i]
  
  true.positive.nano <- candidate.true.positive.nano[candidate.true.positive.nano$`-log10(q)`>=cutoff.nano[i], ]
  ROC.nano$tp[i] <- nrow(true.positive.nano)
  
  false.positive.nano <- candidate.false.positive.nano[candidate.false.positive.nano$`-log10(q)`>=cutoff.nano[i], ]
  ROC.nano$fp[i] <- nrow(false.positive.nano)
}

ROC.nano$tpr <- ROC.nano$tp/nrow(condition.positive)
ROC.nano$fpr <- ROC.nano$fp/nrow(condition.negative)

####### AUC
tpr.nano <- c(0,sort(ROC.nano$tpr, decreasing=FALSE),1)
fpr.nano <- c(0,sort(ROC.nano$fpr, decreasing=FALSE),1)
AUC.nano <- c()
for(i in 2:length(fpr.nano)){
  S = ( fpr.nano[i] - fpr.nano[i-1] ) * ( tpr.nano[i-1] + (tpr.nano[i] - tpr.nano[i-1])/2 )
  AUC.nano <- sum(AUC.nano, S)
}

######

pdf("./fig_ROC.pdf")

plot(fpr.100, tpr.100, type="l", xlim=c(0,1), ylim=c(0,1), col=1,lwd=4, 
     xlab="False Positive Rate",ylab="True Positive Rate", main="ROC of TSS±5kb")#,axes=FALSE, ann=FALSE, frame.plot=TRUE)
par(new=TRUE)
plot(fpr.500, tpr.500, type="l", xlim=c(0,1), ylim=c(0,1), col=2, lwd=4, axes=FALSE, ann=FALSE)
par(new=TRUE)
plot(fpr.10E4, tpr.10E4, type="l", xlim=c(0,1), ylim=c(0,1), col=3, lwd=4, axes=FALSE, ann=FALSE)
par(new=TRUE)
plot(fpr.1M, tpr.1M, type="l", xlim=c(0,1), ylim=c(0,1), col=4, lwd=4, axes=FALSE, ann=FALSE)
par(new=TRUE)
plot(fpr.nano, tpr.nano, type="l", xlim=c(0,1), ylim=c(0,1), col=5, lwd=4, axes=FALSE, ann=FALSE)
points(c(0,1),c(0,1),type="l",col="black", lty=2, lwd=4)


dev.off()




