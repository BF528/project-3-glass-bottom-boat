#load in control samples
control_counts = read.csv("control_counts.csv")
control_counts = as.data.frame(control_counts)
row.names(control_counts) = control_counts[,1]
control_counts = control_counts[,-1]

all_counts = cbind(counts, control_counts)

#splitting up groups
mdata = read.csv('toxgroup_2_rna_info.csv')

ahr = rbind(mdata[mdata$mode_of_action=="AhR",], mdata[mdata$mode_of_action=="Control" & mdata$vehicle=="CMC_.5_%",])
car = rbind(mdata[mdata$mode_of_action=="CAR/PXR",], mdata[mdata$mode_of_action=="Control" & mdata$vehicle=="CORN_OIL_100_%",])
cyto = rbind(mdata[mdata$mode_of_action=="Cytotoxic",], mdata[mdata$mode_of_action=="Control" & mdata$vehicle=="CORN_OIL_100_%",])

ahr_counts = cbind(all_counts$SRR1177998, all_counts$SRR1178001, all_counts$SRR1178003, all_counts$SRR1178030, all_counts$SRR1178040, all_counts$SRR1178056)
row.names(ahr_counts) <- rownames(all_counts)
ahr_counts = as.data.frame(ahr_counts)
AHR_counts <- transform(ahr_counts, SRR1177998 = as.numeric(SRR1177998), SRR1178001 = as.numeric(SRR1178001), SRR1178003 = as.numeric(SRR1178003), SRR1178030 = as.numeric(SRR1178030), SRR1178040 = as.numeric(SRR1178040), SRR1178056 = as.numeric(SRR1178056))
colnames(AHR_counts) <- ahr$Run

car_counts = cbind(all_counts$SRR1177993, all_counts$SRR1177994, all_counts$SRR1177995, all_counts$SRR1178024, all_counts$SRR1178035, all_counts$SRR1178045)
row.names(car_counts) <- rownames(all_counts)
car_counts = as.data.frame(car_counts)
#CAR_counts <- transform(car_counts, SRR1177993 = as.numeric(SRR1177993), SRR1177994 = as.numeric(SRR1177994), SRR1177995 = as.numeric(SRR1177995), SRR1178024 = as.numeric(SRR1178024), SRR1178035 = as.numeric(SRR1178035), SRR1178045 = as.numeric(SRR1178045))
#colnames(car_counts) <- car$Run

cyto_counts = cbind(all_counts$SRR1177966, all_counts$SRR1177969, all_counts$SRR1177970, all_counts$SRR1178024, all_counts$SRR1178035, all_counts$SRR1178045)
row.names(cyto_counts) <- rownames(all_counts)
cyto_counts = as.data.frame(cyto_counts)
#CYTO_counts <- transform(cyto_counts, SRR1177966 = as.numeric(SRR1177966), SRR1177969 = as.numeric(SRR1177969), SRR1177970 = as.numeric(SRR1177970), SRR1178024 = as.numeric(SRR1178024), SRR1178035 = as.numeric(SRR1178035), SRR1178045 = as.numeric(SRR1178045))
#colnames(cyto_counts) <- cyto$Run

# filter out rows that have any zeros for funzies
AHR_counts <- subset(AHR_counts,rowSums(AHR_counts==0)==0)
car_counts <- subset(car_counts,rowSums(car_counts==0)==0)
cyto_counts <- subset(cyto_counts,rowSums(cyto_counts==0)==0)


#AHR
# create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = AHR_counts,
  colData = ahr,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','AhR','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'ahr_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'ahr_deseq_norm_counts.csv')

#CAR
CAR_counts <- transform(car_counts, SRR1177993 = as.numeric(SRR1177993), SRR1177994 = as.numeric(SRR1177994), SRR1177995 = as.numeric(SRR1177995), SRR1178024 = as.numeric(SRR1178024), SRR1178035 = as.numeric(SRR1178035), SRR1178045 = as.numeric(SRR1178045))
colnames(CAR_counts) <- car$Run
#create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = CAR_counts,
  colData = car,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','CAR/PXR','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'car_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'car_deseq_norm_counts.csv')

##CYTO
CYTO_counts <- transform(cyto_counts, SRR1177966 = as.numeric(SRR1177966), SRR1177969 = as.numeric(SRR1177969), SRR1177970 = as.numeric(SRR1177970), SRR1178024 = as.numeric(SRR1178024), SRR1178035 = as.numeric(SRR1178035), SRR1178045 = as.numeric(SRR1178045))
colnames(CYTO_counts) <- cyto$Run
#create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = CYTO_counts,
  colData = cyto,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','Cytotoxic','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'cyto_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'cyto_deseq_norm_counts.csv')

#padjust
#reload results
ahr_deseq <- read.csv("ahr_deseq_results.csv")
ahr_deseq <- as.data.frame(ahr_deseq)
row.names(ahr_deseq) <- ahr_deseq[, 1]
ahr_deseq <- ahr_deseq[,-1]

car_deseq <- read.csv("car_deseq_results.csv")
car_deseq <- as.data.frame(car_deseq)
row.names(car_deseq) <- car_deseq[, 1]
car_deseq <- car_deseq[,-1]

cyto_deseq <- read.csv("cyto_deseq_results.csv")
cyto_deseq <- as.data.frame(cyto_deseq)
row.names(cyto_deseq) <- cyto_deseq[, 1]
cyto_deseq <- cyto_deseq[,-1]

#sort p-values
ahr_deseq_padj <- ahr_deseq[order(ahr_deseq$padj),]
ahr_deseq_padj <- ahr_deseq_padj[ahr_deseq_padj$padj < 0.05,]
ahr_deseq_padj <- ahr_deseq_padj[complete.cases(ahr_deseq_padj), ]

car_deseq_padj <- car_deseq[order(car_deseq$padj),]
car_deseq_padj <- car_deseq_padj[car_deseq_padj$padj < 0.05,]
car_deseq_padj <- car_deseq_padj[complete.cases(car_deseq_padj), ]

cyto_deseq_padj <- cyto_deseq[order(cyto_deseq$padj),]
cyto_deseq_padj <- cyto_deseq_padj[cyto_deseq_padj$padj < 0.05,]
cyto_deseq_padj <- cyto_deseq_padj[complete.cases(cyto_deseq_padj), ]

#create padj filtered csv files
write.csv(ahr_deseq_padj, file="ahr_deseq_results_padj.csv")
write.csv(car_deseq_padj, file="car_deseq_results_padj.csv")
write.csv(cyto_deseq_padj, file="cyto_deseq_results_padj.csv")

#top10
ahr_top10 <- ahr_deseq_padj[1:10,]
write.csv(ahr_top10, file="ahr_top10.csv")

car_top10 <- car_deseq_padj[1:10,]
write.csv(car_top10, file="car_top10.csv")

cyto_top10 <- cyto_deseq_padj[1:10,]
write.csv(cyto_top10, file="cyto_top10.csv")
