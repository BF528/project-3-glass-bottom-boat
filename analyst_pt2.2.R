## part 2 of the analyst code
## analyst: Marina Natividad Avila
## marinant

setwd("/projectnb2/bf528/users/glass_bottom_boat/project_3/analyst")

library(AnnotationDbi)
library(rat2302.db)
library(dplyr)
library(ggplot2)
#columns(rat2302.db)

# 1 ahr, 2 car, 3 cyto
# 1 beta, 2 eco, 3 thio
## microarray data
ds.limma1 <- read.csv("/projectnb/bf528/users/glass_bottom_boat/project_3/analyst/beta-naphthoflavone_limma_padj.csv",
                      header=TRUE)
ds.limma2 <- read.csv("/projectnb/bf528/users/glass_bottom_boat/project_3/analyst/econazole_limma_padj.csv",
                      header=TRUE)
ds.limma3 <- read.csv("/projectnb/bf528/users/glass_bottom_boat/project_3/analyst/thioacetamide_limma_padj.csv",
                      header=TRUE)
## rnaseq data
ds.deseq1 <- read.csv("/projectnb/bf528/users/glass_bottom_boat/project_3/programmer/programmer_deliverables/ahr_deseq_results_padj.csv",
                      header = TRUE)
ds.deseq2 <- read.csv("/projectnb/bf528/users/glass_bottom_boat/project_3/programmer/programmer_deliverables/car_deseq_results_padj.csv",
                      header = TRUE)
ds.deseq3 <- read.csv("/projectnb/bf528/users/glass_bottom_boat/project_3/programmer/programmer_deliverables/cyto_deseq_results_padj.csv",
                      header = TRUE)

ds.limma1 <- as.data.frame(ds.limma1)
ds.limma2 <- as.data.frame(ds.limma2)
ds.limma3 <- as.data.frame(ds.limma3)
ds.deseq1 <- as.data.frame(ds.deseq1)
ds.deseq2 <- as.data.frame(ds.deseq2)
ds.deseq3 <- as.data.frame(ds.deseq3)

# |logFC| < 1.5 and p-value < 0.05
ds.limma1.cutoff <- ds.limma1[which(abs(ds.limma1$logFC) >= 1.5 &
                                      ds.limma1$P.Value < 0.05), ] 
ds.limma2.cutoff <- ds.limma2[which(abs(ds.limma2$logFC) >= 1.5 &
                                      ds.limma2$P.Value < 0.05), ] 
ds.limma3.cutoff <- ds.limma3[which(abs(ds.limma3$logFC) >= 1.5 &
                                      ds.limma3$P.Value < 0.05), ] 
ds.deseq1.cutoff <- ds.deseq1[which(abs(ds.deseq1$log2FoldChange) >= 1.5 &
                                      ds.deseq1$pvalue < 0.05), ]
ds.deseq2.cutoff <- ds.deseq2[which(abs(ds.deseq2$log2FoldChange) >= 1.5 &
                                      ds.deseq2$pvalue < 0.05), ]
ds.deseq3.cutoff <- ds.deseq3[which(abs(ds.deseq3$log2FoldChange) >= 1.5 &
                                      ds.deseq3$pvalue < 0.05), ]

k <- keys(rat2302.db,keytype="PROBEID")
k.col <- select(rat2302.db, keys=k, columns=c("SYMBOL"), keytype = "PROBEID")
l <- keys(rat2302.db, keytype="REFSEQ")
l.col <- select(rat2302.db, keys=l, columns=c("SYMBOL"), keytype="REFSEQ")

## mapIds
#m <- keys(rat2302.db,keytype="PROBEID")
#m.map <- mapIds(rat2302.db, keys=m, column=c("SYMBOL"), keytype="PROBEID")
#r <- keys(rat2302.db,keytype="REFSEQ")
#r.map <- mapIds(rat2302.db, keys=r, column=c("SYMBOL"), keytype="REFSEQ")

## complete.limma
comp.limma1 <- merge(ds.limma1.cutoff, k.col, by.x="X", by.y ="PROBEID", 
                     all.x=TRUE) # left join
comp.limma1 <- comp.limma1[, c(1,8,2:7)]
comp.limma2 <- merge(ds.limma2.cutoff, k.col, by.x="X", by.y ="PROBEID", 
                     all.x=TRUE)
comp.limma2 <- comp.limma2[, c(1,8,2:7)]
comp.limma3 <- merge(ds.limma3.cutoff, k.col, by.x="X", by.y ="PROBEID", 
                     all.x=TRUE)
comp.limma3 <- comp.limma3[, c(1,8,2:7)]

# add column for FC direction positive=0 negative =1
#comp.limma$limmadir <- (ifelse(comp.limma$logFC>0, 0, 1))

comp.deseq1 <- merge(ds.deseq1.cutoff, l.col, by.x="X", by.y="REFSEQ",
                     all.x=TRUE)
comp.deseq1 <- comp.deseq1[, c(1,8,2:7)]
comp.deseq2 <- merge(ds.deseq2.cutoff, l.col, by.x="X", by.y="REFSEQ",
                     all.x=TRUE)
comp.deseq2 <- comp.deseq2[, c(1,8,2:7)]
comp.deseq3 <- merge(ds.deseq3.cutoff, l.col, by.x="X", by.y="REFSEQ",
                     all.x=TRUE)
comp.deseq3 <- comp.deseq3[, c(1,8,2:7)]

###
# Concordance Calculation
###

N <- 13079 # as stated by the authors
n_lim1 <- nrow(comp.limma1)
n_des1 <- nrow(comp.deseq1)
n_01 <- length(intersect(comp.deseq1$SYMBOL, comp.limma1$SYMBOL))
# background corrected intersection x
x1 <- (n_01 * N - n_lim1 * n_des1)/(n_01 + N - n_lim1 - n_des1)
# cross-platform concordance analysis
conc1 <- 2*x1/(n_lim1 + n_des1)

n_lim2 <- nrow(comp.limma2)
n_des2 <- nrow(comp.deseq2)
n_02 <- length(intersect(comp.deseq2$SYMBOL, comp.limma2$SYMBOL))
x2 <- (n_02 * N - n_lim2 * n_des2)/(n_02 + N - n_lim2 - n_des2)
conc2 <- 2*x2/(n_lim2 + n_des2)

n_lim3 <- nrow(comp.limma3)
n_des3 <- nrow(comp.deseq3)
n_03 <- length(intersect(comp.deseq3$SYMBOL, comp.limma3$SYMBOL))
x3 <- (n_03 * N - n_lim3 * n_des3)/(n_03 + N - n_lim3 - n_des3)
conc3 <- 2*x3/(n_lim3 + n_des3)

###
# PLOTS
###
plot.y <- (c(conc1, conc2, conc3)*100)
mplot.x <- c(n_lim1, n_lim2, n_lim3)
rplot.x <- c(n_des1, n_des2, n_des3)

png("micro_concord1.png")
micro <- plot(mplot.x, plot.y, xlab="Treatment Effect", ylab="Concordance of DEG (%)",
              main="Concordance vs Microarray")
dev.off()

png("rna_concord1.png")
rna <- plot(rplot.x, plot.y, xlab="Treatment Effect", ylab="Concordance of DEG (%)",
              main="Concordance vs RNA-seq")
dev.off()
#Between-platform concordance of DEGs against the number of DEGs identified by microarray

median1 <- median(ds.deseq1.cutoff$baseMean)
median2 <- median(ds.deseq2.cutoff$baseMean)
median3 <- median(ds.deseq3.cutoff$baseMean)

#ok so get the names of the genes and then check if they're also in the microarray data
# above mean _.1
med.deseq1 <- comp.deseq1[ which( comp.deseq1$baseMean > median1),]
med.deseq2 <- comp.deseq2[ which( comp.deseq2$baseMean > median2),]
med.deseq3 <- comp.deseq3[ which( comp.deseq3$baseMean > median3),]

n_des1.1 <- nrow(med.deseq1)
n_01.1 <- length(intersect(med.deseq1$SYMBOL, comp.limma1$SYMBOL))
x1.1 <- (n_01.1 * N - n_lim1 * n_des1.1)/(n_01.1 + N - n_lim1 - n_des1.1)
conc1.1 <- 2*x1.1/(n_lim1 + n_des1.1)

n_des2.1 <- nrow(med.deseq2)
n_02.1 <- length(intersect(med.deseq2$SYMBOL, comp.limma2$SYMBOL))
x2.1 <- (n_02.1 * N - n_lim2 * n_des2.1)/(n_02.1 + N - n_lim2 - n_des2.1)
conc2.1 <- 2*x2.1/(n_lim2 + n_des2.1)

n_des3.1 <- nrow(med.deseq3)
n_03.1 <- length(intersect(med.deseq3$SYMBOL, comp.limma3$SYMBOL))
x3.1 <- (n_03.1 * N - n_lim3 * n_des3.1)/(n_03.1 + N - n_lim3 - n_des3.1)
conc3.1 <- 2*x3.1/(n_lim3 + n_des3.1)

# below mean _.2
med.deseq1.2 <- comp.deseq1[ which( comp.deseq1$baseMean < median1),]
med.deseq2.2 <- comp.deseq2[ which( comp.deseq2$baseMean < median2),]
med.deseq3.2 <- comp.deseq3[ which( comp.deseq3$baseMean < median3),]

n_des1.2 <- nrow(med.deseq1.2)
n_01.2 <- length(intersect(med.deseq1.2$SYMBOL, comp.limma1$SYMBOL))
x1.2 <- (n_01.2 * N - n_lim1 * n_des1.2)/(n_01.2 + N - n_lim1 - n_des1.2)
conc1.2 <- 2*x1.2/(n_lim1 + n_des1.2)

n_des2.2 <- nrow(med.deseq2.2)
n_02.2 <- length(intersect(med.deseq2.2$SYMBOL, comp.limma2$SYMBOL))
x2.2 <- (n_02.2 * N - n_lim2 * n_des2.2)/(n_02.2 + N - n_lim2 - n_des2.2)
conc2.2 <- 2*x2.2/(n_lim2 + n_des2.2)

n_des3.2 <- nrow(med.deseq3.2)
n_03.2 <- length(intersect(med.deseq3.2$SYMBOL, comp.limma3$SYMBOL))
x3.2 <- (n_03.2 * N - n_lim3 * n_des3.2)/(n_03.2 + N - n_lim3 - n_des3.2)
conc3.2 <- 2*x3.2/(n_lim3 + n_des3.2)

## barplots
specie <- c(rep("AhR", 3), rep("CAR", 3), rep("Cyto", 3))
condition <- rep(c("overall", "above", "below"))
value <- c(conc1, conc2, conc3, conc1.1, conc2.1, conc3.1, conc1.2,
           conc2.2, conc3.2)
png("concplot.png")
plotdata <- data.frame(specie,condition,value)
f.plot <- ggplot(plotdata, aes(fill=condition, y=value, x=specie)) +
  geom_bar(position="dodge", stat="identity")
dev.off()


