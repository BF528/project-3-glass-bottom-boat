### part 1 of the analyst code
## analyst: Marina Natividad Avila
## marinant

#BiocManager::install(c("rat2302.db"))
library(limma)
library(rat2302.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(ggplot2)

samples <- read.csv('/project/bf528/project_3/groups/group_2_mic_info.csv',as.is=TRUE)

rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',as.is=TRUE,header=TRUE,row.names=1,)

# BETA-NAPHTHOFLAVONE
# ahr
rma.subset_beta <- rma[paste0('X',samples$array_id
                      [samples$vehicle=="CMC_.5_%"])]
design_beta <- model.matrix(~factor
              (samples$chemical[samples$vehicle=="CMC_.5_%"],
              levels=c('Control','BETA-NAPHTHOFLAVONE')))
colnames(design_beta) <- c('Intercept','BETA-NAPHTHOFLAVONE')

# ECONAZOLE
# car/pxr
rma.subset_eco <- rma[paste0('X',samples$array_id
                     [samples$vehicle=="CORN_OIL_100_%"])]
design_eco <- model.matrix(~factor
              (samples$chemical[samples$vehicle=="CORN_OIL_100_%"],
              levels=c('Control','ECONAZOLE')))
colnames(design_eco) <- c('Intercept','ECONAZOLE')

# THIOACETAMIDE
# cyto
rma.subset_thio <- rma[paste0('X',samples$array_id
                             [samples$vehicle=="SALINE_100_%"])]
design_thio <- model.matrix(~factor
                           (samples$chemical[samples$vehicle=="SALINE_100_%"],
                             levels=c('Control','THIOACETAMIDE')))
colnames(design_thio) <- c('Intercept','THIOACETAMIDE')


###
# limma
###

# BETA-NAPHTHOFLAVONE
fit_beta <- lmFit(rma.subset_beta, design_beta)
fit_beta <- eBayes(fit_beta)
t_beta <- topTable(fit_beta, coef=2, n=nrow(rma.subset_beta), adjust='BH')
t_beta <- t_beta[order(t_beta$P.Value),]
write.csv(t_beta,'ahr_limma_results.csv')

# ECONAZOLE
fit_eco <- lmFit(rma.subset_eco, design_eco)
fit_eco <- eBayes(fit_eco)
t_eco <- topTable(fit_eco, coef=2, n=nrow(rma.subset_eco), adjust='BH')
t_eco <- t_eco[order(t_eco$P.Value),]
write.csv(t_eco,'car_limma_results.csv')

# THIOACETAMIDE
fit_thio <- lmFit(rma.subset_thio, design_thio)
fit_thio <- eBayes(fit_thio)
t_thio <- topTable(fit_thio, coef=2, n=nrow(rma.subset_thio), adjust='BH')
t_thio <- t_thio[order(t_thio$P.Value),]
write.csv(t_thio,'cyto_limma_results.csv')

###
# Filtering and Sorting
###

t_beta_padj <- subset(t_beta, adj.P.Val < 0.05)
write.csv(t_beta_padj,'ahr_limma_padj.csv')
t_eco_padj <- subset(t_eco, adj.P.Val < 0.05)
write.csv(t_eco_padj,'car_limma_padj.csv')
t_thio_padj <- subset(t_thio, adj.P.Val < 0.05)
write.csv(t_thio_padj,'cyto_limma_padj.csv')

###
# Top 10 DE genes
###

top_beta <- as.data.frame(t_beta[1:10,]) 
top_beta <- top_beta %>% 
  rownames_to_column(var = "probeid")

top_eco <- as.data.frame(t_eco[1:10,]) 
top_eco <- top_eco %>% 
  rownames_to_column(var = "probeid")

top_thio <- as.data.frame(t_thio[1:10,]) 
top_thio <- top_thio %>% 
  rownames_to_column(var = "probeid")

m <- keys(rat2302.db,keytype="PROBEID")
m.col <- select(rat2302.db, keys=m, columns=c("SYMBOL", "GENENAME"), 
                keytype="PROBEID")

top_beta_names <- merge(top_beta, m.col, by.x="probeid", by.y ="PROBEID")
top_beta_names <- top_beta_names[, c(1,8,9,5:7,2:4)]
write.csv(top_beta_names, "top_ahr.csv")
top_eco_names <- merge(top_eco, m.col, by.x="probeid", by.y ="PROBEID")
top_eco_names <- top_eco_names[, c(1,8,9,5:7,2:4)]
write.csv(top_eco_names, "top_car.csv")
top_thio_names <- merge(top_thio, m.col, by.x="probeid", by.y ="PROBEID")
top_thio_names <- top_thio_names[, c(1,8,9,5:7,2:4)]
write.csv(top_thio_names, "top_cyto.csv")

###
# histograms and scatter plots
###

t_beta_padj <- t_beta_padj %>% 
  rownames_to_column(var = "probeid")
t_eco_padj <- t_eco_padj %>% 
  rownames_to_column(var = "probeid")
t_thio_padj <- t_thio_padj %>% 
  rownames_to_column(var = "probeid")

pdf("combinedplots2.pdf")
par(mfrow=c(2,3))

s1 <- ggplot(t_beta_padj, aes(y=P.Value, x=logFC)) + geom_point()
s1 <- s1 + labs(x = "log FC", y = "P Value", title = "AhR")
print(s1)

s2 <- ggplot(t_eco_padj, aes(y=P.Value, x=logFC)) + geom_point()
s2 <- s2 + labs(x = "log FC", y = "P Value", title = "CAR/PXR")
print(s2)

s3 <- ggplot(t_thio_padj, aes(y=P.Value, x=logFC)) + geom_point()
s3 <- s3 + labs(x = "log FC", y = "P Value", title = "Cytotoxic")
print(s3)

h1 <- ggplot(t_beta_padj, aes(y=logFC, x=probeid)) + geom_col()
h1 <- h1 + labs(x = "Genes", y = "log FC", title = "AhR")  
print(h1)

h2 <- ggplot(t_eco_padj, aes(y=logFC, x=probeid)) + geom_col()
h2 <- h2 + labs(x = "Genes", y = "log FC", title = "CAR/PXR")  
print(h2)

h3 <- ggplot(t_beta_padj, aes(y=logFC, x=probeid)) + geom_col()
h3 <- h3 + labs(x = "Genes", y = "log FC", title = "Cytotoxic")  
print(h3)

dev.off()
