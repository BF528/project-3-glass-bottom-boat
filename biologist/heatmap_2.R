source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")

library("ComplexHeatmap")

ahr <- read.csv(file = "C:/Users/b42na/Documents/BF528/project_3/heatmap/ahr_deseq_norm_counts.csv", stringsAsFactors = FALSE)
cyto <- read.csv(file = "C:/Users/b42na/Documents/BF528/project_3/heatmap/cyto_deseq_norm_counts.csv", stringsAsFactors = FALSE)
car <- read.csv(file = "C:/Users/b42na/Documents/BF528/project_3/heatmap/car_deseq_norm_counts.csv", stringsAsFactors = FALSE)
colnames(ahr)[2:7] <- paste(colnames(ahr)[2:7], "AhR")
colnames(cyto)[2:7] <- paste(colnames(cyto)[2:7], "Cytotoxic")
colnames(car)[2:7] <- paste(colnames(car)[2:7], "CAR/PXR")

cv <- function(row){sd(row)/mean(row)}
cv.cutoff <-1.2



m <- merge(ahr, cyto, by="Genes", all = TRUE)
m <- merge(m, car, by="Genes", all=TRUE)
rownames(m) <- m$Genes
m <- m[,-1]
m[is.na(m)] <- 0

sd(as.matrix(m))/mean(as.matrix(m))


m.m <- median(apply(m, 1,FUN=mean, na.rm=T))
m.mean <- m[apply(m, 1,FUN=mean, na.rm=T)>m.m,]

m.cov <- m.mean[apply(m.mean, 1, FUN=cv)>cv.cutoff,]

my_matrix <- as.matrix(m.cov)
my_matrix <- t(my_matrix)
fontsize <- .8

Heatmap(my_matrix,column_names_side = "top", row_names_side = "left", column_dend_height = unit(4, "cm"), row_dend_width  = unit(4, "cm"), column_names_gp = gpar(cex=fontsize))

