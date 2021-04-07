#histograms

ggplot(ahr_deseq_padj, aes(ahr_deseq_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color="black") +
  xlab("log2 Fold Change") +
  theme_classic() +
  labs(title="Histogram of log2 Fold Change in Expression for AhR Samples")

ggplot(car_deseq_padj, aes(car_deseq_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color="black") +
  xlab("log2 Fold Change") +
  theme_classic() +
  labs(title="Histogram of log2 Fold Change in Expression for CAR/PXR Samples")

ggplot(cyto_deseq_padj, aes(cyto_deseq_padj$log2FoldChange)) +
  geom_histogram(bins = 50, binwidth = 0.5, color="black") +
  xlab("log2 Fold Change") +
  theme_classic() +
  labs(title="Histogram of log2 Fold Change in Expression for Cytotoxic Samples")

#scatter plots
ahr_scatter = ahr_deseq[order(ahr_deseq$log2FoldChange),]
ahr_scatter = data.frame(ahr_scatter)
ahrlabels = rbind(head(ahr_scatter), tail(ahr_scatter))
ahrlabels$gene = rownames(ahrlabels)

ggplot(ahr_deseq, aes(log2FoldChange, -log(pvalue))) + 
  geom_point(color="darkgray") +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=ahrlabels, col = "red")+
  geom_text_repel(data=ahrlabels, aes(label=gene))+
  theme_classic() +
  labs(title="AhR Samples Expression log2 fold change against p values")

car_scatter = car_deseq[order(car_deseq$log2FoldChange),]
car_scatter = data.frame(car_scatter)
carlabels = rbind(head(car_scatter), tail(car_scatter))
carlabels$gene = rownames(carlabels)

ggplot(car_deseq, aes(log2FoldChange, -log(pvalue))) + 
  geom_point(color="darkgray") +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=carlabels, col = "red")+
  geom_text_repel(data=carlabels, aes(label=gene))+
  theme_classic() +
  labs(title="CAR/PXR Samples Expression log2 fold change against p values")

cyto_scatter = cyto_deseq[order(cyto_deseq$log2FoldChange),]
cyto_scatter = data.frame(cyto_scatter)
cytolabels = rbind(head(cyto_scatter), tail(cyto_scatter))
cytolabels$gene = rownames(cytolabels)

ggplot(cyto_deseq, aes(log2FoldChange, -log(pvalue))) + 
  geom_point(color="darkgray") +
  geom_hline(yintercept = -log(0.05), colour="red", linetype="dashed") +
  geom_point(data=ahrlabels, col = "red")+
  geom_text_repel(data=cytolabels, aes(label=gene))+
  theme_classic() +
  labs(title="Cytotoxic Samples Expression log2 fold change against p values")

