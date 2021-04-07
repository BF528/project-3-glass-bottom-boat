SRR1177966 = read.table("feature_counts_SRR1177966.txt", sep = "")
SRR1177966 = SRR1177966[-1,]
SRR1177966 = SRR1177966[, 7]

SRR1177969 = read.table("feature_counts_SRR1177969.txt", sep = "")
SRR1177969 = SRR1177969[-1,]
SRR1177969 = SRR1177969[, 7]

SRR1177970 = read.table("feature_counts_SRR1177970.txt", sep = "")
SRR1177970 = SRR1177970[-1,]
SRR1177970 = SRR1177970[, 7]

SRR1177993 = read.table("feature_counts_SRR1177993.txt", sep = "")
SRR1177993 = SRR1177993[-1,]
SRR1177993 = SRR1177993[, 7]

SRR1177994 = read.table("feature_counts_SRR1177994.txt", sep = "")
SRR1177994 = SRR1177994[-1,]
SRR1177994 = SRR1177994[, 7]

SRR1177995 = read.table("feature_counts_SRR1177995.txt", sep = "")
SRR1177995 = SRR1177995[-1,]
SRR1177995 = SRR1177995[, 7]

SRR1177998 = read.table("feature_counts_SRR1177998.txt", sep = "")
SRR1177998 = SRR1177998[-1,]
SRR1177998 = SRR1177998[, 7]

SRR1178001 = read.table("feature_counts_SRR1178001.txt", sep = "")
SRR1178001 = SRR1178001[-1,]
SRR1178001 = SRR1178001[, 7]

SRR1178003 = read.table("feature_counts_SRR1178003.txt", sep = "")
SRR1178003 = SRR1178003[-1,]
SRR1178003 = SRR1178003[, 7]

geneid = read.table("feature_counts_SRR1177966.txt", sep = "")
geneid = geneid[, 1]
geneid = geneid[-1]

counts = data.frame(geneid, SRR1177966, SRR1177969, SRR1177970, SRR1177993, SRR1177994, SRR1177995, SRR1177998, SRR1178001, SRR1178003)
write.csv(counts, "counts.csv")

boxplts = counts
boxplts[boxplts==0] = 1
boxplot(as.integer(paste(boxplts$SRR1177966)), 
        as.integer(paste(boxplts$SRR1177969)), 
        as.integer(paste(boxplts$SRR1177970)), 
        as.integer(paste(boxplts$SRR1177993)),
        as.integer(paste(boxplts$SRR1177994)),
        as.integer(paste(boxplts$SRR1177995)),
        as.integer(paste(boxplts$SRR1177998)),
        as.integer(paste(boxplts$SRR1178001)),
        as.integer(paste(boxplts$SRR1178003)),
        main = "Count Distributions", names = names(counts)[-1], log = "y", xlab = "Samples", ylab = "Feature Counts")