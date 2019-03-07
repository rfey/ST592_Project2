if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pasilla", version = "3.8")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm", version = "3.8")


library("DESeq2"); packageVersion("DESeq2")
library(dplyr)
library(data.table)
library("BiocParallel")
library('purrr')
library('tibble')
library("pasilla")
library(apeglm)

data <- read.table("~/Documents/ST 592/project 2/input.txt", head = T)
head(data, 2)

filteredCounts <- data.frame()
for (row in 1:nrow(data)){
  if (median(as.numeric(data[row, 1:12])) >1 || median(as.numeric(data[row,13:24]))>1){
    filteredCounts <- rbind(filteredCounts, data[row,])
  }
}

dim(filteredCounts)
cts <- as.matrix(filteredCounts)
# filteredCounts <- data.frame()
# for (row in 1:nrow(cts)){
#    if (median(as.numeric(cts[row, 1:12])) >1 || median(as.numeric(cts[row,13:24]))>1){
#      filteredCounts <- rbind(filteredCounts, cts[row,])
#    }
#}

coldata <- read.csv("~/Documents/ST 592/project 2/fruitfly.csv", row.names=1)
head(coldata, 2)
all(rownames(coldata) == colnames(cts))
countData<- round(cts)
dds <- DESeqDataSetFromMatrix(countData, colData = coldata, design = ~age) 
dds

dds$age <- factor(dds$age, levels = c("young", "old"))
dds <- DESeq(dds)
res <- results(dds)
res

dds$time <- factor(dds$time, levels = c("0", "4", "8", "12", "16", "20"))
design(dds) <-  formula(~age*time)
dds <- DESeq(dds)
res <- results(dds)
res

resultsNames(dds)


alpha = 0.05
res1 <- results(dds, name = "age_old_vs_young", cooksCutoff = FALSE)
res1


if (FC == "NA"){
  print("No fold change threshold chosen.")
  qlf <- glmQLFTest(res1, coef=2) # QL F test with no log fold change threshold imposed
  BH <- decideTestsDGE(qlf, adjust.method='BH',p.value=as.numeric(alpha),lfc=0)
  toWrite <- cbind(rownames(qlf$table),qlf$table)
  colnames(toWrite) <- append(colnames(qlf$table),'mirID',after=0)
  toWrite <- cbind(toWrite,BH[,1])
  colnames(toWrite)[6] <- "DE"
  write.table(toWrite, file = outfile, sep='\t', quote=FALSE, row.names=FALSE)
} else {
  print(paste("Using fold change threshold",FC))
  FC <- as.numeric(FC)
  tr <- glmTreat(fit, coef=2, lfc=log2(FC)) # as above, but with log fold change threshold imposed
  BH <- decideTestsDGE(tr,adjust.method='BH',p.value=as.numeric(alpha),lfc=0)
  toWrite <- cbind(rownames(tr$table),tr$table)
  colnames(toWrite) <- c('mirID',colnames(tr$table))
  toWrite <- cbind(toWrite,BH[,1])
  colnames(toWrite)[6] <- "DE"
  write.table(toWrite, file=outfile, sep='\t', quote=FALSE, row.names=FALSE)
}

res1 = res1[which(res1$padj < alpha),]
res1 = cbind(as(res1, "data.frame"))
res1$comparison <- "age_old_vs_young"
write.csv(res1, file = "~/Documents/ST 592/project 2/DESeq1.csv")





