if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pasilla", version = "3.8")

library("DESeq2"); packageVersion("DESeq2")
library(dplyr)
library(data.table)
library("BiocParallel")
library('purrr')
library('tibble')
library("pasilla")

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
design(dds) <-  ~age*time
dds <- DESeq(dds)
res <- results(dds)
res


resultsNames(dds)
resLFC <- lfcShrink(dds, coef="age_old_vs_young", type="apeglm")
resLFC