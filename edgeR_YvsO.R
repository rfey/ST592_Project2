# rmf 2.20.2019

library(edgeR)

# creates a string vector with all arguments
args = commandArgs(trailingOnly=TRUE)

usage <- 'Rscript edgeR.R <counts table> <log fold change threshold>'

if (length(args)!=2){
   stop(usage)
}

countsTable <- args[1]
FC <- args[2]

outfile <- paste('miRsDE_youngVsOld',FC,'FC.txt',sep='')

rawCounts <- read.delim(countsTable,row.names=1)

# remove miRs with mean expression under 1 for either age (under 1 ok for one but not both)
print(paste('Number of miRs before filtering:',nrow(rawCounts),sep=" "))
filteredCounts <- data.frame()
for (row in 1:nrow(rawCounts)){
    if (mean(as.numeric(rawCounts[row,1:12]))>1 || mean(as.numeric(rawCounts[row,13:24]))>1){
       filteredCounts <- rbind(filteredCounts, rawCounts[row,])
    }
}
print(paste('Number of miRs after filtering:',nrow(filteredCounts),sep=" "))

# not sure how much of this is actually necessary
sampleNames <- gsub("_AS","",colnames(data.frame(filteredCounts)))
libGroup <- gsub("Rep\\d","",sampleNames) # ex sample name: youngRep1_ZT0
ZT <- gsub(".*_","",libGroup)
age <- gsub("_ZT.*","",libGroup)
sampleInfo <- cbind(sampleNames,ZT,age)
sampleInfo <- as.data.frame(sampleInfo)
filteredCountsDGE <- DGEList(counts=filteredCounts,group=libGroup)

#group <- factor(paste(sampleInfo$age,sampleInfo$ZT,sep='_'))
#group <- cbind(sampleInfo,group=group)

# not sure if this is correct
design <- model.matrix(~ age + ZT + age:ZT)
#colnames(design) <- levels(age)

counts.norm.TMM <- calcNormFactors(filteredCountsDGE,method='TMM')
counts.norm <- cpm(counts.norm.TMM, normalized.lib.sizes=TRUE)

counts.disp <- estimateDisp(counts.norm.TMM, design, robust=TRUE)

fit <- glmQLFit(counts.disp, design, robust=TRUE)

# this weird thing is to make the 'makeContrasts' function take variables
#cmd <- paste("contrast <- makeContrasts(", condition1, "-", condition2, ", levels=design)")
#eval(parse(text = cmd))

contrast <- makeContrasts(age$young-age$old,levels=design)

if (FC == "NA"){
   print("No fold change threshold chosen.")
   qlf <- glmQLFTest(fit, contrast=contrast) # QL F test with no log fold change threshold imposed
   toWrite <- cbind(rownames(qlf$table),qlf$table)
   colnames(toWrite) <- cbind('mirID',colnames(qlf$table))
   write.table(toWrite, file = outfile, sep='\t', quote=FALSE, row.names=FALSE)
} else {
   print(paste("Using fold change threshold",FC))
   FC <- as.numeric(FC)
   tr <- glmTreat(fit, contrast=contrast, lfc=log2(FC)) # as above, but with log fold change threshold imposed
   toWrite <- cbind(rownames(tr$table),tr$table)
   colnames(toWrite) <- c('mirID',colnames(tr$table))
   write.table(toWrite, file=outfile, sep='\t', quote=FALSE, row.names=FALSE)
}