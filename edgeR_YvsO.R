# rmf 2.20.2019, last modified 3.5.2019

library(edgeR)

# creates a string vector with all arguments
args = commandArgs(trailingOnly=TRUE)

usage <- 'Rscript edgeR.R <counts table> <log fold change threshold> <alpha>\nLFC threshold optional, use "NA" for no threshold'

if (length(args)!=3){
   stop(usage)
}

countsTable <- args[1]
FC <- args[2]
alpha <- args[3]

outfile <- paste('miRsDE_youngVsOld',alpha,'alpha',FC,'FC.txt',sep='')

rawCounts <- read.delim(countsTable,row.names=1)

# remove miRs with median expression under 1 for either age (under 1 ok for one but not both)
# NOTE: this is filtering on raw counts!
print(paste('Number of miRs before filtering:',nrow(rawCounts),sep=" "))
filteredCounts <- data.frame()
for (row in 1:nrow(rawCounts)){
    if (median(as.numeric(rawCounts[row,1:12]))>1 || median(as.numeric(rawCounts[row,13:24]))>1){
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
group <- factor(paste(sampleInfo$age,sampleInfo$ZT,sep='_'))
group <- cbind(sampleInfo,group=group)

sampleInfo$age <- relevel(sampleInfo$age, ref="young")

# not sure if this is correct
design <- model.matrix(~age + ZT + age:ZT, data=sampleInfo)

#design <- model.matrix(~0+libGroup)
#colnames(design) <- levels()

counts.norm.TMM <- calcNormFactors(filteredCountsDGE,method='TMM')
counts.norm <- cpm(counts.norm.TMM, normalized.lib.sizes=TRUE)

counts.disp <- estimateDisp(counts.norm.TMM, design, robust=TRUE)

fit <- glmQLFit(counts.disp, design, robust=TRUE)

# this weird thing is to make the 'makeContrasts' function take variables
#cmd <- paste("contrast <- makeContrasts(", condition1, "-", condition2, ", levels=design)")
#eval(parse(text = cmd))

#contrast <- makeContrasts(age$young-age$old,levels=design)

if (FC == "NA"){
   print("No fold change threshold chosen.")
   qlf <- glmQLFTest(fit, coef=2) # QL F test with no log fold change threshold imposed
   BH <- decideTestsDGE(qlf,adjust.method='BH',p.value=as.numeric(alpha),lfc=0)
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