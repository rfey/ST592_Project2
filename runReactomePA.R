# rmf 3.7.2019

## for getting error: In library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE : there is no package called ‘org.Dm.eg.db’
## use following code in R terminal
# library(BiocInstaller)
# biocinstallRepos()
# install.packages("org.Dm.eg.db",repos="https://bioconductor.org/packages/3.4/data/annotation")

# USAGE
args = commandArgs(trailingOnly=TRUE)
usage <- 'Rscript runReactomePA.R <gene list> <alpha> <background list>'
if (length(args) != 3) {
   stop(usage)
}

library(ReactomePA)

# ARGUMENTS
geneFile <- args[1]
alpha <- args[2]
bgFile <- args[3]

genes <- readLines(geneFile)
background <- readLines(bgFile)

# MAIN
x <- enrichPathway(gene=genes,organism="fly",pvalueCutoff=as.numeric(alpha),universe=background)
head(as.data.frame(x))