# rmf 3.4.2019, last modified 3.6.2019

# step 1: write reference file
# input: mir family file, conserved target file

./writeMirTargetTable.sh

# step 2: find DE mirs
# using edgeR script: edgeR_YvsO.R
# NOTE that the edgeR script pre-filters mirs by RAW counts (this analysis, raw count of 1)

./compareYoungVsOld.sh ## # output is a list of mirs with raw pvals and DE calls (miRsDE_youngVsOld0.05alphaNAFC.txt)
./getDEmirs.sh  ## splits DE mirs into up and downregulated lists (downDE_youngVsOld0.05alphaNAFC.txt, upDE_youngVsOld0.05alphaNAFC.txt)

# step 3: write background gene list for pathway analysis
# input: mir-target table, list of DE mirs (downDE_youngVsOld0.05alphaNAFC.txt, upDE_youngVsOld0.05alphaNAFC.txt)

./getBGlist.sh

# step 4: write input gene list for pathway analysis
# NOTE: must convert FBIDs to Entrez Gene IDs

./writeInputList.sh
wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz
gunzip gene2ensembl.gz 
./fbid2entrez.sh

# step 5: run reactomePA

./analyzePathways.sh

#### FINAL FILES TO USE ####
# NOTE: all on github

# input for DE analysis
inputTableEdgeR_bestReps_AS.txt

# list of all mirs used in DE analysis-- result of DE analysis
miRsDE_youngVsOld0.05alphaNAFC.txt

# lists of up and down regulated mirs-- after parsing above file
# NOTE: these are input for getting gene lists for pathway analysis
upDE_youngVsOld0.05alphaNAFC.txt
downDE_youngVsOld0.05alphaNAFC.txt

# input gene list files for pathway analysis
inputGenes_upDE_DESeq0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt
inputGenes_downDE_DESeq0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt
inputGenes_upDE_0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt
inputGenes_downDE_0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt

# background gene list for pathway analysis
background_1.0FPKM_Entrez.txt
