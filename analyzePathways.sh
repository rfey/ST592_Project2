#!/bin/bash

# rmf 3.8.2019, last modified 3.13.2019

alpha=0.05

# edgeR
Rscript runReactomePA.R inputGenes_upDE_0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt $alpha background_1.0FPKM_Entrez.txt
mv Rplots.png sigPwys_upDE_edgeR.png # rename default plot name
Rscript runReactomePA.R inputGenes_downDE_0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt $alpha background_1.0FPKM_Entrez.txt
mv Rplots.png sigPwys_downDE_edgeR.png

# DESeq
Rscript runReactomePA.R inputGenes_upDE_DESeq0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt $alpha background_1.0FPKM_Entrez.txt
mv Rplots.png sigPwys_upDE_DESeq.png
Rscript runReactomePA.R inputGenes_downDE_DESeq0.05alphaNAFC_1.0FPKM_0.80PCT_Entrez.txt $alpha background_1.0FPKM_Entrez.txt
mv Rplots.png sigPwys_downDE_DESeq.png
