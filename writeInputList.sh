# rmf 3.4.2019

# final choice:
python writePathwayInputList.py mirTargetTable_0.80PCT.txt background_1.0FPKM.txt downDE_youngVsOld0.05alphaNAFC.txt
python writePathwayInputList.py mirTargetTable_0.80PCT.txt background_1.0FPKM.txt upDE_youngVsOld0.05alphaNAFC.txt


# checking til we get a reasonable number of genes for pathway analysis
# least and most stringent for each alpha

# downregulated alpha 0.05
python writePathwayInputList.py mirTargetTable_0.00PCT.txt background_0.0FPKM.txt downDE_youngVsOld0.05alphaNAFC.txt 
python writePathwayInputList.py mirTargetTable_0.80PCT.txt background_1.0FPKM.txt downDE_youngVsOld0.05alpha1.5FC.txt

# upregulated alpha 0.05
python writePathwayInputList.py mirTargetTable_0.00PCT.txt background_0.0FPKM.txt upDE_youngVsOld0.05alphaNAFC.txt
python writePathwayInputList.py mirTargetTable_0.80PCT.txt background_1.0FPKM.txt upDE_youngVsOld0.05alpha1.5FC.txt

# downregulated alpha 0.01
python writePathwayInputList.py mirTargetTable_0.00PCT.txt background_0.0FPKM.txt downDE_youngVsOld0.01alphaNAFC.txt
python writePathwayInputList.py mirTargetTable_0.80PCT.txt background_1.0FPKM.txt downDE_youngVsOld0.01alpha1.5FC.txt

# upregulated alpha 0.01
python writePathwayInputList.py mirTargetTable_0.00PCT.txt background_0.0FPKM.txt upDE_youngVsOld0.01alphaNAFC.txt
python writePathwayInputList.py mirTargetTable_0.80PCT.txt background_1.0FPKM.txt upDE_youngVsOld0.01alpha1.5FC.txt
