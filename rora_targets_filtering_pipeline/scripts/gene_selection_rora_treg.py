#!/usr/bin/env python3
# coding: utf-8

# In[101]:


import pandas as pd
import numpy as np
import os


# In[102]:


ROOT_PATH="/Users/artemii/Desktop/scripts/find_rora_ko_treg_tfs"
PROCESSED_DATA_PATH=os.path.join(ROOT_PATH, "processed_data")
DATA_PATH=os.path.join(ROOT_PATH, "data")


# In[103]:


promoters_with_db_sites_file = os.path.join(PROCESSED_DATA_PATH,
    "merged_promoter_db_table", "merged_promoter_db_table_all_tissues.tsv")
df = pd.read_csv(promoters_with_db_sites_file, index_col=False, header=None)


# In[104]:

print("\n")
print("-------------------------------------------------------")
print("gene_selection_rora_treg.py started\n\n")
print("Parsing file with promoters containing differential binding...")

# Removing rows containing file names
skip_rows = []
row_num = 0
for i in df.iloc[:,0]:
    if i.startswith("db"):
        skip_rows.append(row_num)
    row_num += 1


# In[105]:


miragaia_promoters_db = pd.read_csv(promoters_with_db_sites_file,
    sep="\t", index_col=False, header=None, skiprows=skip_rows)


# In[106]:


print("Read total " + str(miragaia_promoters_db.shape[0]) + " differential binding sites within marker gene promoters (There are duplicates because of different tissues).")


# In[107]:


gene_list = []
for i in miragaia_promoters_db.iloc[:,12]:
    gene_list.append(i.split("_")[0])
gene_list = set(gene_list)


# In[108]:

print("\n")
print("There are " + str(len(gene_list)) + " tissue Treg and Tmem marker genes containing differential histone mark changes in Rora KO Treg")


# In[109]:


# Find transcription factors in the marker genes with promoters

#tf symbols based on RIKEN (see data/riken_tf/tf-db_riken.txt)
riken_tfs = ["1700003F12RIK","1700014N06RIK","1700020N01RIK","1700090G07RIK","2410141K09RIK","2810021J22RIK","3632451O06RIK","4930430A15RIK","4931423N10RIK","5830417I10RIK","6720456H20RIK","6720489N17RIK","9130019O22RIK","A630089N07RIK","AATF","ABT1","ABTB1","ACYP2","ADNP","AEBP1","AEBP2","AEBP2","AEBP2","AES","AHR","AHRR","AIP","AIRE","ALX3","ALX4","ANK1","ANK1","ANK1","ANK2","ANK2","ANKFY1","ANKHD1","ANKRA2","ANKRD1","ANKRD10","ANKRD10","ANKRD10","ANKRD10","ANKRD10","ANKRD2","ANKRD5","ANKRD6","ANKRD6","ANKRD6","AR","ARHGAP17","ARHGAP17","ARHGAP17","ARHGAP17","ARHGAP17","ARID1A","ARID3A","ARID3B","ARID4A","ARID5B","ARNT","ARNT","ARNT2","ARNTL","ARNTL","ARX","ASB1","ASB1","ASB10","ASB11","ASB12","ASB15","ASB2","ASB3","ASB4","ASB5","ASB6","ASB8","ASB8","ASB8","ASB9","ASCL1","ASCL2","ASCL3","ASH1L","ASH2L","ASH2L","ASXL1","ATF1","ATF2","ATF2","ATF3","ATF4","ATF5","ATF5","ATF6","ATF7","ATF7IP","ATM","ATOH1","AW146020","B3GAT3","B930041F14RIK","BACH1","BACH2","BARD1","BARHL1","BARHL1","BARX1","BARX2","BASP1","BATF","BAZ1A","BAZ1B","BAZ2A","BBX","BC024139","BC024139","BCL11A","BCL11A","BCL11A","BCL11A","BCL3","BCL6","BCL6B","BCLAF1","BCLAF1","BCLAF1","BCOR","BCOR","BCOR","BCOR","BCOR","BCOR","BDP1","BLZF1","BLZF1","BLZF1","BMI1","BMYC","BNC1","BPNT1","BRCA1","BRD8","BRDT","BRDT","BRPF1","BSX","BTBD11","BTBD11","BTF3","BTF3","C1D","CAMK4","CARF","CARHSP1","CARM1","CARM1","CASKIN1","CBFB","CBFB","CBFB","CBFB","CBX2","CBX3","CBX4","CBX8","CCNC","CCNC","CCNH","CCNK","CCNT1","CCNT2","CDC5L","CDK2","CDK2","CDK4","CDK7","CDK9","CDKN2A","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDX1","CDX2","CDX4","CEBPA","CEBPB","CEBPD","CEBPG","CEBPZ","CECR6","CHD4","CHGB","CIC","CIC","CIC","CITED1","CITED2","CITED4","CLOCK","CLP1","CML3","CML3","CNOT7","CNOT8","COPS2","COPS5","CREB1","CREB1","CREB1","CREB3","CREB3L1","CREB3L3","CREB3L4","CREB5","CREBBP","CREM","CREM","CREM","CREM","CREM","CREM","CREM","CREM","CREM","CREM","CREM","CRX","CRX","CSDA","CSDA","CTBP2","CTBP2","CTCF","CTNNBL1","CXXC1","DACH1","DACH1","DACH2","DACH2","DAXX","DAXX","DAZAP2","DBP","DBX1","DBX2","DDIT3","DDX20","DDX54","DDX58","DEAF1","DEAF1","DEAF1","DEDD","DEDD","DEDD2","DEK","DIDO1","DIDO1","DIDO1","DLX1","DLX2","DLX3","DLX4","DLX5","DLX5","DLX6","DMAP1","DMRT1","DMRT2","DMRT3","DMRTA1","DMRTA2","DMRTC2","DMTF1","DMTF1","DNMT1","DNMT1","DNMT1","DNMT1","DNMT3A","DNMT3A","DOT1L","DPF2","DPF3","DPF3","DPF3","DR1","DRAP1","DRG1","E2F1","E2F2","E2F3","E2F4","E2F5","E2F6","E2F6","E4F1","EBF1","EBF2","EBF3","EBF3","EBF3","EBF4","EED","EGR1","EGR2","EGR3","EGR4","EHF","ELAVL2","ELAVL2","ELAVL2","ELAVL2","ELF1","ELF2","ELF3","ELF3","ELF4","ELF5","ELF5","ELK1","ELK3","ELK4","ELL","ELL2","EMX1","EMX2","EN1","EN2","ENPP2","ENPP2","EOMES","EOMES","EP300","EPAS1","ERCC2","ERCC2","ERCC3","ERF","ERG","ESR1","ESR2","ESR2","ESRRA","ESRRB","ESRRB","ESRRG","ESRRG","ESRRG","ESX1","ETS1","ETS1","ETS2","ETV1","ETV1","ETV3","ETV3","ETV4","ETV6","EVX1","EVX2","EZH1","EZH2","EZH2","FAH","FANK1","FBXO24","FBXW7","FBXW7","FBXW7","FEM1A","FEM1B","FEM1C","FEV","FHL2","FHL5","FLI1","FMNL2","FOS","FOSB","FOSL1","FOSL2","FOXA1","FOXA2","FOXA3","FOXB1","FOXB2","FOXC1","FOXC2","FOXD1","FOXD2","FOXD3","FOXD4","FOXE1","FOXE3","FOXF1A","FOXF2","FOXG1","FOXG1","FOXH1","FOXI1","FOXJ1","FOXJ2","FOXK1","FOXL1","FOXL2","FOXM1","FOXM1","FOXN1","FOXN2","FOXN4","FOXO1","FOXO3","FOXP1","FOXP1","FOXP1","FOXP2","FOXP2","FOXP3","FOXP3","FOXP3","FOXP3","FOXP4","FOXP4","FOXP4","FOXQ1","FUS","GABPA","GABPB1","GABPB1","GAS7","GAS7","GATA1","GATA2","GATA3","GATA4","GATA5","GATA6","GBX1","GBX2","GCDH","GCDH","GCM1","GCM2","GFI1","GFI1","GFI1B","GFI1B","GLI2","GLI3","GLIS1","GLIS2","GLRP1","GLS2","GMEB1","GMEB1","GMEB2","GSC","GTF2A1","GTF2A1","GTF2B","GTF2E1","GTF2E2","GTF2E2","GTF2E2","GTF2F2","GTF2H1","GTF2H2","GTF2H4","GTF2I","GTF2I","GTF2I","GTF2I","GTF2I","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF2IRD1","GTF3A","GTF3C1","GTF3C4","GTF3C4","HAND1","HAND2","HBP1","HBP1","HCFC1","HCLS1","HDAC1","HDAC10","HDAC10","HDAC10","HDAC10","HDAC11","HDAC2","HDAC3","HDAC5","HDAC5","HDAC6","HDAC6","HDAC8","HDAC9","HDGFRP2","HELZ","HES1","HES2","HES3","HES5","HES6","HES7","HESX1","HEY1","HEY2","HEYL","HHEX","HIC1","HIC1","HIC2","HIF1A","HIF3A","HIF3A","HILS1","HIPK2","HIPK2","HIRA","HIVEP1","HIVEP2","HIVEP3","HLF","HLX","HMG20A","HMG20B","HMG20B","HMG20B","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA1","HMGA2","HMGB2","HMGB3","HMGCS1","HMGN2","HMX1","HMX2","HMX3","HNF4G","HOMEZ","HOMEZ","HOXA1","HOXA10","HOXA10","HOXA11","HOXA13","HOXA2","HOXA3","HOXA3","HOXA3","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXB1","HOXB13","HOXB2","HOXB3","HOXB3","HOXB4","HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXC10","HOXC11","HOXC12","HOXC13","HOXC4","HOXC5","HOXC5","HOXC6","HOXC8","HOXC9","HOXD1","HOXD10","HOXD11","HOXD12","HOXD13","HOXD3","HOXD3","HOXD3","HOXD4","HOXD4","HOXD8","HOXD9","HR","HSBP1","HSF1","HSF2","HSF2BP","HSF4","HSF4","HSF4","HSF4","HSPB9","HTATIP2","HTATIP2","HTATIP2","HTATIP2","HTATIP2","HTR5A","IFNAR2","IFNAR2","IGHMBP2","IKBKB","IKBKB","IKBKG","IKBKG","IKBKG","IKBKG","IKBKG","IKBKG","IKBKG","ILF3","ILF3","ILF3","ILF3","ING1","ING3","ING4","INSM1","INSM2","INVS","IRF1","IRF1","IRF1","IRF2","IRF2BP1","IRF3","IRF3","IRF4","IRF5","IRF5","IRF6","IRF7","IRF7","IRF7","IRX1","IRX2","IRX3","IRX3","IRX4","IRX5","IRX6","ISL1","ISL2","IVNS1ABP","IVNS1ABP","IVNS1ABP","JMY","JUN","JUNB","KEAP1","KEAP1","KEAP1","KEAP1","KHDRBS1","KHDRBS1","KLF1","KLF12","KLF13","KLF15","KLF16","KLF2","KLF3","KLF4","KLF5","KLF7","L3MBTL2","LASS2","LASS4","LASS5","LASS6","LBH","LDB1","LDB1","LDB2","LDB2","LEF1","LHX1","LHX2","LHX3","LHX4","LHX5","LHX6","LHX6","LHX6","LHX6","LHX8","LHX9","LHX9","LHX9","LIMD1","LMO1","LMO2","LMO2","LMO2","LMO2","LMO4","LMO4","LMO4","LMX1A","LMX1B","LRRC6","LSM11","LSM4","LYL1","MAF","MAFA","MAFB","MAFF","MAFG","MAFK","MAGED1","MAML1","MAP3K12","MAP3K12","MAPK1","MAPK1","MAPK3","MAPK8","MAX","MAX","MAZ","MBD2","MBD3","MBD3L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM7","MCM7","MCM7","MCM8","MECP2","MECP2","MEF2A","MEF2B","MEF2B","MEF2C","MEF2C","MEF2D","MEIS1","MEIS1","MEN1","MEN1","MEN1","MEN1","MEOX1","MEOX2","MGA","MGA","MID1","MID1","MITF","MITF","MITF","MIXL1","MKI67IP","MKL1","MKL1","MLL3","MLL5","MLLT1","MLLT10","MLLT10","MLLT10","MLLT10","MLLT10","MLLT3","MLLT3","MLLT6","MNAT1","MNT","MPL","MPL","MSC","MSX1","MSX2","MSX3","MTA1","MTA2","MTA3","MTA3","MTA3","MTA3","MTF1","MTF2","MTF2","MTF2","MTF2","MTF2","MTF2","MTPN","MUSK","MUSK","MUSK","MUSK","MUSK","MUSK","MXI1","MXI1","MXI1","MYB","MYB","MYBBP1A","MYBL1","MYBL2","MYC","MYC","MYC","MYC","MYCBP","MYCS","MYEF2","MYEF2","MYEF2","MYF5","MYF6","MYNN","MYOCD","MYOCD","MYOD1","MYOG","MYST2","MYST2","MYST2","MYST3","MYST4","MYST4","MYT1","MYT1","MYT1","MYT1","MYT1L","MYT1L","MYT1L","MYT1L","NAB1","NAB2","NAB2","NACA","NACA","NCOA1","NCOA2","NCOA2","NCOA3","NCOA4","NCOA4","NCOA5","NCOA6","NCOA6","NCOR1","NCOR1","NCOR1","NCOR2","NCOR2","NCOR2","NDN","NEDD8","NEUROD1","NEUROD2","NEUROD4","NEUROD6","NEUROG1","NEUROG2","NEUROG3","NFAT5","NFAT5","NFATC1","NFATC1","NFATC1","NFATC1","NFATC1","NFATC1","NFATC2","NFATC2","NFATC2","NFATC2","NFATC2IP","NFATC3","NFATC4","NFATC4","NFE2","NFE2L1","NFE2L1","NFE2L1","NFE2L1","NFE2L1","NFE2L1","NFE2L2","NFE2L3","NFIA","NFIA","NFIA","NFIB","NFIB","NFIB","NFIC","NFIC","NFIL3","NFIX","NFIX","NFIX","NFKB1","NFKB2","NFKB2","NFKB2","NFKBIB","NFKBIE","NFKBIL1","NFX1","NFX1","NFYA","NFYA","NFYB","NFYC","NFYC","NFYC","NFYC","NHLH1","NHLH2","NKD2","NKX1-2","NKX2-2","NKX2-2","NKX2-3","NKX2-4","NKX2-5","NKX2-6","NKX2-9","NKX3-1","NKX6-1","NKX6-2","NKX6-2","NMI","NMI","NMI","NOTCH1","NOTCH2","NOTCH3","NOTCH4","NPAS1","NPAS2","NPAS3","NPTXR","NPTXR","NPTXR","NR0B1","NR0B2","NR1D1","NR1D2","NR1H2","NR1H3","NR1H3","NR1H4","NR1H4","NR1H4","NR1I2","NR1I2","NR1I3","NR1I3","NR1I3","NR2C1","NR2C2","NR2E1","NR2E3","NR2F1","NR2F2","NR2F2","NR2F6","NR3C1","NR3C2","NR4A1","NR4A2","NR4A2","NR4A3","NR5A1","NR5A2","NR5A2","NR6A1","NR6A1","NR6A1","NRARP","NRBF2","NRF1","NRF1","NRF1","NRF1","NRF1","NRF1","NRIP1","NRIP2","NRIP2","NRL","NRL","NSD1","NUDT12","NUPR1","OBOX1","OBOX3","OBOX5","OLIG1","OLIG2","ONECUT1","ONECUT3","ORF63","OSTF1","OTP","OTX1","OTX2","OVOL1","PAPOLA","PAPOLB","PAWR","PAX1","PAX2","PAX3","PAX3","PAX4","PAX4","PAX4","PAX5","PAX6","PAX6","PAX6","PAX6","PAX6","PAX7","PAX8","PAX9","PBX1","PBX1","PBX2","PBX3","PBX4","PBXIP1","PDCD11","PDCD7","PDLIM1","PDLIM4","PEG3","PER1","PER1","PER2","PER3","PFDN1","PGR","PHF1","PHF10","PHF10","PHF12","PHF13","PHF2","PHF5A","PHF7","PHF8","PHF8","PHOX2A","PHOX2B","PHTF1","PHTF1","PHTF1","PHTF1","PIAS1","PITX1","PITX2","PITX2","PITX2","PITX3","PKNOX1","PKNOX1","PKNOX2","PKNOX2","PLA2G6","PLA2G6","PLA2G6","PLA2G6","PLAGL1","PLAGL2","PLRG1","PMFBP1","PML","PML","PMS1","POGK","POGK","POLR2A","POLR2B","POLR2C","POLR2E","POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR2K","POLR2L","POLR3B","POLR3D","POLR3D","POLR3E","POLR3E","POLR3H","POLRMT","POU2AF1","POU2F1","POU2F1","POU2F1","POU2F1","POU2F2","POU2F2","POU2F2","POU2F2","POU2F3","POU3F1","POU3F2","POU3F3","POU3F4","POU4F1","POU4F2","POU4F3","POU5F1","POU5F1","POU6F1","PPARA","PPARA","PPARD","PPARD","PPARG","PPARG","PPARGC1A","PPARGC1A","PPARGC1B","PPP1R10","PPP1R10","PPP1R10","PPP1R12B","PPP1R12B","PPP1R12C","PPP1R13B","PPP1R16B","PPP1R16B","PPP1R1B","PPP2R1A","PPP3CB","PPP5C","PQBP1","PQBP1","PQBP1","PRDM1","PRDM16","PRDM16","PRDM16","PRDM5","PREB","PRKAR1A","PROP1","PROX1","PRR3","PRR3","PRR3","PRR3","PRRX1","PRRX1","PRRX1","PRRX2","PRRXL1","PSMC3","PSMC5","PSMD10","PSMD10","PSMD9","PTF1A","PTGES2","PTMA","PURA","PURB","RAB11A","RAB11B","RAB15","RAB18","RAB25","RAB8A","RAB8B","RABL3","RAI1","RAI1","RAI14","RAI14","RARA","RARA","RARA","RARA","RARB","RARG","RARG","RAX","RB1","RB1","RBAK","RBAK","RBBP7","RBL1","RBL1","RBL2","RBM14","RCOR1","REL","RELA","RELB","REM2","REST","REST","RFX1","RFX1","RFX2","RFX2","RFX3","RFX3","RFX4","RFX4","RFX5","RFXANK","RFXANK","RFXAP","RING1","RNF14","RNF14","RNF14","RNF141","RNF38","RNF38","RNF4","RNPS1","RNPS1","RNPS1","RORA","RORA","RORB","RORB","RORC","RPS6KA4","RPS6KA4","RUNX1","RUNX1","RUNX1","RUNX1","RUNX2","RUNX2","RUNX2","RUNX3","RUVBL1","RUVBL2","RXRA","RXRB","RXRB","RXRB","RXRB","RXRG","RXRG","RYBP","SALL1","SALL2","SALL2","SALL3","SALL4","SALL4","SALL4","SAP18","SAP18","SAP18","SAP30","SATB1","SATB1","SATB1","SATB1","SATB2","SCAND1","SCAP","SCAP","SCMH1","SCMH1","SCRT1","SCX","SEC14L2","SERTAD1","SERTAD2","SERTAD2","SETBP1","SFPI1","SH3D19","SHANK3","SHOX2","SHPRH","SHPRH","SIAH2","SIM1","SIM2","SIN3A","SIN3A","SIN3A","SIN3B","SIN3B","SIRT1","SIRT1","SIRT2","SIRT2","SIRT2","SIRT3","SIRT3","SIRT3","SIRT4","SIRT4","SIRT5","SIRT6","SIRT6","SIRT7","SIX1","SIX2","SIX3","SIX4","SIX5","SIX6","SKI","SMARCA1","SMARCA2","SMARCA2","SMARCA4","SMARCA4","SMARCA4","SMARCA5","SMARCB1","SMARCB1","SMARCE1","SMYD1","SMYD1","SNAI1","SNAI2","SNAI3","SNAPC2","SNAPC4","SNAPC4","SNCAIP","SNCAIP","SNCAIP","SNCAIP","SNIP1","SNRPB","SOLH","SOX1","SOX10","SOX11","SOX12","SOX13","SOX14","SOX15","SOX17","SOX18","SOX2","SOX21","SOX3","SOX30","SOX4","SOX5","SOX5","SOX5","SOX6","SOX6","SOX6","SOX7","SOX8","SOX9","SP1","SP2","SP2","SP3","SP3","SP4","SP4","SP5","SP6","SP7","SPDEF","SPIB","SPIC","SPZ1","SQSTM1","SRA1","SRA1","SRA1","SREBF1","SREBF2","SREBF2","SRF","SRY","SSBP2","SSBP2","SSBP3","SSBP3","SSBP4","SSRP1","SSRP1","SSXB1","ST18","ST18","ST18","ST18","ST18","STAG1","STAT1","STAT1","STAT1","STAT2","STAT3","STAT3","STAT3","STAT4","STAT5A","STAT5A","STAT5B","STAT5B","STAT6","STRAP","SUPT3H","SUPT5H","T","TAF10","TAF11","TAF12","TAF13","TAF1A","TAF1B","TAF1B","TAF1C","TAF3","TAF4A","TAF4B","TAF5","TAF5L","TAF6","TAF6L","TAF6L","TAF6L","TAF6L","TAF6L","TAF7","TAF9","TAF9","TAF9","TAL1","TAL2","TARDBP","TARDBP","TARDBP","TARDBP","TARDBP","TARDBP","TBP","TBPL1","TBR1","TBX1","TBX10","TBX10","TBX15","TBX18","TBX19","TBX2","TBX20","TBX20","TBX20","TBX21","TBX22","TBX22","TBX3","TBX3","TBX4","TBX4","TBX5","TBX6","TCEA1","TCEA1","TCEA1","TCEA2","TCEA3","TCEAL1","TCEB2","TCEB3","TCERG1","TCF12","TCF12","TCF12","TCF12","TCF12","TCF15","TCF19","TCF19","TCF19","TCF20","TCF20","TCF20","TCF21","TCF3","TCF3","TCF3","TCF3","TCF3","TCF3","TCF3","TCF3","TCF4","TCF4","TCF7","TCF7L2","TCF7L2","TCF7L2","TCF7L2","TCF7L2","TCF7L2","TCF7L2","TCF7L2","TCFL5","TEAD1","TEAD1","TEAD1","TEAD2","TEAD3","TEAD3","TEAD3","TEAD4","TEAD4","TEF","TEF","TFAM","TFAP4","TFDP1","TFDP2","TFDP2","TFDP2","TFDP2","TFDP2","TFDP2","TGFB1I1","TGIF2","TH1L","THRA","THRAP3","THRB","THRB","THRSP","TIAL1","TLE1","TLE2","TLE2","TLE3","TLE3","TLE3","TLE4","TLE6","TLX1","TLX2","TLX3","TMPO","TMPO","TMPO","TMPO","TMPO","TMPO","TNFAIP3","TNFAIP3","TOB1","TRERF1","TRERF1","TRIB3","TRIM24","TRIM28","TRIP13","TRIP4","TRIP4","TRIP6","TRP53","TRP53","TRP53BP1","TRP63","TRP63","TRP63","TRP63","TRP63","TRP63","TRP63","TRP63","TRP73","TRP73","TRP73","TRPS1","TRPS1","TRRAP","TSG101","TULP4","TULP4","TWIST1","TWIST2","UBN1","UBP1","UBP1","UBTF","UBTF","UGP2","UGP2","UHRF1","UHRF1","UHRF1","UHRF1","UHRF2","USF1","USF2","USP49","UTF1","VAX1","VAX2","VDR","VGLL2","VSX1","WASL","WASL","WBP7","WHSC2","WT1","XAB2","XBP1","XRCC3","YAF2","YAF2","YY1","ZBTB1","ZBTB3","ZBTB3","ZBTB33","ZBTB33","ZCCHC8","ZDHHC1","ZDHHC15","ZDHHC16","ZDHHC21","ZDHHC5","ZDHHC6","ZDHHC6","ZFA","ZFP1","ZFP1","ZFP105","ZFP108","ZFP109","ZFP110","ZFP111","ZFP112","ZFP113","ZFP128","ZFP13","ZFP131","ZFP142","ZFP143","ZFP146","ZFP148","ZFP161","ZFP180","ZFP189","ZFP191","ZFP192","ZFP192","ZFP192","ZFP2","ZFP2","ZFP2","ZFP2","ZFP202","ZFP207","ZFP207","ZFP207","ZFP207","ZFP207","ZFP219","ZFP219","ZFP219","ZFP219","ZFP235","ZFP238","ZFP238","ZFP239","ZFP239","ZFP248","ZFP263","ZFP27","ZFP27","ZFP275","ZFP275","ZFP277","ZFP277","ZFP281","ZFP281","ZFP282","ZFP287","ZFP292","ZFP295","ZFP295","ZFP295","ZFP319","ZFP334","ZFP35","ZFP354A","ZFP354B","ZFP354C","ZFP36","ZFP36L1","ZFP37","ZFP386","ZFP386","ZFP39","ZFP398","ZFP398","ZFP40","ZFP41","ZFP41","ZFP42","ZFP423","ZFP426","ZFP426","ZFP451","ZFP454","ZFP46","ZFP467","ZFP467","ZFP467","ZFP467","ZFP51","ZFP52","ZFP523","ZFP54","ZFP553","ZFP57","ZFP57","ZFP57","ZFP57","ZFP580","ZFP59","ZFP592","ZFP597","ZFP60","ZFP60","ZFP606","ZFP606","ZFP61","ZFP612","ZFP617","ZFP64","ZFP68","ZFP68","ZFP68","ZFP82","ZFP82","ZFP82","ZFP84","ZFP84","ZFP9","ZFP90","ZFP91","ZFP91","ZFP92","ZFP93","ZFP94","ZFP94","ZFP97","ZFPM1","ZFPM2","ZFX","ZFX","ZFY1","ZFY2","ZHX1","ZHX1","ZHX3","ZIC1","ZIC2","ZIC3","ZIK1","ZIM1","ZMYND11","ZMYND11","ZNRD1","ZSWIM4"]
riken_tfs = set(riken_tfs)

miragaia_set = gene_list

tf_in_db_promoters = {x.casefold() for x in riken_tfs} & {y.casefold() for y in miragaia_set}

print("Of these genes " + str(len(list(tf_in_db_promoters))) + " are transcription factors based on RIKEN database")
print("(Read " + str(len(riken_tfs)) + " from RIKEN DB list)")
print("See data/riken_tf/tf-db_riken.txt")


# In[110]:


# microarray = ["Gm23193", "Gm23540", "Mir103-2", "Gm23482", "Gm13251", "Gm23205", "Cd226", "BE692007", "Muc13", "Mir181b-2", "Gm22005", "Gm22908", "Gm23256", "Gm26140", "Gm25654", "2310040G07Rik", "Gad2", "Gm25321", "Tgtp1", "A630023P12Rik", "Trat1", "Gm23058", "9930111J21Rik2", "Mir3967", "Gm25973", "Gm8369", "Gm25506", "Tlr1", "Parp8", "Gm24431", "Gm24775", "Gm24622", "Npy", "Gm22213", "Csn2", "Gm24927", "Gm23557", "Gm22322", "Slamf6", "Pisd-ps3", "C230066G23Rik", "Sell", "Scarna3b", "Lrrc4", "Gm24569", "2310044K18Rik", "Lncpint", "Atxn7l1os1", "Cxcr1", "Gm22446", "Gm25636", "Scarna3a", "Id3", "Il1rl1", "Mirlet7a-1", "Mir3074-1", "Nt5e", "Tgtp2", "Snord17", "A630089N07Rik", "Tcf7", "Gm22304", "Snord35b", "Frmd5", "Spef2", "Mir17", "Zbtb10", "Trbv23", "Bach2os", "Gm26514", "Tbrg3", "Snord45b", "Gpr83", "Atp10d", "Gm22272", "Gm17218", "Gm25945", "Snord92", "Gm22887", "Gm26419", "Tnfsf4", "Gm25109", "Gm22960", "H2-Aa", "Mir1967", "Zfp874a", "Ikzf2", "Gm22890", "Ercc6l2", "Airn", "Gm26187", "Gm26666", "Slc25a36", "Kcnq1ot1", "Rev3l", "Cpm", "Nrn1", "Zfp831", "Dpp4", "Gm25518", "Gm25515", "Mir17hg", "Gm24142", "2610204G07Rik", "Gdap10", "Gm24878", "Mir1955", "2010016I18Rik", "Tnfrsf26", "Zfp141", "LOC105247533", "Pik3c2a", "Trav6-1", "Gm26255", "Adamtsl4", "Arap2", "Gm26838", "Clec2i", "A430078G23Rik", "Hc", "Gm23513", "Snord37", "Plcxd2", "Gm22305", "Syp", "Olfr99", "Gm26265", "Kdm7a", "Igf1r", "Gm25779", "Adcy4", "Gpnmb", "Gm25913", "Rab6b", "Mir5098", "Arnt2", "Gm25323", "Fam196b", "Gm25196", "Hddc3", "Gm13833", "n-R5s192", "P2rx7", "Gm22572", "Gm16158", "Gm20456", "Plcl1", "2010300F17Rik", "Gm24813", "Gm23054", "Gm24876", "n-R5s33", "Gm22977", "4930595D18Rik", "Ar", "2010111I01Rik", "Smarca5-ps", "n-R5s39", "E030011O05Rik", "St8sia4", "1700054O19Rik", "Itga4", "Neb", "Mir1198", "Gm12000", "Gm26877", "Gm22362", "Gm22200", "1700109H08Rik", "Gimap7", "AI987944", "Gm22925", "Pln", "Gm24445", "Gja1", "Bach2", "Gm16085", "Letm2", "Crebrf", "B230325K18Rik", "Traj61", "Gm20236", "Fbxl4", "Hpgds", "4930403D09Rik", "Smad3", "Fam65b", "Lpin1", "Mtmr7", "Htr2b", "Gm23297", "Tgfbr1", "Ago3", "Gbp4", "Klhl24", "Tanc2", "Gm16523", "Gm15777", "S1pr1", "Zfp329", "4930431P03Rik", "Mplkip", "Simc1", "Pikfyve", "Gm25882", "Gm15573", "Eif3j1", "Fktn", "Rab39b", "Gm23734", "Nf1", "Upf2", "Eml5", "Gm20695", "Gm24053", "Bmp7", "Akap9", "Gm24119", "Gm20139", "Gm23206", "Ankrd26", "Gm3383", "Mreg", "Gm24938", "Zdhhc23", "Mbtd1", "Gm12185", "Gm11672", "Zfp397", "Gm26225", "Trav9-1", "Zfp442", "D630008O14Rik", "Bmpr2", "Ly75", "Klf7", "Snord53", "Map1b", "Trav5-1", "Fnip2", "Gm24411", "Rictor", "2810474O19Rik", "Gm11465", "Rev1", "Gm24916", "Herc1", "Mr1", "A830080D01Rik", "LOC73899", "Gm22738", "Mir18", "n-R5s80", "Stra6", "Lama5", "Macf1", "Plac8", "Entpd7", "Gm22327", "Ptges3l", "Hscb", "Ocel1", "Eci1", "Mir32", "Cyp11a1", "Tmem138", "Gm6594", "Gm25244", "Ube2t", "Gm8973", "Cks2", "Gm12891", "Rpl29-ps5", "Shcbp1", "Mir684-2", "Tsx", "C4a", "Defa25", "Suox", "Ska3", "Prkaca", "Nr4a3", "Gm23500", "S100a1", "Gm22224", "Idi2", "Gm12166", "Gm25865", "Tprkb", "Kifc1", "Adgrg5", "Adck1", "Rgs16", "Igf2bp2", "Gm25990", "mt-Tr", "Sema6d", "Ppp1r2-ps3", "Erdr1", "Emilin2", "Gm8724", "Cdkn3", "Lce1h", "Cdkn2b", "Gm25720", "Cdc20", "Sema7a", "Chst11", "Traj44", "Spats2", "Gm23723", "Galnt3", "Gm6685", "Nudt14", "Olfr820", "Gm14029", "Gm21887", "Clnk", "Hip1r", "Mnd1", "Ndufb3", "Gm19705", "Ighv2-6-8", "Cd8b1", "mt-Tm", "Sgms2", "1700012B09Rik", "Gm25663", "Tspan6", "Traj40", "Gm5327", "Gm25227", "Gm6195", "Dnajb13", "mt-Ty", "mt-Tt", "Gm16494", "Ccl4", "Gm7676", "Gm10320", "Isg15", "Pilrb2", "Lrrk1", "Traj53", "Kcnk7", "Gm22723", "Spata24", "Gm25989", "Lag3", "Ppm1h", "Gm25985", "mt-Tf", "Xdh", "Plk1", "mt-Tc", "Gm23389", "Gm15473", "Gzmk", "Gm25360", "Snora30", "Ccnb2", "Olfr1269", "Batf3", "Vash1", "Fkbp11", "Apoa2", "mt-Tk", "B4galt5", "Gpr141", "Hmga1", "Lax1", "Rab27a", "Gm13275", "Adora2b", "Spry4", "Ifng", "Traj56", "mt-Tq", "Gm7075", "Havcr2", "Cldnd2", "Myl6b", "Sult2b1", "Pgam1", "Apoo-ps", "Ighv14-4", "Rpl7a-ps5", "Gch1", "Pxdc1", "Trav6d-4", "Nkg7", "Spry2", "Dennd3", "Gm23444", "Gm24497", "Gm25939", "Trav13d-4", "Avil", "Trav5-4", "Acbd7", "Pak6", "Pmf1", "Gm22369", "Plek", "Gm24415", "Slamf7", "Zbtb32", "Rpl21-ps7", "Cdc25c", "n-R5s28", "Psmb5-ps", "Tnp2", "Tjp1", "Klrc1", "1500009L16Rik", "Cd40lg", "Zfp683", "Clybl", "Arsb", "Gm6625", "Slc37a2", "Hist1h3e", "Gem", "Gm5566", "Gm3238", "Gm26088", "Nkain1", "Gm15135", "Gm25813", "Gm11937", "Mpzl2", "Rps13-ps1", "S100a5", "Spp1", "Il31", "Gm16224", "mt-Tv", "Trim12a", "AF357399", "Hist1h2ag", "Acpp", "Ly6c2", "Ccl3", "Gm8871", "Ccl5", "Klrd1", "Cd8a", "Pcbd2", "Gpr33", "Khdc1c", "Hid1", "Gm26014", "Hist1h2bk", "Gm10566", "Magohb", "mt-Tg", "Gzmc", "Ifitm1", "Mir1932", "Pdcd1", "Gm8922", "Xcl1", "Gm5152", "Gm4609", "Cenpw", "Gm6460", "Hist1h2bc", "Il1a", "Ly6g5b", "Crabp2", "Gm8906", "Gm6455", "Tmem171", "Insl6", "Fam20a", "Fasl", "Srgap3", "Gm8879", "Il17a", "Il2", "Cd160", "Speer1", "Atp8b4", "Gm2663", "LOC102638993", "Scimp", "Il6", "Gm26527", "Cd244", "Eomes", "Gm12214", "Il23a", "4933402N22Rik", "Crtam", "Il3", "Il5", "Gstm5", "Fcer1g", "Prf1", "Csf2", "Gm1966", "Klra5", "Trim30d", "Akr1c18", "Gzmg", "Hbb-bt", "Gzmf", "Gzmd", "Gzme"]
print("\n")
print("-------------------------------------------------------")
print("Step 3: Filtering based on RORA KO Treg differential gene expression")
print("-------------------------------------------------------\n\n")
print("Reading RORA KO Treg gene expression data...")

rora_treg_invitro_gene_expression_file = os.path.join(DATA_PATH, "rora_treg_microarray", "rora_treg_il4_invitro_expression.txt")

expression_data_df = pd.read_csv(rora_treg_invitro_gene_expression_file)
expression_data_df = expression_data_df.rename(columns={"P.Value": "p_value"})
expression_data_df_filtered = expression_data_df.query("p_value < 0.05").query("abs(logFC) > 0.5")
expression_data_df_filtered = expression_data_df_filtered[~expression_data_df_filtered.loc[:,"SYMBOL"].isna()]
microarray = set(expression_data_df_filtered.loc[:, "SYMBOL"])

db_promotors_diff_expressed_in_rora_treg = microarray.intersection(miragaia_set)

print("\n")
print(f'The following tissue T cell marker genes have changed histone marks and are differentially expressed (p. value < 0.05, FC > ±1.2) between WT and Rora KO Treg \n{db_promotors_diff_expressed_in_rora_treg}')


# In[111]:

print("\n")
print("-------------------------------------------------------")
print("Step 4: Filtering based on the presence of RORA Binding sites")
print("-------------------------------------------------------\n\n")
print("Reading annotation of Fang et al 2014 murine liver RORA ChIP-seq...")

rora_chipseq_annotation_file = os.path.join(PROCESSED_DATA_PATH,
                                            "fang_et_al_rora_chipseq_annotation_homer",
                                            "fang_et_al_rora_chipseq_annotation_homer.tsv")

rora_chipseq_annotation = pd.read_csv(rora_chipseq_annotation_file, sep="\t")


# In[112]:


rora_peaks_fang_et_al_liver = set(rora_chipseq_annotation.iloc[:,15].dropna())
print("Read RORA binding signals in " + str(len(rora_peaks_fang_et_al_liver)) + " genes with gene names")


# In[113]:


selection_with_rora_sites = rora_peaks_fang_et_al_liver.intersection(db_promotors_diff_expressed_in_rora_treg)

print(f'Of the previous list, the following have RORA binding sites \n{selection_with_rora_sites}')


# In[114]:

print("\n")
print("-------------------------------------------------------")
print("Step 5: Filtering based on tissue Treg transcription factor network gene list")
print("-------------------------------------------------------\n\n")
print("To reinforce the focus on the genes, regulating Treg identity and function in non-lymphoid tissues")
print("we further filtered our candidate genes using known transcription factor genes, exerting broad control")
print("of Treg.")
print("Data taken from DiSpirito et al. 2018 Table S2, see data/dispirito_et_al/list_of_tissue_treg_tf_network/aat5861_Table_S2.xlsx (Sheet 2)")
print("and data/dispirito_et_al/source_article/ for publication")
print("Reading gene list of tissue-specific Treg transcription factors from DiSpirito et al. 2018...")

dispirito_tf = ["Gata3","Nr1d1","Jdp2","Crem","Bach2","Jund","Fosb","Junb","Elf1","Ets1","Nfatc1","Rela","Relb","Nfkb1","Rel","Nfatc2","Nfkb2","Jdp2","Crem","Atf4","Bach2","Maff","Jund","Maf","Nfil3","Junb","Xbp1","Mafg","Batf","Fosl2","Atf6","Atf3","Bach1","Fosb","Hlf","Runx3","Runx2","Nr1d2","Rorc","Vdr","Nr3c1","Nr4a1","Rfx1","Gata3","Fli1","Etv3","Jdp2","Crem","Atf4","Bach2","Jun","Maff","Jund","Dbp","Junb","Fosl2","Atf6","Atf3","Fos","Fosb","Hlf","Pparg","Runx2","Rora","Nr3c1","Nr4a1","Fli1","Runx3","Jdp2","Crem","Bach2","Maf","Junb","Rel","Nfkb2"]
dispirito_tf = set(dispirito_tf)


# In[116]:


de_db_marker_gene_promotors_dispirito_tf = selection_with_rora_sites.intersection(dispirito_tf)
print("\n")
print(f'When further filtering the set using DiSpirito et al. TF network list, the result is the following: \n{de_db_marker_gene_promotors_dispirito_tf}')
