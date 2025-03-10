#! /usr/bin/env Rscript
# Miragaia et al 2019 Table S3 - Marker genes for Treg and Tmem subpopulations found per organ. Related to Figure 2
# data/miragaia_et_al/tissue_treg_tmem_genesets/mmc4.xlsx
# ±10kb promoters

ROOT_PATH="/Users/artemii/rorafoxp3_2023/code/rora-foxp3/rora-foxp3/rora_targets_filtering_pipeline/"
OUTPUT_PATH=paste0(ROOT_PATH, "processed_data/miragaia_extended_promoter_regions/")


PROMOTERS_PATH=paste0(ROOT_PATH, "processed_data/miragaia_epd_promoters/")

dir.create(PROMOTERS_PATH)
setwd(PROMOTERS_PATH)
files <- c("bln-tmem-activated-mouse_epdnew_QUuxV.bed",
           "bln-tmem-lymphoid-mouse_epdnew_EOLUT.bed",
           "bln-tmem-th1-mouse_epdnew_Y8WvW.bed",
           "bln-tmem-th17-mouse_epdnew_CLXQF.bed",
           "bln-treg-effector-mouse_epdnew_2ZsD0.bed",
           "bln-treg-lymphoid-mouse_epdnew_G8g1W.bed",
           "bln-treg-nlt-like-mouse_epdnew_chWf4.bed",
           "bln-treg-stat1-mouse_epdnew_0CIeZ.bed",
           "colon-lymphoid-tissue-mouse_epdnew_KKxsm.bed",
           "colon-nlt-mouse_epdnew_4ZTfj.bed",
           "colon-stress-mouse_epdnew_akHQB.bed",
           "colon-suppress-mouse_epdnew_S9B7g.bed",
           "colon-tmem-lt-like-mouse_epdnew_tXEkA.bed",
           "colon-tmem-nlt-mouse_epdnew_vdEx4.bed",
           "colon-tmem-stress-mouse_epdnew_OtRaM.bed",
           "colon-tmem-stress-mouse_epdnew_p2bEm.bed",
           "colon-tmem-th1-mouse_epdnew_tljbV.bed",
           "colon-tmem-th2-mouse_epdnew_91exD.bed",
           "colon-tmem-th17-mouse_epdnew_Hz2gs.bed",
           "mln-tmem-activated-mouse_epdnew_J2K5d.bed",
           "mln-tmem-cxcl10-mouse_epdnew_OyenH.bed",
           "mln-tmem-cycling-mouse_epdnew_y3zhP.bed",
           "mln-tmem-lymphoid-egressing-mouse_epdnew_BBtYE.bed",
           "mln-tmem-lymphoid-resident-mouse_epdnew_vSEI5.bed",
           "mln-tmem-stress-mouse_epdnew_nRJ0b.bed",
           "mln-tmem-th17-mouse_epdnew_17oTr.bed",
           "mln-treg-effector-mouse_epdnew_xXthl.bed",
           "mln-treg-lymphoid-mouse_epdnew_nCHWU.bed",
           "mln-treg-nlt-like-mouse_epdnew_mMFgn.bed",
           "skin-tmem-nlt-mouse_epdnew_ObiLI.bed",
           "skin-treg-nlt-mouse_epdnew_oU87U.bed",
           "skin-treg-tmem-low-quality-mouse_epdnew_eRauJ.bed")

for (i in files){
  setwd(PROMOTERS_PATH)
  print(paste0("Extending ±10kb from promoter coordinate for file: ", i))
  df <- read.csv(i, sep="\t", header = F)
  df <- as.data.frame(df)
  df[,2] <- as.numeric(df[,2])
  df[,3] <- as.numeric(df[,3])
  for (j in c(1:nrow(df))){
    df[j,2] <- df[j,2]-10000
    df[j,3] <- df[j,3]+10000
  }
  name <- paste("extended+-10kb-", as.character(i), sep = "")
  write.table(df, file = name, quote = F, col.names = F, row.names = F, sep = "\t")
}

print("Done extending regions")
