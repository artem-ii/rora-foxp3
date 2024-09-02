CEL_FILES_PATH="/Volumes/artemii-4tb/data/artemii-rorafoxp3-for-geo/total-ewat-microarray/total-ewat-cel"
WD_PATH="/Users/artemii/Yandex.Disk.localized/Inserm/Manuscripts/RORa_Treg/Manuscript/drafts/V11/27-august-2021/papers-for-data/misc"

setwd(WD_PATH)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("oligo", force = TRUE)

library("oligo")
library("affycoretools")
library("pd.mogene.2.0.st");
library("genefilter")
library(RColorBrewer)
library(gplots)

# Reading Data, Correcting Background, annotating
setwd(CEL_FILES_PATH)
celFiles <- list.celfiles(full.names=T)
rawData <- read.celfiles(celFiles)
eset <- rma(rawData, target="core")
eset <- getMainProbes(eset)
eset <- annotateEset(eset, "mogene20sttranscriptcluster.db")
eset <- eset[!is.na(fData(eset)$ENTREZID),]

# Expression table
exp <- exprs(eset)
exp <- as.data.frame(exp)
exp < cbind.data.frame(fData(eset)$ENTREZID, fData(eset)$SYMBOL, exp)
colnames(exp) <- c("ENTREZ_ID","SYMBOL", "WTCdfast","WTCdfast","WTCdfast",
                   "WTCdfast","WTCdfast","KOCdfast","KOCdfast",
                   "KOCdfast","KOCdfast","KOCdfast","WTHFDfast",
                   "WTHFDfast","WTHFDfast","KOHFDfast","KOHFDfast",
                   "KOHFDfast","WTHFDfast_ins","WTHFDfast_ins",
                   "KOHFDfast_ins", "KOHFDfast_ins")
geno <- c("WT","WT","WT","WT","WT","KO","KO","KO",
          "KO","KO","KO","KO","KO","WT","WT","WT","KO","KO",
          "KO","WT","WT","WT","WT","WT","KO","KO","KO","WT","WT","KO",
          "KO")
treatment <- c("Cdfast", "Cdfast", "Cdfast",
               "Cdfast", "Cdfast", "Cdfast", "Cdfast",
               "Cdfast", "Cdfast", "Cdfast", "HFDfast",
               "HFDfast", "HFDfast", "HFDfast", "HFDfast",
               "HFDfast", "HFDfast_ins", "HFDfast_ins",
               "HFDfast_ins", "HFDfast_ins")
insulin <- c("n", "n", "n",
             "n", "n", "n", "n",
             "n", "n", "n", "n",
             "n", "n", "n", "n",
             "n", "y", "y",
             "y", "y")
exp <- rbind.data.frame(c("","",geno), c("","",insulin), c("","",treatment), exp)
write.table(x=exp, file="rorafoxp3hfdmicroarray.txt", sep = "\t", row.names = F, col.names = T)

# Add variables to eSet
eset@phenoData@data$genotype <- factor(c("WT","WT","WT","WT","WT","KO","KO","KO",
                                         "KO","KO","WT","WT","WT","KO","KO",
                                         "KO","WT","WT","KO","KO"))
eset@phenoData@data$treatment <- factor(c("Cdfast", "Cdfast", "Cdfast",
                                          "Cdfast", "Cdfast", "Cdfast", "Cdfast",
                                          "Cdfast", "Cdfast", "Cdfast", "HFDfast",
                                          "HFDfast", "HFDfast", "HFDfast", "HFDfast",
                                          "HFDfast", "HFDfast_ins", "HFDfast_ins",
                                          "HFDfast_ins", "HFDfast_ins"))
eset@phenoData@data$insulin <- factor(c("n", "n", "n",
                                        "n", "n", "n", "n",
                                        "n", "n", "n", "n",
                                        "n", "n", "n", "n",
                                        "n", "y", "y",
                                        "y", "y"))
eset@phenoData@varMetadata <- rbind(eset@phenoData@varMetadata, data.frame(row.names = c("genotype",
                                                                                         "treatment","insulin"),labelDescription = c("genotype", "treatment",
                                                                                                                                     "insulin"),channel = c("_ALL_","_ALL","_ALL_")))
save(eset, file = "rorafoxp3hfdmicroarray.RData")

library(limma)
Conditions <- factor(c("WTCdfast","WTCdfast","WTCdfast",
                       "WTCdfast","WTCdfast","KOCdfast","KOCdfast",
                       "KOCdfast","KOCdfast","KOCdfast","WTHFDfast",
                       "WTHFDfast","WTHFDfast","KOHFDfast","KOHFDfast",
                       "KOHFDfast","WTHFDfast_ins","WTHFDfast_ins",
                       "KOHFDfast_ins", "KOHFDfast_ins"))

design <- model.matrix(~0+Conditions)
colnames(design) <- levels(Conditions)
fit <- lmFit(eset, design)

# CD
contrast.matrix <- makeContrasts(cdfast = KOCdfast - WTCdfast,
                                 levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, coef=1, adjust="BH", number = 20000) #BH  BenjaminiHochberg adjusted P value
write.table(x = table, file = "KOvsWT_CD_total_eWAT_limma.txt", sep =",", col.names = NA, row.names=T)

# HFD
contrast.matrix <- makeContrasts(cdfast = KOHFDfast - WTHFDfast,
                                 levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, coef=1, adjust="BH", number = 20000) #BH  BenjaminiHochberg adjusted P value
write.table(x = table, file = "KOvsWT_HFD_total_eWAT_limma.txt", sep =",", col.names = NA, row.names=T)

volcanoplot(fit2, coef=1, highlight=20, names = fit2$genes$SYMBOL)