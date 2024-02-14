if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

bach2exprs <- read.csv('/Users/artemii/Yandex.Disk.localized/Inserm/Data/chip-seq-analysis/6-may-2021/Bach2 represses effector programmes to stabilize Treg-mediated immune homeostasis/bach2-fpkm', sep = "\t", header = T)
#scatter.smooth(bach2exprs.log2)
keep <- apply(bach2exprs[2:7], 1, function(x) !any(x<5))
bach2exprs <- bach2exprs[keep,]
bach2exprs.log2 <- apply(bach2exprs[,2:7], 2, log2)
rownames(bach2exprs.log2) <- bach2exprs$gene_id
means.wt <- apply(bach2exprs.log2[,1:3],1,mean)
means.ko <- apply(bach2exprs.log2[,4:6],1,mean)
bach2exprs.log2.mean <- as.matrix(cbind(means.wt,means.ko))
#bach2exprs.log2 <- bach2exprs.log2[means > 1,]
roraDE.up <- c("Gm10566", "Ndufb3", "Rpl21-ps7", "Ppp1r11", "Havcr2", "Gm22369", "Gm24272", "Sult2b1", "Trim30d", "Gzme", "C4a", "Gm11893", "Nrbf2", "Cenpm", "Mir1932", "Gm6195", "Tsx", "Gm23669", "Immp1l", "Magohb", "Polr2k", "Cdkn3", "Pcbd2", "Gm23389", "Tyms", "1700010B08Rik", "Gm10696", "Gm22319", "Gm4609", "Gm5566", "Gm25244", "Myl6b", "Gm23128", "Pak6", "Erdr1", "Lax1", "Cenpw", "Gm6594", "Cmc2", "Apoo-ps", "Dnaja3", "Fasl", "Trav13d-4", "2900026A02Rik", "Tprkb", "Srgap3", "Suox", "1810013A23Rik", "Adck1", "Gm25989", "Gm22224", "Gm25201", "Hdgfrp3", "Pgam1", "BC031181", "Gm25342", "Gigyf1", "Pmf1", "Prss57", "Sarnp", "3110040N11Rik", "Fkbp2", "Scarb2", "Gabarapl2", "Gm16494", "Gm9855", "Tmem171", "Cd244", "B4galt5", "Gm8615", "Zbtb8os", "Ndufv1", "1500009L16Rik", "Traf5", "Cdk1", "Csf2", "Rnaseh2c", "Bsdc1", "Ccne1", "n-R5s10", "Ppp1r2-ps3", "Platr2", "Lce1h", "A730063M14Rik", "Gm23500", "Pomk", "Eci1", "Nphp1", "Gm15849", "Gm23444", "Gm24497", "Gm25939", "Cd160", "Cks2", "Gm12891", "n-R5s28", "Arsb", "Spry2", "Gm2663", "Gm8973", "Gm6625", "Gins1", "Grpel1", "Oog1", "Ivd", "Smpd5", "Lrrc8a", "Defa25", "Ppm1h", "Mrpl11", "Khdc1c", "Ido1", "Mrpl15", "Gzmd", "Gm6788", "Bpnt1", "Gm5327", "Mir383", "Cdc20", "H2afv", "Hist1h2ag", "Gm8394", "Eomes", "Gm22956", "Mpzl2", "Gm8702", "Gm6685", "Il6", "Zfp683", "Tyw3", "Gm3934", "Cib1", "Otulin", "Trav5-4", "Trav3d-3", "2810428I15Rik", "1700037H04Rik", "Gm20500", "Igkv5-37", "Hip1r", "Olfr419", "Emc10", "Gm23795", "Cyb5b", "Gpha2", "Slc37a2", "Eif1b", "mt-Tc", "Hist1h2bk", "Scimp", "Tmem138", "Gm3238", "Sct", "mt-Tq", "Lzts1", "Gm9857", "Gm7075", "Cenpa", "mt-Tg", "Prkaca", "Arhgap19", "Gm22681", "Slc25a14", "Gm12614", "Gm9797", "mt-Tf", "Ptges3l", "Cuta", "Gm24502", "Mrpl51", "S100a1", "Gm22937", "Mfsd11", "Gm13620", "I830127L07Rik", "Akr1c18", "Oip5", "Gm15946", "Ddx39", "Spp1", "Nasp", "Churc1", "Nol7", "Olfr820", "Ddt", "Serf1", "Tomm22", "Chsy1", "Il31", "Mphosph6", "Nme3", "Cenpn", "Hist2h2ac", "Rplp0", "Psme2", "1700028B04Rik", "Mir3110", "Zfp30", "Gm26088", "Vmn1r159", "Gm29719", "Ly6c2", "Lsm5", "Akip1", "Grhpr", "LOC105245385", "Clybl", "Gm25360", "1700012B09Rik", "B4gat1", "Sumo1", "LOC102639794", "Lipt2", "Gm23974", "mt-Tv", "Gm8050", "Gzmc", "Trim12a", "Gm7589", "Vmn2r86", "Gm23666", "Mir598", "D930030I03Rik", "Insl6", "Rnf113a1", "Adgrg5", "Fam189b", "Mrps33", "S100a5", "Mir654", "Cracr2a", "Gm22361", "Hist1h1b", "Apoa2", "Samd8", "Dnajc15", "Gzmf", "Spink2", "Ybey", "Flot2", "Mrrf", "Xdh", "2010320O07Rik", "Alg8", "Dpy30", "Nsmf", "Sema7a", "Dnajc19-ps", "Cenpp", "Cox19", "Smg8", "Hmx2", "Pgap2", "Sh3bp2", "Gm22880", "1700039E22Rik", "Mir148b", "Psmb10", "Gm6444", "Gm22370", "Gm8879", "Cyp3a25", "Atp5f1", "Pdk3", "Slc9a5", "Emilin2", "Dut", "Med11", "Ccl28", "Gm22327", "Gramd1b", "Gm24728", "Cxcl16", "Gstp2", "Mgat2", "Mip", "Abracl", "Txnl4b", "Olfr1269", "4930578G10Rik", "mt-Tn", "Sema6d", "Hdgf", "Rpl21-ps14", "Cbwd1", "Cdc25c", "Acbd7", "Tacc3", "Haus5", "Acad8", "Gm10482", "Tlr2", "Gm20305", "Ngdn", "LOC102635566", "Olfr720", "Gm4775", "Col18a1", "Gm23723", "Gm25663", "Prim2", "Gm10320", "Tmem33", "Gm25985", "Fkbp3", "Zbtb32", "Isg15", "Cdkn2b", "Rala", "Gm5900", "Krt1", "Lcn4", "Pbx4", "2310040G24Rik", "Il3", "Gm8906", "Gm5093", "Ggcx", "Rnaseh2b", "Tmed9", "Mex3b", "Ndufb2", "Anxa2", "Gm26896", "Gstm5", "Gm15473", "Gm24415", "9130023H24Rik", "Rps6ka2", "Spc24", "Adgrg3", "Bcl2a1c", "Gm11559", "Mnd1", "Lamtor5", "Gins2", "Islr", "1700063A18Rik", "Lair1", "4933402N22Rik", "Rgs16", "Fam92a", "Tmem194", "Ndufb4", "Speer1", "Gm26737", "Gm3436", "Gm14029", "Gm7762", "Trav13-2", "Avil", "Serpinb6b", "Hilpda", "Natd1", "Sdhd", "Ccl5", "Yif1a", "Prdx4", "Fmnl3", "Gm17103", "Slc36a3os", "Psenen", "Bub1", "Sgms2", "Noct", "Bad", "LOC100043918", "Gtf2e2", "Ocel1", "Arl8a", "Gm25720")
roraDE.down <- c("B230325K18Rik", "Kif13b", "Ptcd3", "4930431P03Rik", "Frmd5", "Scaf11", "B3gnt5", "Naa16", "Gm22260", "Vcpip1", "Rere", "Zfp36l2", "Zyg11b", "Rap1gds1", "Akap9", "Trbj2-7", "Kdm5a", "Rab6b", "Gm22304", "Gm24569", "Mir134", "Rnf2", "Gm22845", "Gm22925", "Trib1", "Phc3", "Id3", "Akap13", "Gm26254", "Gm15879", "Ambra1", "Gm25895", "Phkb", "Iqgap2", "Akt3", "Mettl20", "Apom", "Gm21944", "Elmo1", "Gm24126", "Gm15334", "Marf1", "Gm10851", "Mycbp2", "Mir190b", "Aftph", "B230219D22Rik", "Ppp3cb", "Gm14093", "Prrc2c", "Dqx1", "Kbtbd2", "Gm24938", "Gm22322", "Emp2", "F630111L10Rik", "Tigd2", "Gm22362", "Bach2os", "Gm24878", "Arap2", "Gm15735", "Gad2", "Sms", "Il6ra", "Gm22540", "Sacs", "Gm15777", "Nt5c3", "Gm26515", "Ap4e1", "Mir1967", "Bmp2k", "Hectd1", "Gm26493", "Upf2", "Zfp160", "Enpp1", "Smad3", "Gm20695", "LOC105247533", "Ttll4", "Mbtd1", "H2-M10.6", "Gm3371", "4930558J22Rik", "Qtrtd1", "Sp4", "Avpr2", "Gabbr1", "Gm22005", "2310044K18Rik", "Gm24762", "Gm25817", "Nkiras1", "Rbpj", "Mir181b-2", "Dclre1a", "Cnst", "E030011O05Rik", "Jmy", "Arhgef12", "Pnpla8", "Olfr1157", "Rny1", "Gm22504", "Fam208b", "Muc13", "Wapl", "Fam214a", "Gm12000", "Usp32", "Ifngr2", "Phtf1os", "Gm23391", "Nup153", "Gm23293", "Gm25109", "A630072M18Rik", "Ralgapb", "Amd1", "Klhl2", "Gm25613", "Gm17218", "Adhfe1", "Whsc1l1", "6820431F20Rik", "Plcl1", "Ercc6", "Klrg1", "Srgap2", "Dnajc10", "Huwe1", "Fam65b", "Olfr165", "Gm15346", "Gm26409", "D630008O14Rik", "Mir302c", "Hdac4", "Srrm2", "Pisd-ps3", "Crebrf", "Gm22083", "4930549G23Rik", "Kdm4c", "Foxred2", "S1pr1", "Dnajb7", "Bmp7", "Ago3", "Zfp955b", "Mr1", "Mplkip", "Ly75", "Arhgef7", "2610020H08Rik", "Gm26518", "Mrgpra6", "Rev1", "AU020206", "Zdhhc23", "Zfp750", "Sav1", "Slc19a2", "Gm25945", "Aim2", "Lncpint", "Mir3074-1", "AB010352", "Inppl1", "n-R5s34", "Gm17767", "St6gal1", "Kdm7a", "Lemd3", "Rpl3", "AU041133", "Gm16033", "Zfp442", "Neb", "Klhl24", "Igkv12-41", "Snord92", "Rnf219", "Lrch3", "A830080D01Rik", "Gm22977", "N4bp2l1", "Mfap1a", "Phtf2", "Gm26187", "Tbc1d12", "Gm22593", "Fcho2", "Prkag2os1", "Herc1", "Rnmtl1", "Abcc4", "2900060B14Rik", "Jmjd1c", "Gm24235", "Gm23054", "4930594M22Rik", "Gm26265", "Efr3a", "Atm", "Gpr83", "Atp13a3", "Gm22787", "Gm26419", "Gm23209", "Zfp81", "Tgtp2", "Zfp59", "Gm16085", "Polr3g", "Trim68", "Kmt2a", "Yod1", "Gm22213", "Fam46a", "D730005E14Rik", "Zfp143", "Trav9-1", "2010016I18Rik", "Gm25506", "Gm25654", "Mreg", "AA474408", "Ppp1r12a", "Slamf6", "Syt11", "Gm22311", "Scarna3a", "Peg13", "Gm13251", "Klf7", "Satb1", "Rab39b", "Gm5528", "Gm22621", "Gm25136", "Herc4", "Lrrc4", "Gm24119", "Foxp3", "5430403G16Rik", "Zfp141", "Cyp4b1-ps2", "Smo", "Kdm3b", "Tex15", "Cd226", "Lrrc58", "Gm26406", "Fbxl4", "Gm16158", "Pros1", "Muc3a", "Bckdhb", "Atxn7l1os1", "Zfp329", "Gm26255", "Insr", "Gm23494", "Shprh", "1700109H08Rik", "Gm10477", "Nf1", "Vmn2r68", "Lama5", "Gm10197", "Gm25973", "Ankrd26", "Eif3j1", "Gm26666", "2310040G07Rik", "Gm22200", "A630089N07Rik", "Cpm", "Trbv2", "LOC102634683", "n-R5s192", "Fam120c", "Fbxw7", "Ankrd11", "Gm25196", "Rictor", "Chd6", "Magi3", "Cxcr1", "9930111J21Rik2", "Polr2h", "Dpp4", "Macf1", "Rev3l", "Gimap7", "Gm22765", "Mllt4", "Gm26785", "Gm15182", "A630023P12Rik", "Gm22272", "Nars2", "Gm25882", "P2rx7", "Ighv1-58", "Zfp397", "Hpgds", "Gm24775", "Gm22887", "n-R5s154", "Mir103-2", "Gm23193", "Mtmr7", "Mxi1", "Parp8", "Alg6", "Smarca5-ps", "Zfp451", "Letm2", "Gm25323", "Zfp831", "Gm15601", "Snord17", "Simc1", "Gm25636", "A430078G23Rik", "Mirlet7a-1", "Plcxd2", "Fnip2", "Zfp944", "Gm24142", "Tlr1", "Gm15342", "Arnt2", "Airn", "Mir3967", "Gm22908", "Gm23256", "Gm26140", "Zfp874a", "BE692007", "4930595D18Rik", "Dync2h1", "Tgtp1", "Slc25a36", "Gm23482", "Gm17056", "Ercc6l2", "Gm24927")
bach2exprs.log2.mean.roraDE.up <- bach2exprs.log2.mean[row.names(bach2exprs.log2.mean) %in% roraDE.up,]
bach2exprs.log2.mean.roraDE.down <- bach2exprs.log2.mean[row.names(bach2exprs.log2.mean) %in% roraDE.down,]

bach2exprs.log2.roraDE.up <- bach2exprs.log2[row.names(bach2exprs.log2) %in% roraDE.up,]
bach2exprs.log2.roraDE.down <- bach2exprs.log2[row.names(bach2exprs.log2) %in% roraDE.down,]
#sort the genes instead of clustering them?
#Here I plotted means. I think it's better to plot all the replicates

bach2exprs.log2.mean.roraDE.up <- as.data.frame(bach2exprs.log2.mean.roraDE.up)
bach2exprs.log2.mean.roraDE.up$sort <- (bach2exprs.log2.mean.roraDE.up.ranked$means.wt - bach2exprs.log2.mean.roraDE.up.ranked$means.ko)
bach2exprs.log2.mean.roraDE.up.sorted <- bach2exprs.log2.mean.roraDE.up[order(bach2exprs.log2.mean.roraDE.up$sort),]
bach2exprs.log2.mean.roraDE.up.sorted <- bach2exprs.log2.mean.roraDE.up.sorted[,1:2]
bach2exprs.log2.mean.roraDE.up.sorted <- as.matrix(bach2exprs.log2.mean.roraDE.up.sorted)

bach2exprs.log2.mean.roraDE.down <- as.data.frame(bach2exprs.log2.mean.roraDE.down)
bach2exprs.log2.mean.roraDE.down$sort <- (bach2exprs.log2.mean.roraDE.down.ranked$means.wt - bach2exprs.log2.mean.roraDE.down.ranked$means.ko)
bach2exprs.log2.mean.roraDE.down.sorted <- bach2exprs.log2.mean.roraDE.down[order(bach2exprs.log2.mean.roraDE.down$sort),]
bach2exprs.log2.mean.roraDE.down.sorted <- bach2exprs.log2.mean.roraDE.down.sorted[,1:2]
bach2exprs.log2.mean.roraDE.down.sorted <- as.matrix(bach2exprs.log2.mean.roraDE.down.sorted)

heatmap.2(bach2exprs.log2.mean.roraDE.up.sorted,
          Colv=FALSE, 
          Rowv = FALSE,
          #dendrogram="row",
          scale="row",
          col="bluered",
          trace="none",
          #ColSideColors=rep(c("green","orange"), each=3),
          main="Bach2 WT vs KO - up in Rora KO geneset",
          xlab="Samples")

heatmap.2(bach2exprs.log2.mean.roraDE.down.sorted,
          Colv=FALSE, 
          Rowv = FALSE,
          #dendrogram="row",
          scale="row",
          col="bluered",
          trace="none",
          #ColSideColors=rep(c("green","orange"), each=3),
          main="Bach2 WT vs KO - down in Rora KO geneset",
          xlab="Samples")

#will plot all the replicates
bach2exprs.log2.roraDE.up <- as.data.frame(bach2exprs.log2.roraDE.up)
bach2exprs.log2.roraDE.up$sort <- (bach2exprs.log2.mean.roraDE.up$means.wt - bach2exprs.log2.mean.roraDE.up$means.ko)
bach2exprs.log2.roraDE.up.sorted <- bach2exprs.log2.roraDE.up[order(bach2exprs.log2.roraDE.up$sort),]
bach2exprs.log2.roraDE.up.sorted <- bach2exprs.log2.roraDE.up.sorted[,1:6]
bach2exprs.log2.roraDE.up.sorted <- as.matrix(bach2exprs.log2.roraDE.up.sorted)

heatmap.2(bach2exprs.log2.roraDE.up.sorted,
          Colv=FALSE, 
          Rowv = FALSE,
          #dendrogram="row",
          scale="row",
          col="bluered",
          trace="none",
          #ColSideColors=rep(c("green","orange"), each=3),
          main="Bach2 WT vs KO - up in Rora KO geneset",
          xlab="Samples")


bach2exprs.log2.roraDE.down <- as.data.frame(bach2exprs.log2.roraDE.down)
bach2exprs.log2.roraDE.down$sort <- (bach2exprs.log2.mean.roraDE.down$means.wt - bach2exprs.log2.mean.roraDE.down$means.ko)
bach2exprs.log2.roraDE.down.sorted <- bach2exprs.log2.roraDE.down[order(bach2exprs.log2.roraDE.down$sort),]
bach2exprs.log2.roraDE.down.sorted <- bach2exprs.log2.roraDE.down.sorted[,1:6]
bach2exprs.log2.roraDE.down.sorted <- as.matrix(bach2exprs.log2.roraDE.down.sorted)

heatmap.2(bach2exprs.log2.roraDE.down.sorted,
          Colv=FALSE, 
          Rowv = FALSE,
          #dendrogram="row",
          scale="row",
          col="bluered",
          trace="none",
          #ColSideColors=rep(c("green","orange"), each=3),
          main="Bach2 WT vs KO - up in Rora KO geneset",
          xlab="Samples")


#clusters

#First clusters and then heat map
#That way you get to extract the clusters as they appear in the heatmap

hc_Clust <- hclust(dist(bach2exprs.log2.mean.roraDE.up))
Clust_names <- cutree(hc_Clust, h=20) ###value for h...This is not the number of
## clusters ....play with it
###use k-means k for defining number of clusters
###look at the cluster tree
# pdf("clustering_tree.pdf")
plot(hc_Clust)
Clust_names <- cutree(hc_Clust, h=20)
rect.hclust(hc_Clust, h=20)
# dev.off()
head (Clust_names)
summary (Clust_names)

#Maybe will do a barcode plot?

bach2exprs.log2.mean.roraDE.up.ranked <- as.data.frame(bach2exprs.log2.mean.roraDE.up)
rank <- ifelse((bach2exprs.log2.mean.roraDE.up.ranked$means.wt - bach2exprs.log2.mean.roraDE.up.ranked$means.ko) < 0, -1, 1)
bach2exprs.log2.mean.roraDE.up.ranked$rank <- rank

bach2exprs.log2.mean.roraDE.down.ranked <- as.data.frame(bach2exprs.log2.mean.roraDE.down)
rank <- ifelse((bach2exprs.log2.mean.roraDE.down.ranked$means.wt - bach2exprs.log2.mean.roraDE.down.ranked$means.ko) < 0, -1, 1)
bach2exprs.log2.mean.roraDE.down.ranked$rank <- rank


