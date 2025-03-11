#!/bin/bash
echo 
echo ---------------------------------
echo find_treg_tf_targets.sh run at
date 
echo ---------------------------------
echo 
echo Script requires R 4.4.0, bedtools 2.31.1 and python3.10 with pandas and numpy installed
echo
echo ---------------------------------
R --version
echo 
python3 --version
echo 
bedtools --version
echo ---------------------------------

MY_PATH="/Users/artemii/rorafoxp3_2023/code/rora-foxp3/rora-foxp3/rora_targets_filtering_pipeline"

ROOT_PATH=MY_PATH # !TODO: Assign path to rora_targets_filtering_pipeline directory on your machine here

PROCESSED_DATA_PATH="$ROOT_PATH/processed_data"
SCRIPTS_PATH="$ROOT_PATH/scripts"
PROMOTERS_PATH="$PROCESSED_DATA_PATH/miragaia_epd_promoters"
EXTENDED_PROMOTERS_PATH="$PROCESSED_DATA_PATH/miragaia_extended_promoter_regions"
DATA_PATH="$ROOT_PATH/data"
RORA_TREG_HISTONE_DB_PATH="$DATA_PATH/rora_ko_treg_histone_db_pepr"
H3K27AC_UP_RORA_TREG_BEDFILE="$RORA_TREG_HISTONE_DB_PATH/h3k27ac-db-pepr__PePr_chip2_peaks.bed"
H3K27AC_DOWN_RORA_TREG_BEDFILE="$RORA_TREG_HISTONE_DB_PATH/h3k27ac-db-pepr__PePr_chip1_peaks.bed"
H3K27ME3_UP_RORA_TREG_BEDFILE="$RORA_TREG_HISTONE_DB_PATH/h3k27me3-db-pepr__PePr_chip2_peaks.bed"
H3K27ME3_DOWN_RORA_TREG_BEDFILE="$RORA_TREG_HISTONE_DB_PATH/h3k27me3-db-pepr__PePr_chip1_peaks.bed"

OUTPUT_EXTENDED_PROMOTER_REGIONS="$PROCESSED_DATA_PATH/miragaia_epd_promoters"
OUTPUT_HISTONE_DB_PROMOTERS_PATH="$PROCESSED_DATA_PATH/miragaia_promoters_with_histone_db"
MERGED_TABLE_PATH=$PROCESSED_DATA_PATH/merged_promoter_db_table

echo; echo
echo ---------------------------------
echo Step 1: Finding chromatin changes in RORA KO Treg relevant for Treg identity and function in non-lymphoid tissues
echo ---------------------------------
echo
echo To filter genes potentially regulated by RORA in Treg and important for controlling Treg fate,
echo identity and function in the context of non-lymphoid tissues we used marker gene lists
echo of Treg and Tmem cells in colon, skin and draining lymph nodes.
echo Data were taken from Miragaia et al. 2019. See data/miragaia_et_al/tissue_treg_tmem_genesets/mmc4.xlsx
echo and data/miragaia_et_al/source_article for the publication.
echo
echo Promoter coordinates of tissue Treg and Tmem marker genes were found using EPDnew database
echo ±10kb from promoter coordinates will be checked for changed histone marks between WT and RORA KO Treg
echo
echo Executing extend_promoter_regions.R to extend promoter regions...
echo ---------------------------------
echo

$SCRIPTS_PATH/extend_promoter_regions.R

echo 
echo "Moving files to output directory..."
echo ---------------------------------
echo

mkdir $EXTENDED_PROMOTERS_PATH
cd $PROMOTERS_PATH
mv extended* $EXTENDED_PROMOTERS_PATH

echo
echo ---------------------------------
echo Step 2: Intersecting RORA KO Treg changed histone mark peaks with extended promoter regions...
echo ---------------------------------
echo 
mkdir $OUTPUT_HISTONE_DB_PROMOTERS_PATH
cd $EXTENDED_PROMOTERS_PATH
echo Executing bedtools intersect scripts...
echo Finding differential RORA KO Treg histone mark signals in
echo
for f in *.bed
do
  echo $f
  # intersect Miragaia et al promoters with db sites
  # echo "H3K27Ac..."
  bedtools intersect -wa -wb -a $H3K27AC_UP_RORA_TREG_BEDFILE -b $f > "$OUTPUT_HISTONE_DB_PROMOTERS_PATH/db-h3k27ac-dn-$f-intersect.bed" 
  bedtools intersect -wa -wb -a $H3K27AC_DOWN_RORA_TREG_BEDFILE -b $f > "$OUTPUT_HISTONE_DB_PROMOTERS_PATH/db-h3k27ac-up-$f-intersect.bed" 
  # echo "H3K27me3..."
  bedtools intersect -wa -wb -a $H3K27ME3_UP_RORA_TREG_BEDFILE -b $f > "$OUTPUT_HISTONE_DB_PROMOTERS_PATH/db-h3k27me3-dn-$f-intersect.bed" 
  bedtools intersect -wa -wb -a $H3K27ME3_DOWN_RORA_TREG_BEDFILE -b $f > "$OUTPUT_HISTONE_DB_PROMOTERS_PATH/db-h3k27me3-up-$f-intersect.bed"
done

cd $OUTPUT_HISTONE_DB_PROMOTERS_PATH

for f in *-intersect.bed
#echo "Processing $f"
do
  echo $f >> merged_promoter_db_table_all_tissues.tsv
  cat $f >> merged_promoter_db_table_all_tissues.tsv
done

mv merged_promoter_db_table_all_tissues.tsv $MERGED_TABLE_PATH/merged_promoter_db_table_all_tissues.tsv

echo
echo Done finding differential binding in extended promoter regions
echo
echo For full list of histone mark changes found in Miragaia et al. marker genes with tissue and cell type annotation
echo see processed_data/merged_promoter_db_table/merged_promoter_db_table_all_tissues.tsv

echo
echo ---------------------------------
echo Executing gene_selection_rora_treg.py for further filtering steps...
echo
$SCRIPTS_PATH/gene_selection_rora_treg.py
echo
echo ---------------------------------
echo Done filtering
echo ---------------------------------
echo 
echo find_treg_tf_targets.sh executed
echo
