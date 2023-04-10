#!/bin/bash
mkdir tfpromoters
mkdir promoters
for f in p-*
do
  echo "Processing $f file..."
  # intersect Dispirito promoters with TF promoters
  bedtools intersect -wa -a "/media/serlia/Storage 2/artemii/results-chipseq/29-october-2020/TF-promoters-mouse_epdnew_CvOu1-5000TSS.bed" -b $f > "$f-tf-promoters.bed"
  mv "$f-tf-promoters.bed" "tfpromoters/$f-tf-promoters.bed"
  cd tfpromoters
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip1_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27ac-dn-$f-tf-promoters.bed" 
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip2_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27ac-up-$f-tf-promoters.bed" 
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip1_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27me3-dn-$f-tf-promoters.bed" 
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip2_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27me3-up-$f-tf-promoters.bed"
  cd ..
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip1_peaks.bed' -b $f > "db-h3k27ac-dn-$f-promoters.bed"
  mv "db-h3k27ac-dn-$f-promoters.bed" "promoters/db-h3k27ac-dn-$f-promoters.bed"
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip2_peaks.bed' -b $f > "db-h3k27ac-up-$f-promoters.bed"
  mv "db-h3k27ac-up-$f-promoters.bed" "promoters/db-h3k27ac-up-$f-promoters.bed"
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip1_peaks.bed' -b $f > "db-h3k27me3-dn-$f-promoters.bed"
  mv "db-h3k27me3-dn-$f-promoters.bed" "promoters/db-h3k27me3-dn-$f-promoters.bed"
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip2_peaks.bed' -b $f > "db-h3k27me3-up-$f-promoters.bed"
  mv "db-h3k27me3-up-$f-promoters.bed" "promoters/db-h3k27me3-up-$f-promoters.bed"
done
echo "merging db regions"
cd tfpromoters
for f in db-*
do
  echo "$f" > temp
  cat temp $f >> all-db-tf-promoters-dispirito.txt
  rm temp
done
cd ..
cd promoters
for f in db-*
do
  echo "$f" > temp
  cat temp $f >> all-db-promoters-dispirito.txt
  rm temp
done
