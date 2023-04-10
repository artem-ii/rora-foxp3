#!/bin/bash

mkdir tfpromoters
mkdir tfpromoters/heatmaps
for f in *.bed
do
  echo "Processing $f file..."
  # intersect Teichmann promoters with TF promoters
  bedtools intersect -wa -a "/media/serlia/Storage 2/artemii/results-chipseq/29-october-2020/TF-promoters-mouse_epdnew_CvOu1-5000TSS.bed" -b "$f" > "$f-tf-promoters.bed"
  mv "$f-tf-promoters.bed" "tfpromoters/$f-tf-promoters.bed"
  cd tfpromoters
  computeMatrix scale-regions -p 32 -S '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-1.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-2.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-3.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-1.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-2.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-3.fltrd.bw' -R "$f-tf-promoters.bed" -b 5000 --skipZeros -o "$f-h3k27ac-matrix" --outFileNameMatrix "$f-h3k27ac-matrix.tab"
  computeMatrix scale-regions -p 32 -S '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-k27-1_S7_R1_001.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-k27-2_S8_R1_001.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-k27-3_S9_R1_001.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-k27-1_S10_R1_001.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-k27-2_S11_R1_001.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-ko-k27-3_S12_R1_001.fltrd.bw' -R "$f-tf-promoters.bed" -b 5000 --skipZeros -o "$f-h3k27me3-matrix" --outFileNameMatrix "$f-h3k27me3.tab"
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip1_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27ac-dn-$f-tf-promoters.bed" 
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip2_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27ac-up-$f-tf-promoters.bed" 
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip1_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27me3-dn-$f-tf-promoters.bed" 
  bedtools intersect -wa -wb -a '/media/serlia/Storage 2/artemii/rora-treg-chipseq/B2306_2307/fastq/treg-fastq-cat/filtered-indexed-bams/pepr-results/h3k27ac-db-pepr__PePr_chip2_peaks.bed' -b "$f-tf-promoters.bed" > "db-h3k27me3-up-$f-tf-promoters.bed"
  plotHeatmap -m "$f-h3k27me3-matrix" -o "heatmaps/$f-h3k27me3-heatmap.svg" --kmeans 12 --outFileSortedRegions "heatmaps/$f-h3k27me3-sorted-regions.bed"
  plotHeatmap -m "$f-h3k27ac-matrix" -o "heatmaps/$f-h3k27ac-heatmap.svg" --kmeans 12 --outFileSortedRegions "heatmaps/$f-h3k27ac-sorted-regions.bed"
  cd ..
done
