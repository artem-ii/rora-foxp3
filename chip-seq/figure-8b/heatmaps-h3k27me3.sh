#!/bin/bash
for f in h3k27me3-matrices/*
do
  echo "Processing $f file..."
  plotHeatmap -m $f -o $f.svg --kmeans 12 --outFileSortedRegions $f.bed
done
