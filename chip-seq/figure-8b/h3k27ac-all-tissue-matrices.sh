#!/bin/bash
for f in *.bed
do
  echo "Processing $f file..."
  computeMatrix scale-regions -p 32 -S '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-1.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-2.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-3.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-1.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-2.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-3.fltrd.bw' -R "$f" -b 5000 --skipZeros -o "processed-$f" --outFileNameMatrix "$f.tab"
done
