#! /bin/bash/

wiggletools write h3k27ac-wt-mean.wig mean '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-1.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-2.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/wt-h3k27ac-3.fltrd.bw'

wiggletools write h3k27ac-ko-mean.wig mean  '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-1.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-2.fltrd.bw' '/media/serlia/Storage 2/artemii/results-chipseq/histones/bigwig/ko-h3k27ac-3.fltrd.bw' 

wigToBigWig h3k27ac-wt-mean.wig http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes h3k27ac-wt-mean.bw

wigToBigWig h3k27ac-ko-mean.wig http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes h3k27ac-ko-mean.bw

#it actually subtracts the second file (b2) from the first one
bigwigCompare -p 32 -b1 h3k27ac-ko-mean.bw -b2 h3k27ac-wt-mean.bw --operation "subtract" -o subtract-h3k27ac-signal -of "bigwig"

bigwigCompare -p 32 -b1 h3k27ac-ko-mean.bw -b2 h3k27ac-wt-mean.bw --operation "ratio" -o ratio-h3k27ac-signal -of "bigwig"

computeMatrix reference-point -p 32 -S 'subtract-h3k27ac-signal.bw' -R '/media/serlia/Storage 2/artemii/results-chipseq/2-may-2021/bach2-chipseq-1120740-top-3000.bed' -bs 50 --referencePoint center -b 3000 -a 3000 --skipZeros -o "processed-h3k27ac-bach2-ref-point-3000-1" --outFileNameMatrix "bach2-h3k27ac-ref-point-3000-1"


  
plotHeatmap -m "processed-h3k27ac-bach2-ref-point-3000" -o "bach2-h3k27ac-ref-point-3000-kmeans12.svg" --kmeans 12 --outFileSortedRegions "bach2-h3k27ac.tab-ref-point-3000-kmeans12.bed"

# no clustering, but sort
plotHeatmap -m "processed-h3k27ac-bach2-ref-point-3000" -o "bach2-h3k27ac-ref-point-3000-no-cluster.svg" --outFileSortedRegions "bach2-h3k27ac.tab-ref-point-3000-no-cluster.bed"

# 3 clusters
plotHeatmap -m "processed-h3k27ac-bach2-ref-point-3000" -o "bach2-h3k27ac-ref-point-3000-kmeans3.svg" --kmeans 3 --outFileSortedRegions "bach2-h3k27ac.tab-ref-point-3000-kmeans3.bed"

#three clusters and also will plot BACH2 signal
computeMatrix reference-point -p 32 -S 'subtract-h3k27ac-signal.bw' '/media/serlia/Storage 2/artemii/results-chipseq/2-may-2021/Bach2-gsm1120740-converted-to-mm10.bw' -R '/media/serlia/Storage 2/artemii/results-chipseq/2-may-2021/bach2-chipseq-1120740-top-3000.bed' -bs 50 --referencePoint center -b 3000 -a 3000 --skipZeros -o "processed-h3k27ac-bach2-ref-point-3000" --outFileNameMatrix "bach2-h3k27ac-ref-point-3000"

computeMatrix reference-point -p 32 -S 'ratio-h3k27ac-signal.bw' '/media/serlia/Storage 2/artemii/results-chipseq/2-may-2021/Bach2-gsm1120740-converted-to-mm10.bw' -R '/media/serlia/Storage 2/artemii/results-chipseq/2-may-2021/bach2-chipseq-1120740-top-3000.bed' -bs 50 --referencePoint center -b 3000 -a 3000 --skipZeros -o "processed-h3k27ac-bach2-ref-point-3000-ratio" --outFileNameMatrix "bach2-h3k27ac-ref-point-3000-ratio"

plotHeatmap -m "processed-h3k27ac-bach2-ref-point-3000" -o "bach2-h3k27ac-ref-point-3000-kmeans3-with-bach2.svg" --colorMap Greys Reds --kmeans 3 --outFileSortedRegions "bach2-h3k27ac.tab-ref-point-3000-kmeans3.bed"

plotHeatmap -m "processed-h3k27ac-bach2-ref-point-3000-ratio" -o "bach2-h3k27ac-ref-point-3000-kmeans3-with-bach2-ratio.svg" --colorMap Greys Reds --kmeans 3 --outFileSortedRegions "bach2-h3k27ac.tab-ref-point-3000-kmeans3-ratio.bed"
