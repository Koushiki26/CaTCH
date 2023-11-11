#!/bin/bash

cat hg38_phastCons100way.bed | awk '$1=="chr1"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr1.txt
cat hg38_phastCons100way.bed | awk '$1=="chr2"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr2.txt
cat hg38_phastCons100way.bed | awk '$1=="chr3"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr3.txt
cat hg38_phastCons100way.bed | awk '$1=="chr4"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr4.txt
cat hg38_phastCons100way.bed | awk '$1=="chr5"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr5.txt
cat hg38_phastCons100way.bed | awk '$1=="chr6"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr6.txt
cat hg38_phastCons100way.bed | awk '$1=="chr7"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr7.txt
cat hg38_phastCons100way.bed | awk '$1=="chr8"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr8.txt
cat hg38_phastCons100way.bed | awk '$1=="chr9"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr9.txt
cat hg38_phastCons100way.bed | awk '$1=="chr10"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr10.txt
cat hg38_phastCons100way.bed | awk '$1=="chr11"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr11.txt
cat hg38_phastCons100way.bed | awk '$1=="chr12"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr12.txt
cat hg38_phastCons100way.bed | awk '$1=="chr13"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr13.txt
cat hg38_phastCons100way.bed | awk '$1=="chr14"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr14.txt
cat hg38_phastCons100way.bed | awk '$1=="chr15"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr15.txt
cat hg38_phastCons100way.bed | awk '$1=="chr16"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr16.txt
cat hg38_phastCons100way.bed | awk '$1=="chr17"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr17.txt
cat hg38_phastCons100way.bed | awk '$1=="chr18"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr18.txt
cat hg38_phastCons100way.bed | awk '$1=="chr19"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr19.txt
cat hg38_phastCons100way.bed | awk '$1=="chr20"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr20.txt
cat hg38_phastCons100way.bed | awk '$1=="chr21"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr21.txt
cat hg38_phastCons100way.bed | awk '$1=="chr22"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chr22.txt
cat hg38_phastCons100way.bed | awk '$1=="chrX"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chrX.txt
cat hg38_phastCons100way.bed | awk '$1=="chrY"{print $2, $4}' > ~/Documents/koushiki/project/conservation/chrY.txt
