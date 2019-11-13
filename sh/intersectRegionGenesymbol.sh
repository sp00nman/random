#!/bin/bash

FILES=$1

cat $FILES | while read sample
do
	echo "Processing $sample"
	# careful sometimes "chr" is need, sometimes not
	awk 'BEGIN{OFS=""}{print $2,"\t", $3,"\t",$4,"\t",$5,"\t",$6}' $sample.txt >$sample.bed
	
	#remove first line
	echo "$(tail -n +2 $sample.bed)" > $sample.bed
	
	#liftover
	#liftOver -bedPlus=3 $sample.bed /home/schischlikf2/resources/hg19ToHg38.over.chain.gz $sample-Hg38.bed $sample-Hg38.bed_unmapped
	
	#intersect bed
	bedtools intersect \
	-a /data/Lab_ruppin/TARGET/OS/gencode.v19.annotation.bed \
	-b $sample.bed \
	-wao \
	-f 1.0 >$sample-intersect_cn.bed
	
	#replace -1 values with 0  $12 & $11
	awk 'BEGIN{OFS="\t"}{if($12==0) {$11=0; print $0} else {$11=$11; print $0}}' $sample-intersect_cn.bed >$sample-repl-intersect_cn.bed
	#extract column $11 or $19
	awk '{print $11}' $sample-repl-intersect_cn.bed >$sample-intvalue-intersect_cn.bed 
	
done
