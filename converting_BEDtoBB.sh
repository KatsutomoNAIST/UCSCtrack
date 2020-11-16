#!/bin/sh

#  converting_BEDtoBB.sh
#  
#
#  Created by okamurak on 2020/11/13.
#

## Input bed file created by WritingFrameBED.R

cd /Users/okamurak/Documents/ExonFrame_BED

## Sort the bed file by coordinate
sort -k1,1 -k2,2n readingframe.bed > readingframe.sorted.bed

##convert bed to bigbed
bedToBigBed readingframe.sorted.bed hg38.chrom.sizes readingframe.bb

## track type=bigBed name="ReadingFrame" description="Reading frame of exon" itemRgb="On" bigDataUrl=https://github.com/KatsutomoNAIST/UCSCtrack/raw/master/readingframe.bb

