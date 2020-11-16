rm(list=ls())
setwd("/Users/okamurak/Documents/ExonFrame_BED")
Refseq_GTF <- read.table("hg38.ncbiRefSeq.gtf", sep = "\t", header = FALSE)
## GTF file downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/

## Select forward strand, frame0 exons, remove first exons
Frame <- subset(subset(Refseq_GTF, V8=="0"), V7=="+")
#x<- NULL
#y<- NULL
#i =1

#while (i < length(Frame[,3])) {
  #if (Frame[i,3]=="start_codon") {
   # x <- Frame[i,]
   # y <- rbind(y,x)
   # i=i+1} else {
     # i=i+1}}

#Frame$V10 <- paste(Frame$V1, "_", Frame$V4)
#y$V10 <- paste(y$V1, "_", y$V4)

#NoFirstExon <- Frame[ ! Frame$V10 %in% y$V10, ]
CDS <- subset(Frame, V3=="CDS")

BEDw0 <- CDS[,c(1,4,5)] 
BEDw0[,2] <- BEDw0[,2]-1
BEDw0[,4] <- "Frame0"
BEDw0[,5] <- 0
BEDw0[,6] <- "+"
BEDw0[,7] <- BEDw0[,2]
BEDw0[,8] <- BEDw0[,3]
BEDw0[,9] <- "255,0,0"
BEDw0[,10] <- paste(BEDw0[,1], "_", BEDw0[,2], "_", BEDw0[,3])


## Select reverse strand, frame0 exons, remove first exons
Frame <- subset(subset(Refseq_GTF, V8=="0"), V7=="-")
#x<- NULL
#y<- NULL
#i =1

#while (i < length(Frame[,3])) {
  #if (Frame[i,3]=="start_codon") {
   # x <- Frame[i,]
    #y <- rbind(y,x)
    #i=i+1} else {
     # i=i+1}}

#Frame$V10 <- paste(Frame$V1, "_", Frame$V5)
#y$V10 <- paste(y$V1, "_", y$V5)

#NoFirstExon <- Frame[ ! Frame$V10 %in% y$V10, ]
CDS <- subset(Frame, V3=="CDS")

BEDc0 <- CDS[,c(1,4,5)] 
BEDc0[,2] <- BEDc0[,2]-1
BEDc0[,4] <- "Frame0"
BEDc0[,5] <- 0
BEDc0[,6] <- "-"
BEDc0[,7] <- BEDc0[,2]
BEDc0[,8] <- BEDc0[,3]
BEDc0[,9] <- "200,0,0"
BEDc0[,10] <- paste(BEDc0[,1], "_", BEDc0[,2], "_", BEDc0[,3])

## Select forward strand, frame1 exons, remove first exons
Frame <- subset(subset(Refseq_GTF, V8=="1"), V7=="+")

CDS <- subset(Frame, V3=="CDS")


BEDw1 <- CDS[,c(1,4,5)] 
BEDw1[,2] <- BEDw1[,2]-1
BEDw1[,4] <- "Frame1"
BEDw1[,5] <- 0
BEDw1[,6] <- "+"
BEDw1[,7] <- BEDw1[,2]
BEDw1[,8] <- BEDw1[,3]
BEDw1[,9] <- "0,255,0"
BEDw1[,10] <- paste(BEDw1[,1], "_", BEDw1[,2], "_", BEDw1[,3])

## Select reverse strand, frame1 exons, remove first exons
Frame <- subset(subset(Refseq_GTF, V8=="1"), V7=="-")
CDS <- subset(Frame, V3=="CDS")

BEDc1 <- CDS[,c(1,4,5)] 
BEDc1[,2] <- BEDc1[,2]-1
BEDc1[,4] <- "Frame1"
BEDc1[,5] <- 0
BEDc1[,6] <- "-"
BEDc1[,7] <- BEDc1[,2]
BEDc1[,8] <- BEDc1[,3]
BEDc1[,9] <- "0,200,0"
BEDc1[,10] <- paste(BEDc1[,1], "_", BEDc1[,2], "_", BEDc1[,3])


## Select forward strand, frame1 exons, remove first exons
Frame <- subset(subset(Refseq_GTF, V8=="2"), V7=="+")

CDS <- subset(Frame, V3=="CDS")


BEDw2 <- CDS[,c(1,4,5)] 
BEDw2[,2] <- BEDw2[,2]-1
BEDw2[,4] <- "Frame2"
BEDw2[,5] <- 0
BEDw2[,6] <- "+"
BEDw2[,7] <- BEDw2[,2]
BEDw2[,8] <- BEDw2[,3]
BEDw2[,9] <- "0,0,255"
BEDw2[,10] <- paste(BEDw2[,1], "_", BEDw2[,2], "_", BEDw2[,3])

## Select reverse strand, frame2 exons, remove first exons
Frame <- subset(subset(Refseq_GTF, V8=="2"), V7=="-")
CDS <- subset(Frame, V3=="CDS")

BEDc2 <- CDS[,c(1,4,5)] 
BEDc2[,2] <- BEDc2[,2]-1
BEDc2[,4] <- "Frame2"
BEDc2[,5] <- 0
BEDc2[,6] <- "-"
BEDc2[,7] <- BEDc2[,2]
BEDc2[,8] <- BEDc2[,3]
BEDc2[,9] <- "0,0,200"
BEDc2[,10] <- paste(BEDc2[,1], "_", BEDc2[,2], "_", BEDc2[,3])



BED <- NULL
BED <- rbind(BEDw0, BEDc0, BEDw1, BEDc1, BEDw2, BEDc2)
BED <- BED[ !grepl("_fix", BED$V10), ]
BED <- BED[ !grepl("_alt", BED$V10), ]
BED2 <- BED[,1:9]




write.table(BED2, sep="\t", file="readingframe.bed", row.names=F, col.names=F)

## after this, remove double quotes, and add the following line at the beginning
## track type=bigBed name="ReadingFrame" description="Reading frame of exon" itemRgb="On" bigDataUrl=https://github.com/KatsutomoNAIST/UCSCtrack/raw/master/readingframe.bb







