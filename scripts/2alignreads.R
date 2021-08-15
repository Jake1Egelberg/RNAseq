#***********************************************************
#*************************RNA SEQ **************************
#***********************************************************

#---------------------LOADING PARMS----------------------

#Choose primary workflow file path
file.path<-"D:\\RNAseq"

#Load user-set parms file
parms<-read.delim(paste(file.path,"\\parms.txt",sep=""),sep=":")

#Redefine parms for R
index.file<-trimws(parms[which(parms$RNA_SEQ_PARAMETERS=="index.file"),2])
paired.end.status<-as.logical(trimws(parms[which(parms$RNA_SEQ_PARAMETERS=="paired.end.status"),2]))
ref.genome<-trimws(parms[which(parms$RNA_SEQ_PARAMETERS=="ref.genome"),2])

#Load packages
library(BiocManager)
library(Rsubread)
library(stringr)
set.seed(42)

#Create text file to update user
update<-data.frame(Update="Status")

#Remove existing progress files
progress.files<-list.files(path=paste(file.path,"\\progress",sep=""),full.names = TRUE)
file.remove(progress.files)

setwd(paste(file.path,"\\progress",sep=""))
write.table(update,"ALIGNING READS.txt")

#---------------------ALIGN READS----------------------

#Identify fastq files in 1fastqfiles folder
fastq.files <- list.files(path = paste(file.path,"\\1fastqfiles\\",sep=""), 
                          pattern = ".fastq.gz$", 
                          full.names = TRUE)

if(paired.end.status==TRUE){
  fastq.files.pair <- list.files(path = paste(file.path,"\\1fastqfiles\\pair",sep=""), 
                                 pattern = ".fastq.gz$", 
                                 full.names = TRUE)
}

#Align reads to index
if(paired.end.status==FALSE){
  align(index = paste(file.path,"\\buildindex\\",index.file,sep=""),
        readfile1 = fastq.files)
} else if(paired.end.status==TRUE){
  align(index = paste(file.path,"\\buildindex\\",index.file,sep=""),
        readfile1 = fastq.files,
        readfile2 = fastq.files.pair)
}

setwd(paste(file.path,"\\progress",sep=""))
write.table(update,"ALIGNMENT COMPLETE.txt")

