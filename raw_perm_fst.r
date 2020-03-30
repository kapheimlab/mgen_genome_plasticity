#!/usr/bin/Rscript


#libraries --------------------------------------------------------------------
library(IRanges)
library(plyr)
library(plotrix)
require(ggplot2)
require(ggthemes)
require(plyr)
require(wesanderson)
require(reshape)
require(scales)
require(grid)
library(tidyr)



#functions --------------------------------------------------------------------
fst.trim <- function(fst){
  names(fst)[5] = "fst"
  temp = fst$fst
  temp[temp<0] = 0
  fst$fst = temp
  return(fst)
}



##########################
# Permute Fst 
size1 <- 9
size2 <- 9
for(i in c(1:100)){
targ <- paste(sample(c(0,1,2,3,4,5,10,14,16,6,7,8,9,11,12,13,15,17),size1),collapse = ",") #sample positions to be taken from VCF file
back <- paste(sample(c(0,1,2,3,4,5,10,14,16,6,7,8,9,11,12,13,15,17),size2),collapse = ",")

fil <- paste("pFSTperm", i, sep="")
system(paste("./wcFst --target", targ, "--background", back, "--deltaaf 0.05 --file filtered_excl_snps_r4_v2_5_4_05.recode.vcf  --type PL  >", fil, sep = " "))
}

#add snp column
for(i in c(1:100)){
	system(paste("sed -i 's/\t/_/'", fil))

}


##########################
# Run Fst 
samps1 <- c()
samps2 <- c()
targ <- paste(samps1,collapse = ",") 
back <- paste(samps2,collapse = ",") 
fil <- 'wcFST0NOMISST.out'
system(paste("./pFst --target", targ, "--background", back, "--deltaaf 0.05 --file filtered_excl_snps_r4_v2_5_4_05.recode.vcf  --type PL  >", fil, sep = " "))



samps1 <- c()
samps2 <- c()
targ <- paste(samps1,collapse = ",") 
back <- paste(samps2,collapse = ",") 
fil <- 'pFST0NOMISST.out'
system(paste("./wcFst --target", targ, "--background", back, "--deltaaf 0.05 --file filtered_excl_snps_r4_v2_5_4_05.recode.vcf  --type PL  >", fil, sep = " "))



















#  #######################
#  #!/bin/bash
#  directory="/media/data1/megalop/mgen_pop/New"
#  if [ ! -d $directory ]; then
#    echo "Error: Directory doesn't exist"
#    exit 1
#  fi
#
#  for file in $directory/pFSTperm*
#  do
#    if [ -f $file ]; then
#        cat $file | sort > $file.tmp
#        mv -f  $file.tmp $file
#    fi
#  done
#
#  directory="/media/data1/megalop/mgen_pop/New"
#  data=
#  for f in "${directory}"/pFSTperm*
#  do
#    if [ ! "$data" ]
#    then
#      data="$(sort "$f")"
#      continue
#    fi
#    data="$(join <(sort "$f") /dev/stdin <<< "$data")"
#  done
#  echo "$data"





#  #concatenate files 
#  sort pFSTperm1 > pFSTperm1s
#  sort pFSTperm2 > pFSTperm2s

#  join pFSTperm1s pFSTperm2s > test








#  i = 1
# fil = paste("/media/data1/megalop/mgen_pop/New/pFSTperm", i, sep="")
#  fil = read.table(fil, header=F)
#  fil = fil[c(1,4)]
#  temp = fil$V4
#  temp[which(temp<0)] = 0
#  fil$V4 = temp
#  df = fil
#
#
#  for(i in c(3:10)){
#  	fil = paste("/media/data1/megalop/mgen_pop/New/pFSTperm", i, sep="")
#  	fil = read.table(fil, header=F, colClasses = rep("character",4))
#  	fil = fil[c(1,4)]
#  	temp = fil$V4
#  	temp[which(as.numeric(temp)<0)] = 0
#  	fil$V4 = temp
#  	df = join(df, fil, by = "V1")
#  	print(i)
#
#  }





