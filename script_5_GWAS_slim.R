rm(list=ls(all=TRUE))
#Part 1 load

#load packages

library(ggplot2)
library(ggpubr)
library(GWASpoly)

#load scripts
source("inputs/deduplicateTWimproved.R")
source("inputs/GWASpoly_manhattan.R")
source("inputs/GWASpoly_QQ.R")
source("inputs/GWASpoly_trait_hist.R")
        
#set inputs
        
geno.input<-"outputs/dat.geno.with.pheno.csv"
pheno.input<-"outputs/dat.pheno.with.geno.csv"
map.input <- "inputs/map_TWedit.csv"
        
#read data
geno <-read.csv(geno.input,header=T,row.names = 1)
pheno<-read.csv(pheno.input,header=T,row.names = 1)
map  <-read.csv(map.input,header=T)   

#Part 2 - data prep

#form traits vector
traits<-names(pheno)

#traits trait.86 and trait.88 have too low phenotype diversity so are excluded
traits<-traits[!traits %in% c("trait.86", "trait.88")]

#remove unmapped markers - removes 122 markers
unmapped<-geno[ ,!names(geno) %in% map$Marker ]
geno<-geno[ ,!names(geno) %in% names(unmapped)]

#run some initial QC to reduce marker set to useful markers
#ensure no low MAF in the starting dataframe (some genotypes removed without phenotype data)
geno  <- geno[ , sapply(geno, function(x) length(which(x == 0)) >= 12) ]
geno  <- geno[ , sapply(geno, function(x) length(which(x == 1)) >= 12) ]

#reduce duplicates by skimming on 100%
geno_filt.1 <- deduplicate(geno ,1)

#take full marker set for storing score
gwas.scores<- data.frame(Marker=names(geno_filt.1))
gwas.scores  <- merge(gwas.scores,map,by="Marker")
gwas.scores  <- gwas.scores[with(gwas.scores, order(Chr, Pos)), ]

#form empty major hits table
hits <- data.frame(Trait = as.character(),
           Model=as.character(),
           Threshold=as.numeric(), 
           Marker=as.character(), 
           Chrom = as.character(), 
           Position = as.numeric(),
           Ref=as.numeric(), 
           Alt = as.numeric(), 
           Score=as.numeric(), 
           Effect=as.numeric())

#form table for storing sig thresholds and number of markers

thresholds<-data.frame(Trait = as.character(),
           Threshold=as.numeric(),
           n.marker=as.integer()) 
           
#Part 3 - GWAS for each trait

for(i in traits){

  
  # check order of genotypes and phenotypes is the same
    if (!identical(row.names(geno_filt.1), row.names(pheno))) {
      stop("genotypes and phenotypes do not match")
    }  
  
    #reduce genotype data down to those with phenotypes for each trait
    geno_file  <- geno_filt.1[!is.na(pheno[,i]),]
    pheno_file <- pheno[!is.na(pheno[,i]),]
    
    #if NA in phenos has reduced genotype number, repeat QC
    
    if (!nrow(geno_file)==nrow(geno_filt.1)){
      #ensure no low MAF in the starting dataframe (some genotypes removed without phenotype data)
      geno_file  <- geno_file [ , sapply(geno_file , function(x) length(which(x == 0)) >= 12) ]
      geno_file  <- geno_file [ , sapply(geno_file , function(x) length(which(x == 1)) >= 12) ]
      
      #reduce duplicates by skimming on 100%
      geno_file<- deduplicate(geno_file, 1)
      
    }
    
    # check order of genotypes and phenotypes is the same
    if (!identical(row.names(geno_file), row.names(pheno_file))) {
      stop("genotypes and phenotypes do not match")
    } 
    
    #form id col
    geno.names <- row.names(geno_file)
    
    #form geno file
    geno_file  <- as.data.frame(t(geno_file))
    geno_file  <- merge(map,geno_file, by.x="Marker", by.y ="row.names")
    geno_file  <- geno_file[with(geno_file, order(Chr, Pos)), ]
    
    # Check if it exists, and create it if not
    if (!dir.exists("outputs/GWASouts")) {
      dir.create("outputs/GWASouts")
    }
    #write genotype file
    write.csv(geno_file,paste0("outputs/GWASouts/gwaspoly.",i,".geno.csv"), row.names=F)
      
    #form phenotype file
    pheno_file <- data.frame(id =geno.names ,pheno_file[,i])
    names(pheno_file)[2] <- i
    write.csv(pheno_file,paste0("outputs/GWASouts/gwaspoly.",i,".phen.csv"), row.names=F)
      
    #read GWAS object
    #read data
    
    gpoly <- read.GWASpoly(ploidy=1, pheno.file=paste0("outputs/GWASouts/gwaspoly.",i,".phen.csv"), 
                           geno.file=paste0("outputs/GWASouts/gwaspoly.",i,".geno.csv"),
                            format="numeric", n.traits=1, delim=",")
      
    #set K
    gpoly <- set.K(gpoly,LOCO=F)
      
    #run GWAS
    gpoly.scan <- GWASpoly(data=gpoly,
                                     models="additive", params = NULL )
    #set thresh
    gpoly.scan.bonf <- set.threshold(gpoly.scan,method="Bonferroni",level=0.05)
    
    #control inf values which are common for some DUS traits
    gpoly.scan.bonf@scores[[i]]$additive[
      is.infinite(gpoly.scan.bonf@scores[[i]]$additive)] <- 307.6527#min p
    
    #run plotting functions
    GWASpoly_manhattan(gpoly.scan.bonf , i, "outputs/GWASouts")
    GWASpoly_QQ(gpoly.scan.bonf , i, "outputs/GWASouts")
    GWASpoly_trait_hist(gpoly.scan.bonf , i, "outputs/GWASouts")
      
    #save scores
    scores<-as.data.frame(gpoly.scan.bonf@scores) 
    names(scores) <- i
    gwas.scores<-merge(gwas.scores,scores,
                       by.x="Marker", by.y=0, all=T )
    gwas.scores  <- gwas.scores[with(gwas.scores, order(Chr, Pos)), ]
    
    #save thresholds
    thresholds <- rbind(thresholds, data.frame(
      Trait=i,
      Threshold=gpoly.scan.bonf@threshold[1],
      n.marker=length(gpoly.scan.bonf@scores[[1]]$additive)
      ))
    
    
    #save hits across window of 15mb and above score of 10
    hits.temp<-get.QTL(data=gpoly.scan.bonf,
                       traits=i,
                       models="additive",bp.window=15e6)
    hits.temp<-hits.temp[hits.temp$Score>10,]
    hits <- rbind(hits,hits.temp)
}

#write final outputs
write.csv(hits,"outputs/GWASouts/gwas.major.loci.csv")
write.csv(gwas.scores,"outputs/GWASouts/gwas.all.scores.csv")
write.csv(thresholds ,"outputs/GWASouts/gwas.thresholds.csv")
 


