rm(list=ls(all=TRUE))

#Part 1 - load

#load packages
library(dplyr)
library(doParallel)
library(foreach)

#load scripts
source("inputs/deduplicateTWimproved.R")

#index data
geno.input  <-"outputs/dat.geno.with.pheno.csv"
pheno.input <-"outputs/dat.pheno.with.geno.csv"
year.input  <-"inputs/application_year.csv"

#read data
geno <-read.csv(geno.input,header=T,row.names = 1)
pheno<-read.csv(pheno.input,header=T,row.names = 1)
year<-read.csv(year.input, header = T)


# Check if it exists, and create it if not
if (!dir.exists("outputs/predict_candidates")) {
  dir.create("outputs/predict_candidates")
}

#Part 2 - prepare data and settings

# check order of genotypes and phenotypes is the same
if (!identical(row.names(geno), row.names(pheno))) {
  stop("genotypes and phenotypes do not match")
} 

#perform initial marker QC
geno  <- geno[ , sapply(geno, function(x) length(which(x == 0)) >= 12) ]
geno  <- geno[ , sapply(geno, function(x) length(which(x == 1)) >= 12) ]

#reduce duplicates by skimming on 100%
geno <- deduplicate(geno ,1)

#recode markers for prediction
geno <- geno %>% 
  mutate(across(everything(), ~ifelse(.==0,-1,.)))

#save trait list
traits<-names(pheno)

#traits trait.86 and trait.88 have too low phenotype diversity so are excluded
traits<-traits[!traits %in% c("trait.86", "trait.88")]

#set year vector
entry_year<-unique(year$ap.year)
          

#set up parallel
cl <- makeCluster(length(traits))
registerDoParallel(cl)

# Part 3 - run prediction

all_traits <- foreach(i = traits, .packages = c("randomForestSRC")) %dopar% {
  
  geno_filt  <- geno[!is.na(pheno[,i]),]
  pheno_filt <- pheno[!is.na(pheno[,i]),]
  
  #run further marker QC if lost lines for each trait
  if(!nrow(geno_filt)==nrow(geno)){
    #run some initial QC to reduce marker set to useful markers
    #ensure no low MAF in the starting dataframe (some genotypes removed without phenotype data)
    geno_filt  <- geno_filt[ , sapply(geno_filt, function(x) length(which(x == -1)) >= 12) ]
    geno_filt  <- geno_filt[ , sapply(geno_filt, function(x) length(which(x ==  1)) >= 12) ]
    
    #reduce duplicates by skimming on 100%
    geno_filt<- deduplicate(geno_filt ,1)
  }
  
  #save input files for possible reruns
 write.csv(geno_filt,paste0("outputs/predict_candidates/geno_pred_input_",i, ".csv"))
 write.csv(pheno_filt,paste0("outputs/predict_candidates/pheno_pred_input_",i, ".csv"))
  
  #create list for results
  
  results <- list()

  
  for (j in 1:length(entry_year)){
    
  candidates <- year$AFP[year$ap.year==entry_year[j]]
  
  pheno_train <- pheno_filt[!row.names(pheno_filt) %in% candidates, ]
  pheno_test <- pheno_filt[row.names(pheno_filt) %in% candidates, ]
  geno_train <- geno_filt[!row.names(geno_filt) %in% candidates, ]
  geno_test <- geno_filt[row.names(geno_filt) %in% candidates, ]
  
  #control for cases where all test data is NA for a trait
  #will run as long as one genotype has trait data
  if(nrow(pheno_test)==0){
    
    results[[j]] <-"no test genotypes with trait data"
    
  }else{
  
  #set data frames for test and train
  
  rf.train <- data.frame(y = as.factor(pheno_train[,i]),geno_train)
  rf.test  <- data.frame(y = as.factor(pheno_test[,i]),geno_test)
  
  #run rf model as default
  model <- rfsrc(y~.,data= rf.train)
  
  #run predictions
  predicted <- predict(model, rf.test)$class
  
  #assign results to list 
  results[[j]]<-data.frame(Line= row.names(rf.test), 
                         observed = rf.test$y, 
                         predicted = predicted)
  #assign list element a name
  names(results)[j] <- paste0(entry_year[j])
  }
  }
  results
}
stopCluster(cl)  

#set names of all_traits
names(all_traits) <- traits

#save image
save.image(file="outputs/predict_candidates/script_7_data.RData")
library(randomForestSRC)
writeLines(capture.output(sessionInfo()), "outputs/predict_candidates/sessionInfo_script_7.txt")





  