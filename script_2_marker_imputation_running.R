#Part1 - set up
rm(list = ls())

#load packages

library(missForest)
library(doParallel)
library(dplyr)

#set input

to_impute.loc <- "outputs/dat.clean.1.csv"

#read input

to_impute<-read.csv(to_impute.loc, header=T, row.names = 1)

#Part2 - prepare data

#count number of values before edits
table(unlist(to_impute),useNA="always")

#prepare data for imputation

to_impute <- to_impute %>%
  mutate(across(everything(), ~na_if(., 1))) %>%
  mutate(across(everything(), ~ ifelse(. == 2, 1, .))) %>%
  mutate(across(everything(), as.factor))

#counts after edits, 1s to NA and 2s to 1s
table(unlist(to_impute),useNA="always")


#Part3 - impute

#set cores

registerDoParallel(cores = 60)

#run missForest

print("imputation started")
date()

missing.imp.forest<-missForest(to_impute,ntree=150,parallelize = "variables", maxiter = 5, verbose = TRUE)

print("imputation finished")
date()

#finish parallel
stopImplicitCluster()

#Part4 inspect results

#save results to ouput
write.csv(as.data.frame(missing.imp.forest$ximp),"outputs/dat.clean.1.imputed.csv")

#print OOBerror. 0 = good, 1 = bad 
missing.imp.forest$OOBerror

