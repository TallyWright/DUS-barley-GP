#Part1 - set up
rm(list=ls(all=TRUE))

#load packages

library(tidyverse)

#set inputs

geno.input  <- "outputs/dat.clean.1.imputed.csv" 
pheno.input <- "inputs/pheno_input_TWedit.csv" 

#form name conversion (done in excel)
#write.csv(data.frame(geno.names=row.names(geno)),
          #"geno_to_pheno_name_convert.csv", row.names = F)

name.convert <- "inputs/geno_to_pheno_name_convert.csv"


#read inputs

pheno<-read.csv(pheno.input,header=T, row.names = 1)
geno<-read.csv(geno.input,header=T, row.names = 1)
name.convert<-read.csv(name.convert, header=T)

#Part2 - fix names

#add AFP to pheno names

row.names(pheno) <- paste("AFP",row.names(pheno),sep="")

#convert genotype names so they would match phenotype names

if(identical(row.names(geno),name.convert$geno.names)){
  row.names(geno) <- name.convert$geno.name.pheno
  print("names converted")
}

#Part3 - run edits to phenotype data

pheno <- pheno %>%
  #mark 0 as NA 
  mutate(across(everything(), ~ ifelse(. == 0, NA, .))) %>% 
  #These traits are scored on 1-9 but should be treated as 1-2 to be fair in the analysis 
  mutate(across(c(trait.2, trait.86, trait.76), ~ ifelse(. == 9, 2, .))) %>%
  #change '2' alternative season habit scores to NA
  mutate(trait.90= ifelse(trait.90 == 2, NA, trait.90)) %>%
  #change '3' for trait.90 scores to '2' to be fair in the analysis
  mutate(trait.90= ifelse(trait.90 == 3, 2, trait.90)) %>%
  #Remove the extra trait that is mostly NAs
  dplyr::select(-trait.95)%>%
  #mark the six rows with non "255" calls as "255" for trait 22, NA control added
  mutate(trait.22.edit = ifelse(!is.na(trait.39) & trait.39==2 & trait.22.edit < 255, 255, trait.22.edit))%>%
  #mark the six rows with non "255" calls as "255" for trait 20, NA control added
  mutate(trait.20 = ifelse(!is.na(trait.39) & trait.39 == 2 & trait.20 < 255, 255, trait.20)) %>%
  #then mark all 255 as NA, e.g. trait.20 has 255 calls as it's not relevant to six row
  mutate(across(everything(), ~ ifelse(. == 255, NA, .)))%>%
  #change name of 22,30,41 to keep names simple
  rename_with(~ c("trait.22", "trait.30", "trait.41"), 
              all_of(c("trait.22.edit", "trait.30.edit", "trait.41.edit"))) 

#Part4 geno/pheno crossover 


#reduce pheno to lines in geno and sort pheno by row names
pheno.slim <- pheno %>%
  #First reduce set down to those with genotype data
  filter(row.names(pheno) %in% row.names(geno)) %>%
  #sort by row names
  rownames_to_column(var = "rowname") %>%
  arrange(rowname) %>%
  #set the row names back
  column_to_rownames(var = "rowname")

#inspect missing characters per genotype

na.per.geno <- apply(pheno.slim, 1, function(x){sum(is.na(x))})

#plot

jpeg("outputs/Missing_DUS_characters_per_genotype.jpeg", res = 300, height = 5, width = 5, units = "in")
hist(na.per.geno, breaks=30, xlab = "Missing data points", main = "Missing DUS characters per genotype")
abline(v=5,col=2)
dev.off()

#note those with more than 5 missing

high.na.pheno <- names(na.per.geno[na.per.geno >5])

#remove with over 5 missing DUS characters

pheno.slim <- pheno.slim  %>%
  filter(!row.names(pheno.slim) %in% high.na.pheno)

#fixing a couple of errors/uncertains 
fixes<-read.csv("inputs/checked_issues.csv", header=T)
#remove uncertain genotype
pheno.slim<-pheno.slim[
  -which(row.names(pheno.slim)== fixes[1,1]),]
#recall line to spring
pheno.slim[fixes[2,1], "trait.90"] <- 2
#recall line to NA
pheno.slim[fixes[3,1], "trait.90"] <- NA
pheno.slim[fixes[4,1], "trait.90"] <- NA
#change six row to two row 
pheno.slim[fixes[5,1], "trait.39"] <- 1
pheno.slim[fixes[6,1], "trait.39"] <- 1
pheno.slim[fixes[7,1], "trait.39"] <- 1

#reduce geno to lines in pheno

geno.slim <- geno[row.names(geno) %in% row.names(pheno.slim),]

#check that the pheno.slim names match the geno.slim names match

identical(row.names(pheno.slim), row.names(geno.slim))

#Part5 - save table info on pheno in geno

#Produce table showing counts per trait
tables<-pheno.slim %>%
  map(~ table(.x, useNA = "always") %>% as.data.frame())
#combine the list of data frames into a single data frame
combined.tables <- bind_rows(tables, .id = "trait")
#pivot to wide format
combined.tables <- combined.tables %>%
  pivot_wider(names_from = .x, values_from = Freq, values_fill = list(Freq = 0))




#Part6 - save final outputs

write.csv(geno.slim,"outputs/dat.geno.with.pheno.csv")
write.csv(pheno.slim,"outputs/dat.pheno.with.geno.csv")
write.csv(combined.tables,"outputs/pheno.with.geno.trait.table.csv",row.names = F)
