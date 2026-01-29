###Part 1 workspace load and inputs
rm(list = ls())

#load packages

library(dplyr)
library(tibble)

#read inputs

combined.data.loc <- "inputs/geno_input.txt"
name.convert.loc  <- "inputs/name_convert_and_duplicates.csv" 
marker.converts.loc   <-  "inputs/marker_conversions.csv"
allele.miss.match.loc <- "inputs/allele_missmatches_in_convert.csv"
genotypes.to.remove.loc <- "inputs/genotypes_removed_from_combined.csv" 

#read in inputs

dat<-read.table(combined.data.loc,header = T,sep = "\t", colClasses = "character", check.names=FALSE)
name.convert <- read.csv(name.convert.loc,header=T, check.names = FALSE)
marker.converts <- read.csv(marker.converts.loc, header=T, row.names = 1)
allele.miss.match <- read.csv(allele.miss.match.loc, header=T, row.names = 1)
genotypes.to.remove <- read.csv(genotypes.to.remove.loc, header=T)

#form output folder

dir.create("outputs", showWarnings = F)

###Part 2 data prep

#merge and convert names in data

dat <- name.convert %>%
  left_join(dat, by="Line/Marker") %>%
  column_to_rownames(var = "AFP") %>%
  {assign("dat.name.info", select(., 1:3), envir = .GlobalEnv); .} %>%
  select(-(1:3)) 
rm(name.convert)

#perform initial minor data edits

dat <- dat %>%
  mutate(across(everything(), ~na_if(., "--"))) %>%  #Replace "--" string with NA
  select(where(~!any(. == "N/A", na.rm = TRUE))) %>%  # Remove columns containing "N/A" (only 5 cases)
  select(-`JHI-Hv50k-2016-513784`) #marker had incorrect allele het calls
  
#save and remove monomorphic markers (need this for allele convert)

rm_monomorphic_m <- dat %>%
  select(where(~n_distinct(na.omit(.)) == 1 | all(is.na(.))))

dat <- dat %>%
  select(where(~n_distinct(na.omit(.)) > 1 & any(!is.na(.))))

###Part 3 marker allele conversion

#create function for recalling

replace_values <- function(column, allele_a, allele_b){
  column <- ifelse(column == allele_a, 0, column)
  column <- ifelse(column == allele_b, 2, column)
  return(column)
}

#run recalling script

dat.num <- dat %>%
  select(-all_of(row.names(allele.miss.match)))%>% # remove several markers where calls and conversion made little sense (30 markers)
  mutate(across(everything(), ~ifelse(. == "AA", "A", .)))%>%
  mutate(across(everything(), ~ifelse(. == "CC", "C", .)))%>%
  mutate(across(everything(), ~ifelse(. == "GG", "G", .)))%>%
  mutate(across(everything(), ~ifelse(. == "TT", "T", .)))%>%
  mutate(across(everything(), ~ifelse(. %in% c("AC", "AG", "AT", "CG", "CT", "GT"), 1, .)))%>%
  mutate(across(everything(), ~replace_values(.,
                                               allele_a = marker.converts[match(cur_column(), rownames(marker.converts)), "Allele.A"],
                                               allele_b = marker.converts[match(cur_column(), rownames(marker.converts)), "Allele.B"])))%>%
  mutate(across(everything(), as.numeric))

###Part 4 data quality control

#Reduce genotypes based on previous genotype investigation

dat.num <- dat.num %>%
  rownames_to_column(var="row_names") %>% 
  filter(!row_names %in% genotypes.to.remove$Line) %>% #should leave 1171
  column_to_rownames(var="row_names")

#plot missing and hets %

jpeg("outputs/marker_missing_call_rates.jpeg", units="in", height=6, width=12, res=300)
par(mfrow=(c(1,2)))

na_count <-apply(dat.num, 2, function(x) {sum(is.na(x))} )
na_perc <- (na_count/nrow(dat.num))*100
hist(na_perc, breaks = 50,main="Percentage missing calls in markers", xlab="%")
hist(na_perc, breaks = 100,ylim=c(0,1000), xlim = c(0,20),main="Limited axis graph", xlab="%")
abline(v=3,col=2)
dev.off()

rm(na_count,na_perc)

jpeg("outputs/marker_heterozgosity_rates.jpeg", units="in", height=6, width=12, res=300)
par(mfrow=(c(1,2)))

hets_count <- apply(dat.num, 2, function(x) {length(which(x==1))} )
hets_perc <- (hets_count/nrow(dat.num))*100
hist(hets_perc, main = "Percentage of heterozgosity in markers", xlab = "%",breaks=50)
hist(hets_perc, main = "Limited axis graph", xlab = "%",breaks=100,xlim = c(0,20),ylim=c(0,1000))
abline(v=3,col=2)
dev.off()

rm(hets_count,hets_perc)

#create function for cleaning 

marker.clean<-function(data, na.thresh = 10,
                       het.thresh = 10, rare.hom.thresh = 15){
  
  # Check if na.thresh are numeric and valid
  if (!is.numeric(na.thresh) || na.thresh < 0 || na.thresh> 100) {
    stop("The na.thresh has to be numeric and between 0 and 100.")
  }
  
  if (!is.numeric(het.thresh) || het.thresh< 0 || het.thresh > 100) {
    stop("The het.thresh has to be numeric and between 0 and 100.")
  }
  
  if (!is.numeric(rare.hom.thresh) || rare.hom.thresh< 0 || rare.hom.thresh > 50) {
    stop("The rare.hom.thresh has to be numeric and between 0 and 50.")
  }
  
  #print NA removal number 
  excluded_columns_na <- data %>%
    select(where(~ mean(is.na(.)) * 100 > na.thresh))%>%
    ncol()
  print(paste("Number of markers excluded with higher than", na.thresh,
              "% NA:", excluded_columns_na))
  
  #clean on NA
  data <- data %>%
    select(where(~ mean(is.na(.)) * 100 <= na.thresh ))
  
  #print hets removal number 
  excluded_columns_hets <- data %>%
    select(where(~ mean(. == 1, na.rm = TRUE)   * 100 > het.thresh ))%>%
    ncol()
  print(paste("Number of markers excluded with higher than", het.thresh,
              "% hets:", excluded_columns_hets))
  
  #clean on hets
  data <- data %>%
    select(where(~ mean(. == 1, na.rm = TRUE)   * 100 <=   het.thresh )) #Remove columns with het percentage above the threshold, note this will estimate % hets with excluding NA. So % hets of present data. 
  
  #print 0 count removal number 
  excluded_columns_0_count <- data %>%
    select(where(~ sum(. == 0, na.rm = TRUE) < rare.hom.thresh  ))%>% 
  ncol()
  print(paste("Number of markers excluded with fewer than", rare.hom.thresh,
              "'0' hom calls:", excluded_columns_0_count))
  
  #clean on 0 count
  data <- data %>%
    select(where(~ sum(. == 0, na.rm = TRUE)  >=  rare.hom.thresh  ))
  
  #print 2 count removal number 
  excluded_columns_2_count <- data %>%
    select(where(~ sum(. == 2, na.rm = TRUE) < rare.hom.thresh  ))%>%
    ncol()
  print(paste("Number of markers excluded with fewer than", rare.hom.thresh,
              "'2' hom calls:", excluded_columns_2_count))
  
  #clean on 2 count
  data <- data %>%
    select(where(~ sum(. == 2, na.rm = TRUE)  >=  rare.hom.thresh  ))
  
  return(data)
  
}

#Form cleaned data, remove markers with high hets, high NA, monomorphic and low MAF

print("number before clean")
dim(dat.num)

dat.clean<-marker.clean(dat.num, na.thresh = 0, het.thresh = 3, rare.hom.thresh = 12)

print("number after clean")
dim(dat.clean)

write.csv(dat.clean,"outputs/dat.clean.1.csv")
