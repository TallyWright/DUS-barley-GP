rm(list=ls(all=TRUE))

#Part 1 - load and prepare

#load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(correctR)
library(caret)
library(multcompView)
library(ggplot2)
library(tidyr)

#load workspace from script 7 
load("outputs/predict_candidates/script_7_data.RData")

#load inputs
class.widths.input<-"inputs/barley.class.widths.edit.csv"

#read data
class.widths<- read.csv(class.widths.input, header=T)

#fix edit to trait.90 class width 
class.widths$class.width[class.widths$trait=="trait.90"] <- 1

#set entry year
entry_year<-unique(year$ap.year)

#make output folder for saving results 

# Check if it exists, and create it if not
if (!dir.exists("outputs/predict_final_plots")) {
  dir.create("outputs/predict_final_plots")
}

# Part 2 - now compare subsets via the two approaches

#form cleaner version
clean_all_traits <- all_traits

#form empty entry to fill empty cases (where all data for candidates were missing for a trait in a year)
filler <-data.frame(Line=NA, observed = NA,predicted = NA)

#in three cases, there was no data for a certain trait within a year

#fix missing for 2010 trait.9
clean_all_traits$trait.9[[1]] <- filler
names(clean_all_traits$trait.9)[1] <- "2010"


#now do the same for trait20 2018 and 2019
clean_all_traits$trait.20[[9]] <- filler
clean_all_traits$trait.20[[10]] <- filler
names(clean_all_traits$trait.20)[c(9,10)] <-c("2018", "2019")



#create data frame for results

subset_results <- data.frame()

for (i in c(1:length(entry_year))){

    #first bind the results for each year.  
    
    bound <- bind_rows(
      lapply(names(clean_all_traits), function(x) {
        clean_all_traits[[x]][[i]] %>% #extract df of each sublist
          mutate(source = x)    #move trait name to col
      })
    )
    
    #remove NAs if present, artifact of missing trait for particular year
    bound <- bound[!is.na(bound$Line),] 
    
    #for the observations and the predictions, convert to wide from long 
    result_obs <- bound %>%         
      pivot_wider(
        id_cols = source,        
        names_from = Line,                
        values_from = observed) %>%
      as.data.frame()
    
    #convert source to rownames
    rownames(result_obs) <- result_obs$source
    result_obs$source <- NULL
    
    result_pred <- bound %>%         
      pivot_wider(
        id_cols = source,        
        names_from = Line,                
        values_from = predicted) %>%
      as.data.frame()
    
    #convert source to rownames
    rownames(result_pred) <- result_pred$source
    result_pred$source <- NULL
    
    #convert to numeric from factor
    result_obs[] <- lapply(result_obs, function(x) as.numeric(as.character(x)))
    result_pred[] <- lapply(result_pred, function(x) as.numeric(as.character(x)))
    
    #make sure traits match predicted and observed data to original pheno
    pheno.full <- pheno[,names(pheno) %in% row.names(result_obs)]
    pheno.full<-as.data.frame(t(pheno.full))
    
    # check order of traits is the same as predictions (pheno.ful and results_obs/results_preds should have same order)
    if (!identical(row.names(pheno.full), row.names(result_obs))) {
      stop("pheno.full and result_obs order is wrong")
    } 
    
    #run same for class widths
    cw.df <- class.widths[class.widths$trait %in% row.names(result_obs),]
    
    # check order of cw is the same as predictions
    if (!identical(cw.df$trait, row.names(result_obs))) {
      stop("cw vs. phenotypes do not match order")
    } 
    
    #remove candidates from pheno file names(result_obs)
    ref.col<- pheno.full[,!names(pheno.full) %in% names(result_obs)]
    
    # check number of entries are the same in result_obs and result_pred
    if (!ncol(result_obs) == ncol(result_pred)){
      stop("phenotype numbers do not match")
    }
    
    ent.no <- ncol(result_obs)
    
    #then compare
    
    for (j in c(1:ent.no)){
      
      #obs
      
      obs_subset <-apply(ref.col, 2, function(x){
        
        ifelse(is.na(result_obs[,j]) | is.na(x), NA,
             abs(result_obs[,j] - x) >= cw.df$class.width)})
      
        #pull cols with all FALSE (ignores NA)
        obs_subset  <- colnames(obs_subset)[
          apply(obs_subset, 2, function(x) {
            x_clean<- na.omit(x)
            length(x_clean) > 0 && all(x_clean == FALSE)
          })
        ]
        
      #preds
        
        pred_subset <-apply(ref.col, 2, function(x){
          
          ifelse(is.na(result_pred[,j]) | is.na(x), NA,
                 abs(result_pred[,j] - x) >= cw.df$class.width)})
        
        #pull cols with all FALSE (ignores NA)
        pred_subset  <- colnames(pred_subset)[
          apply(pred_subset, 2, function(x) {
            x_clean<- na.omit(x)
            length(x_clean) > 0 && all(x_clean == FALSE)
          })
        ] 
          
       #save
       
          subset_results <- rbind (subset_results,
                                   
                                   data.frame( testing.year = entry_year[i],  
                                               entry = names(result_obs)[j],
                                               obs_subset_size = length(obs_subset),
                                               pred_subset_size = length(pred_subset),
                                               matches = length(intersect(obs_subset, pred_subset)),
                                               prop_of_obs_found = mean(obs_subset %in% pred_subset) * 100,
                                               obs_subset_colapse = paste(obs_subset, collapse = ";"),
                                               pred_subset_colapse = paste(pred_subset, collapse = ";")
                                               ) )
    }

}

# Part 3 - Summary plots for subsets

#convert to long format for plotting
res_table_long <- pivot_longer(subset_results, cols = c("obs_subset_size", "pred_subset_size"), names_to = "indepen", values_to = "depen")

plot_subset_compare<-ggplot(res_table_long, aes(x = indepen, y = depen)) +
  geom_boxplot(outlier.shape = NA, fill = "aliceblue", color = "black") +  
  geom_jitter(width = 0.1,color = "darkblue",alpha = 0.5) +   
  geom_line(aes(group = entry),color = "coral",alpha = 0.5) +  
  labs(x = "Phenotypes", y = "Genotypes required for field testing") +
  scale_x_discrete(labels = c("obs_subset_size" = "Observed", "pred_subset_size" = "Predicted"))+
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = paste0("Mean: ",round(..y..))),
    vjust = -22,
    color = "black"
  ) +
  theme_minimal(base_size = 14)+
  coord_cartesian(ylim = c(0, 125))

#reduce observations with one or less in the obs subset
reduced_subset_results<- subset_results[subset_results$obs_subset_size > 1,]

plot_prop_subsets<- ggplot(reduced_subset_results, aes(x = prop_of_obs_found)) +
  geom_histogram(binwidth = 5, fill = "darkblue", color = "aliceblue", alpha = 0.5) +
  labs(x = "Percentage of observed found in predicted (%)", y = "Frequency") +
  theme_minimal(base_size = 14) +
  annotate("text", x = 60, y = 17, 
           label = paste("Mean = ", round(mean(reduced_subset_results$prop_of_obs_found, na.rm=T)), "%",sep=""), 
           colour = "black")+
  annotate("text", x = 60, y = 16, 
           label = paste("Median = ", round(median(reduced_subset_results$prop_of_obs_found, na.rm=T)), "%",sep=""), 
           colour = "black")
#save
ggsave("outputs/predict_final_plots/fig_subsets.jpeg", ggarrange(plot_subset_compare,plot_prop_subsets, ncol=2,labels = c("a)", "b)") ), height = 5, width = 10)

#Run updated tests for significance test between using prediction accuracies

# Part 4 - sig test for cor

#fix NA issue with trait.25
dat <-read.csv("outputs/predict_outs/results_all_folds_cor_trait.25.csv", header=T)

write.csv(dat,"outputs/predict_outs/orginal_before_impute_trait.25.csv")

#estimate % missing values = rrBLUP_int=5, rf/rf_tune/rf_weight2 = 1.  
colSums(is.na(dat))

#fill impacted cells with mean
dat$rrBLUP_int[is.na(dat$rrBLUP_int)]<-mean(dat$rrBLUP_int,na.rm=T)
dat$rf[is.na(dat$rf)]<-mean(dat$rf,na.rm=T)
dat$rf_tune[is.na(dat$rf_tune)]<-mean(dat$rf_tune,na.rm=T)
dat$rf_weight2[is.na(dat$rf_weight2)]<-mean(dat$rf_weight2,na.rm=T)

#now save imputed over original
write.csv(dat,"outputs/predict_outs/results_all_folds_cor_trait.25.csv", row.names = F)

#clean enviro
rm(list=ls(all=TRUE))
      

#read in file list and only keep corr
files_to_read <- list.files("outputs/predict_outs")
files_to_read <- files_to_read [grepl("^results_all_folds_cor_trait\\.", files_to_read )]

#save just traits
traits<-gsub("results_all_folds_cor_|\\.csv", "", files_to_read)

#create letter results table to save to

results_cor <- data.frame(method=c("rrBLUP", 
                                   "rrBLUP_int", 
                                   "rf", 
                                   "rf_tune", 
                                   "rf_weight2", 
                                   "rf_weight10"))

#run loop for plots and table

for (i in 1:length(traits)){
  
  #read in data
  dat <- read.csv(paste0("outputs/predict_outs/results_all_folds_cor_",traits[i],".csv"), header = T)
  dat_pheno <- read.csv(paste0("outputs/predict_outs/pheno_data_",traits[i],"_.csv"), 
                        header = T,row.names = 1)
  
  #form test set examples
  outer_folds<- createFolds(dat_pheno[,traits[i]], k = 10, list = TRUE, returnTrain = FALSE)
  
  #average size of test
  test.set.size<-round(sum(sapply(outer_folds, length))/10)
  
  #form training set examples and average size
  train.set.size<- round(sum(sapply(outer_folds, function(x){
    length(setdiff(1:length(dat_pheno[,traits[i]]),x))}))/10)
  
  #remove NA methods
  dat <- dat[,colSums(is.na(dat))<nrow(dat)]
  
  #form pairs
  pairs <- combn(names(dat)[-c(1:3)], 2, simplify = FALSE)
  
  #run p-values on each test  
  test.results<- lapply(pairs, function(z){
    
    temp <- dat[, c("iteration", "fold", z)] %>%
      pivot_longer(cols = all_of(z), names_to = "model", values_to = "values")
    
    names(temp)[names(temp) == "iteration"] <- "r"
    names(temp)[names(temp) == "fold"]      <- "k"
    
    temp <- as.data.frame(temp)
    
    repkfold_ttest(temp, 
                   n1 =train.set.size , 
                   n2 = test.set.size, 
                   k=10, r=20) 
  })
  
  #unlist p-values
  p_values<-unlist(sapply(test.results,function(x){x[2]}))
  
  #adjust for multiple corrections
  p_values.A<-p.adjust(p_values, method = "holm")
  
  #set names
  names(p_values.A)<-sapply(pairs, function(x){paste(x[1], x[2], 
                                                     sep="-")})  
  #save letters
  letters <-as.data.frame(multcompLetters(p_values.A)[1])
  
  #convert to wide format
  long_data <- pivot_longer(
    dat,
    cols = names(dat)[-c(1:3)],
    names_to = "method",
    values_to = "values"
  )
  
  #change factor order
  long_data$method <- factor(long_data$method, levels = c("rrBLUP", "rrBLUP_int", "rf", "rf_tune", "rf_weight2", "rf_weight10"))
  
  #drop empty levels 
  long_data$method<-droplevels(long_data$method)
  
  #aggregate step 1
  plot_data <- aggregate(values~method+iteration, data=long_data, mean)
  
  #aggregate step 2
  plot_data <- aggregate(values~method, data=plot_data, 
                         FUN = function(x){c(mean=mean(x),
                                             se=sd(x)/sqrt(length(x)))})
  
  plot_data <- do.call(data.frame, plot_data)
  names(plot_data) <- c("method", "mean", "se")
  
  #run merge for data and letters
  plot_data<-merge(plot_data,letters, by.x="method", by.y=0)
  
  #form table to save as main results
  temp<-setNames(data.frame(plot_data$method, 
                            paste0(sprintf("%.3f", plot_data$mean),plot_data$Letters)), 
                 c("method", traits[i]))
  
  #save to main results
  results_cor <- merge(results_cor,temp, by="method", all = T)
  
  #form boxplots with letter labels overlaid
  
  label_offset<- diff(range(plot_data$mean)) * 0.05 
  
  plot_point  <- ggplot(plot_data, aes(x = method, y = mean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_text(aes(label = Letters, y = mean + label_offset + se)) +
    theme_minimal()+
    labs(x="Method",
         y = expression(italic("r")),
         title =sub("trait\\.", "Characteristic no. ", traits[i]) )
  
  ggsave(paste0("outputs/predict_final_plots/cor_point_plot_with_letters_",
                traits[i],".jpeg"), plot_point, 
         height = 4, 
         width = 5)
  
  
}
#save results

#save correlation letters 

#sort on of methods 
results_cor <- results_cor[match(c("rrBLUP", 
                                   "rrBLUP_int", 
                                   "rf", 
                                   "rf_tune", 
                                   "rf_weight2", 
                                   "rf_weight10"), results_cor$method), ]
#move methods to row names
row.names(results_cor) <- results_cor$method
results_cor <-results_cor[,-1]

#sort on order of traits
trait_num <- as.numeric(gsub("trait\\.", "", names(results_cor)))
results_cor<-results_cor[,names(results_cor)[order(trait_num)]]

#save final data
write.csv(results_cor,"outputs/predict_final_plots/cor_letters_table.csv")

#############
# Part 5 - scatter for cor
#############

#clean enviro
rm(list=ls(all=TRUE))

#read in file list and only keep corr
files_to_read <- list.files("outputs/predict_outs")
files_to_read <- files_to_read [grepl("^results_cor_trait\\.", files_to_read )]

#read in h2 results for plotting 
h2<-read.csv("outputs/table.h2.results.csv", header=T)

#save trait name vector
trait_name<-sub("^results_cor_(trait\\.[0-9]+)\\.csv$", "\\1", files_to_read)

#save results object
results_cor <- data.frame()

#run loop for bind

for (i in 1:length(files_to_read)){
  
  #read in data
  dat <- read.csv(paste0("outputs/predict_outs/",files_to_read[i]), header = T)
  
  #convert to long format
  n <- names(dat)
  dat_long<-dat
  
  dat_long <- dat_long %>%         
    pivot_longer(
      cols = all_of(n),        
      names_to = "model",                
      values_to = "accuracy") %>%
    as.data.frame()
  
  #fix factor order
  dat_long$model <-factor(dat_long$model, levels= c("rrBLUP",
                                                    "rrBLUP_int",
                                                    "rf",
                                                    "rf_tune",
                                                    "rf_weight2",
                                                    "rf_weight10"))
  #clean 
  dat_long <- dat_long[ !is.na(dat_long$accuracy),]
  dat_long <- droplevels(dat_long)
  
  #save means for cor results and anova resuts
  
  results_cor <-rbind(results_cor,
                      data.frame(traits = trait_name[i],  
                                 as.data.frame(t(colMeans(dat, na.rm=T)))))
}
#merge and plot for cor plot
results_cor_h2<-merge(results_cor, h2, by="traits")

dat_long_corr <- tidyr::pivot_longer(results_cor_h2, cols = c(rrBLUP,rf), names_to = "Variable", values_to = "Value")

scat_corr<-ggplot(dat_long_corr, aes(x = h2, y = Value, color = Variable)) +
  geom_point() +
  coord_cartesian(ylim = c(0.2, 1))+
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = expression(italic("h"^2)),
    y = expression(italic("r"))) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.1))+
  guides(color = guide_legend(title = NULL))

ggsave("outputs/predict_final_plots/scatter_cor.jpeg", scat_corr, width = 6, height = 5)

write.csv(results_cor,
          "outputs/predict_final_plots/results_cor.csv", row.names = F)


# Part 6 - sig test for PCCC


#clean enviro
rm(list=ls(all=TRUE))


#read in file list and only keep pcccr
files_to_read <- list.files("outputs/predict_outs")
files_to_read <- files_to_read [grepl("^results_all_folds_pccc_trait\\.", files_to_read )]

#save just traits
traits<-gsub("results_all_folds_pccc_|\\.csv", "", files_to_read)

#create letter results table to save to

results_pccc <- data.frame(method=c("rrBLUP", 
                                    "rrBLUP_int", 
                                    "rf", 
                                    "rf_tune", 
                                    "rf_weight2", 
                                    "rf_weight10"))

#run loop for plots and table

for (i in 1:length(traits)){
  
  #read in data
  dat <- read.csv(paste0("outputs/predict_outs/results_all_folds_pccc_",traits[i],".csv"), header = T)
  dat_pheno <- read.csv(paste0("outputs/predict_outs/pheno_data_",traits[i],"_.csv"), 
                        header = T,row.names = 1)
  
  #form test set examples
  outer_folds<- createFolds(dat_pheno[,traits[i]], k = 10, list = TRUE, returnTrain = FALSE)
  
  #average size of test
  test.set.size<-round(sum(sapply(outer_folds, length))/10)
  
  #form training set examples and average size
  train.set.size<- round(sum(sapply(outer_folds, function(x){
    length(setdiff(1:length(dat_pheno[,traits[i]]),x))}))/10)
  
  #remove NA methods
  dat <- dat[,colSums(is.na(dat))<nrow(dat)]
  
  #form pairs
  pairs <- combn(names(dat)[-c(1:3)], 2, simplify = FALSE)
  
  #run p-values on each test  
  test.results<- lapply(pairs, function(z){
    
    temp <- dat[, c("iteration", "fold", z)] %>%
      pivot_longer(cols = all_of(z), names_to = "model", values_to = "values")
    
    names(temp)[names(temp) == "iteration"] <- "r"
    names(temp)[names(temp) == "fold"]      <- "k"
    
    temp <- as.data.frame(temp)
    
    repkfold_ttest(temp, 
                   n1 =train.set.size , 
                   n2 = test.set.size, 
                   k=10, r=20) 
  })
  
  #unlist p-values
  p_values<-unlist(sapply(test.results,function(x){x[2]}))
  
  #adjust for multiple corections
  p_values.A<-p.adjust(p_values, method = "holm")
  
  #set names
  names(p_values.A)<-sapply(pairs, function(x){paste(x[1], x[2], 
                                                     sep="-")})  
  #save letters
  letters <-as.data.frame(multcompLetters(p_values.A)[1])
  
  #convert to wide format
  long_data <- pivot_longer(
    dat,
    cols = names(dat)[-c(1:3)],
    names_to = "method",
    values_to = "values"
  )
  
  #change factor order
  long_data$method <- factor(long_data$method, levels = c("rrBLUP", "rrBLUP_int", "rf", "rf_tune", "rf_weight2", "rf_weight10"))
  
  #drop empty levels 
  long_data$method<-droplevels(long_data$method)
  
  #aggregate step 1
  plot_data <- aggregate(values~method+iteration, data=long_data, mean)
  
  #aggregate step 2
  plot_data <- aggregate(values~method, data=plot_data, 
                         FUN = function(x){c(mean=mean(x),
                                             se=sd(x)/sqrt(length(x)))})
  
  plot_data <- do.call(data.frame, plot_data)
  names(plot_data) <- c("method", "mean", "se")
  
  #run merge for data and letters
  plot_data<-merge(plot_data,letters, by.x="method", by.y=0)
  
  #form table to save as main results
  temp<-setNames(data.frame(plot_data$method, 
                            paste0(sprintf("%.3f", plot_data$mean),plot_data$Letters)), 
                 c("method", traits[i]))
  
  #save to main results
  results_pccc <- merge(results_pccc,temp, by="method", all = T)
  
  #form boxplots with letter labels overlaid
  
  if(traits[i]=="trait.90"){ 
    label_offset<- diff(range(plot_data$mean)) * 0.15 
  }else{
    label_offset<- diff(range(plot_data$mean)) * 0.05 }
  
  plot_point  <- ggplot(plot_data, aes(x = method, y = mean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_text(aes(label = Letters, y = mean + label_offset + se)) +
    theme_minimal()+
    labs(x="Method",
         y = "PCCC",
         title =sub("trait\\.", "Characteristic no. ", traits[i]) )
  
  ggsave(paste0("outputs/predict_final_plots/pccc_point_plot_with_letters_",
                traits[i],".jpeg"), plot_point, 
         height = 4, 
         width = 5)
  
  
}
#save results

#save PCCC letters 

#sort on of methods 
results_pccc <- results_pccc[match(c("rrBLUP", 
                                     "rrBLUP_int", 
                                     "rf", 
                                     "rf_tune", 
                                     "rf_weight2", 
                                     "rf_weight10"), results_pccc$method), ]
#move methods to row names
row.names(results_pccc) <- results_pccc$method
results_pccc <-results_pccc[,-1]

#sort on order of traits
trait_num <- as.numeric(gsub("trait\\.", "", names(results_pccc)))
results_pccc<-results_pccc[,names(results_pccc)[order(trait_num)]]

#save final data
write.csv(results_pccc,"outputs/predict_final_plots/pccc_letters_table.csv")

# Part 7 - scatter for PCCC

#clean enviro
rm(list=ls(all=TRUE))

#read in file list and only keep pccc
files_to_read <- list.files("outputs/predict_outs")
files_to_read <- files_to_read [grepl("^results_pccc_trait\\.", files_to_read )]

#read in h2 results for plotting 
h2<-read.csv("outputs/table.h2.results.csv", header=T)

#save trait name vector
trait_name<-sub("^results_pccc_(trait\\.[0-9]+)\\.csv$", "\\1", files_to_read)

#save results object
results_pccc <- data.frame()

#run loop for bind

for (i in 1:length(files_to_read)){
  
  #read in data
  dat <- read.csv(paste0("outputs/predict_outs/",files_to_read[i]), header = T)
  
  #first col was not used, so can be removed
  dat <- dat[,-1]
  
  #convert to wide format
  n <- names(dat)
  dat_long<-dat
  
  dat_long <- dat_long %>%         
    pivot_longer(
      cols = all_of(n),        
      names_to = "model",                
      values_to = "accuracy") %>%
    as.data.frame()
  
  #fix factor order
  dat_long$model <-factor(dat_long$model, levels= c("rrBLUP_int",
                                                    "rf",
                                                    "rf_tune",
                                                    "rf_weight2",
                                                    "rf_weight10"))
  #clean 
  dat_long <- dat_long[ !is.na(dat_long$accuracy),]
  dat_long <- droplevels(dat_long)
  
  results_pccc <-rbind(results_pccc,
                       data.frame(traits = trait_name[i],  
                                  as.data.frame(t(colMeans(dat, na.rm=T)))))
}
#merge and plot for pccc plot
results_pccc_h2<-merge(results_pccc, h2, by="traits")

dat_long_pcccr_h2 <- tidyr::pivot_longer(results_pccc_h2, cols = c(rrBLUP_int,rf), names_to = "Variable", values_to = "Value")

scatter_pccc <- ggplot(dat_long_pcccr_h2, aes(x = h2, y = Value, color = Variable)) +
  geom_point() +
  coord_cartesian(ylim = c(0.2, 1))+
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = expression(italic("h"^2)),
    y = "PCCC") +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.1))+
  guides(color = guide_legend(title = NULL))

ggsave("outputs/predict_final_plots/scatter_pccc.jpeg", scatter_pccc, width = 6, height = 5)

#save output
write.csv(results_pccc,
          "outputs/predict_final_plots/results_pccc.csv", row.names = F)


writeLines(capture.output(sessionInfo()), "script.8.sessionInfo.txt")