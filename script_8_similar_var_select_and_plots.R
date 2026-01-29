rm(list=ls(all=TRUE))

#Part 1 - load and prepare

#load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(see)

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

# Part 4 - Summary plots for predictions

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

#add check for NA
results_NA_cor <- data.frame()

#run loop for bind

for (i in 1:length(files_to_read)){
  
  #read in data
  dat <- read.csv(paste0("outputs/predict_outs/",files_to_read[i]), header = T)
  
  #run significance test
  #iteration treated as a block
  
  #convert to long format
  
  n <- names(dat)
  dat_long<-dat
  dat_long$iteration <- as.factor(c(1:20))
  
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
  
  #run mixed model
  
  m <- lmer(accuracy ~ model + (1|iteration), data =dat_long)
  
  #save means for cor results and anova resuts
  
  results_cor <-rbind(results_cor,
                      data.frame(traits = trait_name[i],  
                                 as.data.frame(t(colMeans(dat, na.rm=T))),
                                 as.data.frame(anova(m))))
  #count na
  
  results_NA_cor <- rbind(results_NA_cor,
                          data.frame(traits =trait_name[i], 
                                     as.data.frame(t(colSums(is.na(dat))))))
  
  #form box plots
  
  #extract BLUEs/emm object
  emm <- emmeans(m, ~ model)

  #pair-wise comparisons with tukey adjustment
  res.cld<-multcomp::cld(emm, adjust = "tukey", Letters=LETTERS, sort = F)
  
  #form boxplots with letter labels overlaid
  
  p <- ggplot(res.cld, aes(x = model, y = emmean)) +
        geom_point(size = 3) +
         geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
         geom_text(aes(label = .group, y = emmean + 0.01)) +
         theme_minimal()+
        labs(x="Method",
            y = expression(italic("r")),
            title =trait_name[i] )
  
  #produce same plot without sig group labels
  
  q <- ggplot(res.cld, aes(x = model, y = emmean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
    theme_minimal()+
    labs(x="Method",
         y = expression(italic("r")),
         title =trait_name[i] )
  
  ggsave(paste0("outputs/predict_final_plots/boxplot.cor.pred.sig.groups.method",
                trait_name[i],
                 ".jpeg"),
                p, width = 5, height = 4)
  ggsave(paste0("outputs/predict_final_plots/boxplot.cor.pred.no.groups.method",
                trait_name[i],
                ".jpeg"),
         q, width = 5, height = 4)
}

#save output

write.csv(results_cor,
          "outputs/predict_final_plots/results_cor.csv", row.names = F)

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

########
#PCCCC #
########

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

#add check for NA
results_NA_pccc <- data.frame()

#run loop for bind

for (i in 1:length(files_to_read)){
  
  #read in data
  dat <- read.csv(paste0("outputs/predict_outs/",files_to_read[i]), header = T)
  
  #first col was not used, so can be removed
  
  dat <- dat[,-1]
  
  #run significance test
  #iteration treated as a block
  
  #convert to wide format
  n <- names(dat)
  dat_long<-dat
  dat_long$iteration <- as.factor(c(1:20))
  
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
  
  #run mixed model
  
  m <- lmer(accuracy ~ model + (1|iteration), data =dat_long)
  
  #save means for pccc results and anova resuts
  
  results_pccc <-rbind(results_pccc,
                      data.frame(traits = trait_name[i],  
                                 as.data.frame(t(colMeans(dat, na.rm=T))),
                                 as.data.frame(anova(m))))
  #count na
  
  results_NA_pccc <- rbind(results_NA_pccc,
                          data.frame(traits =trait_name[i], 
                                     as.data.frame(t(colSums(is.na(dat))))))
  
  #form box plots
  
  #extract BLUEs/emm object
  emm <- emmeans(m, ~ model)
  
  #pair-wise comparisons with tukey adjustment
  res.cld<-multcomp::cld(emm, adjust = "tukey", Letters=LETTERS, sort = F)
  
  #form boxplots with letter labels overlaid
  p <- ggplot(res.cld, aes(x = model, y = emmean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
    geom_text(aes(label = .group, y = emmean + 0.01)) +
    theme_minimal()+
    labs(x="Method",
         y = expression(italic("PCCC")),
         title =trait_name[i] )
  
  #form boxplots with no letter labels overlaid
  q <- ggplot(res.cld, aes(x = model, y = emmean)) +
    geom_point(size = 3)+ 
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
    theme_minimal()+
    labs(x="Method",
         y = expression(italic("PCCC")),
         title =trait_name[i] )
  
  ggsave(paste0("outputs/predict_final_plots/boxplot.PCCC.pred.sig.groups.method",
                trait_name[i],
                ".jpeg"),
         p, width = 5, height = 4)
  
  ggsave(paste0("outputs/predict_final_plots/boxplot.PCCC.pred.no.groups.method",
                trait_name[i],
                ".jpeg"),
         q, width = 5, height = 4)
}

#save output
write.csv(results_pccc,
          "outputs/predict_final_plots/results_pccc.csv", row.names = F)

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

writeLines(capture.output(sessionInfo()), "script.8.sessionInfo.txt")