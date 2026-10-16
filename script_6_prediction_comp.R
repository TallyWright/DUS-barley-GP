rm(list=ls(all=TRUE))

#Part 1 - load

#load packages
library(dplyr)
library(doParallel)
library(foreach)
library(purrr)

#load scripts
source("inputs/deduplicateTWimproved.R")

#added logging info function
log_progress <- function(trait, iteration, fold, method, logfile) {
  cat(sprintf("%s | trait=%s | iter=%s | fold=%s | method=%s | completed\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), trait, iteration, fold, method),
      file = logfile, append = TRUE)
}

#index data
geno.input<-"outputs/dat.geno.with.pheno.csv"
pheno.input<-"outputs/dat.pheno.with.geno.csv"

#weighting
major.loci.input <- "outputs/GWASouts/gwas.major.loci.csv"

#read data
geno <-read.csv(geno.input,header=T,row.names = 1)
pheno<-read.csv(pheno.input,header=T,row.names = 1)
major.loci<-read.csv(major.loci.input,header=T,row.names = 1)

# Check if it exists, and create it if not
if (!dir.exists("outputs/predict_outs")) {
  dir.create("outputs/predict_outs")
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

#recode GWAS for rrBLUP
geno <- geno %>% 
  mutate(across(everything(), ~ifelse(.==0,-1,.)))

#set number of iterations
iterations <- 20

#save trait list
traits<-names(pheno)

#traits trait.86 and trait.88 have too low phenotype diversity so are excluded
traits<-traits[!traits %in% c("trait.86", "trait.88")]

#set tuning parameters 
tuning_values <- list(ntrees=c(250, 500, 750),
                    mtry=c(50, 150, 250, 350),
                      nodesize=c(1, 3, 6)) 

all_combinations <- cross(tuning_values)

#set up parallel
cl <- makeCluster(length(traits))
registerDoParallel(cl)

# Part 3 - run prediction

foreach(i = traits, .packages = c("rrBLUP", "randomForestSRC","caret")) %dopar% {
    
    geno_filt  <- geno[!is.na(pheno[,i]),]
    pheno_filt <- pheno[!is.na(pheno[,i]),]
    id <- row.names(geno_filt)
    
    
	logfile <- paste0("outputs/predict_outs/progress_trait_", i, ".txt")
	cat(sprintf("=== Trait %s started: %s ===\n", i, Sys.time()), file = logfile, append = TRUE)
	
	#run further marker QC if lost lines for each trait
    if(!nrow(geno_filt)==nrow(geno)){
      #run some initial QC to reduce marker set to useful markers
      #ensure no low MAF in the starting dataframe (some genotypes removed without phenotype data)
      geno_filt  <- geno_filt[ , sapply(geno_filt, function(x) length(which(x == -1)) >= 12) ]
      geno_filt  <- geno_filt[ , sapply(geno_filt, function(x) length(which(x ==  1)) >= 12) ]
      
      #reduce duplicates by skimming on 100%
      geno_filt<- deduplicate(geno_filt ,1)
    }
    
    #save final data used for each trait
	 write.csv(geno_filt, 
              paste0("outputs/predict_outs/geno_data_", i, "_.csv"))
     write.csv(pheno_filt, 
              paste0("outputs/predict_outs/pheno_data_", i, "_.csv"))
	
	#form results data frame for cor
    results_cor <- as.data.frame(matrix(NA, nrow = iterations, ncol =  6)) #the number of prediction methods
    # fix names
    names(results_cor) <- c("rrBLUP","rrBLUP_int","rf","rf_tune", "rf_weight2","rf_weight10" )
    
    #replicate results table for pccc
    results_pccc <- results_cor
    
    #set results for best params for each fold in tuning step
    results_tuning <- data.frame(trait=as.character(),
                                 iteration=as.numeric(),
                                 outer_fold=as.numeric(),
                                 ntree=as.numeric(),
                                 mtry=as.numeric(), 
                                 nodesize=as.numeric(),
                                 pccc = as.numeric())
	
	#create object to save original folds 
	results_all_folds_cor <- data.frame(trait=character(), 
				iteration=numeric(),
				fold=numeric(),
				rrBLUP = numeric(),
				rrBLUP_int = numeric(),
				rf = numeric(),
				rf_tune = numeric(),
				rf_weight2 =  numeric(),
				rf_weight10 = numeric())
				
	results_all_folds_pccc <- data.frame(trait=character(), 
				iteration=numeric(),
				fold=numeric(),
				rrBLUP = numeric(),
				rrBLUP_int = numeric(),
				rf = numeric(),
				rf_tune = numeric(),
				rf_weight2 =  numeric(),
				rf_weight10 = numeric())											 
    	
	
	#r-replicated k-fold cross validation
    
	#set seed here for reproducibility (does not correlate fold iterations)
    set.seed(42) 
	
    #this part is the r-replicated
    for (j in c(1:iterations)){
      
      #r-replicated k-fold cross validation
      outer_folds<- createFolds(pheno_filt[,i], k = 10, list = TRUE, returnTrain = FALSE)
      
      #create results dataframe for average of outer_folds
      results_folds_cor<- as.data.frame(matrix(NA, nrow = 10, ncol = 6)) #number of prediction methods
      results_folds_pccc<- as.data.frame(matrix(NA, nrow = 10, ncol = 6)) #number of prediction methods
      
      #counter
      counter = 1
      
        for (k in outer_folds){
          
          #define test and train index
          test <- k
          train <- setdiff(1:nrow(pheno_filt),k)
          
          #set test
          pheno_test  <- pheno_filt[test, i]
          geno_test  <- as.matrix(geno_filt[test,])
          
          #set train 
          pheno_train <- pheno_filt[train,i]
          geno_train <- as.matrix(geno_filt[train,])
		  
		  #log progress
		  log_progress(i, j, counter, "prediction start", logfile)
      
          ######################
          #run ridge regression
          ######################
          
          mod <- mixed.solve(pheno_train,
                             Z=geno_train,
                             K=NULL,
                             SE=FALSE)
          
          #extract marker effects
          e <- mod$u
          e <- as.matrix(e)
          #G_test * e = geno_test * the marker effects
          pred_test <- geno_test %*% e
          #predicted value based on marker effects of the training population with the grand mean added in
          pred_test <- as.vector(pred_test) + as.vector(mod$beta) 
          #determine correlation accuracy
          results_folds_cor[counter, 1] <- cor(pred_test, pheno_test, use="complete")
          #save results with rounding to integer
          results_folds_cor[counter, 2]  <- cor(as.integer(round(pred_test)), pheno_test, use="complete")
          results_folds_pccc[counter, 2] <- mean(as.integer(round(pred_test))==pheno_test, na.rm=TRUE)
		  
		  #log progress
		  log_progress(i, j, counter, "rrBLUP", logfile)
          
		  ######################
          #random forest default
          ######################
          
          #set data frames for test and train
          
          rf.train <- data.frame(y = as.factor(pheno_train),geno_train)
          rf.test  <- data.frame(y = as.factor(pheno_test),geno_test)
          
          #run rf model as default
          model <- rfsrc(y~.,data= rf.train)
          
          #run predictions
          predicted <- predict(model, rf.test)$class
          
          #save cor
          results_folds_cor[counter, 3] <- cor(as.numeric(as.character(predicted)),as.numeric(as.character(rf.test$y)))
          results_folds_pccc[counter, 3] <- mean(as.character(predicted)==as.character(rf.test$y),
                                                 na.rm=TRUE)
          #log progress
		  log_progress(i, j, counter, "rf_default", logfile)
          
		  ##########################
          #random forest with tuning
          ##########################
          
          # Variable for best combo of hyperparameters and PCCC produced.
          best_params <- list(pccc=-Inf)
          
          #set inner folds
          inner_folds <- createFolds(rf.train$y, k = 5, list = TRUE, returnTrain = FALSE)
      
          #for each combination of hyperparams
          for (l in 1:length(all_combinations)) {
            
            #saves hyperparams for this run
            flags <- all_combinations[[l]]
            
            #for each inner fold
            for (m in inner_folds) {
            
            #set inner training vs. inner test just using outfold training 
            DataInnerTraining <- rf.train[setdiff(1:nrow(rf.train),m), ]
            DataInnerTesting <- rf.train[m, ] 
            
            #run model on inner fold with this combo of hyperprams
            tuning_model <- rfsrc(y ~ ., data=DataInnerTraining,
                                  ntree=flags$ntree,
                                  mtry=flags$mtry, nodesize=flags$nodesize)
            tuning_predictions <- predict(tuning_model, newdata=DataInnerTesting)$
              class
            
            #Compute PCCC for this combo of hyperparameters
            current_pccc <- mean(as.character(DataInnerTesting$y) == 
                                   as.character(tuning_predictions), na.rm=TRUE)
            
            ##overwrite best params if pccc is higher than current
            if (current_pccc > best_params$pccc) {
              best_params <- flags
              best_params$pccc <- current_pccc
            }
            }
          }
          
          #run final model for this outer fold
          model <- rfsrc(y ~ ., data= rf.train, 
                         ntree=best_params$ntree,
                         mtry=best_params$mtry, 
                         nodesize=best_params$nodesize)
          
          predicted <- predict(model, newdata=rf.test)$class
        
          #save cor for this outerfold
          results_folds_cor[counter, 4] <- cor(as.numeric(as.character(predicted)),as.numeric(as.character(rf.test$y)))  
          results_folds_pccc[counter, 4] <- mean(as.character(predicted)==
                                                   as.character(rf.test$y), na.rm=TRUE)
          #save best tuning results for summary
          results_tuning <- rbind(results_tuning, data.frame(
                                      trait=as.character(i),
                                       iteration=j,
                                       outer_fold=counter,
                                       ntree=best_params$ntree,
                                       mtry=best_params$mtry, 
                                       nodesize=best_params$nodesize,
                                       pccc = best_params$pccc))
          #log progress
		  log_progress(i, j, counter, "rf_tuning", logfile)  
          
		  ##########################
          #random forest with weighting
          ##########################
          
          #check if trait has major loci
          if(i %in% major.loci$Trait){
            
            #pull weights for each trait
            markers.to.wt <-major.loci$Marker[major.loci$Trait==i]
            
            #this will break if major loci are not in data
            if (!all(markers.to.wt %in% names(rf.train))){
              cat("major loci missing ")
              break
            }
            
            #set all weights as 1   
            wts<-rep(1,length(rf.train[,2:ncol(rf.train)]))
            
            #set major loci markers to 2
            wts[which(names(rf.train[,2:ncol(rf.train)])%in% markers.to.wt)] <- 2
            
            #Run GS with both weights
            model_1 <- rfsrc(y ~ ., 
                             data=rf.train, 
                             xvar.wt =wts,
                             split.wt=wts )
            #run predictions
            predicted <- predict(model_1, rf.test)$class
            
            #save cor
            results_folds_cor[counter, 5] <- cor(as.numeric(as.character(predicted)),as.numeric(as.character(rf.test$y)))
            results_folds_pccc[counter, 5] <- mean(as.character(predicted)==
                                                     as.character(rf.test$y), na.rm=TRUE)
            
            #log progress
		    log_progress(i, j, counter, "rf_weight2", logfile)
            
			#set all weights as 1   
            wts<-rep(1,length(rf.train[,2:ncol(rf.train)]))
            
            #set major loci markers to 10
            wts[which(names(rf.train[,2:ncol(rf.train)])%in% markers.to.wt)] <- 10
            
            #Run GS with both weights
            model_2 <- rfsrc(y ~ ., 
                             data=rf.train, 
                             xvar.wt =wts,
                             split.wt=wts )
            #run predictions
            predicted <- predict(model_2, rf.test)$class
            
            #save cor
            results_folds_cor[counter, 6] <- cor(as.numeric(as.character(predicted)),as.numeric(as.character(rf.test$y)))
            results_folds_pccc[counter, 6] <- mean(as.character(predicted)==
                                                     as.character(rf.test$y), na.rm=TRUE)
            #log progress
		    log_progress(i, j, counter, "rf_weight10", logfile)
            
          }#if for major loci
          
          counter <- counter + 1
        }
      
      #update main results with means of outfolds 
      results_cor[j , "rrBLUP"     ] <- mean(results_folds_cor$V1, na.rm=TRUE)
      results_cor[j , "rrBLUP_int" ] <- mean(results_folds_cor$V2, na.rm=TRUE)
      results_cor[j , "rf"         ] <- mean(results_folds_cor$V3, na.rm=TRUE)
      results_cor[j , "rf_tune"    ] <- mean(results_folds_cor$V4, na.rm=TRUE)
      results_cor[j , "rf_weight2"  ] <- mean(results_folds_cor$V5, na.rm=TRUE)
      results_cor[j , "rf_weight10"  ] <- mean(results_folds_cor$V6, na.rm=TRUE)
      
      results_pccc[j , "rrBLUP"     ] <- mean(results_folds_pccc$V1)
      results_pccc[j , "rrBLUP_int" ] <- mean(results_folds_pccc$V2)
      results_pccc[j , "rf"         ] <- mean(results_folds_pccc$V3)
      results_pccc[j , "rf_tune"    ] <- mean(results_folds_pccc$V4)
      results_pccc[j , "rf_weight2"  ] <- mean(results_folds_pccc$V5)
      results_pccc[j , "rf_weight10"  ] <- mean(results_folds_pccc$V6)

	  #save all fold results
	  
	  #rename before merge
      names(results_folds_cor)  <- c("rrBLUP","rrBLUP_int","rf","rf_tune","rf_weight2","rf_weight10")
      names(results_folds_pccc) <- c("rrBLUP","rrBLUP_int","rf","rf_tune","rf_weight2","rf_weight10")
  
	  results_all_folds_cor <- rbind(results_all_folds_cor,
		cbind(data.frame(trait = i, iteration = j, fold = 1:10), results_folds_cor))

	  results_all_folds_pccc <- rbind(results_all_folds_pccc,
		cbind(data.frame(trait = i, iteration = j, fold = 1:10), results_folds_pccc))

      }
    
    #save mean results for each trait
    
  write.csv(results_cor, 
              paste0("outputs/predict_outs/results_cor_", i, ".csv"),
              row.names = F)
    
  write.csv(results_pccc, 
              paste0("outputs/predict_outs/results_pccc_", i, ".csv"),
              row.names=F)
      
  write.csv(results_tuning, 
              paste0("outputs/predict_outs/results_tuning_", i, ".csv"),
              row.names=F)
			  
  write.csv(results_all_folds_cor,
                paste0("outputs/predict_outs/results_all_folds_cor_", i, ".csv"),
              row.names=F)
			  
  write.csv(results_all_folds_pccc,
                paste0("outputs/predict_outs/results_all_folds_pccc_", i, ".csv"),
              row.names=F)
  
}

stopCluster(cl)    
