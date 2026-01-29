rm(list=ls(all=TRUE))
#Part 1 - load

#load packages

library(ape)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(sommer)

#load
source("inputs/deduplicateTWimproved.R")

#set inputs

geno.input<-"outputs/dat.geno.with.pheno.csv"
pheno.input<-"outputs/dat.pheno.with.geno.csv"


#read data

geno<-read.csv(geno.input,header=T,row.names = 1)
pheno<-read.csv(pheno.input,header=T,row.names = 1)

#will make general reduction in duplicate markers and low allele counts
geno  <- geno[ , sapply(geno, function(x) length(which(x == 0)) >= 12) ]
geno  <- geno[ , sapply(geno, function(x) length(which(x == 1)) >= 12) ]

#reduce duplicates by skimming on 100%
geno <- deduplicate(geno ,1)

#Part 2 - distance matrices

#geno euclidean with scaling
geno.scaled <- scale(geno)
geno.eucl.dist<-dist(geno.scaled, method = "euclidean")

#geno manhattan
geno.man.dist<-dist(geno,method = "manhattan")
geno.man.dist<-geno.man.dist/ncol(geno)

#pheno euclidean with scaling
pheno.scaled <- scale(pheno)
pheno.eucl.dist<-dist(pheno.scaled, method = "euclidean")

#pheno manhattan with min-max normalisation

pheno.norm<-as.data.frame(lapply(pheno, function(x){
  mn <- min(x, na.rm = T)
  mx <- max(x, na.rm = T)
  (x - mn)/(mx - mn)
}))

rownames(pheno.norm) <- rownames(pheno)

pheno.man.dist.norm<-dist(pheno.norm, method = "manhattan")
pheno.man.dist.norm<-pheno.man.dist.norm/ncol(pheno.norm)

#without normalisation
pheno.man.dist<-dist(pheno, method = "manhattan")
pheno.man.dist<-pheno.man.dist/ncol(pheno)

#correlate matrices
man_man_cor<-mantel(geno.man.dist, pheno.man.dist, 
                    method="pearson")

man_mannorm_cor<-mantel(geno.man.dist, pheno.man.dist.norm, 
                        method="pearson")


#save mantel results for annotation
mantel.m.r.1<-round(man_man_cor$statistic,2)
mantel.m.p.1 <- man_man_cor$signif
mantel.txt.1<-as.expression(bquote("Mantel test:" ~ italic(r) == .(mantel.m.r.1) * ","
                   ~ italic(P) == .(mantel.m.p.1)))

mantel.m.r.2<-round(man_mannorm_cor$statistic,2)
mantel.m.p.2 <- man_mannorm_cor$signif
mantel.txt.2<-as.expression(bquote("Mantel test:" ~ italic(r) == .(mantel.m.r.2) * ","
                                   ~ italic(P) == .(mantel.m.p.2)))

#form scatter plot between matrices
man.dist.vecs.1<-data.frame(genetic=c(geno.man.dist), phenotype=c(pheno.man.dist))
man.dist.vecs.2<-data.frame(genetic=c(geno.man.dist), phenotype=c(pheno.man.dist.norm))

man.scatter.plot.1<-ggscatter(man.dist.vecs.1, x = "genetic", y = "phenotype", shape = 20,alpha = 0.2,size = 0.3,
          cor.coef = F,
          xlab="Genotype Manhattan Distance", 
          ylab = "Phenotype Manhattan Distance")+
  annotate("text", x = 0.12, y = 2.8, label = mantel.txt.1,  size = 4, parse=T)+
  stat_density_2d(#for smoothing
    aes(fill = after_stat(level)),
    geom = "polygon"
  )+
  scale_fill_viridis_c(name = "Density")+
  theme(legend.position = c(0.95, 0.15))

man.scatter.plot.2<-ggscatter(man.dist.vecs.2, x = "genetic", y = "phenotype", shape = 20,alpha = 0.2,size = 0.3,
                              cor.coef = F,
                              xlab="Genotype Manhattan Distance", 
                              ylab = "Min-Max Scaled Phenotype Manhattan Distance")+
  annotate("text", x = 0.12, y = 0.5, label = mantel.txt.2,  size = 4, parse=T) +
  stat_density_2d(#for smoothing
    aes(fill = after_stat(level)),
    geom = "polygon"
  )+
  scale_fill_viridis_c(name = "Density")+
  theme(legend.position = c(0.95, 0.15))


#Part3 - PCoA

#Geno plot

#calculate PCoA for geno
pcoa.geno.eucl <- pcoa(geno.eucl.dist, correction="none", rn=NULL)

#extract first two PCoA gor geno
pcoa.geno.vec <- as.data.frame(pcoa.geno.eucl$vectors[,1:2])

#merge PCoA with pheno for plotting for geno
pcoa.geno.vec <- merge(pcoa.geno.vec, pheno, by="row.names")

pcoa.geno.plot<-ggplot(pcoa.geno.vec, aes(x = Axis.1, y = Axis.2, color = as.factor(trait.90), shape = as.factor(trait.39))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("1"="darkcyan", "2"="purple4"),#Customize color scale
                     labels = c("1" = "Winter type", "2" = "Spring type"),
                     guide = guide_legend(order = 1)) +  
  scale_shape_manual(values = c(16, 17),#Customize shape scale
                     labels = c("1" = "Two-row", "2" = "Six-row"),
                     guide = guide_legend(order = 2)) +  
  labs(
       x = paste("PCoA_1 (",round(100*pcoa.geno.eucl$values$Relative_eig[1]),"%)", sep="") ,
       y = paste("PCoA_2 (",round(100*pcoa.geno.eucl$values$Relative_eig[2]),"%)", sep=""),
       color = "Seasonal type",
       shape = "No. of ear rows")+
  theme_minimal()+
  theme(
    legend.title = element_text(margin = margin(r = 10))
  )

#Pheno plot

#calculate PCoA for pheno
pcoa.pheno.eucl <- pcoa(pheno.eucl.dist, correction="none", rn=NULL)

#extract first two PCoA gor pheno
pcoa.pheno.vec <- as.data.frame(pcoa.pheno.eucl$vectors[,1:2])

#merge PCoA with pheno for plotting for pheno
pcoa.pheno.vec <- merge(pcoa.pheno.vec, pheno, by="row.names")

pcoa.pheno.plot<-ggplot(pcoa.pheno.vec, aes(x = Axis.1, y = Axis.2, color = as.factor(trait.90), shape = as.factor(trait.39))) +
  geom_point(size = 2) +
  scale_color_manual(values = c("1"="darkcyan", "2"="purple4"),#Customize color scale
                     labels = c("1" = "Winter type", "2" = "Spring type"),
                     guide = guide_legend(order = 1)) +  
  scale_shape_manual(values = c(16, 17),#Customize shape scale
                     labels = c("1" = "Two-row", "2" = "Six-row"),
                     guide = guide_legend(order = 2)) + 
  labs(
    x = paste("PCoA_1 (",round(100*pcoa.pheno.eucl$values$Relative_eig[1]),"%)", sep="") ,
    y = paste("PCoA_2 (",round(100*pcoa.pheno.eucl$values$Relative_eig[2]),"%)", sep=""),
    color = "Seasonal type",
    shape = "No. of ear rows")+
  theme_minimal()+
  theme(
    legend.title = element_text(margin = margin(r = 10))
  )

#form merged plot
pcoa.merged.plot.1<-ggarrange(pcoa.geno.plot, pcoa.pheno.plot,
          ncol = 2, nrow =1,
          common.legend = TRUE, legend = "right", labels = c("a)","b)"),  label.x = -0.01)

combined.scatter.plot.2 <- ggarrange(man.scatter.plot.1,man.scatter.plot.2, ncol = 2, nrow =1,labels = c("c)","d)"))
combined.pcoa.scatter<-ggarrange(pcoa.merged.plot.1,combined.scatter.plot.2, ncol = 1, nrow = 2) 

#save plot
ggsave("outputs/combined.pcoa.scatter.jpeg",combined.pcoa.scatter,width = 15,height = 13)

#Part4 - H2 Calc

#recode =  coded as {-1,0,1}.

geno.rr <- geno %>% 
  mutate(across(everything(), ~ifelse(.==0,-1,.)))

#form traits vector
traits<-names(pheno)

#form table to save h2
h2.results<-data.frame(traits = traits, h2 = NA, se=NA)

#should be true
identical(row.names(pheno), row.names(geno.rr))

#estimate h2 per trait

for (i in traits){
  
geno_filt  <- geno.rr[!is.na(pheno[,i]),]
pheno_filt <- pheno[!is.na(pheno[,i]),]
id <- row.names(geno_filt)

#run further marker QC if lost lines for each trait
if(!nrow(geno_filt)==nrow(geno.rr)){
  #run some initial QC to reduce marker set to useful markers
  #ensure no low MAF in the starting dataframe (some genotypes removed without phenotype data)
  geno_filt  <- geno_filt[ , sapply(geno_filt, function(x) length(which(x == -1)) >= 12) ]
  geno_filt  <- geno_filt[ , sapply(geno_filt, function(x) length(which(x ==  1)) >= 12) ]
  
  #reduce duplicates by skimming on 100%
  geno_filt<- deduplicate(geno_filt ,1)
}

# check order of genotypes and phenotypes is the same
if (!identical(row.names(geno_filt), row.names(pheno_filt))) {
  stop("genotypes and phenotypes do not match")
}  

#convert geno to matrix
geno_filt_mat<-as.matrix(geno_filt)

#Calculate realized additive relationship matrix
grm <-sommer::A.mat(geno_filt_mat)

#form formula for model
form <- as.formula(paste(i, "~ 1"))

#run model, no dominance to estimate (inbreds)
ans.A <- mmes(form, 
random=~vsm(ism(id),Gu=grm), 
rcov=~units, 
data=pheno_filt,verbose = F)

#save h2
h2.results[which(h2.results$traits==i),2]<-
  vpredict(ans.A, h2 ~ (V1) / ( V1+V2) )[1] 

#save se
h2.results[which(h2.results$traits==i),3]<-
  vpredict(ans.A, h2 ~ (V1) / ( V1+V2) )[2] 

}

#save h2 results 

write.csv(h2.results, "outputs/table.h2.results.csv", row.names = F)
