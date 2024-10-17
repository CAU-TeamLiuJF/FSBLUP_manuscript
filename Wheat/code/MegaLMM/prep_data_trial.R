setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/MegaLMM")

library(data.table)
library(rrBLUP)
library(tidyr)

trial = as.numeric(commandArgs(t=T)[1])
if(is.na(trial)) trial = 1
foldid = as.numeric(commandArgs(t=T)[2])
if(is.na(foldid)) foldid = 1

print(sprintf(' 当前prep_data_trial.R  Trial 为 Trial %d, Fold(SLURM_ARRAY_TASK_ID) %d',trial,foldid))

#MD2014
#trial = 6 

folder = sprintf('Trial_%02d',trial)
try(dir.create(folder))
try(dir.create(sprintf('%s/Results',folder)))


# load full data

BLUEs = fread('data/Krause_et_al_2018_Yield_BLUEs.csv',data.table=F)
BLUPs = fread('data/Krause_et_al_2018_Yield_iid_BLUPs.csv',data.table=F)
BLUEs$BLUP = BLUPs$Grain_Yield_iid_BLUP
BLUEs$Trial = paste(BLUEs$`Breeding Cycle`,BLUEs$Managed_Treatment,sep='_')
trialName = unique(BLUEs$Trial)[trial]
BLUEs = subset(BLUEs,Trial == trialName)

## r$> unique(BLUEs$Trial)
##  [1] "2013-14_Moderate Drought" "2013-14_Optimal Bed"      "2013-14_Heat"             "2013-14_Severe Drought"   "2013-14_Optimal Flat"     ## "2014-15_Moderate Drought" "2014-15_Optimal Bed"     
##  [8] "2014-15_Heat"             "2014-15_Severe Drought"   "2014-15_Optimal Flat"     "2015-16_Moderate Drought" "2015-16_Optimal Bed"      ## "2015-16_Heat"             "2015-16_Severe Drought"  
## [15] "2015-16_Optimal Flat"     "2016-17_Moderate Drought" "2016-17_Optimal Bed"      "2016-17_Heat"             "2016-17_Severe Drought"   ## "2016-17_Optimal Flat"


# BLUPs in trials 16:20 look bad. The correlation with the BLUEs is <.26, and the prediction accuracy is pretty bad. I'll use BLUEs here instead.
if(trial %in% 16:20) BLUEs$BLUP = BLUEs$Grain_Yield_BLUE

geno = fread('data/Krause_et_al_2018_Genotypes.csv',data.table=F)
rownames(geno) = geno[,1]
geno = as.matrix(geno[,-1])
geno = geno[rownames(geno) %in% BLUEs$GID,]

HTP = fread('data/Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv',data.table = F)
HTP$Trial = paste(HTP$Breeding_Cycle,HTP$Managed_Treatment,sep='_')
HTP = subset(HTP,Trial == trialName)
HTP_tall = pivot_longer(HTP,cols = Hyper_BLUE_398nm:Hyper_BLUE_847nm,values_to = 'Hyper')
HTP_tall$Phenotype = paste(HTP_tall$name,HTP_tall$Phenotyping_Date,sep='::')

HTP_wide = as.data.frame(pivot_wider(HTP_tall,names_from = Phenotype,values_from = Hyper,id_cols = GID))
geno = geno[rownames(geno) %in% HTP_wide$GID,]
BLUEs = BLUEs[match(rownames(geno),BLUEs$GID),]
HTP_wide = HTP_wide[match(BLUEs$GID,HTP_wide$GID),]
rownames(HTP_wide) = HTP_wide$GID
HTP_wide = as.matrix(HTP_wide[,-1])


# clean geno data

# drop markers with > 50% missing, or MAF < 0.05
geno = geno[,colMeans(!is.na(geno)) >  0.5 & (0.5-abs(0.5-colMeans(geno,na.rm = T)) > 0.05)]
geno[is.na(geno)] = matrix(colMeans(geno,na.rm=T),nr = nrow(geno),nc = ncol(geno),byrow=T)[is.na(geno)]

# calculate K and D
K = A.mat(2*geno-1)
D = as.matrix(dist(geno))

# save data
write.csv(BLUEs,file = sprintf('%s/BLUES.csv',folder),row.names=F)
write.csv(geno,file = sprintf('%s/X.csv',folder),row.names=T)
write.csv(HTP_wide,file = sprintf('%s/HTP.csv',folder),row.names=T)
write.csv(K,file = sprintf('%s/K.csv',folder),row.names=T)
write.csv(D,file = sprintf('%s/D.csv',folder),row.names=T)

