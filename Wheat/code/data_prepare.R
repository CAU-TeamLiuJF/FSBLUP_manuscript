setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix")

#source("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/code/calT.R")

data_wd = "/public/home/liujf/workspace/xueyh/public_data/wheat"

trials = 1
trials = commandArgs(trailingOnly = TRUE)[1]
trials = as.numeric(trials)

print(paste0("传入的Trails值为: ", trials))
options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
library(rrBLUP)


extract_genetic <- function(x, gmat) {
    temp_df <- matrix(NA, nrow = nrow(x), ncol = ncol(x)) %>% 
        as.data.frame()
    rownames(temp_df) <- rownames(x)

    same_id = intersect(rownames(x), colnames(gmat))
    x = x[same_id,]
    G = gmat[same_id, same_id]

    for (i in seq_len(ncol(x))) {

        pheno = x[, c(1, i)] %>%
            rename(id = 1, y = 2)

        temp_df[same_id, i] <- mixed.solve(pheno$y, K = G)$u

        if( i %% 100 == 0) { cat("i:", i, " /", ncol(x), "\n") }

    }

    return(temp_df)
}


# load full data
for (trials in 17:20) {
print(paste0("Trails: ", trials, " /20"))
BLUEs = fread(paste0(data_wd, '/Krause_et_al_2018_Yield_BLUEs.csv'),data.table=F)
BLUPs = fread(paste0(data_wd, '/Krause_et_al_2018_Yield_iid_BLUPs.csv'),data.table=F)
pheno <- dplyr::full_join(BLUEs, BLUPs, by = join_by(GID, `Breeding Cycle`, Managed_Treatment)) %>% 
  dplyr::mutate(Trial = paste(`Breeding Cycle`, Managed_Treatment, sep='_')) %>% 
  tibble::as_tibble()

# pheno %>% count(Managed_Treatment, `Breeding Cycle`)

trialName = unique(pheno$Trial)[trials]
# trialName = "2014-15_Moderate Drought"
pheno = dplyr::filter(pheno, Trial == trialName) %>% select(-`Breeding Cycle`, -Managed_Treatment)
print(paste0("当前为trialName为: ", trialName))
#print(paste0("当前为pheno为: ", paste(dim(pheno), collapse = "_")))

# kinship matrix subset with pheno id

A <- readr::read_csv(paste0(data_wd, "/Krause_et_al_2018_Pedigree.csv"), col_names = T, show_col_types = F) %>% 
  dplyr::rename(ID = 1) %>% 
  dplyr::filter(ID %in% pheno$GID) %>%
  tibble::column_to_rownames("ID") %>% 
  dplyr::select(!!paste0(pheno$GID)) %>% 
  as.matrix()

geno = fread(paste0(data_wd, '/Krause_et_al_2018_Genotypes.csv'), data.table = F)
rownames(geno) = geno[,1]
geno = as.matrix(geno[,-1])
geno = geno[rownames(geno) %in% pheno$GID,]
#print(paste0("当前为geno为: ", paste(dim(geno), collapse = "_")))
geno = geno[,colMeans(!is.na(geno)) >  0.5 & (0.5-abs(0.5-colMeans(geno,na.rm = T)) > 0.05)]
geno[is.na(geno)] = matrix(colMeans(geno,na.rm=T),nr = nrow(geno),nc = ncol(geno),byrow=T)[is.na(geno)]
G = rrBLUP::A.mat(2*geno-1)
D = as.matrix(dist(geno))

HTP <- readr::read_csv(paste0(data_wd, "/Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv"), col_names = T, show_col_types = F) %>% 
  dplyr::mutate(Trail = paste(Breeding_Cycle, Managed_Treatment, sep = '_')) %>% 
  dplyr::filter(Trail == trialName) %>% 
  tidyr::pivot_longer(cols = Hyper_BLUE_398nm:Hyper_BLUE_847nm, values_to = "Hyper") %>% 
  dplyr::mutate(Phenotype = paste(name, Phenotyping_Date, sep='::')) %>% 
  tidyr::pivot_wider(names_from = Phenotype, values_from = Hyper, id_cols = GID) %>%
  dplyr::filter(GID %in% pheno$GID) %>%
  #.[match(pheno$GID, .$GID), ] %>% 
  tibble::column_to_rownames("GID") %>% 
  #extract_genetic(G)%>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
  #calT(2)


pheno <- pheno %>%
	filter(GID %in% colnames(HTP), GID %in% colnames(G))
HTP <- HTP[as.character(pheno$GID), as.character(pheno$GID)]
G <- G[as.character(pheno$GID), as.character(pheno$GID)]
# save data
Trial = paste0(str_sub(trialName, start = 9, end = 9), str_extract(str_split(trialName, " ")[[1]][2], "\\w"), str_sub(trialName, start = 1, end = 4))
if(str_detect(Trial, "NA")) {Trial = str_remove(Trial, "NA")}
if(!dir.exists(paste0("full_data/", Trial))) {dir.create(paste0("full_data/", Trial), recursive = T)}
#if(!dir.exists(paste0(Trial, "/result"))) {dir.create(paste0(Trial, "/result"), recursive = T)}

if(trials %in% 16:20) pheno$Grain_Yield_iid_BLUP = pheno$Grain_Yield_BLUE

wheat.pheno = pheno
save(A, G, D, HTP, geno, wheat.pheno, file = paste0("full_data/", Trial, "/data_adj.rdata"))
}






# classfication data

classification <- data.frame(
  date = c("140110", "140117", "140130", "140207", "140214", "140219", "140227", "140311", "140317"),
  label = c("VEG", "VEG", "VEG", "VEG", "VEG", "HEAD", "HEAD", "GF", "GF")
)

HTP_VEG <- readr::read_csv(paste0(data_wd, "/Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv"), col_names = T, show_col_types = F) %>% 
  dplyr::mutate(Trail = paste(Breeding_Cycle, Managed_Treatment, sep = '_')) %>% 
  dplyr::filter(Trail == trialName) %>% 
  tidyr::pivot_longer(cols = Hyper_BLUE_398nm:Hyper_BLUE_847nm, values_to = "Hyper") %>% 
  dplyr::filter(Phenotyping_Date %in% classification$date[classification$label == "VEG"]) %>%
  dplyr::mutate(Phenotype = paste(name, Phenotyping_Date, sep='::')) %>% 
  tidyr::pivot_wider(names_from = Phenotype, values_from = Hyper, id_cols = GID) %>%
  dplyr::filter(GID %in% pheno$GID) %>%
  #.[match(pheno$GID, .$GID), ] %>% 
  tibble::column_to_rownames("GID") %>% 
  #extract_genetic(G)%>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
cat("\n HTP_VEG done !!! \n")
HTP_HEAD <- readr::read_csv(paste0(data_wd, "/Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv"), col_names = T, show_col_types = F) %>% 
  dplyr::mutate(Trail = paste(Breeding_Cycle, Managed_Treatment, sep = '_')) %>% 
  dplyr::filter(Trail == trialName) %>% 
  tidyr::pivot_longer(cols = Hyper_BLUE_398nm:Hyper_BLUE_847nm, values_to = "Hyper") %>% 
  dplyr::filter(Phenotyping_Date %in% classification$date[classification$label == "HEAD"]) %>%
  dplyr::mutate(Phenotype = paste(name, Phenotyping_Date, sep='::')) %>% 
  tidyr::pivot_wider(names_from = Phenotype, values_from = Hyper, id_cols = GID) %>%
  dplyr::filter(GID %in% pheno$GID) %>%
  #.[match(pheno$GID, .$GID), ] %>% 
  tibble::column_to_rownames("GID") %>% 
  #extract_genetic(G)%>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
cat("\n HTP_HEAD done !!! \n")
HTP_GF <- readr::read_csv(paste0(data_wd, "/Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv"), col_names = T, show_col_types = F) %>% 
  dplyr::mutate(Trail = paste(Breeding_Cycle, Managed_Treatment, sep = '_')) %>% 
  dplyr::filter(Trail == trialName) %>% 
  tidyr::pivot_longer(cols = Hyper_BLUE_398nm:Hyper_BLUE_847nm, values_to = "Hyper") %>% 
  dplyr::filter(Phenotyping_Date %in% classification$date[classification$label == "GF"]) %>%
  dplyr::mutate(Phenotype = paste(name, Phenotyping_Date, sep='::')) %>% 
  tidyr::pivot_wider(names_from = Phenotype, values_from = Hyper, id_cols = GID) %>%
  dplyr::filter(GID %in% pheno$GID) %>%
  #.[match(pheno$GID, .$GID), ] %>% 
  tibble::column_to_rownames("GID") %>% 
  #extract_genetic(G)%>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }
cat("\n HTP_GF done !!! \n")

HTP_VEG <- HTP_VEG[as.character(pheno$GID), as.character(pheno$GID)]
HTP_HEAD <- HTP_HEAD[as.character(pheno$GID), as.character(pheno$GID)]
HTP_GF <- HTP_GF[as.character(pheno$GID), as.character(pheno$GID)]

save(A, G, D, HTP, HTP_VEG, HTP_HEAD, HTP_GF, geno, wheat.pheno, file = paste0("full_data/OF2013/data_adj_classfication_ori.rdata"))