setwd("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2")
suppressPackageStartupMessages(library(tidyverse))
library(rrBLUP)
suppressPackageStartupMessages(library(data.table))
library(pedigree)

sce <- as.numeric(commandArgs(t=T)[1])
if(is.na(sce)) sce = 1

arg <- as.numeric(commandArgs(t=T)[2])
if(is.na(arg)) arg = 1

print(paste0("Rscript || Scenario: ", sce, " Rep: ", arg))


tbv <- read_delim("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/tbv.txt", show_col_types = F, col_names = F) %>% 
  mutate(id = row_number())

pheno_all <- read_delim("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/my_bv.txt", show_col_types = F) %>% 
  left_join(tbv[,3:4], by = "id") %>% 
  rename(y = 2, tbv = 3) %>%
  mutate(ngen = c(rep(1,1100), rep(2:11, each = 2000)))

set.seed(arg)
pheno <- pheno_all %>% 
  group_by(ngen) %>%
  slice_sample(prop = 0.1) %>%
  ungroup() %>% 
  arrange(ngen, id)
  #%>% mutate(y = case_when(ngen == 11 ~ NA, TRUE ~ y))

pedi <- read_delim("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/ped", col_names = F, show_col_types = F) %>% 
  rename(id = 1, sire = 2, dam = 3) %>% 
  as.data.frame()
pedi <- pedi[order(orderPed(pedi)),]
makeA(pedi, which = pedi[, 1] %in% pheno$id)
amat <- fread("A.txt", header = F) %>% 
  bind_rows(data.frame(V1 = .$V2, V2 = .$V1, V3 = .$V3)) %>% 
  unique() %>%  
  tidyr::pivot_wider(names_from = V2, values_from = V3, values_fill = 0) %>% 
  column_to_rownames("V1") %>% 
  as.matrix  

if(sce == 1) {

    gmat <- fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/G", data.table = F) %>%
      filter(V1 %in% pheno$id, V2 %in% pheno$id) %>% #pheno$id[pheno$ngen %in% 6:11]
      mutate(V1 = as.character(V1), V2 = as.character(V2)) %>%
      bind_rows(data.frame(V1 = .$V2, V2 = .$V1, V3 = .$V3)) %>% 
      unique() %>%  
      tidyr::pivot_wider(names_from = V2, values_from = V3, values_fill = 0) %>% 
      column_to_rownames("V1") %>% 
      as.matrix


    tmat <- fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/M.txt", data.table = F) %>%
      mutate(id = row_number()) %>% 
      filter(row_number() %in% as.character(pheno$id)) %>%  #pheno$id[pheno$ngen %in% 10:11]
      column_to_rownames("id") %>% 
      as.matrix() %>% 
      scale() %>% 
      { tcrossprod(.) / ncol(.) }

    if(!fs::dir_exists(paste0("scenario_", sce, "/rep_", arg))) fs::dir_create(paste0("scenario_", sce, "/rep_", arg))
    save(pheno, gmat, tmat, amat, file = paste0("scenario_", sce, "/rep_", arg, "/data.rdata"))

} else if (sce == 2) {

    gmat <- fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/G", data.table = F) %>%
      filter(V1 %in% pheno$id, V2 %in% pheno$id) %>% #pheno$id[pheno$ngen %in% 6:11]
      mutate(V1 = as.character(V1), V2 = as.character(V2)) %>%
      bind_rows(data.frame(V1 = .$V2, V2 = .$V1, V3 = .$V3)) %>% 
      unique() %>%  
      tidyr::pivot_wider(names_from = V2, values_from = V3, values_fill = 0) %>% 
      column_to_rownames("V1") %>% 
      as.matrix

    tmat <- fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/M.txt", data.table = F) %>%
      mutate(id = row_number()) %>% 
      filter(row_number() %in% as.character(pheno$id[pheno$ngen %in% 10:11])) %>%  #pheno$id[pheno$ngen %in% 10:11]
      column_to_rownames("id") %>% 
      as.matrix() %>% 
      scale() %>% 
      { tcrossprod(.) / ncol(.) }

    if(!fs::dir_exists(paste0("scenario_", sce, "/rep_", arg))) fs::dir_create(paste0("scenario_", sce, "/rep_", arg))
    save(pheno, gmat, tmat, amat, file = paste0("scenario_", sce, "/rep_", arg, "/data.rdata"))

} else if (sce == 3) {

    gmat <- fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/G", data.table = F) %>% 
      filter(V1 %in% pheno$id[pheno$ngen %in% 6:11], V2 %in% pheno$id[pheno$ngen %in% 6:11]) %>% #pheno$id[pheno$ngen %in% 6:11]
      mutate(V1 = as.character(V1), V2 = as.character(V2)) %>%
      bind_rows(data.frame(V1 = .$V2, V2 = .$V1, V3 = .$V3)) %>% 
      unique() %>%  
      tidyr::pivot_wider(names_from = V2, values_from = V3, values_fill = 0) %>% 
      column_to_rownames("V1") %>% 
      as.matrix
    
    tmat <- fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/M.txt", data.table = F) %>%
      mutate(id = row_number()) %>% 
      filter(row_number() %in% as.character(pheno$id[pheno$ngen %in% 10:11])) %>%  #pheno$id[pheno$ngen %in% 10:11]
      column_to_rownames("id") %>% 
      as.matrix() %>% 
      scale() %>% 
      { tcrossprod(.) / ncol(.) }

    if(!fs::dir_exists(paste0("scenario_", sce, "/rep_", arg))) fs::dir_create(paste0("scenario_", sce, "/rep_", arg))
    save(pheno, gmat, tmat, amat, file = paste0("scenario_", sce, "/rep_", arg, "/data.rdata"))

} else {
    stop("Invalid scenario")
}


