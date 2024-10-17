setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig")


options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
library(rrBLUP)

load("data/pig_data2.rdata")
id.code <- read_delim("data/id.code.txt")


ge <- readRDS("data/yy_gene_fpkm_clean.RDS") %>% 
  rownames_to_column("id") %>% 
  left_join(id.code, by = "id") %>% 
  select(-id) %>% 
  filter(num %in% rownames(TT)) %>% 
  column_to_rownames("num") %>% 
  as.data.frame()

ge[1:10,1:10]

gmat <- GG[rownames(TT), colnames(TT)]
ge <- ge[rownames(TT), ]

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

ge_add <- extract_genetic(ge, gmat)
saveRDS(ge_add, "data/yy_gene_fpkm_clean_add.RDS")

TT2 <- ge_add %>% 
  as.matrix() %>% 
  scale() %>% 
  { tcrossprod(.) / ncol(.) }


geno <- fread("data/geno.raw", data.table = F) %>% 
  select(-c(1, 3:6), id = 2) 

pheno <- phe %>% rename(id = 1)

save(AA, GG, TT, TT2, pheno, geno, file = paste0("data/data_adj.rdata"))

load("data/data_adj.rdata")
