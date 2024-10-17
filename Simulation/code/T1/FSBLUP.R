wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rrBLUP))
suppressPackageStartupMessages(library(FSBLUP))

samnum <- 20 ##  repeat times
cvnum <- 5 ##  cross validation times


for (s in 1:3) {
    for (r in 1:20) {
        # s=1;r=1
        cat("scenario_", s, "rep_", r, "\n")
        results <- data.frame()
        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/scenario_", s,"/rep_", r, "/data_adj.rdata"))

        pheno_all <- pheno %>%
            select(GID = 1, yc = 2, tbv, ngen) %>%
            mutate(GID = as.character(GID), yc = as.numeric(yc), GID2 = GID, GID3 = GID) %>% 
            as.data.frame()

        fs <- FSBLUP(phe = pheno_all, trait_col = 2, M1 = amat, M2 = gmat, M3 = tmat, po.ngen = 10, 
            po.gs.point_num = 25, po.bi.max_iter = 10, po.bi.threshold = 1e-4, 
            stas.phe.col = 3, stas.type = "cor", return.matrix = T)

        gc()

        pheno_all <- pheno_all %>% mutate(partition = case_when(ngen %in% 1:10 ~ "train", TRUE ~ "test"))
        pheno_all$yNA <- pheno_all$yc
        nas <- pheno_all$partition == "test"
        pheno_all$yNA[nas] <- NA

        ## rrBLUP
        res_k <- tryCatch(mixed.solve(pheno_all$yNA, K = fs), error = {fs=adj_pos(fs);mixed.solve(pheno_all$yNA, K = fs)})
        res_k <- res_k$u[nas]

        results <- bind_rows(
            results,
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "FSBLUP",
                pearson = cor(pheno_all$tbv[nas], res_k, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ res_k)$coefficients[2]
            )
        )

        if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/replicate"))) {
            dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/replicate"), recursive = T)
        }
        write_csv(results, paste0(wd, "/scenario_", s, "/rep_", r, "/replicate/FSBLUP_r_scenario_", s, "_rep_", r,"_fsmat.csv"))

    }
}
