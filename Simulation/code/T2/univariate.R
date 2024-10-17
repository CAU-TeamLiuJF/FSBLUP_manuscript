wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
library(rrBLUP)
suppressPackageStartupMessages(library(lme4))
library(lme4qtl)
library(BGLR)
library(FSBLUP)

adj_pos <- function(x) {
    if (!is.null(x)) {
        diag(x) <- diag(x) + 1e-06
        B <- try(chol(x), silent = TRUE)
            if (inherits(B, what = "try-error")) {
                x <- Matrix::nearPD(x)$mat %>% as.matrix()
                rm(B)
            }
    }
    return(x)
}


geno <- data.table::fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_c2021/r_datasim/p0_mrk_001.txt", data.table = F) %>% 
    column_to_rownames("ID") %>% 
    as.matrix()

geno[is.na(geno)] = matrix(colMeans(geno,na.rm=T),nr = nrow(geno),nc = ncol(geno),byrow=T)[is.na(geno)]


for (s in 1) {

    for (r in 1:20) {
        cat("scenario_", s, "rep_", r, "\n")
        results <- data.frame()
        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/scenario_", s,"/rep_", r, "/data.rdata"))

        pheno_all <- pheno %>%
        select(GID = 1, yc = 2, tbv, ngen) %>%
        mutate(GID = as.character(GID), yc = as.numeric(yc), GID2 = GID, GID3 = GID) %>% 
        as.data.frame()

        amat <- adj_pos(amat)
        gmat <- adj_pos(gmat)
        tmat <- adj_pos(tmat)

        pheno_all <- pheno_all %>% mutate(partition = case_when(ngen %in% 1:10 ~ "train", TRUE ~ "test"))
        pheno_all$yNA <- pheno_all$yc
        nas <- pheno_all$partition == "test"
        pheno_all$yNA[nas] <- NA

        ## rrBLUP - A
        res_a <- mixed.solve(pheno_all$yNA, K = amat)
        res_a <- res_a$u[nas]

        ## rrBLUP - G
        res_k <- mixed.solve(pheno_all$yNA, K = gmat)
        h2 <- res_k$Vu / (res_k$Vu + res_k$Ve)
        res_k <- res_k$u[nas]
        
        ## rrBLUP - T
        res_t <- mixed.solve(pheno_all$yNA, K = tmat)
        res_t <- res_t$u[nas]

        ## lme4 - G + T
        goo <- gmat[!nas, !nas]
        gno <- gmat[nas, !nas]
        too <- tmat[!nas, !nas]
        tno <- tmat[nas, !nas]
        
        res_gt <- relmatLmer(yc ~ (1 | GID) + (1 | GID2), data = pheno_all[!nas, ], relmat = list(GID = gmat, GID2 = tmat))
        u1 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(res_gt)[[1]]))[rownames(goo), 1]
        u2 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(res_gt)[[2]]))[rownames(too), 1]
        pred_lme4qtl_gt <- gno %*% MASS::ginv(goo) %*% u1 + tno %*% MASS::ginv(too) %*% u2

        ## lme4 - G + T + A
        goo <- gmat[!nas, !nas]
        gno <- gmat[nas, !nas]
        too <- tmat[!nas, !nas]
        tno <- tmat[nas, !nas]
        aoo <- amat[!nas, !nas]
        ano <- amat[nas, !nas]

        res_gta <- relmatLmer(yc ~ (1 | GID) + (1 | GID2)+ (1 | GID3), data = pheno_all[!nas, ], relmat = list(GID = gmat, GID2 = tmat, GID3 = amat))
        u1 <- as.matrix(t(res_gta@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(res_gta)[[1]]))[rownames(goo), 1]
        u2 <- as.matrix(t(res_gta@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(res_gta)[[2]]))[rownames(too), 1]
        u3 <- as.matrix(t(res_gta@optinfo$relmat$relfac$GID3) %*% as.matrix(ranef(res_gta)[[3]]))[rownames(aoo), 1]
        pred_lme4qtl_gta <- gno %*% MASS::ginv(goo) %*% u1 + tno %*% MASS::ginv(too) %*% u2 + ano %*% MASS::ginv(aoo) %*% u3

        ##BGLR - BL - default priors
        X <- geno[pheno_all$GID,]
        bglr_BL = BGLR(y = pheno_all$yNA,ETA = list(list(X = X,model = 'BL')),burnIn = 5000,nIter = 50000,verbose=F)
        bglr_BLs = X[nas,] %*% bglr_BL$ETA[[1]]$b


        # result
        results <- bind_rows(
            results,
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "A",
                pearson = cor(pheno_all$tbv[nas], res_a, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ res_a)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "G",
                pearson = cor(pheno_all$tbv[nas], res_k, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ res_k)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "T",
                pearson = cor(pheno_all$tbv[nas], res_t, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ res_t)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "G + T",
                pearson = cor(pheno_all$tbv[nas], pred_lme4qtl_gt, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ pred_lme4qtl_gt)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "G + T + A",
                pearson = cor(pheno_all$tbv[nas], pred_lme4qtl_gta, use = 'pairwise.complete.obs'),
                bias = lm(pheno_all$tbv[nas] ~ pred_lme4qtl_gta)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "BL",
                pearson = cor(pheno_all$tbv[nas], bglr_BLs, use = 'pairwise.complete.obs'),
                bias = lm(pheno_all$tbv[nas] ~ bglr_BLs)$coefficients[2]
            )
        )

    if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/univariate"))) {
        dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/univariate"), recursive = T)
    }   
        write_csv(results, paste0(wd, "/scenario_", s, "/rep_", r, "/univariate/GT_scenario_", s, "_rep_", r,".csv"))
    }
    
}


for (s in 2) {
    
    for (r in 1:20) {
        results <- data.frame()
        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/scenario_", s,"/rep_", r, "/data.rdata"))

        pheno_all <- pheno %>%
        select(GID = 1, yc = 2, tbv, ngen) %>%
        mutate(GID = as.character(GID), yc = as.numeric(yc), GID2 = GID) %>% 
        as.data.frame()

        amat <- adj_pos(amat)
        gmat <- adj_pos(gmat)
        tmat <- adj_pos(tmat)

        
        pheno_all <- pheno_all %>% mutate(partition = case_when(ngen %in% 1:10 ~ "train", TRUE ~ "test"))
        pheno_all$yNA <- pheno_all$yc
        nas <- pheno_all$partition == "test"
        pheno_all$yNA[nas] <- NA

        ## rrBLUP - A
        res_a <- mixed.solve(pheno_all$yNA, K = amat)
        res_a <- res_a$u[nas]

        ## rrBLUP - G
        res_k <- mixed.solve(pheno_all$yNA, K = gmat)
        h2 <- res_k$Vu / (res_k$Vu + res_k$Ve)
        res_k <- res_k$u[nas]


        pheno_all_2 <- pheno_all %>% 
            filter(ngen %in% 10:11)

        pheno_all_2$yNA <- pheno_all_2$yc
        nas2 <- pheno_all_2$partition == "test"
        pheno_all_2$yNA[nas2] <- NA
        
        ## rrBLUP - T
        res_t <- mixed.solve(pheno_all_2$yNA, K = tmat)
        res_t <- res_t$u[nas2]

        ## lme4 - G + T
        gmat <- gmat[pheno_all_2$GID, pheno_all_2$GID]
        goo <- gmat[!nas2, !nas2]
        gno <- gmat[nas2, !nas2]
        too <- tmat[!nas2, !nas2]
        tno <- tmat[nas2, !nas2]
        
        res_gt <- relmatLmer(yc ~ (1 | GID) + (1 | GID2), data = pheno_all_2[!nas, ], relmat = list(GID = gmat, GID2 = tmat))
        u1 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(res_gt)[[1]]))[rownames(goo), 1]
        u2 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(res_gt)[[2]]))[rownames(too), 1]
        pred_lme4qtl_gt <- gno %*% MASS::ginv(goo) %*% u1 + tno %*% MASS::ginv(too) %*% u2
        

        # result
        results <- rbind(
            results,
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "A",
                pearson = cor(pheno_all$tbv[nas], res_a, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ res_a)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "G",
                pearson = cor(pheno_all$tbv[nas], res_k, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ res_k)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "T",
                pearson = cor(pheno_all_2$tbv[nas2], res_t, use = "pairwise.complete.obs"),
                bias = lm(pheno_all_2$tbv[nas2] ~ res_t)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "G + T",
                pearson = cor(pheno_all_2$tbv[nas2], pred_lme4qtl_gt, use = "pairwise.complete.obs"),
                bias = lm(pheno_all_2$tbv[nas2] ~ pred_lme4qtl_gt)$coefficients[2]
            )
        )

    if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/univariate"))) {
        dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/univariate"), recursive = T)
    }
        write_csv(results, paste0(wd, "/scenario_", s, "/rep_", r, "/univariate/univariate_scenario_", s, "_rep_", r,".csv"))
    }

}


for (s in 3) {
    
    for (r in 1:20) {
        results <- data.frame()
        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/scenario_", s,"/rep_", r, "/data.rdata"))

        pheno_all <- pheno %>%
        select(GID = 1, yc = 2, tbv, ngen) %>%
        mutate(GID = as.character(GID), yc = as.numeric(yc), GID2 = GID) %>% 
        as.data.frame()

        amat <- adj_pos(amat)
        gmat <- adj_pos(gmat)
        tmat <- adj_pos(tmat)

        
        pheno_all <- pheno_all %>% mutate(partition = case_when(ngen %in% 1:10 ~ "train", TRUE ~ "test"))
        pheno_all$yNA <- pheno_all$yc
        nas <- pheno_all$partition == "test"
        pheno_all$yNA[nas] <- NA

        ## rrBLUP - A
        res_a <- mixed.solve(pheno_all$yNA, K = amat)
        res_a <- res_a$u[nas]


        pheno_all_3 <- pheno_all %>% 
            filter(ngen %in% 6:11)

        pheno_all_3$yNA <- pheno_all_3$yc
        nas3 <- pheno_all_3$partition == "test"
        pheno_all_3$yNA[nas3] <- NA

        gamt <- gmat[pheno_all_3$GID, pheno_all_3$GID]
        ## rrBLUP - G
        res_k <- mixed.solve(pheno_all_3$yNA, K = gmat)
        h2 <- res_k$Vu / (res_k$Vu + res_k$Ve)
        res_k <- res_k$u[nas3]


        pheno_all_2 <- pheno_all %>% 
            filter(ngen %in% 10:11)

        pheno_all_2$yNA <- pheno_all_2$yc
        nas2 <- pheno_all_2$partition == "test"
        pheno_all_2$yNA[nas2] <- NA
        
        ## rrBLUP - T
        res_t <- mixed.solve(pheno_all_2$yNA, K = tmat)
        res_t <- res_t$u[nas2]

        ## lme4 - G + T
        gmat <- gmat[pheno_all_2$GID, pheno_all_2$GID]
        goo <- gmat[!nas2, !nas2]
        gno <- gmat[nas2, !nas2]
        too <- tmat[!nas2, !nas2]
        tno <- tmat[nas2, !nas2]
        
        res_gt <- relmatLmer(yc ~ (1 | GID) + (1 | GID2), data = pheno_all_2[!nas, ], relmat = list(GID = gmat, GID2 = tmat))
        u1 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(res_gt)[[1]]))[rownames(goo), 1]
        u2 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(res_gt)[[2]]))[rownames(too), 1]
        pred_lme4qtl_gt <- gno %*% MASS::ginv(goo) %*% u1 + tno %*% MASS::ginv(too) %*% u2
        

        # result
        results <- rbind(
            results,
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "A",
                pearson = cor(pheno_all$tbv[nas], res_a, use = "pairwise.complete.obs"),
                bias = lm(pheno_all$tbv[nas] ~ res_a)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "G",
                pearson = cor(pheno_all_3$tbv[nas3], res_k, use = "pairwise.complete.obs"),
                bias = lm(pheno_all_3$tbv[nas3] ~ res_k)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "T",
                pearson = cor(pheno_all_2$tbv[nas2], res_t, use = "pairwise.complete.obs"),
                bias = lm(pheno_all_2$tbv[nas2] ~ res_t)$coefficients[2]
            ),
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "G + T",
                pearson = cor(pheno_all_2$tbv[nas2], pred_lme4qtl_gt, use = "pairwise.complete.obs"),
                bias = lm(pheno_all_2$tbv[nas2] ~ pred_lme4qtl_gt)$coefficients[2]
            )
        )

    if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/univariate"))) {
        dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/univariate"), recursive = T)
    } 
        write_csv(results, paste0(wd, "/scenario_", s, "/rep_", r, "/univariate/univariate_scenario_", s, "_rep_", r,".csv"))
    }

}