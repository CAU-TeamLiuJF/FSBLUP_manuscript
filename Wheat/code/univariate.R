wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rrBLUP))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lme4qtl))
suppressPackageStartupMessages(library(BGLR))
suppressPackageStartupMessages(library(FSBLUP))

suppressPackageStartupMessages(source("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/data/Estimate_gcor_prediction.R"))

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

# read genotype data
geno <- data.table::fread("MegaLMM/data/Krause_et_al_2018_Genotypes.csv", data.table = F)
rownames(geno) = geno[,1]
geno = as.matrix(geno[,-1])



samnum <- 20 ##  repeat times 
cvnum <- 2 ##  cross validation times
Trial <- 'OF2013'

times <- data.frame() ## record time

## GBLUP
start_time <- Sys.time()
#load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/full_data/", Trial, "/data_adj.rdata"))
pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID), GID2 = GID) %>% 
    as.data.frame()


G <- G[pheno_all$GID, pheno_all$GID]

G_c <- G
G <- adj_pos(G)
results <- data.frame()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
        #cv_s=1; j=1
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$BLUE
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA

        Knn = G_c[nas,nas]
        sKnn = svd(Knn)
        ## rrBLUP - G
        res_k <- mixed.solve(pheno_all$yNA, K = G)
        res_k <- res_k$u[nas]

        # result
        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                Method = "GBLUP",
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = res_k),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']],
                bias = lm(pheno_all$BLUP[nas] ~ res_k)$coefficients[2]
            )
        )
        
    }
    end_time <- Sys.time()
    times <- rbind(times, data.frame(Trial = Trial, method = "GBLUP", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | GBLUP | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | GBLUP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
}
if(!dir.exists(paste0(getwd(), "/", Trial, "/results"))) { dir.create(paste0(getwd(), "/", Trial, "/results"),recursive = T) }
readr::write_csv(results, paste0(getwd(), "/", Trial, "/results/GBLUP_2x5__rep_20_20240808_adj.csv"))
if(!dir.exists(paste0(getwd(), "/", Trial, "/test_time"))) { dir.create(paste0(getwd(), "/", Trial, "/test_time"),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/", Trial, "/test_time/GBLUP_2x5__rep_20_times_20240808_adj.csv"))



## GBLUP_HTP
times <- data.frame() ## record time
start_time <- Sys.time()
#load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/full_data/", Trial, "/data_adj.rdata"))
pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID), GID2 = GID) %>% 
    as.data.frame()

G <- G[pheno_all$GID, pheno_all$GID]
HTP <- HTP[pheno_all$GID, pheno_all$GID]
G_c <- G
HTP <- adj_pos(HTP)
results <- data.frame()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$BLUE
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA

        Knn = G_c[nas,nas]
        sKnn = svd(Knn)

        ## rrBLUP - HTP kernel
        res_t <- mixed.solve(pheno_all$yNA, K = HTP)
        res_t <- res_t$u[nas]

        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                Method = "GBLUP_HTP",
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = res_t),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']],
                bias = lm(pheno_all$BLUP[nas] ~ res_t)$coefficients[2]
            )
        )
        
    }
    end_time <- Sys.time()
    times <- rbind(times, data.frame(Trial = Trial, method = "GBLUP_HTP", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | GBLUP_HTP | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | GBLUP_HTP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
}
if(!dir.exists(paste0(getwd(), "/", Trial, "/results"))) { dir.create(paste0(getwd(), "/", Trial, "/results"),recursive = T) }
readr::write_csv(results, paste0(getwd(), "/", Trial, "/results/GBLUP_HTP_2x5__summum_20_20240808_adj.csv"))
if(!dir.exists(paste0(getwd(), "/", Trial, "/test_time"))) { dir.create(paste0(getwd(), "/", Trial, "/test_time"),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/", Trial, "/test_time/GBLUP_HTP_2x5__summum_20_times_20240808_adj.csv"))



## GBLUP_HTP+G
times <- data.frame() ## record time
start_time <- Sys.time()
#load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/full_data/", Trial, "/data_adj.rdata"))
pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID), GID2 = GID) %>% 
    as.data.frame()

G <- G[pheno_all$GID, pheno_all$GID]
HTP <- HTP[pheno_all$GID, pheno_all$GID]
G_c <- G
HTP <- adj_pos(HTP)
results <- data.frame()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$BLUE
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA

        Knn = G_c[nas,nas]
        sKnn = svd(Knn)

        goo <- G[!nas, !nas]
        gno <- G[nas, !nas]
        too <- HTP[!nas, !nas]
        tno <- HTP[nas, !nas]

        ## lme4 - G + HTP
        res_gt <- relmatLmer(BLUP ~ (1 | GID) + (1 | GID2), data = pheno_all[!nas, ], relmat = list(GID = G, GID2 = HTP))
        u1 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(res_gt)[[1]]))[rownames(goo), 1]
        u2 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(res_gt)[[2]]))[rownames(too), 1]
        pred_lme4qtl_gt <- gno %*% MASS::ginv(goo) %*% u1 + tno %*% MASS::ginv(too) %*% u2


        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                Method = "G+HTP",
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = pred_lme4qtl_gt),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']],
                bias = lm(pheno_all$BLUP[nas] ~ pred_lme4qtl_gt)$coefficients[2]
            )
        )
        
    }
    end_time <- Sys.time()
    times <- rbind(times, data.frame(Trial = Trial, method = "G+HTP", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | G+HTP | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | G+HTP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
}
if(!dir.exists(paste0(getwd(), "/", Trial, "/results"))) { dir.create(paste0(getwd(), "/", Trial, "/results"),recursive = T) }
readr::write_csv(results, paste0(getwd(), "/", Trial, "/results/G+HTP_2x5__summum_20_20240808_adj.csv"))
if(!dir.exists(paste0(getwd(), "/", Trial, "/test_time"))) { dir.create(paste0(getwd(), "/", Trial, "/test_time"),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/", Trial, "/test_time/G+HTP_2x5__summum_20_times_20240808_adj.csv"))



## RKHS
times <- data.frame() ## record time
start_time <- Sys.time()
#load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/full_data/", Trial, "/data_adj.rdata"))
pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID), GID2 = GID) %>% 
    as.data.frame()

G <- G[pheno_all$GID, pheno_all$GID]
D <- D[pheno_all$GID, pheno_all$GID]
G_c <- G
results <- data.frame()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$BLUE
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA
        
        Knn = G_c[nas,nas]
        sKnn = svd(Knn)

        ## rrBLUP - RKHS
        res_RKHS = kin.blup(pheno_all, geno = 'GID',pheno = 'yNA',GAUSS = T,K = D)$pred
        res_RKHS = res_RKHS[nas]

        # result
        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                Method = "RKHS",
                pearson_blup = cor(pheno_all$BLUP[nas], res_RKHS, use = 'pairwise.complete.obs'),
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = res_RKHS),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']],
                bias = lm(pheno_all$BLUP[nas] ~ res_RKHS)$coefficients[2]
            )
        )
        
    }
    end_time <- Sys.time()
    times <- rbind(times, data.frame(Trial = Trial, method = "RKHS", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | RKHS | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | RKHS | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
}
if(!dir.exists(paste0(getwd(), "/", Trial, "/results"))) { dir.create(paste0(getwd(), "/", Trial, "/results"),recursive = T) }
readr::write_csv(results, paste0(getwd(), "/", Trial, "/results/RKHS_2x5__summum_20_20240808_adj.csv"))
if(!dir.exists(paste0(getwd(), "/", Trial, "/test_time"))) { dir.create(paste0(getwd(), "/", Trial, "/test_time"),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/", Trial, "/test_time/RKHS_2x5__summum_20_times_20240808_adj.csv"))



## BL
times <- data.frame() ## record time
start_time <- Sys.time()
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))

pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID), GID2 = GID) %>% 
    as.data.frame()

G <- G[pheno_all$GID, pheno_all$GID]
X = geno[rownames(geno) %in% pheno_all$GID,]
X = X[,colMeans(!is.na(X)) >  0.5 & (0.5-abs(0.5-colMeans(X,na.rm = T)) > 0.05)]
X[is.na(X)] = matrix(colMeans(X,na.rm=T),nr = nrow(X),nc = ncol(X),byrow=T)[is.na(X)]
G_c <- G
results <- data.frame()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$BLUE
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA

        Knn = G_c[nas,nas]
        sKnn = svd(Knn)

        ##BGLR - BL - default priors
        bglr_BL = BGLR(y = pheno_all$yNA,ETA = list(list(X = X,model = 'BL')),burnIn = 5000,nIter = 50000,verbose=F)
        bglr_BLs = X[nas,] %*% bglr_BL$ETA[[1]]$b

        # result
        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                Method = "BL",
                pearson_blup = cor(pheno_all$BLUP[nas], bglr_BLs, use = 'pairwise.complete.obs'),
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = bglr_BLs),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]
            )
        )
        
    }
    end_time <- Sys.time()
    times <- rbind(times, data.frame(Trial = Trial, method = "BL", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | BL | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | BL | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
}
if(!dir.exists(paste0(getwd(), "/", Trial, "/results"))) { dir.create(paste0(getwd(), "/", Trial, "/results"),recursive = T) }
readr::write_csv(results, paste0(getwd(), "/", Trial, "/results/BL_2x5__summum_20.csv"))
if(!dir.exists(paste0(getwd(), "/", Trial, "/test_time"))) { dir.create(paste0(getwd(), "/", Trial, "/test_time"),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/", Trial, "/test_time/BL_2x5__summum_20_times.csv"))

## GBLUP_HTP+G+A
times <- data.frame() ## record time
start_time <- Sys.time()
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))
#load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/full_data/", Trial, "/data_adj.rdata"))
pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID), GID2 = GID, GID3 = GID) %>% 
    as.data.frame()

A <- A[pheno_all$GID, pheno_all$GID]
G <- G[pheno_all$GID, pheno_all$GID]
HTP <- HTP[pheno_all$GID, pheno_all$GID]
G_c <- G
HTP <- adj_pos(HTP)
results <- data.frame()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
        # cv_s <- 1; j <- 1
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$BLUE
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA

        Knn = G_c[nas,nas]
        sKnn = svd(Knn)

        goo <- G[!nas, !nas]
        gno <- G[nas, !nas]
        too <- HTP[!nas, !nas]
        tno <- HTP[nas, !nas]
        aoo <- A[!nas, !nas]
        ano <- A[nas, !nas]

        ## lme4 - G + HTP
        res_gt <- relmatLmer(BLUP ~ (1 | GID) + (1 | GID2) + (1 | GID3), data = pheno_all[!nas, ], relmat = list(GID = G, GID2 = HTP, GID3 = A))
        u1 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(res_gt)[[1]]))[rownames(goo), 1]
        u2 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(res_gt)[[2]]))[rownames(too), 1]
        u3 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID3) %*% as.matrix(ranef(res_gt)[[3]]))[rownames(aoo), 1]
        pred_lme4qtl_gta <- gno %*% MASS::ginv(goo) %*% u1 + tno %*% MASS::ginv(too) %*% u2 + ano %*% MASS::ginv(aoo) %*% u3


        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                Method = "G+HTP+A",
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = pred_lme4qtl_gta),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']],
                bias = lm(pheno_all$BLUP[nas] ~ pred_lme4qtl_gta)$coefficients[2]
            )
        )
    }
    end_time <- Sys.time()
    times <- rbind(times, data.frame(Trial = Trial, method = "G+HTP+A", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | G+HTP+A | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | G+HTP+A | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
}
if(!dir.exists(paste0(getwd(), "/", Trial, "/results"))) { dir.create(paste0(getwd(), "/", Trial, "/results"),recursive = T) }
readr::write_csv(results, paste0(getwd(), "/", Trial, "/results/G+HTP+A_5x20.csv"))
if(!dir.exists(paste0(getwd(), "/", Trial, "/test_time"))) { dir.create(paste0(getwd(), "/", Trial, "/test_time"),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/", Trial, "/test_time/G+HTP+A_5x20.csv"))



## FSBLUP

acc <- function(x1, x2, G = G_c) {
  Knn <- G[names(x2), names(x2)]
  sKnn <- svd(Knn)
  val <- estimate_gcor(data.frame(ID=names(x2),obs = x1,pred = x2),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]
}

Trial <- c("OF2013")

times <- data.frame() ## record time
start_time <- Sys.time()

load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/full_data/", Trial, "/data_adj.rdata"))
#load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))
G_c <- G

pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID))

A <- A[pheno_all$GID, pheno_all$GID]
G <- G[pheno_all$GID, pheno_all$GID]
HTP <- HTP[pheno_all$GID, pheno_all$GID]
G_C <- G_c[pheno_all$GID, pheno_all$GID]

## FS matrix
fs = FSBLUP(phe = pheno_all, trait_col = "BLUE", M1 = A, M2 = G, M3 = HTP,
				po.crv.num = 2, po.crv.rep.num = 2, po.gs.point_num = 25, 
				po.bi.max_iter = 10, po.bi.threshold = 1e-4, 
				stas.type = "other", stas.phe.col = 3, stas.fn = acc, return.matrix = T)
gc()

results <- data.frame()
for (cv_s in seq_len(samnum)) {
        start_time1 <- Sys.time()
        for(j in seq_len(cvnum)) {
            set.seed(cv_s)
            pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1/cvnum), times = cvnum)))
            pheno_all$yNA <- pheno_all$BLUE
            nas = pheno_all$partition == j
            pheno_all$yNA[nas] <- NA

            Knn = G_c[nas,nas]
            sKnn = svd(Knn)

            ## rrBLUP 
            res_k <- tryCatch(mixed.solve(pheno_all$yNA, K = fs), error = {fs=adj_pos(fs);mixed.solve(pheno_all$yNA, K = fs)})
            res_k <- res_k$u[nas]

            # result
            results = rbind(
                results,
                data.frame(
                    cv_seed = cv_s,
                    cvnum = j,
                    Method = 'FSBLUP',
                    g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = res_k),Knn,sKnn,method = 'MCMCglmm', normalize = T)[['g_cor']],
                    bias = lm(pheno_all$BLUE[nas] ~ res_k)$coefficients[2]
                )
            )
        }
	end_time <- Sys.time()
    times <- rbind(times, data.frame(Trial = Trial, method = "FSBLUP", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | FSBLUP | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | FSBLUP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))

    }


if(!dir.exists(paste0(getwd(), "/full_data/", Trial, "/FSBLUP"))) { 
    dir.create(paste0(getwd(), "/full_data/", Trial, "/FSBLUP"), recursive = T) 
}

readr::write_csv(results, paste0(getwd(), "/full_data/", Trial, "/FSBLUP/FSBLUP_2x20.csv"))
