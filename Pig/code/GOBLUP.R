wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rrBLUP))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lme4qtl))


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

Traits <- c("age", "bf", "czs", "chzs")
cvnum = 5

times <- data.frame()
results <- data.frame()
start_time <- Sys.time()
for (tr in 1:4) {

    start_time1 <- Sys.time()

    # tr <- 1
    load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata"))

    pheno_all <- pheno %>% 
        dplyr::select(ID = id, y = Traits[tr]) %>% 
        dplyr::filter(ID %in% colnames(GG)) %>%
        dplyr::mutate(ID = as.character(ID), ID2 = ID)

    ## HTP and G
    A <- GG
    G <- TT
    id_1 <- colnames(A)
    id_2 <- colnames(G)
    id_1 <- id_1[!id_1 %in% id_2]
    A11 <- A[id_1, id_1]
    A12 <- A[id_1, id_2]
    A21 <- A[id_2, id_1]
    A22 <- A[id_2, id_2]
    if(!is.null(A22)) {
        diag(A22) <- diag(A22) + 1e-06
        B <- try(chol(A22), silent = TRUE)
        if (inherits(B, what = "try-error")) {
            A22 <- Matrix::nearPD(A22)$mat %>% as.matrix()
            rm(B)
        }
    }
    iA22 <- solve(A22)
    G <- G[colnames(A22), colnames(A22)]
    ## make the mathematical expectation of two matrix equal
    #LHS <- matrix(c(mean(diag(G)), 1, mean(diag(G[, ncol(G):1])), 1), nrow = 2)
    #RHS <- c(mean(diag(A22)), mean(diag(A22[, ncol(A22):1])))
    #G <- solve(LHS, RHS)[1] * G + solve(LHS, RHS)[2]
    nind = nrow(G)
    avg_sum_g = sum(G) / (nind * nind)
    avg_sum_a = sum(A22) / (nind * nind)
    avg_diag_g = sum(diag(G)) / nind
    avg_diag_a = sum(diag(A22)) / nind
    sw = (avg_sum_a - avg_diag_a) / (avg_sum_g - avg_diag_g)
    sm = avg_diag_a - sw * avg_diag_g
    G = sw * G + sm
    G = 0.95 * G + 0.05 * A22
    H11 <- A11 + A12 %*% iA22 %*% (G - A22) %*% iA22 %*% A21
    H12 <- A12 %*% iA22 %*% G
    H21 <- G %*% iA22 %*% A21
    H <- cbind(rbind(H11, H21), rbind(H12, G))
    H <- H / mean(diag(H))
    H <- H[pheno_all$ID, pheno_all$ID]
    G <- GG[pheno_all$ID, pheno_all$ID]



    for (i in 1:20) {
        start_time2 <- Sys.time()
        for (j in 1:5) {
            cat("trait", colnames(pheno)[tr + 1], " | cv_seed", i, " | cvnum", j, "\n")
            # i=1; j=1

            set.seed(i)
            pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1/cvnum), times = cvnum)))
            pheno_all$yNA <- pheno_all$y
            nas = pheno_all$partition == j
            pheno_all$yNA[nas] <- NA


            hoo <- H[!nas, !nas]
            hno <- H[nas, !nas]
            goo <- G[!nas, !nas]
            gno <- G[nas, !nas]


            res_gt <- lme4qtl::relmatLmer(yNA ~ (1 | ID) + (1 | ID2), data = pheno_all, relmat = list(ID = H, ID2 = G))
            u1 <- Matrix::as.matrix(t(as.matrix(res_gt@optinfo$relmat$relfac$ID)) %*% Matrix::as.matrix(ranef(res_gt)[[1]]))[rownames(hoo), 1]
            u2 <- Matrix::as.matrix(t(as.matrix(res_gt@optinfo$relmat$relfac$ID2)) %*% Matrix::as.matrix(ranef(res_gt)[[2]]))[rownames(goo), 1]

            res_k <- mixed.solve(c(u1, hno %*% MASS::ginv(hoo) %*% u1), K = G[c(colnames(hno), rownames(hno)), c(colnames(hno), rownames(hno))])
            pred_lme4qtl_gt <- c(res_k$u[rownames(hno)]) + c(gno %*% MASS::ginv(goo) %*% u2)

            # result
            results <- rbind(
                results,
                data.frame(
                    trait = colnames(pheno)[tr + 1],
                    cv_seed = i,
                    cvnum = j,
                    method = "GOBLUP",
                    pearson = cor(pheno_all$y[nas], pred_lme4qtl_gt, use = 'pairwise.complete.obs'),
                    bias = lm(pheno_all$y[nas] ~ pred_lme4qtl_gt)$coefficients[2]
                )
            )

        }
        end_time <- Sys.time()
        times <- rbind(times, data.frame(
            trait = colnames(pheno)[tr + 1], 
            method = "GOBLUP", 
            rep = i, 
            trait_start_time = start_time1,
            start_time = start_time2, 
            end_time = end_time, 
            spend_time = difftime(end_time, start_time2, units = "secs")
        ))
    }

    readr::write_csv(results, paste0(getwd(), "/results/", colnames(pheno)[tr + 1], "/GOBLUP_5x20.csv"))
    readr::write_csv(times, paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1], "/GOBLUP_5x20_times.csv"))
}

cat("GOBLUP done !!! spend time:", difftime(end_time, start_time, units = "secs"), "\n")