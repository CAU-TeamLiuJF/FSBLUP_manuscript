wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rrBLUP))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lme4qtl))


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

Trials = c('MD2013', 'OB2013', 'H2013', 'SD2013', 'OF2013', 'MD2014', 'OB2014', 'H2014', 'SD2014',
    'OF2014', 'MD2015', 'OB2015', 'H2015', 'SD2015', 'OF2015', 'MD2016', 'OB2016', 'H2016', 'SD2016', 'OF2016')

samnum <- 20 ##  repeat times 
cvnum <- 2 ##  cross validation times
Trials <- 'OF2013'



for(Trial in Trials){
times <- data.frame() ## record time
start_time <- Sys.time()

#load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/full_data/", Trial, "/data_adj.rdata"))
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))

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

        goo <- G[!nas, !nas]
        gno <- G[nas, !nas]
        too <- HTP[!nas, !nas]
        tno <- HTP[nas, !nas]

        Knn = G_c[nas,nas]
        sKnn = svd(Knn)

        ## lme4 - G + HTP
        res_gt <- relmatLmer(BLUP ~ (1 | GID) + (1 | GID2), data = pheno_all[!nas, ], relmat = list(GID = G, GID2 = HTP))
        u1 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID) %*% as.matrix(ranef(res_gt)[[1]]))[rownames(goo), 1]
        u2 <- as.matrix(t(res_gt@optinfo$relmat$relfac$GID2) %*% as.matrix(ranef(res_gt)[[2]]))[rownames(too), 1]


        ## rrBLUP - G
        res_k <- mixed.solve(c(u2, tno %*% MASS::ginv(too) %*% u2), K = G[c(colnames(gno), rownames(gno)), c(colnames(gno), rownames(gno))])

        pred_lme4qtl_gt <- c(gno %*% MASS::ginv(goo) %*% u1) + res_k$u[rownames(gno)]


        # result
        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                method = "GOBLUP",
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = pred_lme4qtl_gt),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']],
                bias = lm(pheno_all$BLUP[nas] ~ pred_lme4qtl_gt)$coefficients[2]
            )
        )

    }
    end_time <- Sys.time()

    times <- rbind(times, data.frame(Trial = Trial, method = "GOBLUP", rep = cv_s, start_time = start_time1, end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")))
    cat(paste0("\n\n [", Trial, " | GOBLUP | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | GOBLUP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n"))

}

if(!dir.exists(paste0(getwd(), "/", Trial, "/results"))) { dir.create(paste0(getwd(), "/", Trial, "/results"),recursive = T) }
readr::write_csv(results, paste0(getwd(), "/", Trial, "/results/GOBLUP_2x5__rep_20_20241007_ori.csv"))
if(!dir.exists(paste0(getwd(), "/", Trial, "/test_time"))) { dir.create(paste0(getwd(), "/", Trial, "/test_time"),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/", Trial, "/test_time/GOBLUP_2x5__rep_20_times_20241007.csv"))

}




