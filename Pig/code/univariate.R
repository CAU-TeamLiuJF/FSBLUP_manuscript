wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rrBLUP))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lme4qtl))
suppressPackageStartupMessages(library(BGLR))
suppressPackageStartupMessages(library(FSBLUP))

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

tr = commandArgs(trailingOnly = TRUE)[1]
print(paste0("Pig || Trait: ",tr))
tr = as.numeric(tr)

for(tr in 1:4) {

    samnum <- 20 ##  repeat times
    cvnum <- 5 ##  cross validation times

    times <- data.frame() ## record time
    start_time <- Sys.time()
    load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata"))


    pheno_all <- pheno %>% 
        dplyr::select(ID = 1, y = (tr + 1)) %>% 
        dplyr::filter(ID %in% colnames(AA)) %>%
        dplyr::mutate(ID = as.character(ID), ID2 = ID) %>% 
        as.data.frame()

    GG <- GG[pheno_all$ID, pheno_all$ID]

    results <- data.frame()
    for (cv_s in seq_len(samnum)) {
        start_time1 <- Sys.time()
        for (j in seq_len(cvnum)) {
            #cv_s=1; j=1
            set.seed(cv_s)
            pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
            pheno_all$yNA <- pheno_all$y
            nas <- pheno_all$partition == j
            pheno_all$yNA[nas] <- NA


            ## rrBLUP - G
            res_k <- mixed.solve(pheno_all$yNA, K = GG)
            res_k <- res_k$u[nas]

            # result
            results <- rbind(
                results,
                data.frame(
                    trait = colnames(pheno)[tr + 1],
                    cv_seed = cv_s,
                    cvnum = j,
                    Method = "GBLUP",
                    pearson = cor(pheno_all$y[nas], res_k, use = 'pairwise.complete.obs'),
                    bias = lm(pheno_all$y[nas] ~ res_k)$coefficients[2]
                )
            )
        }
        end_time <- Sys.time()
        times <- rbind(times, data.frame(
            trait = colnames(pheno)[tr + 1], 
            method = "GBLUP", 
            rep = cv_s, 
            start_time = start_time1, 
            end_time = end_time, 
            spend_time = difftime(end_time, start_time1, units = "secs")
        ))
        cat(paste0("\n\n [ Pig | ", colnames(pheno)[tr + 1]," | GBLUP | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
        cat(paste0(" [ Pig | ", colnames(pheno)[tr + 1]," | GBLUP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
    }
    if(!dir.exists(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]), recursive = T) }
    readr::write_csv(results, paste0(getwd(), "/results/", colnames(pheno)[tr + 1], "/GBLUP_5x20.csv"))
    if(!dir.exists(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]),recursive = T) }
    readr::write_csv(times, paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1], "/GBLUP_5x20_times.csv"))

}

 ABLUP

for(tr in 1:4) {
    # tr=1
    samnum <- 20 ##  repeat times
    cvnum <- 5 ##  cross validation times

    times <- data.frame() ## record time
    start_time <- Sys.time()
    load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata"))


    pheno_all <- pheno %>% 
        dplyr::select(ID = 1, y = (tr + 1)) %>% 
        dplyr::filter(ID %in% colnames(AA)) %>%
        dplyr::mutate(ID = as.character(ID), ID2 = ID) %>% 
        as.data.frame()

    AA <- AA[pheno_all$ID, pheno_all$ID]

    results <- data.frame()
    for (cv_s in seq_len(samnum)) {
        start_time1 <- Sys.time()
        for (j in seq_len(cvnum)) {
            #cv_s=1; j=1
            set.seed(cv_s)
            pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
            pheno_all$yNA <- pheno_all$y
            nas <- pheno_all$partition == j
            pheno_all$yNA[nas] <- NA


            ## rrBLUP - G
            res_k <- mixed.solve(pheno_all$yNA, K = AA)
            res_k <- res_k$u[nas]

            # result
            results <- rbind(
                results,
                data.frame(
                    trait = colnames(pheno)[tr + 1],
                    cv_seed = cv_s,
                    cvnum = j,
                    Method = "ABLUP",
                    pearson = cor(pheno_all$y[nas], res_k, use = 'pairwise.complete.obs'),
                    bias = lm(pheno_all$y[nas] ~ res_k)$coefficients[2]
                )
            )
        }
        end_time <- Sys.time()
        times <- rbind(times, data.frame(
            trait = colnames(pheno)[tr + 1], 
            method = "ABLUP", 
            rep = cv_s, 
            start_time = start_time1, 
            end_time = end_time, 
            spend_time = difftime(end_time, start_time1, units = "secs")
        ))
        cat(paste0("\n\n [ Pig | ", colnames(pheno)[tr + 1]," | ABLUP | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
        cat(paste0(" [ Pig | ", colnames(pheno)[tr + 1]," | ABLUP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
    }
    if(!dir.exists(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]), recursive = T) }
    readr::write_csv(results, paste0(getwd(), "/results/", colnames(pheno)[tr + 1], "/ABLUP_5x20.csv"))
    if(!dir.exists(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]),recursive = T) }
    readr::write_csv(times, paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1], "/ABLUP_5x20_times.csv"))

}



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
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata"))
pheno_all <- pheno %>% 
        dplyr::select(ID = 1, y = 2) %>% 
        dplyr::filter(ID %in% colnames(AA)) %>%
        dplyr::mutate(ID = as.character(ID), ID2 = ID) %>% 
        as.data.frame()
X <- geno %>% filter(as.character(id) %in% pheno_all$ID) %>% column_to_rownames("id") %>% as.matrix()

X = X[,colMeans(!is.na(X)) >  0.5 & (0.5-abs(0.5-colMeans(X,na.rm = T)) > 0.05)]
X[is.na(X)] = matrix(colMeans(X,na.rm=T),nr = nrow(X),nc = ncol(X),byrow=T)[is.na(X)]
X = X[pheno_all$ID, ]

DD <- as.matrix(dist(X))

for(tr in 1:4) {

    samnum <- 20 ##  repeat times
    cvnum <- 5 ##  cross validation times

    times <- data.frame() ## record time
    start_time <- Sys.time()
    load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata"))


    pheno_all <- pheno %>% 
        dplyr::select(ID = 1, y = (tr + 1)) %>% 
        dplyr::filter(ID %in% colnames(AA)) %>%
        dplyr::mutate(ID = as.character(ID), ID2 = ID) %>% 
        as.data.frame()

    D <- DD[pheno_all$ID, pheno_all$ID]

    results <- data.frame()
    for (cv_s in seq_len(samnum)) {
        start_time1 <- Sys.time()
        for (j in seq_len(cvnum)) {
            #cv_s=1; j=1
            set.seed(cv_s)
            pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
            pheno_all$yNA <- pheno_all$y
            nas <- pheno_all$partition == j
            pheno_all$yNA[nas] <- NA


            ## rrBLUP - G
            res_k <- kin.blup(pheno_all, geno = "ID", pheno = "yNA", K = D, GAUSS = T)$pred
            res_k <- res_k[nas]

            # result
            results <- rbind(
                results,
                data.frame(
                    trait = colnames(pheno)[tr + 1],
                    cv_seed = cv_s,
                    cvnum = j,
                    Method = "RKHS",
                    pearson = cor(pheno_all$y[nas], res_k, use = 'pairwise.complete.obs'),
                    bias = lm(pheno_all$y[nas] ~ res_k)$coefficients[2]
                )
            )
        }
        end_time <- Sys.time()
        times <- rbind(times, data.frame(
            trait = colnames(pheno)[tr + 1], 
            method = "RKHS", 
            rep = cv_s, 
            start_time = start_time1, 
            end_time = end_time, 
            spend_time = difftime(end_time, start_time1, units = "secs")
        ))
        cat(paste0("\n\n [ Pig | ", colnames(pheno)[tr + 1]," | RKHS | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
        cat(paste0(" [ Pig | ", colnames(pheno)[tr + 1]," | RKHS | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
    }
    if(!dir.exists(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]), recursive = T) }
    readr::write_csv(results, paste0(getwd(), "/results/", colnames(pheno)[tr + 1], "/RKHS_5x20.csv"))
    if(!dir.exists(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]),recursive = T) }
    readr::write_csv(times, paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1], "/RKHS_5x20_times.csv"))

}




## BL
times <- data.frame() ## record time
start_time <- Sys.time()
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata"))


pheno_all <- pheno %>% 
    dplyr::select(ID = 1, y = (tr + 1)) %>% 
    dplyr::filter(ID %in% colnames(AA)) %>%
    dplyr::mutate(ID = as.character(ID), ID2 = ID) %>% 
    as.data.frame()

GG <- GG[pheno_all$ID, pheno_all$ID]

X <- geno %>% filter(as.character(id) %in% pheno_all$ID) %>% column_to_rownames("id") %>% as.matrix()

X = X[,colMeans(!is.na(X)) >  0.5 & (0.5-abs(0.5-colMeans(X,na.rm = T)) > 0.05)]
X[is.na(X)] = matrix(colMeans(X,na.rm=T),nr = nrow(X),nc = ncol(X),byrow=T)[is.na(X)]
X = X[pheno_all$ID, ]
#G_c <- G
results <- data.frame()
start_time <- Sys.time()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
        #cv_s = 1; j = 1
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$y
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA

        #Knn = G_c[nas,nas]
        #sKnn = svd(Knn)

        ##BGLR - BL - default priors
        suppressWarnings(try(dir.create(paste0('BGLR_dir/', colnames(pheno)[tr + 1]))))
        bglr_BL = BGLR(y = pheno_all$yNA,ETA = list(list(X = X,model = 'BL')),burnIn = 5000,nIter = 10000,verbose=F,
            saveAt = paste0('BGLR_dir/',colnames(pheno)[tr + 1],'/fold_', cv_s,'_', j))
        bglr_BLs = X[nas,] %*% bglr_BL$ETA[[1]]$b

        # result
        results <- rbind(
            results,
            data.frame(
                trait = colnames(pheno)[tr + 1],
                cv_seed = cv_s,
                cvnum = j,
                Method = "BL",
                pearson = cor(pheno_all$y[nas], bglr_BLs, use = 'pairwise.complete.obs'),
                bias = lm(pheno_all$y[nas] ~ bglr_BLs)$coefficients[2]
                #g_cor = estimate_gcor(data.frame(ID=pheno_all$ID[nas],obs = pheno_all$BLUP[nas],pred = bglr_BLs),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]
            )
        )
    }
    end_time <- Sys.time()
    times <- rbind(times, data.frame(
        trait = colnames(pheno)[tr + 1], 
        method = "BL", 
        rep = cv_s, 
        start_time = start_time1, 
        end_time = end_time, 
        spend_time = difftime(end_time, start_time1, units = "secs")
    ))
    cat(paste0("\n\n [ Pig | ", colnames(pheno)[tr + 1]," | BL | ", cv_s, "/", samnum,"] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
    cat(paste0(" [ Pig | ", colnames(pheno)[tr + 1]," | BL | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
    if(!dir.exists(paste0(getwd(), "/process/", colnames(pheno)[tr + 1], "/BL"))) { dir.create(paste0(getwd(), "/process/", colnames(pheno)[tr + 1], "/BL"), recursive = T) }
    readr::write_csv(results, paste0(getwd(), "/process/", colnames(pheno)[tr + 1], "/BL/BL_5x20_rep_", cv_s, ".csv"))
    #readr::write_csv(times, paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1], "/BL/BL_5x20_times_rep_", cv_s, ".csv"))
}
if(!dir.exists(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/results/", colnames(pheno)[tr + 1]), recursive = T) }
readr::write_csv(results, paste0(getwd(), "/results/", colnames(pheno)[tr + 1], "/BL_5x20.csv"))
if(!dir.exists(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]))) { dir.create(paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1]),recursive = T) }
readr::write_csv(times, paste0(getwd(), "/test_time/", colnames(pheno)[tr + 1], "/BL_5x20_times.csv"))


## FSBLUP

Traits <- c("age", "bf", "czs")

times <- data.frame() ## record time
for (Trait in Traits) {
    load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata"))

    start_time <- Sys.time()

    pheno_all <- pheno %>%
    dplyr::select(ID = id, y = Trait) %>%
    dplyr::filter(ID %in% colnames(AA)) %>%
    dplyr::mutate(ID = as.character(ID)) %>%
    as.data.frame()

    G <- GG[pheno_all$ID, pheno_all$ID]
    A <- AA[pheno_all$ID, pheno_all$ID]

    # fs matrix
    fs <- FSBLUP(
        phe = pheno_all, trait_col = 2, M1 = A, M2 = G, M3 = TT,
        po.crv.num = 5, po.crv.rep.num = 2, po.gs.point_num = 25,
        po.bi.max_iter = 10, po.bi.threshold = 1e-4, stas.type = "cor", return.matrix = T
    )

    gc()

    results <- data.frame()
    for (cv_s in seq_len(samnum)) {
        start_time1 <- Sys.time()
        # cv_s=1;j=1;Trait="age"
        for (j in seq_len(cvnum)) {
            set.seed(cv_s)
            # set.seed(1)
            pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
            pheno_all$yNA <- pheno_all$y
            nas <- pheno_all$partition == j
            pheno_all$yNA[nas] <- NA

            ## rrBLUP
            res_k <- tryCatch(mixed.solve(pheno_all$yNA, K = fs), error = {fs=adj_pos(fs);mixed.solve(pheno_all$yNA, K = fs)})
            res_k <- res_k$u[nas]

            # result
            results <- rbind(
                results,
                data.frame(
                    trait = Trait,
                    Method = "FSBLUP",
                    cv_seed = cv_s,
                    cvnum = j,
                    pearson = cor(pheno_all$y[nas], res_k, use = "pairwise.complete.obs"),
                    bias = lm(pheno_all$y[nas] ~ res_k)$coefficients[2]
                )
            )
        }

        end_time <- Sys.time()
        times <- rbind(times, data.frame(
            Trait = Trait, method = "FSBLUP", rep = cv_s, start_time = start_time1, end_time = end_time,
            spend_time = difftime(end_time, start_time1, units = "secs")
        ))

        cat(paste0("\n\n [", Trait, " | FSBLUP | ", cv_s, "/20] spend_time: ", difftime(end_time, start_time1, units = "secs"), " \n"))
        cat(paste0(" [", Trait, " | FSBLUP | TotalConsuming ] spend_time: ", difftime(end_time, start_time, units = "secs"), " \n\n"))
    }


    readr::write_csv(results, paste0(wd, "/results/", Trait, "/FSBLUP_5x20_fsmat.csv"))
}
