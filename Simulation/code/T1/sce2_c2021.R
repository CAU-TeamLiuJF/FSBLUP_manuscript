wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))


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


result <- data.frame()
for (s in 2) {
    
    for (r in 1:20) {
        cat("\nScenario: ", s, "Rep: ", r, "\n total: ", r,"/20\n")
        #s=2;r=1

        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/scenario_", s,"/rep_", r, "/data.rdata"))
        

        

        n.feat <- 1200

        M <- fread("data/M.txt", data.table = F) %>% 
            slice_tail(n = 400) %>%
            as.matrix()

        M <- t((t(M)-colMeans(M))/apply(M,2,sd))/sqrt(n.feat)
        noom <- 2110-400

        MMpsinv <- solve(tcrossprod(M)+diag(400)*0.0001)
        rm(M)
        gc()
        

        yy <- pheno %>%
            filter(ngen %in% 1:10) %>%
            pull(y)


        sigma2m <- 2
        h2m <- 0.61
        h2 <- 0.47
        h2ar <- 0.20
        c2m <- (h2-h2ar)/h2m
        s2y <- h2m*sigma2m/(h2-h2ar)
        var_ar <- h2ar*s2y
        var_epsilon <-(1-c2m-h2ar)*s2y

        eta1 <- var_epsilon/(c2m*s2y)
        ## eta1 also equals var_epsilon/sigma2m
        eta2 <- var_epsilon/var_ar

        #eta1 <- (1-c2m-h2ar)/c2m
        #eta2 <- (1-c2m-h2ar)/h2ar

        zeta <- (1-h2m)/h2m

        G <- gmat %>% as.matrix()
        GGzeta <- h2m*G
        diag(GGzeta) <- diag(GGzeta)+(1-h2m)


        Omega.inv <- solve(GGzeta)
        GGzeta22.inv <- solve(GGzeta[-(1:noom),-(1:noom)])
        rm(GGzeta)
        gc()

        Omega.inv[-(1:noom),-(1:noom)] <- Omega.inv[-(1:noom),-(1:noom)] + MMpsinv - GGzeta22.inv 
        rm(MMpsinv)
        rm(GGzeta22.inv)
        gc() 


        Ginv <- gmat %>% as.matrix() %>% solve()


        ## X contains separate intercept for those with/without omics

        RHS <- c(sum(yy[1:1710]),sum(yy[1711:1910]), yy[1:1910], rep(0,200), yy[1:1910], rep(0,200))

        LHS11 <- matrix(c(1710,0,0,200),2,2)
        LHS12 <- rbind(t(c(rep(1,1710),rep(0,400))),t(c(rep(0,1710),rep(1,200),rep(0,200))))
        ## LHS13 the same as LHS12 
        LHS22 <- cbind(rbind(diag(1910), matrix(0,200,1910)), matrix(0,2110,200)) + eta1*Omega.inv
        LHS23 <- cbind(rbind(diag(1910), matrix(0,200,1910)), matrix(0,2110,200))
        LHS33 <- cbind(rbind(diag(1910), matrix(0,200,1910)), matrix(0,2110,200)) + eta2*Ginv


        rm(Omega.inv)
        gc()

        LHS <- rbind(cbind(LHS11, LHS12, LHS12), cbind(t(LHS12), LHS22, LHS23), cbind(t(LHS12),t(LHS23), LHS33))

        rm(LHS11, LHS12, LHS22, LHS23, LHS33)
        gc()

        SOL1 <- solve(LHS, RHS)

        write(SOL1, file="SOL1")

        ## SOL1 is phenotype prediction from omics M%*%alpha


        ## now we aim for the breeding values

        RHS <- c(sum(SOL1[2+(1:2110)]), SOL1[2+(1:2110)])

        LHS11 <- 2110
        LHS12 <- t(rep(1,2110))
        LHS22 <- diag(2110) + zeta*Ginv

        LHS <- rbind(cbind(LHS11, LHS12), cbind(t(LHS12), LHS22))

        rm(LHS11,LHS12,LHS22)
        gc()

        SOL2 <- solve(LHS, RHS)

        write(SOL2, file="SOL2")


        hat.a.r <- SOL1[2+2110+(1:2110)]
        hat.a.m <- SOL2[2:2111]


        #write.table(cbind(hat.a.m,hat.a.r), file="EBVss2", quote=FALSE, row.names=FALSE)

        ## the last 2000 individuals of those are used for validation


        EBV1sum <- hat.a.m + hat.a.r

        tbv <- pheno %>% pull(tbv)

        ebv <- cbind(id=1:2110,hat.a.m,hat.a.r,EBV1sum, tbv)
    

        # result
        results <- rbind(
            data.frame(
                Scenario = s, 
                Rep = r,
                method = "GOBLUP",
                pearson = cor(EBV1sum[1911:2110], tbv[1911:2110], use = "pairwise.complete.obs"),
                bias = lm(tbv[1911:2110] ~ EBV1sum[1911:2110])$coefficients[2]
            )
        )
        result <- rbind(result, results)

    if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/goblup"))) {
        dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/goblup"), recursive = T)
    }   
        write.csv(ebv, paste0(wd, "/scenario_", s, "/rep_", r, "/goblup/goblup_ebv_scenario_", s, "_rep_", r,".csv"), row.names = F)
        write.csv(results, paste0(wd, "/scenario_", s, "/rep_", r, "/goblup/goblup_results_scenario_", s, "_rep_", r,".csv"), row.names = F)
    }
    
}

write_csv(result, paste0(wd, "/result/goblup_scenario2.csv"))
