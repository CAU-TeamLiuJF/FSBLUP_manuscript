wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
library(rrBLUP)
suppressPackageStartupMessages(library(lme4))
library(lme4qtl)

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
for (s in 1) {
    
    for (r in 1:20) {
        cat("\nScenario: ", s, "Rep: ", r, "\n total: ", (s-1)*20+r,"/20\n")
        #s=1;r=1

        #print(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/scenario_", s,"/rep_", r, "/data.rdata"))
        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/scenario_", s,"/rep_", r, "/data.rdata"))
        

        Ginv <- gmat %>% as.matrix() %>% solve()

        n.feat <- 1200

        M <- data.table::fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_2/M.txt", data.table = F) %>% 
            filter(row_number() %in% colnames(Ginv)) %>%
            as.matrix()

        M <- t((t(M)-colMeans(M))/apply(M,2,sd))/sqrt(n.feat)
        

        yy <- pheno %>%
            filter(ngen %in% 1:10) %>%
            pull(y)



        h2m <- 0.61
        h2 <- 0.47
        h2ar <- 0.20
        c2m <- (h2-h2ar)/h2m

        eta1 <- (1-c2m-h2ar)/c2m
        eta2 <- (1-c2m-h2ar)/h2ar

        zeta <- (1-h2m)/h2m

        ### computing G.hat:


        LHS11 <- 2110
        LHS12 <- t(rep(1,2110))
        LHS22 <- diag(2110) + zeta*Ginv

        LHS <- rbind(cbind(LHS11, LHS12), cbind(t(LHS12), LHS22))

        rm(LHS11,LHS12,LHS22)
        gc()

        G.hat <- solve(LHS, rbind(colSums(M), M))[-1,]

        rm(LHS)
        gc()


        LHS11 <- 1910
        LHS12 <- t(colSums(M[1:1910,]))
        LHS13 <- t(c(rep(1,1910),rep(0,200)))
        LHS22 <- crossprod(M[1:1910,]) + eta1*diag(n.feat)
        LHS23 <- cbind(t(M[1:1910,]),matrix(0,n.feat,200))
        LHS33 <- rbind(cbind(diag(1910),matrix(0,1910,200)),matrix(0,200,2110)) + eta2*Ginv

        LHS <- rbind(cbind(LHS11, LHS12, LHS13), cbind(t(LHS12), LHS22, LHS23),cbind(t(LHS13),t(LHS23), LHS33))

        rm(LHS11,LHS12,LHS13,LHS22,LHS23,LHS33)
        gc()

        SOL <- solve(LHS, c(sum(yy[1:1910]), t(t(M[1:1910,])%*%yy[1:1910]), t(c(yy[1:1910],rep(0,200)))))

        hat.a.m <- as.vector(G.hat%*%SOL[1+(1:n.feat)])
        hat.a.r <- SOL[-c(1,1+(1:n.feat))] 

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

write_csv(result, paste0(wd, "/result/goblup_scenario1.csv"))
