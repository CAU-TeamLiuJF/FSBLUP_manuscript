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
for (s in 3) {
    
    for (r in 1:20) {
        cat("\nScenario: ", s, "Rep: ", r, "\n total: ", r,"/20\n")
        #s=3;r=1

        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/scenario_", s,"/rep_", r, "/data.rdata"))

        ## A and G
        A <- amat
        G <- gmat
        id_1 <- colnames(A)
        id_2 <- colnames(G)
        id_1 <- id_1[!id_1 %in% id_2]
    
        A11 <- A[id_1, id_1]
        A12 <- A[id_1, id_2]
        A21 <- A[id_2, id_1]
        A22 <- A[id_2, id_2]
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
        a = (avg_sum_a - avg_diag_a) / (avg_sum_g - avg_diag_g)
        b = avg_diag_a - a * avg_diag_g
        G = a * G + b
        G = 0.95 * G + 0.05 * A22
            
        G <- b * G + (1 - b) * A22
        H11 <- A11 + A12 %*% iA22 %*% (G - A22) %*% iA22 %*% A21
        H12 <- A12 %*% iA22 %*% G
        H21 <- G %*% iA22 %*% A21
        H <- cbind(rbind(H11, H21), rbind(H12, G))
        H <- H / mean(diag(H))
        rm(A11, A12, A21, A22, iA22, G, H11, H12, H21, B, A, avg_sum_g, avg_sum_a, avg_diag_g, avg_diag_a, sw, sm, nind)
        

        G <- H[rownames(tmat), colnames(tmat)] %>% as.matrix()


        n.feat <- 1200

        M <- fread("data/M.txt", data.table = F) %>% 
            filter(row_number() %in% colnames(G)) %>%
            as.matrix()

        M <- t((t(M)-colMeans(M))/apply(M,2,sd))/sqrt(n.feat)
        

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

        ### computing G.hat:


        LHS11 <- 400
        LHS12 <- t(rep(1,400))
        LHS22 <- diag(400) + zeta*solve(G)

        LHS <- rbind(cbind(LHS11, LHS12), cbind(t(LHS12), LHS22))

        rm(LHS11,LHS12,LHS22)
        gc()

        G.hat <- solve(LHS, rbind(colSums(M), M))[-1,]

        rm(LHS)
        gc()

        E.hat <- M - G.hat
        rm(M)
        gc()

        Q <- rbind(cbind(tcrossprod(G.hat), tcrossprod(G.hat,E.hat)),cbind(tcrossprod(E.hat,G.hat), tcrossprod(E.hat))) 

        rm(G.hat, E.hat)
        gc()

        Qinv <- solve(Q+diag(nrow(Q))*0.0001*mean(diag(Q)))
        ## since Q is singular
        rm(Q)
        gc()


        Ginv <- H %>% as.matrix() %>% solve()

        Qext.inv.HH <- Ginv/h2m 
        Qext.inv.HH[-(1:1710),-(1:1710)] <- Ginv[-(1:1710),-(1:1710)]/h2m + Qinv[1:400,1:400] - solve(G)/h2m 
        rm(G)
        gc()

        Qext.inv.HI <- matrix(0, 2110,400)
        Qext.inv.HI[-(1:1710),] <- Qinv[1:400,-(1:400)]

        Qext.inv.II <- Qinv[-(1:400),-(1:400)]

        rm(Qinv)
        gc()

        Qext.inv <- rbind(cbind(Qext.inv.HH, Qext.inv.HI),cbind(t(Qext.inv.HI), Qext.inv.II))

        rm(Qext.inv.HH, Qext.inv.HI, Qext.inv.II)
        gc()

        r.no <- var_epsilon/(var_epsilon+ (c2m*s2y)*(1-h2m))
        ## since the residual variance for individuals without omics is different from the ones with omics 

        LHS11 <- matrix(c(1710*r.no,0,0,200),2,2)
        XZ <- rbind(t(c(rep(1,1710)*r.no,rep(0,400))),t(c(rep(0,1710),rep(1,200),rep(0,200))))
        LHS12 <- cbind(XZ,XZ[,-(1:1710)])
        LHS13 <- XZ
        rm(XZ)
        ZZ <- cbind(rbind(diag(1910), matrix(0,200,1910)), matrix(0,2110,200))
        ZZ[(1:1710),(1:1710)] <- ZZ[(1:1710),(1:1710)]*r.no
        LHS22 <- rbind(cbind(ZZ, ZZ[,-(1:1710)]), cbind(ZZ[-(1:1710),], ZZ[-(1:1710),-(1:1710)])) + eta1*Qext.inv
        LHS23 <- rbind(ZZ,ZZ[-(1:1710),])
        LHS33 <- ZZ + eta2*Ginv

        rm(ZZ, Ginv, Qext.inv)
        gc()


        LHS <- rbind(cbind(LHS11, LHS12, LHS13), cbind(t(LHS12), LHS22, LHS23), cbind(t(LHS13),t(LHS23), LHS33))

        rm(LHS11, LHS12, LHS13, LHS22, LHS23, LHS33)
        gc()

        RHS <- c(
            sum(yy[1:1710])*r.no,sum(yy[1711:1910]),c(yy[1:1710]*r.no,yy[1711:1910],rep(0,200)),
            c(yy[1711:1910],rep(0,200)),c(yy[1:1710]*r.no,yy[1711:1910],rep(0,200))
        )

        SOL <- solve(LHS, RHS)

        hat.a.r <- SOL[2+2110+400+(1:2110)]
        hat.a.m <- SOL[2+(1:2110)]


        EBV1sum <- hat.a.m + hat.a.r

        tbv <- pheno %>% pull(tbv)

        ebv <- cbind(id=1:2110,hat.a.m,hat.a.r,EBV1sum, tbv)
    

        # result
        results <- rbind(
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "GOBLUP",
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

write_csv(result, paste0(wd, "/result/goblup_scenario3.csv"))
