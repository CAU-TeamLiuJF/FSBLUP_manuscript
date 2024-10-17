wd <- "/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/"
setwd(paste0(wd))

options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
library(BGLR)

geno <- data.table::fread("/public/home/liujf/workspace/xueyh/TempWork/simulate_c2021/r_datasim/p0_mrk_001.txt", data.table = F) %>% 
    column_to_rownames("ID") %>% 
    as.matrix()

geno[is.na(geno)] = matrix(colMeans(geno,na.rm=T),nr = nrow(geno),nc = ncol(geno),byrow=T)[is.na(geno)]

times <- data.frame() ##  record time
#start_time <- Sys.time()
#for (s in 1) {
#    
#    for (r in 1:20) {
#		start_time1 <- Sys.time()
#        results <- data.frame()
#        #print(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/scenario_", s,"/rep_", r, "/data.rdata"))
#        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/scenario_", s,"/rep_", r, "/data.rdata"))
#
#        pheno_all <- pheno %>%
#        select(GID = 1, yc = 2, tbv, ngen) %>%
#        mutate(GID = as.character(GID), yc = as.numeric(yc), GID2 = GID) %>% 
#        as.data.frame()
#        
#        pheno_all <- pheno_all %>% mutate(partition = case_when(ngen %in% 1:10 ~ "train", TRUE ~ "test"))
#        pheno_all$yNA <- pheno_all$yc
#        nas <- pheno_all$partition == "test"
#        pheno_all$yNA[nas] <- NA
#
#		
#		##BGLR - BB - default priors
#		X <- geno[pheno_all$GID,]
#        bglr_BB = BGLR(y = pheno_all$yNA,ETA = list(list(X = X,model = 'BayesB')),burnIn = 5000,nIter = 50000,verbose=F)
#        bglr_BBs = X[nas,] %*% bglr_BB$ETA[[1]]$b
#        
#
#        # result
#        results <- rbind(
#            results,
#			data.frame(
#                Scenario = s, 
#                Rep = r,
#                Method = "BB",
#                pearson = cor(pheno_all$tbv[nas], bglr_BBs, use = 'pairwise.complete.obs'),
#                bias = lm(pheno_all$tbv[nas] ~ bglr_BBs)$coefficients[2]
#            )
#        )
#	end_time <- Sys.time()
#	print(paste0("\n\nBayesB calTime No.", r, " /20: ", difftime(end_time, start_time1, units = "secs")))
#	print(paste0("TotalTime:", difftime(end_time, start_time, units = "secs"), "\n\n"))
#			
#    if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/bayes"))) {
#        dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/bayes"), recursive = T)
#    }   
#        write_csv(results, paste0(wd, "/scenario_", s, "/rep_", r, "/bayes/bayesB_scenario_", s, "_rep_", r,".csv"))
#    }
#    
#}
#	times <- rbind(times, data.frame(Trait = "t2", method = "BayesB", start_time = start_time, end_time = end_time, 
#            spend_time = difftime(end_time, start_time, units = "secs")))


start_time <- Sys.time()
for (s in 1) {
    
    for (r in 1:20) {
		start_time1 <- Sys.time()
        results <- data.frame()
        #print(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_c2021/scenario_", s,"/rep_", r, "/data.rdata"))
        load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2/scenario_", s,"/rep_", r, "/data.rdata"))

        pheno_all <- pheno %>%
        select(GID = 1, yc = 2, tbv, ngen) %>%
        mutate(GID = as.character(GID), yc = as.numeric(yc), GID2 = GID) %>% 
        as.data.frame()
        
        pheno_all <- pheno_all %>% mutate(partition = case_when(ngen %in% 1:10 ~ "train", TRUE ~ "test"))
        pheno_all$yNA <- pheno_all$yc
        nas <- pheno_all$partition == "test"
        pheno_all$yNA[nas] <- NA


        if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/bayes"))) {
            dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/bayes"), recursive = T)
        }
        setwd(paste0(wd, "/scenario_", s, "/rep_", r, "/bayes"))
       
        ##BGLR - BC - default priors
        suppressWarnings(try(dir.create('BGLR_dir')))
        X <- geno[pheno_all$GID,]
        X <- X[pheno_all$GID, ]
        bglr_BC = BGLR(y = pheno_all$yNA,ETA = list(list(X = X,model = 'BL')),burnIn = 5000,nIter = 50000,verbose=F,saveAt = sprintf('BGLR_dir/fold_%d_%s', s, r))
        bglr_BCs = X[nas,] %*% bglr_BC$ETA[[1]]$b
		        

        # result
        results <- rbind(
            results,
            data.frame(
                Scenario = s, 
                Rep = r,
                Method = "BL",
                pearson = cor(pheno_all$tbv[nas], bglr_BCs, use = 'pairwise.complete.obs'),
                bias = lm(pheno_all$tbv[nas] ~ bglr_BCs)$coefficients[2]
            )
        )
	end_time <- Sys.time()
	print(paste0("\n\nBL calTime NO.", r, " /20: ", difftime(end_time, start_time1, units = "secs")))
	print(paste0("TotalTime:", difftime(end_time, start_time, units = "secs"), "\n\n"))
    if (!dir.exists(paste0(wd, "/scenario_", s, "/rep_", r, "/bayes"))) {
        dir.create(paste0(wd, "/scenario_", s, "/rep_", r, "/bayes"), recursive = T)
    }   
        write_csv(results, paste0(wd, "/scenario_", s, "/rep_", r, "/bayes/BL_scenario_", s, "_rep_", r,".csv"))
    }
    
}
	times <- rbind(times, data.frame(Trait = "t2", method = "BL", start_time = start_time, end_time = end_time, 
            spend_time = difftime(end_time, start_time, units = "secs")))
write.table(times, file = paste0("BL_method_spendTimes_", format(Sys.Date(),format = "%Y_%m_%d"), ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)			