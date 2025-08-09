setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_revision")
library(tidyverse)
library(cli)
library(KAML)
library(data.table)
source("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/data/Estimate_gcor_prediction.R")

arg <- as.numeric(commandArgs(t=T)[1])
if(is.na(arg)) arg = 1

# read genotype data
geno <- data.table::fread("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix/MegaLMM/data/Krause_et_al_2018_Genotypes.csv", data.table = F)
rownames(geno) = geno[,1]
geno = as.matrix(geno[,-1])

samnum <- 20 ##  repeat times 
cvnum <- 2 ##  cross validation times
Trial <- 'OF2013'

Trial <- c("SD2013", "MD2013", "OB2013", "H2013", "OF2013",
        "MD2014", "OB2014", "H2014", "SD2014", "OF2014",
        "MD2015", "OB2015", "H2015", "SD2015", "OF2015",
        "MD2016", "OB2016", "H2016", "SD2016", "OF2016")  #'OF2013'
Trial <- Trial[arg]		
cli_h2("KAML | {Trial}")		

times <- data.frame()
start_time <- Sys.time()
load(paste0("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat/percent_sample/", Trial, "/data.rdata"))

pheno_all <- wheat.pheno %>% 
    dplyr::select(GID, BLUE = 2, BLUP = 3) %>% 
    dplyr::filter(GID %in% colnames(A)) %>%
    dplyr::mutate(GID = as.character(GID), GID2 = GID) %>% 
    as.data.frame()

G <- G[pheno_all$GID, pheno_all$GID]
X = geno[pheno_all$GID,]
X = X[,colMeans(!is.na(X)) >  0.5 & (0.5-abs(0.5-colMeans(X,na.rm = T)) > 0.05)]
X[is.na(X)] = matrix(colMeans(X,na.rm=T),nr = nrow(X),nc = ncol(X),byrow=T)[is.na(X)]
G_c <- G

pth=paste0(getwd(), "/runs/kaml/wheat/", Trial)
if(!dir.exists(pth)) dir.create(pth, recursive = T)
setwd(pth)
fwrite(t(X * 2), "geno.txt",quote = F, sep = "\t", col.names = F, row.names = F)
map=tibble(V1 = colnames(X), 
	V2 =  str_remove(str_extract(V1, "S[^_]+"), "S"),
	V3 = str_extract(V1, "\\d+$")
	)
fwrite(map, "map.txt",quote = F, sep = "\t", col.names = F, row.names = F)
KAML.Data(numfile="geno.txt", mapfile="map.txt", out="wheat", sep = "\t")

results <- data.frame()
for (cv_s in seq_len(samnum)) {
    start_time1 <- Sys.time()
    for (j in seq_len(cvnum)) {
		#cv_s = 1; j = 1
        set.seed(cv_s)
        pheno_all$partition <- sample(seq_len(cvnum), size = nrow(pheno_all), replace = TRUE, prob = c(rep((1 / cvnum), times = cvnum)))
        pheno_all$yNA <- pheno_all$BLUE
        nas <- pheno_all$partition == j
        pheno_all$yNA[nas] <- NA

        Knn = G_c[nas,nas]
        sKnn = svd(Knn)

		fwrite(pheno_all, "pheno.txt", quote = F, sep = "\t", col.names = T, row.names = F, na = "NA")
	
		kaml <- KAML(pfile = "pheno.txt", pheno = 5, gfile = "wheat", cpu = 1)
        #bglr_BL = BGLR(y = pheno_all$yNA,ETA = list(list(X = X,model = 'BL')),burnIn = 5000,nIter = 50000,verbose=F)
        #bglr_BLs = X[nas,] %*% bglr_BL$ETA[[1]]$b
		
        # result
        results <- rbind(
            results,
            data.frame(
                cv_seed = cv_s,
                cvnum = j,
                Method = "kaml",
                pearson_blup = cor(pheno_all$BLUP[nas], kaml$gebv[nas], use = 'pairwise.complete.obs'),
                g_cor = estimate_gcor(data.frame(ID=pheno_all$GID[nas],obs = pheno_all$BLUP[nas],pred = kaml$gebv[nas]),Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']],
				time =  difftime(Sys.time(), start_time, units = "secs")
            )
        )
        
    }
    
	cli_alert_info("Progress: {cv_s} /{samnum} | {j} /{cvnum}")
    cat(paste0("\n\n [", Trial, " | KAML | ", cv_s, "/", samnum,"] spend_time: ", difftime(Sys.time(), start_time1, units = "secs"), " \n"))
    cat(paste0(" [", Trial, " | KAML | TotalConsuming ] spend_time: ", difftime(Sys.time(), start_time, units = "secs"), " \n\n"))
}

write_csv(results, paste0(pth, "/kaml_2*5.csv"))



#
#setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_revision")
Trial <- c("SD2013", "MD2013", "OB2013", "H2013", "OF2013",
        "MD2014", "OB2014", "H2014", "SD2014", "OF2014",
        "MD2015", "OB2015", "H2015", "SD2015", "OF2015",
        "MD2016", "OB2016", "H2016", "SD2016", "OF2016")
wresult = expand_grid(t = Trial) %>%
	rowwise %>%
	mutate(
		pth = paste0(getwd(), "/runs/kaml/wheat/", t, "/kaml_2*5.csv"),
		#file = fs::file_exists(pth)
		cvr = list(read_csv(pth, show_col_types = F)),
		pth = NULL
	) %>%
	ungroup %>%
	unnest(cvr) %>%
	group_by(t) %>%
	summarise(cor = mean(g_cor), cor_sd = sd(g_cor), time = sum(time)) %>%
	ungroup
	