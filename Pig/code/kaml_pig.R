setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_revision")
library(tidyverse)
library(cli)
library(KAML)
library(data.table)

arg <- as.numeric(commandArgs(t=T)[1])
if(is.na(arg)) arg = 1

result = expand_grid(t = c("age", "bf", "czs"), r = 1:20)

result = result[arg, ]

cli_h1("KAML | Pig")
cli_h2("{result$t} | r{result$r}")

results <- data.frame()

start_time <- Sys.time()

load("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/data_adj.rdata")

# 去掉geno列名最后两个字符 (e.g. _A, _T)
colnames(geno) <- c(colnames(geno)[1], str_sub(colnames(geno)[-1], end = -3))
X <- geno %>% column_to_rownames("id") %>% as.matrix()
X = X[,colMeans(!is.na(X)) >  0.5 & (0.5-abs(0.5-colMeans(X,na.rm = T)) > 0.05)]
X[is.na(X)] = matrix(colMeans(X,na.rm=T),nr = nrow(X),nc = ncol(X),byrow=T)[is.na(X)]
X = X[as.character(pheno$id), -1]

#
#file.copy("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/data/pig.map.txt", 
#	"/public/home/liujf/workspace/xueyh/TempWork/h_matrix_revision/runs/kaml/pig/map.txt")


wkp <- paste0(getwd(), "/runs/kaml/pig/", result$t, "/rep_", result$r)

if(!dir.exists(wkp)) { dir.create(wkp, recursive = T) }

setwd(wkp)
	  
fwrite(t(X), "geno.txt",quote = F, sep = "\t", col.names = F, row.names = F)
file.copy("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_revision/runs/kaml/pig/map.txt", wkp)
KAML.Data(numfile="geno.txt", mapfile="map.txt", out="pig", sep = "\t")

for (r in 1:5) {
	set.seed(result$r)
    phenos <- pheno %>% 
		dplyr::select(ID = 1, y = result$t) %>% 
		dplyr::filter(ID %in% colnames(AA)) %>%
		dplyr::mutate(
			ID = as.character(ID), 
			ID2 = ID,
			partition = sample(seq_len(5), size = nrow(pheno), replace = TRUE, prob = c(rep((1 / 5), times = 5))),
			yNA = y
		) %>% 
		as.data.frame()

    nas <- phenos$partition == r
    phenos$yNA[nas] <- NA

	fwrite(phenos, "pheno.txt", quote = F, sep = "\t", col.names = T, row.names = F, na = "NA")
	
	kaml <- KAML(pfile = "pheno.txt", pheno = 5, gfile = paste0("pig"), cpu = 1)
	
	results <- rbind(
            results,
            data.frame(
				Trait = result$t, 
                Rep = result$r,
				cvnum = r,
                Method = "KAML",
                pearson = cor(phenos$y[nas], kaml$gebv[nas], use = "pairwise.complete.obs"),
                bias = lm(phenos$y[nas] ~ kaml$gebv[nas])$coefficients[2],
				time = difftime(Sys.time(), start_time, units = "secs")
            )
        )
	 cli_alert_info("Progress: {r} /5")
}

write_csv(results, paste0("kaml_", arg, ".csv"))
		
cli_alert_success("Task < {result$t} | r{result$r} > done !")

setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_revision")
presult = expand_grid(t = c("age", "bf", "czs"), r = 1:20) %>%
	mutate(rn = row_number()) %>%
	rowwise %>%
	mutate(
		pth = paste0(getwd(), "/runs/kaml/pig/", t, "/rep_", r, "/kaml_", rn, ".csv"),
		#file = fs::file_exists(pth)
		cvr = list(read_csv(pth, show_col_types = F)),
		pth = NULL
	) %>%
	ungroup %>%
	unnest(cvr) %>%
	group_by(t) %>%
	summarise(cor = mean(pearson), cor_sd = sd(pearson)) %>%
	ungroup
	