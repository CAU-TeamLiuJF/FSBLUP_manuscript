setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_revision")
library(tidyverse)
library(cli)
library(KAML)
library(data.table)


arg <- as.numeric(commandArgs(t=T)[1])
if(is.na(arg)) arg = 1

result = expand_grid(t = c("t1", "t2"), s = 1:3, r = 1:20)

result = result[arg, ]

cli_h2("{result$t} | s{result$s} | r{result$r}")

results <- data.frame()

map <- fread(paste0(getwd(), "/simulation/", ifelse(result$t == "t1", "c2021", "c2"), "/mrk_map.txt"))

wkp <- paste0(getwd(), "/simulation/", result$t, "/scenario_", result$s, "/rep_", result$r)
load(file.path(wkp, "data.rdata"))

nas <- pheno$partition == "test"

pheno$yNA[nas] <- NA  

setwd(wkp)
	  
fwrite(t(geno), "geno.txt", quote = F, sep = "\t", col.names = F, row.names = F)
fwrite(pheno, "pheno.txt", quote = F, sep = "\t", col.names = T, row.names = F, na = "NA")
fwrite(map, "map.txt", quote = F, sep = "\t", col.names = F, row.names = F)

start_time <- Sys.time()

KAML.Data(numfile="geno.txt", mapfile="map.txt", out="t1", sep = "\t")

kaml <- KAML(pfile = "pheno.txt", pheno = 5, gfile = paste0("t1"), cpu = 1)

results <- rbind(
            results,
            data.frame(
				Trait = result$t,
                Scenario = result$s, 
                Rep = result$r,
                Method = "KAML",
                pearson = cor(pheno$tbv[nas], kaml$gebv[nas], use = "pairwise.complete.obs"),
                bias = lm(pheno$tbv[nas] ~ kaml$gebv[nas])$coefficients[2],
				time = difftime(Sys.time(), start_time, units = "secs")
            )
        )

write_csv(results, paste0("kaml_", arg, ".csv"))
		
cli_alert_success("Task < {result$t} | s{result$s} | r{result$r} > done !")
