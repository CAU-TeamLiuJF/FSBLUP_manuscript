setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix")
library(tidyverse)
library(ggsci)
library(showtext)
font_add("Arial", "arial.ttf")
showtext::showtext_auto()
library(Cairo)
library(patchwork)


Trials <- tibble(Trials = c("2013-14_Moderate Drought", "2013-14_Optimal Bed", "2013-14_Heat", "2013-14_Severe Drought", "2013-14_Optimal Flat",
                           "2014-15_Moderate Drought", "2014-15_Optimal Bed", "2014-15_Heat", "2014-15_Severe Drought", "2014-15_Optimal Flat",
                           "2015-16_Moderate Drought", "2015-16_Optimal Bed", "2015-16_Heat", "2015-16_Severe Drought", "2015-16_Optimal Flat",
                           "2016-17_Moderate Drought", "2016-17_Optimal Bed", "2016-17_Heat", "2016-17_Severe Drought", "2016-17_Optimal Flat"),
                    Ab = c("MD2013", "OB2013", "H2013", "SD2013", "OF2013",
                           "MD2014", "OB2014", "H2014", "SD2014", "OF2014",
                           "MD2015", "OB2015", "H2015", "SD2015", "OF2015",
                           "MD2016", "OB2016", "H2016", "SD2016", "OF2016")) %>% 
        mutate(
            Trial = if_else(row_number() < 10, paste0("Trial_0", row_number()), paste0("Trial_", row_number())),
        ) %>% 
        rowwise %>% 
        mutate(
            Managed_Treatment = str_split(Trials, "_")[[1]][2],
            Breeding_Cycle = str_split(Trials, "_")[[1]][1]) %>% 
        ungroup


trials_20 <- read_csv("results/Wheat_final_results_20_trials.csv", show_col_types = F) %>%
	bind_rows(wresult[,-4])
	
of2013_results <- read_csv("results/Wheat_OF2013_final_results.csv")

mix_band_of2013_result <- read_csv("results/Wheat_OF2013_G_P_COR_BLUE.csv") %>%
  mutate(
    wavelength = str_extract(ph_data, "(?<=BLUE_)\\d+(?=nm)") %>% as.numeric,
    test_date = case_when(
      str_trim(str_sub(ph_data, -4, -1)) == "0110" ~ "10-Jan",
      str_trim(str_sub(ph_data, -4, -1)) == "0117" ~ "17-Jan",
      str_trim(str_sub(ph_data, -4, -1)) == "0130" ~ "30-Jan",
      str_trim(str_sub(ph_data, -4, -1)) == "0207" ~ "07-Feb",
      str_trim(str_sub(ph_data, -4, -1)) == "0214" ~ "14-Feb",
      str_trim(str_sub(ph_data, -4, -1)) == "0219" ~ "19-Feb",
      str_trim(str_sub(ph_data, -4, -1)) == "0227" ~ "27-Feb",
      str_trim(str_sub(ph_data, -4, -1)) == "0311" ~ "11-Mar",
      str_trim(str_sub(ph_data, -4, -1)) == "0317" ~ "17-Mar"
    ),
    test_date = fct_relevel(test_date, "10-Jan", "17-Jan", "30-Jan", "07-Feb",
                            "14-Feb", "19-Feb", "27-Feb", "11-Mar", "17-Mar"
                            ),
    grow_stage = case_when(
      test_date %in% c("10-Jan", "17-Jan", "30-Jan", "07-Feb", "14-Feb") ~ "VEG",
      test_date %in% c("19-Feb", "27-Feb") ~ "HEAD",
      test_date %in% c("11-Mar", "17-Mar") ~ "GF",
    ),
    grow_stage = fct_relevel(grow_stage, "VEG", "HEAD", "GF")
  )

calTime <- read_csv("results/Wheat_OF2013_calTime.csv", show_col_types = F) %>% 
  mutate(
    method = fct_relevel(method, "PBLUP", "GBLUP", "BayesianLasso", "KAML", "RKHS", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP")
  )

  
p1=of2013_results %>% 
    filter(
        method %in% c("PBLUP", "GBLUP", "BayesianLasso", "RKHS", "KAML", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP", "FSBLUP_CV1")
    ) %>%
    mutate(
        method = fct_relevel(method, "PBLUP", "GBLUP", "BayesianLasso", "RKHS", "KAML", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP", "FSBLUP_CV1")
    ) %>%
    ggplot(aes(x = method, y = g_cor_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = g_cor_mean, ymax = g_cor_mean + g_cor_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
	geom_text(aes(y = g_cor_mean + g_cor_sd + 0.005,  # 在误差条上方添加标签
      label = sprintf("%.3f", g_cor_mean)  # 格式化数值 
	 ), position = position_dodge(0.1), vjust = -0.5,  size = 2.3,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
        "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial"),
		line = element_line(linewidth = 0.5)
	) +
    coord_cartesian(ylim = c(0, 0.65))
	
p2=mix_band_of2013_result %>% 
    ggplot(aes(x = wavelength)) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_line(aes(y = g_cor, color = "Genetic")) +
    geom_line(aes(y = p_cor, color = "Phenotypic")) +
    facet_wrap(test_date ~ grow_stage, scales = "fixed", ncol = 4) + 
    labs(x = "Wavelength", y = "Correlation", title = NULL) +  #
    scale_color_manual(values=c('Genetic' = 'red','Phenotypic' = 'black'), name = 'Correlation') +
    theme_bw() +
    theme(
        legend.position = c(0.75,0.09),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 8, family = "Arial", angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 12, face = "bold", family = "Arial",
            margin = margin(t = -10, r = 0, b = 10, l = 0, unit = "pt")),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.background = element_rect(color = "black", linewidth = 0.5),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(face = "bold", size = 10, family = "Arial"),
		line = element_line(linewidth = 0.5)
	)

p3=calTime %>%
    mutate(
        method = fct_relevel(method, "PBLUP", "GBLUP", "BayesianLasso", "RKHS", "KAML", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP")
    ) %>%
    ggplot(aes(x = method)) +
    geom_col(aes(y = log_time, fill = method), position = "identity", color = "black") +
	geom_text(aes(y = log_time + 0.005,  
      label = sprintf("%.3f", log_time) 
	 ), position = position_dodge(0.1), vjust = -0.5,  size = 2.3,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
        "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
	labs(x = "", y = expression("Calculating Time ("*Log["10"]*", s)")) +
    theme_bw() + 
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial"),
		line = element_line(linewidth = 0.5)
	)


p4=of2013_results %>% 
    mutate(
        method = fct_relevel(method, "FSBLUP_VEG", "FSBLUP_HEAD", "FSBLUP_GF", "FSBLUP")
    ) %>%
    filter(method %in% c("FSBLUP", "FSBLUP_VEG", "FSBLUP_HEAD", "FSBLUP_GF")) %>%
    ggplot(aes(x = method, y = g_cor_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = g_cor_mean, ymax = g_cor_mean + g_cor_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
	geom_text(aes(y = g_cor_mean + g_cor_sd + 0.005,
      label = sprintf("%.3f", g_cor_mean) 
	 ), position = position_dodge(0.1), vjust = -0.5,  size = 2.3,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    scale_fill_npg(alpha = 0.85) +
    coord_cartesian(ylim = c(0.4, 0.65)) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial"),
		line = element_line(linewidth = 0.5)
	)

design = "AAAABBBB
          AAAABBBB
          AAAABBBB
          CCCCDDDD
          CCCCDDDD
          CCCCDDDD"

opw = wrap_plots(p1, p2, p3, p4, design = design) +
    plot_annotation(tag_levels = "A") +
    plot_layout(guides = "auto") &
    theme(
        axis.title.x = element_text(size = 12,margin = margin(t = -20, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(size = 12, family = "Arial", color = "black"),
        plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        plot.tag = element_text(size = 18, face = "bold", family = "Arial", hjust = 0, vjust = 0),
    )
ggsave("plot/Wheat_Prediction_Accuracy.pdf", opw, width = 180, height = 225, units = "mm", device = CairoPDF)



## 20 trials 
op20 = trials_20 %>% 
    mutate(
        method = fct_relevel(method, "PBLUP", "GBLUP", "BayesianLasso", "RKHS", "KAML", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP")
    ) %>% 
    ggplot(aes(x = method, y = g_cor_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = g_cor_mean, ymax = g_cor_mean + g_cor_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    geom_text(aes(y = g_cor_mean + g_cor_sd + 0.005, label = round(g_cor_mean, 3)), vjust = -0.5, size = 2.5) +
    facet_grid(Managed_Treatment ~ Breeding_Cycle, scales = "fixed") +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
        "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial")) +
    coord_cartesian(ylim = c(0, 1))

ggsave("plot/Wheat_Prediction_Accuracy_20trials.pdf", op20, width = 320, height = 225, units = "mm", device = CairoPDF)



## all trials

univariate_acc_all <- read_csv("results/Wheat_final_results_20_trials.csv")

Cairo::CairoPDF("Prediction_Accuracy_20trials.pdf", width = 12, height = 8)
univariate_acc_all %>% 
    mutate(
        method = fct_relevel(method, "ABLUP", "GBLUP", "BayesianLasso", "RKHS", "HBLUP", "GBLUP+H","GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP")
    ) %>% 
    ggplot(aes(x = method, y = g_cor_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = g_cor_mean, ymax = g_cor_mean + g_cor_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    facet_grid(Managed_Treatment ~ Breeding_Cycle, scales = "fixed") +
    labs(x = NULL, y = NULL, title = NULL) +
    scale_fill_npg(alpha = 0.85) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial")) +
    coord_cartesian(ylim = c(0, 1))
dev.off()
