setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/")
library(tidyverse)
library(ggsci)
library(showtext)
font_add("Arial", "arial.ttf")
showtext::showtext_auto()
library(Cairo)
library(patchwork)

# accuracy
result <- read_csv("results/pig_results.csv", show_col_types = F) %>%
    mutate(
		method = case_when(method == "HAT" ~ "ssGBLUP_AT", method == "HGT" ~ "ssGBLUP_GT", TRUE ~ method),
        method = fct_relevel(method, "PBLUP", "GBLUP", "BayesianLasso", "RKHS", "KAML",  
                             "ssGBLUP_AT", "ssGBLUP_GT", "GOBLUP", "FSBLUP")
    )
# time
resultt <- read_csv("results/pig_times.csv", show_col_types = F) %>%
    mutate(
		method = case_when(method == "HAT" ~ "ssGBLUP_AT", method == "HGT" ~ "ssGBLUP_GT", TRUE ~ method),
        method = fct_relevel(method, "PBLUP", "GBLUP", "BayesianLasso", "RKHS", "KAML",  
                             "ssGBLUP_AT", "ssGBLUP_GT", "GOBLUP", "FSBLUP")
    )

p1=result %>% 
    ggplot(aes(x = method, y = pearson_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
	geom_text(aes(y = pearson_mean + pearson_sd + 0.005, 
      label = sprintf("%.3f", pearson_mean)   
	 ), position = position_dodge(0.1), vjust = -0.5,  size = 4,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    facet_wrap(~ trait, scales = "fixed") +
    labs(x = NULL, y = "Prediction Performance", title = NULL, tag = "A") +
    scale_fill_npg(alpha = 0.85) +
    theme_bw(base_size = 36) +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.title = element_text(size = 16, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 16, face = "bold", family = "Arial"),
        title = element_text(size = 16, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 14, family = "Arial")) +
    coord_cartesian(ylim = c(0.1, 0.37))

p2=resultt %>% 
    ggplot(aes(x = method, y = logtime, fill = method, group = method)) +
    geom_bar(stat = "identity", color = "black") +
	geom_text(aes(y = logtime + 0.01,  
      label = sprintf("%.3f", logtime)   
	 ), position = position_dodge(0.1), vjust = -0.5,  size = 4,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    scale_fill_npg() +
    facet_wrap(~trait, scales = "fixed") +
    labs(x = NULL, y = expression("Calculating Time ("*Log["10"]*", s)"), title = NULL, tag = "B") +
    coord_cartesian(ylim = c(0.1, 5.5)) +
    theme_bw(base_size = 36) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 10, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
        axis.title = element_text(size = 16, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 16, face = "bold", family = "Arial"),
        title = element_text(size = 16, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 14, family = "Arial"))

op = p1/p2 +
    plot_layout(guides = 'collect', axis_title = "collect") & 
    theme(
        axis.title.x = element_text(size = 14,margin = margin(t = -20, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
        plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        plot.tag = element_text(size = 18, face = "bold", family = "Arial", hjust = 0, vjust = 0))
		
ggsave("plot/Pig_Prediction_Accuracy.pdf", op, width = 300, height = 225, units = "mm", device = CairoPDF)
