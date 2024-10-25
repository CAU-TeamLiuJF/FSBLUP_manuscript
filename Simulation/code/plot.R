setwd("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2")
library(tidyverse)
library(ggsci)
library(showtext)
font_add("Arial", "arial.ttf")
showtext::showtext_auto()
library(Cairo)
library(patchwork)


final_result <- read_csv("results/simulate_final_result.csv") %>% 
    mutate(
        Method = fct_relevel(Method, "ABLUP", "GBLUP", "TBLUP", "GBLUP+T", "GBLUP+T+A", 
                             "ssGBLUP_AG", "ssGBLUP_AT", "ssGBLUP_GT", "GOBLUP", "FSBLUP")
    )

p1 = final_result %>%
    filter(Trait == "T1") %>%
    ggplot(aes(x = Method, y = pearson_mean, fill = Method, group = Method)) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    facet_grid(Trait ~ Scenario, scales = "free_x", space = "free") +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    scale_fill_npg(alpha = 0.85) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
        #axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 12, family = "Arial")
     ) +
    coord_cartesian(ylim = c(0.1, 0.61))
p2 = final_result %>%
    filter(Trait == "T2") %>%
    ggplot(aes(x = Method, y = pearson_mean, fill = Method, group = Method)) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    facet_grid(Trait ~ Scenario, scales = "free_x", space = "free") +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    scale_fill_npg(alpha = 0.85) +
    theme_bw() +
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 12, family = "Arial")
     ) +
    coord_cartesian(ylim = c(0, 0.4))

Cairo::CairoPDF("Prediction_Accuracy.pdf", width = 12, height = 7)
p1 / p2 +
    plot_layout(guides = 'collect', axis_title = "collect") & #& theme(plot.margin = margin(0.1, 0.1, 0, 0.1, unit = "cm")) 
    theme(
        axis.title.x = element_text(size = 14),  # margin = margin(t = -20, r = 0, b = 10, l = 0, unit = "pt")
        axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
        #plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        plot.tag = element_text(size = 18, face = "bold", family = "Arial", hjust = 0, vjust = 0),
    )
dev.off()






