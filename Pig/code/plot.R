setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_pig/")
library(tidyverse)
library(ggsci)
library(showtext)
font_add("Arial", "arial.ttf")
showtext::showtext_auto()
library(Cairo)
library(patchwork)

final_results <- read_csv("results/Prediction_Accuracy.csv", show_col_types = F) %>%
    group_by(trait, Method) %>% 
    summarise(
        pearson_sd = sd(pearson_mean, na.rm = TRUE),
        pearson_mean = mean(pearson_mean, na.rm = TRUE),
        bias_sd = sd(bias_mean, na.rm = TRUE),
        bias_mean = mean(bias_mean, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    mutate(
        Method = fct_relevel(Method, "ABLUP", "GBLUP", "BayesLasso", "GOBLUP", "FSBLUP")
    )

##  histgram plot
p1=final_results %>% 
    ggplot(aes(x = Method, y = pearson_mean, fill = Method, group = Method)) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    facet_wrap(~ trait, scales = "fixed") +
    labs(x = NULL, y = "Prediction Performance", title = NULL, tag = "A") +
    scale_fill_npg(alpha = 0.85) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
        axis.title = element_text(size = 16, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 16, face = "bold", family = "Arial"),
        title = element_text(size = 16, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 14, family = "Arial"))
    #guides(fill = guide_legend(ncol = 6)) +
    #coord_cartesian(ylim = c(0, 1))



## cal time
times_result <- read_csv("results/times.csv") %>%  
    mutate(
        method = fct_relevel(method, "ABLUP", "GBLUP", "BayesLasso", "GOBLUP", "FSBLUP")
    )

p2=times_result %>% 
    ggplot(aes(x = method, y = log_time, fill = method, group = method)) +
    #geom_col(position = "identity", color = "black") +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_npg() +
    facet_wrap(~trait, scales = "fixed") +
    labs(x = NULL, y = "Calculating Times(Log10, S)", title = NULL, tag = "B") +
    coord_cartesian(ylim = c(0.1, 5.5)) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
        axis.title = element_text(size = 16, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 16, face = "bold", family = "Arial"),
        title = element_text(size = 16, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 14, family = "Arial"))
dev.off()


Cairo::CairoPDF("Prediction_Accuracy.pdf", width = 12, height = 8)
p1/p2 +
    plot_layout(guides = 'collect', axis_title = "collect") & #& theme(plot.margin = margin(0.1, 0.1, 0, 0.1, unit = "cm")) 
    theme(
        axis.title.x = element_text(size = 14,margin = margin(t = -20, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
        plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        plot.tag = element_text(size = 18, face = "bold", family = "Arial", hjust = 0, vjust = 0))
dev.off()

