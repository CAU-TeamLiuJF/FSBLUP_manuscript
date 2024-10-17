setwd("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2")
library(tidyverse)
library(ggsci)
library(showtext)
font_add("Arial", "/public/home/liujf/workspace/xueyh/font/arial.ttf")
showtext::showtext_auto()
library(Cairo)
library(patchwork)

wd=getwd()

para <- expand.grid(x=1:3, y=1:20)

## 高低遗传力结果整合
final_result <- bind_rows(
        read_csv("result/simulate_result_20rep_high_h2_20241008.csv", show_col_types = F) %>% 
            mutate(Trait = "T1"), 
        read_csv("result/simulate_result_20rep_low_h2_20241008.csv", show_col_types = F) %>% 
            mutate(Trait = "T2")
    ) %>% 
    arrange(Scenario, Trait, Method) %>% 
    filter(
        !(Scenario == 2 & Method %in% c("G + T", "T")) & 
        !(Scenario == 3 & Method %in% c("G + T", "G", "T")),
        !Method %in% c("BL")
    ) %>%
    mutate(
        Method = fct_relevel(Method, "A", "G", "T", "G + T", "G + T + A", "H_AG", "H_AT", "H_GT", "GOBLUP", "FS"),
        Method = case_when(Method == "Mix" ~ "FS", , .default = Method),
        Scenario = case_when(Scenario == 1 ~ "Scenario 1", Scenario == 2 ~ "Scenario 2", Scenario == 3 ~ "Scenario 3")
    ) #%>% write_csv("result/simulate_final_result_20241008.csv")

final_result <- read_csv("result/simulate_final_result_20241008.csv") %>% 
    mutate(
        Method = fct_relevel(Method, "ABLUP", "GBLUP", "TBLUP", "GBLUP+T", "GBLUP+T+A", "ssGBLUP_AG", "ssGBLUP_AT", "ssGBLUP_GT", "GOBLUP", "FSBLUP")
    ) %>% 
    filter(!Method %in% "GBLUP+T+A")

p1 = final_result %>%
    filter(Trait == "T1") %>%
    ggplot(aes(x = Method, y = pearson_mean, fill = Method, group = Method)) +
    #geom_text(aes(label = round(pearson_mean, 3)), vjust = -0.5, size = 5) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    facet_grid(Trait ~ Scenario, scales = "free_x", space = "free") +
    #geom_text(aes(label = round(pearson_mean, 2)), vjust = -0.5, size = 3) + 
    #scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
    #    "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
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
    #guides(fill = guide_legend(ncol = 6)) +
    coord_cartesian(ylim = c(0.1, 0.61))
p2 = final_result %>%
    filter(Trait == "T2") %>%
    ggplot(aes(x = Method, y = pearson_mean, fill = Method, group = Method)) +
    #geom_text(aes(label = round(pearson_mean, 3)), vjust = -0.5, size = 5) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    facet_grid(Trait ~ Scenario, scales = "free_x", space = "free") +
    #geom_text(aes(label = round(pearson_mean, 2)), vjust = -0.5, size = 3) + 
    #scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
    #    "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
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
    #guides(fill = guide_legend(ncol = 6)) +
    coord_cartesian(ylim = c(0, 0.4))

Cairo::CairoPDF("plot/Prediction_Accuracy_all_20241008.pdf", width = 12, height = 7)
p1 / p2 +
    plot_layout(guides = 'collect', axis_title = "collect") & #& theme(plot.margin = margin(0.1, 0.1, 0, 0.1, unit = "cm")) 
    theme(
        axis.title.x = element_text(size = 14),  # margin = margin(t = -20, r = 0, b = 10, l = 0, unit = "pt")
        axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
        #plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        plot.tag = element_text(size = 18, face = "bold", family = "Arial", hjust = 0, vjust = 0),
    )
dev.off()

Cairo::CairoPDF("plot/Prediction_Bias_all_20241008.pdf", width = 12, height = 7)
final_result %>%
    ggplot(aes(x = Method, y = bias_mean, fill = Method, group = Method)) +
    #geom_text(aes(label = round(pearson_mean, 3)), vjust = -0.5, size = 5) +
    geom_errorbar(aes(ymin = bias_mean, ymax = bias_mean + bias_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    facet_grid(Trait ~ Scenario, scales = "free_x", space = "free") +
    #geom_text(aes(label = round(pearson_mean, 2)), vjust = -0.5, size = 3) + 
    #scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
    #    "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
    labs(x = "", y = "unbiasedness", title = "") +
    scale_fill_npg(alpha = 0.85) +
    theme_bw() +
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 12, family = "Arial")) #+
    #guides(fill = guide_legend(ncol = 6)) +
    #coord_cartesian(ylim = c(0, 0.7))
dev.off()






