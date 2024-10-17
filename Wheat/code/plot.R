setwd("/public/home/liujf/workspace/xueyh/TempWork/h_matrix_wheat_mix")
library(tidyverse)
library(ggsci)
library(showtext)
font_add("Arial", "/public/home/liujf/workspace/xueyh/font/arial.ttf")
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


final_results <- read_csv("result_full_data/Wheat_OF2013_final_results_20241009.csv", show_col_types = F)
Cairo::CairoPDF("plot/OF2013_Prediction_Accuracy_20240923.pdf", width = 8, height = 6)
p1=final_results %>% 
    filter(
        Method %in% c("GBLUP", "BayesLasso", "RKHS", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP", "FSBLUP_CV1")
    ) %>%
    mutate(
        Method = fct_relevel(Method, "GBLUP", "BayesLasso", "RKHS", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP", "FSBLUP_CV1")
    ) %>%
    ggplot(aes(x = Method, y = g_cor_mean, fill = Method, group = Method)) +
    geom_errorbar(aes(ymin = g_cor_mean, ymax = g_cor_mean + g_cor_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    #geom_text(aes(label = round(g_cor_mean, 3)), vjust = -0.5, size = 5) +
    #scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
    #    "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    scale_fill_npg(alpha = 0.85) +
    theme_bw() +
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial")) +
    #guides(fill = guide_legend(ncol = 6)) +
    coord_cartesian(ylim = c(0, 0.65))
dev.off()


mix_band_of2013_result <- read_csv("results/Wheat_OF2013_G_P_COR_BLUE.csv")
p2=mix_band_of2013_result %>%
    mutate(
        test_date = fct_relevel(test_date, "10-Jan", "17-Jan", "30-Jan", "07-Feb", "14-Feb", "19-Feb", "27-Feb", "11-Mar", "17-Mar"),
        chpt = fct_relevel(chpt, "VEG", "HEAD", "GF")
    ) %>% 
    ggplot(aes(x = wavelength)) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_line(aes(y = g_cor, color = "Genetic")) +
    geom_line(aes(y = p_cor, color = "Phenotypic")) +
    #geom_ribbon(aes(ymin = g_cor_low, ymax = g_cor_high), alpha = 0.3) +
    #facet_wrap(~test_date, ncol = 3)  +
    facet_wrap(test_date~chpt, scales = "fixed", ncol = 4) + 
    labs(x = "Wavelength", y = "Correlation", title = NULL) +  #
    scale_color_manual(values=c('Genetic' = 'red','Phenotypic' = 'black'), name = 'Correlation') +
    #theme_minimal() +
    theme_bw() +
    theme(
        legend.position = c(0.75,0.09),   #c(0.94, 0.38)
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 8, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title.y = element_text(size = 12, face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 12, face = "bold", family = "Arial",
            margin = margin(t = -10, r = 0, b = 10, l = 0, unit = "pt")),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.background = element_rect(color = "black", linewidth = 0.5),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(face = "bold", size = 10, family = "Arial"))
dev.off()





calTime <- read_csv("result_full_data/Wheat_OF2013_calTime_20240923.csv", show_col_types = F)
p3=calTime %>% 
    mutate(
        method = fct_relevel(method, "GBLUP", "BayesLasso", "RKHS", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP")
    ) %>%
    filter(method %in% c("GBLUP", "BayesLasso", "RKHS", "HBLUP", "GBLUP+H", "GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP")) %>% 
    ggplot(aes(x = method)) +
    geom_col(aes(y = log_time, fill = method), position = "identity", color = "black") + 
    #geom_boxplot(aes(y = log_time, fill = method), position = "identity") +
    scale_fill_npg(alpha = 0.85) +
    labs(x = "", y = "Calculating Time(Log10, S)") +  #expression("Calculating Time ("*Log["10"]*", s)")
    #theme_minimal() +
    theme_bw() + 
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial"))
        #+coord_cartesian(ylim = c(0, 7))


p4=final_results %>% 
    mutate(
        Method = fct_relevel(Method, "FSBLUP_VEG", "FSBLUP_HEAD", "FSBLUP_GF", "FSBLUP")
    ) %>%
    filter(Method %in% c("FSBLUP", "FSBLUP_VEG", "FSBLUP_HEAD", "FSBLUP_GF")) %>%
    ggplot(aes(x = Method, y = g_cor_mean, fill = Method, group = Method)) +
    geom_errorbar(aes(ymin = g_cor_mean, ymax = g_cor_mean + g_cor_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    #geom_text(aes(label = round(g_cor_mean, 3)), vjust = -0.5, size = 5) +
    #scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
    #    "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    scale_fill_npg(alpha = 0.85) +
    coord_cartesian(ylim = c(0.4, 0.62)) +
    theme_bw() +
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial"))
dev.off()



design = "AAAABBBB
          AAAABBBB
          AAAABBBB
          CCCCDDDD
          CCCCDDDD
          CCCCDDDD"

CairoPDF("plot/OF2013_figures_20241009.pdf", width = 8, height = 9)
wrap_plots(p1, p2, p3, p4, design = design) +
    plot_annotation(tag_levels = "A") +
    plot_layout(guides = "auto") &
    theme(
        axis.title.x = element_text(size = 12,margin = margin(t = -20, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(size = 12, family = "Arial", color = "black"),
        plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        plot.tag = element_text(size = 18, face = "bold", family = "Arial", hjust = 0, vjust = 0),
    )
dev.off()




## all trials

univariate_acc_all <- read_csv("result_full_data/final_results_all_trials_20241009..csv")

Cairo::CairoPDF("plot/Prediction_Accuracy_20241009.pdf", width = 12, height = 8)
univariate_acc_all %>% 
    mutate(
        Method = fct_relevel(Method, "ABLUP", "GBLUP", "BayesLasso", "RKHS", "HBLUP", "GBLUP+H","GBLUP+H+A", "MegaLMM", "GOBLUP", "FSBLUP")
    ) %>% 
    ggplot(aes(x = Method, y = g_cor_mean, fill = Method, group = Method)) +
    geom_errorbar(aes(ymin = g_cor_mean, ymax = g_cor_mean + g_cor_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
    #geom_text(aes(label = round(g_cor_mean, 3)), vjust = -0.5, size = 5) +
    facet_grid(Managed_Treatment ~ Breeding_Cycle, scales = "fixed") +
    #geom_text(aes(label = round(pearson_mean, 2)), vjust = -0.5, size = 3) + 
    #scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
    #    "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#F5A623FF", "#4099A7FF")) +
    labs(x = NULL, y = NULL, title = NULL) +
    scale_fill_npg(alpha = 0.85) +
    theme_bw() +
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 30, hjust = 1, color = "black"),
        axis.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 10, family = "Arial")) +
    #guides(fill = guide_legend(ncol = 6)) +
    coord_cartesian(ylim = c(0, 1))
dev.off()