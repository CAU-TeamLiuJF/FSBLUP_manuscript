setwd("/public/home/liujf/workspace/xueyh/TempWork/h_mat_sim_2")
library(tidyverse)
library(ggsci)
library(showtext)
font_add("Arial", "arial.ttf")
showtext::showtext_auto()
library(Cairo)
library(patchwork)
library(ggh4x)


final_result <- read_csv("results/simulate_final_result.csv") %>% 
    mutate(
        method = fct_relevel(method, "PBLUP", "GBLUP", "BayesianLasso", "RKHS", "KAML", "TBLUP", "GBLUP+T", "GBLUP+T+A", 
                             "ssGBLUP_AG", "ssGBLUP_AT", "ssGBLUP_GT", "MegaLMM", "GOBLUP", "FSBLUP")
    )

p1 = final_result %>%
    filter(trait == "T1") %>%
    ggplot(aes(x = method, y = pearson_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
	geom_text(aes(y = pearson_mean + pearson_sd + 0.005,  # 在误差条上方添加标签
      label = sprintf("%.3f", pearson_mean)  # 格式化数值 
	 ), position = position_dodge(0.1), vjust = -0.5,  size = 4,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    facet_grid(trait ~ scenario, scales = "free_x", space = "free") +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    #scale_fill_npg(alpha = 0.85) +
	scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
        "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#EE3377FF", "#DDAA33FF", "#33BBEEFF", "#BB5566FF", "#CC3311FF", "#332288FF")) +
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
    filter(trait == "T2") %>%
    ggplot(aes(x = method, y = pearson_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
	geom_text(aes(y = pearson_mean + pearson_sd + 0.005,  # 在误差条上方添加标签
      label = sprintf("%.3f", pearson_mean)  # 格式化数值 
	 ), position = position_dodge(0.1), vjust = -0.5,  size = 4,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    facet_grid(trait ~ scenario, scales = "free_x", space = "free") +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    #scale_fill_npg(alpha = 0.85) +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
        "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#EE3377FF", "#DDAA33FF", "#33BBEEFF", "#BB5566FF" , "#CC3311FF", "#332288FF")) +
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

#Cairo::CairoPDF("plot/Simulation_Prediction_Accuracy.pdf", width = 12, height = 7)
ops = p1 / p2 +
    plot_layout(guides = 'collect', axis_title = "collect") & #& theme(plot.margin = margin(0.1, 0.1, 0, 0.1, unit = "cm")) 
    theme(
        axis.title.x = element_text(size = 14),  # margin = margin(t = -20, r = 0, b = 10, l = 0, unit = "pt")
        axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
        #plot.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        plot.tag = element_text(size = 18, face = "bold", family = "Arial", hjust = 0, vjust = 0),
    )
#dev.off()

ops = final_result %>%
    #filter(trait == "T1") %>%
    ggplot(aes(x = method, y = pearson_mean, fill = method, group = method)) +
    geom_errorbar(aes(ymin = pearson_mean, ymax = pearson_mean + pearson_sd), width = 0.15, position = position_dodge(0.1)) +
    geom_col(position = "identity", color = "black") +
	geom_text(aes(y = pearson_mean + pearson_sd + 0.005,  # 在误差条上方添加标签
      label = sprintf("%.3f", pearson_mean)  # 格式化数值 
	 ), position = position_dodge(0.1), vjust = 0,  size = 3,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    facet_grid(trait ~ scenario, scales = "free", space = "free") +
    labs(x = NULL, y = "Prediction Performance", title = NULL) +
    #scale_fill_npg(alpha = 0.85) +
	scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
        "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#EE3377FF", "#DDAA33FF", "#33BBEEFF", "#BB5566FF", "#CC3311FF", "#332288FF")) +
    theme_bw() +
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 12, family = "Arial")
     ) +
	 facetted_pos_scales(
        y = list(
            trait == "T1" ~ scale_y_continuous(limits = c(0.1, 0.61)),
            trait == "T2" ~ scale_y_continuous(limits = c(0.1, 0.4))
        )
    )
ggsave("plot/Simulation_Prediction_Accuracy2.pdf", ops, width = 250, height = 200, units = "mm", device = CairoPDF)


opst = final_result %>%
    mutate(logtime = log10(time)) %>%
    ggplot(aes(x = method, y = logtime, fill = method, group = method)) +
	geom_col(position = "identity", color = "black") +
	geom_text(aes(y = logtime + 0.01,  # 在误差条上方添加标签
      label = sprintf("%.3f", logtime)  # 格式化数值 
	 ), position = position_dodge(0.1), vjust = -0.3,  size = 3,
		color = "black",
		fontface = "bold",
		family = "Arial"
	) +
    facet_grid(trait ~ scenario, scales = "free", space = "free") +
    labs(x = NULL, y = expression("Calculating Time ("*Log["10"]*", s)")) +
    #scale_fill_npg(alpha = 0.85) +
	scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
        "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#EE3377FF", "#DDAA33FF", "#33BBEEFF", "#BB5566FF", "#CC3311FF", "#332288FF")) +
    theme_bw() +
    theme(
        legend.position = "none",
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 12, family = "Arial", angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 14, face = "bold", family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        title = element_text(size = 12, face = "bold", family = "Arial"),
        strip.background = element_rect(fill = alpha("grey80", 0.85)),
        strip.text = element_text(face = "bold", size = 12, family = "Arial")
     )
ggsave("plot/Simulation_Time.pdf", opst, width = 250, height = 200, units = "mm", device = CairoPDF)
