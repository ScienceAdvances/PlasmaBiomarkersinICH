outdir <- "结果/两个模型交集基因"
dir.create(outdir, showWarnings = F,recursive = TRUE)
Canton::using(tidyverse,ggvenn, glmnet, pROC, data.table, magrittr)

a=fread("结果/ICH_CI模型基因筛选/marker.csv") %>% pull(1)
b=fread("结果/ICH_HTN模型基因筛选/marker.csv") %>% pull(1)

marker_list <- list(ICH_CI = a, ICH_HTN = b)

p <- ggvenn::ggvenn(marker_list,
    show_percentage = F, text_size = 8,
    fill_color = c("#1B9E77", "#D95F02", "#7570B3"), fill_alpha = 0.7,
    stroke_size = 0.5, stroke_color = c("white")
)
p
ggsave(str_glue("{outdir}/marker_venn_plot.pdf"), width = 6, height = 6)
