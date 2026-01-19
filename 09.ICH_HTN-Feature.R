source('/amax/home/wangcy/Alex/CODEHUB/RCode/unique_exprs.R')

output_dir <- "结果/ICH_HTN模型基因筛选"
dir.create(output_dir, showWarnings = F,recursive = TRUE)
Canton::using(ggvenn, glmnet, pROC, data.table, magrittr)
cancidate <- fread('结果/marker筛选/candidates.csv') %>% pull(1)

fdata <- data.table::fread("Data/Set1_mRNA_tpm.xls", data.table = F) %>% dplyr::rename(Feature=Symbol) %>% unique_exprs  %>% 
    tibble::column_to_rownames("Feature") %>% t() %>% as.data.frame()

pdata <- data.table::fread("Data/Set1/meta.tsv", data.table = F)

cancidate <- intersect(colnames(fdata),cancidate)
X <- fdata[,cancidate]

intersect(rownames(fdata), pdata$Sample2)
label <- pdata %>%
    select(Sample2, Group) %>%
    tibble::column_to_rownames("Sample2") %>% dplyr::filter(Group %in% c('HTN',"ICH"))
X = X[rownames(label),]
X <- log2(X+1)
seed <- 0
source('/amax/home/wangcy/Alex/CODEHUB/RCode/LASSO_Feature.R')
source('/amax/home/wangcy/Alex/CODEHUB/RCode/RF_Feature.R')
source('/amax/home/wangcy/Alex/CODEHUB/RCode/SVM_Feature.R')
library(tidyverse)
lasso_res <- LASSO_Feature(as.matrix(X), factor(label$Group), outdir = output_dir, seed = seed)
lasso_res
rf_res <- RF_Feature(feature = as.data.frame(t(X)), signature = cancidate, label = label, outdir = output_dir, seed = seed)
rf_res

svm_res <- SVM_Feature(X, signature = cancidate, y = label, outdir = output_dir, seed = seed)
marker_list <- list(LASSO = lasso_res, RandomForest = rf_res, SVM = svm_res)
marker <- purrr::reduce(marker_list, base::intersect)
marker
fwrite(data.table(GeneID = marker), str_glue("{output_dir}/marker.csv"))

marker_list <- list(LASSO = lasso_res, RandomForest = rf_res, SVM = svm_res)

p <- ggvenn::ggvenn(marker_list,
    show_percentage = F, text_size = 8,
    fill_color = c("#1B9E77", "#D95F02", "#7570B3"), fill_alpha = 0.7,
    stroke_size = 0.5, stroke_color = c("white")
)
p
ggsave(str_glue("{output_dir}/marker_venn_plot.pdf"), width = 6, height = 6)
