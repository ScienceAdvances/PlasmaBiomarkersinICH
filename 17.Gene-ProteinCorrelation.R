source('/amax/home/wangcy/Alex/CODEHUB/RCode/unique_exprs.R')
conflicted::conflicts_prefer(dplyr::select)

dir.create(outdir, showWarnings = F,recursive = TRUE)
Canton::using(tidyverse, ggvenn, glmnet, pROC, data.table, magrittr)
cancidate <- fread('结果/marker筛选/candidates.csv') %>% pull(1)

outdir <- "结果/相关性散点图/ICH-CI-HTN"
pdata <- data.table::fread("Data/Set1/meta.tsv", data.table = F) %>% dplyr::filter(Group %in% c('ICH','CI','HTN'))
outdir <- "结果/相关性散点图/ICH-CI"
pdata <- data.table::fread("Data/Set1/meta.tsv", data.table = F) %>% dplyr::filter(Group %in% c('ICH','CI'))
outdir <- "结果/相关性散点图/ICH-HTN"
pdata <- data.table::fread("Data/Set1/meta.tsv", data.table = F) %>% dplyr::filter(Group %in% c('ICH','HTN'))
outdir <- "结果/相关性散点图/CI-HTN"
pdata <- data.table::fread("Data/Set1/meta.tsv", data.table = F) %>% dplyr::filter(Group %in% c('CI','HTN'))


Canton::mkdir(outdir)
mRNA <- data.table::fread("Data/Set1_mRNA_tpm.xls", data.table = F) %>% dplyr::rename(Feature=Symbol) %>% unique_exprs  %>% 
    tibble::column_to_rownames("Feature") %>% t() %>% as.data.frame()

cancidate <- intersect(colnames(mRNA),cancidate)
mRNA <- mRNA[,cancidate]
mRNA <- log2(mRNA+1)

protein <- data.table::fread("Data/Set1_protein.xls", data.table = F) %>% dplyr::select(-ENTREZID,-proteinID) %>% 
    dplyr::rename(Feature=SYMBOL) %>% unique_exprs  %>% 
    tibble::column_to_rownames("Feature") %>% t() %>% as.data.frame()
protein <- log2(protein+1)
protein <- protein[,cancidate]
protein=merge(protein,select(pdata,Sample,Sample2),by.x=0,by.y='Sample') %>% select(-Row.names) %>% column_to_rownames('Sample2')

common_sample = intersect(rownames(mRNA), rownames(protein))

mRNA=mRNA[common_sample,]
protein=protein[common_sample,]

for(i in cancidate){ 
df <- data.frame(x = mRNA[,i,drop=T], y = protein[,i,drop=T])
# 计算Spearman相关性（关键一步！）
spearman_test <- cor.test(df$x, df$y, method = "spearman")
rho <- sprintf("%.2f", spearman_test$estimate)  # 保留两位小数
pval <- sprintf("%.3f", spearman_test$p.value)  # 保留三位小数
# 绘制散点图 + 标注结果
p=ggplot(df, aes(x = x, y = y) ) +
  geom_point(alpha = 0.7, size = 2,color=Canton::hue('NPG')[2]) +  # 点有点透明，不刺眼
  labs(title = "", x = glue::glue("{i} mRNA"), y = glue::glue("{i} Protein")) +
  geom_smooth(method = "lm", se = FALSE, color = Canton::hue('NPG')[1], linetype = "solid",alpha=0.6) + 
  theme_bw(base_size = 16) +
  ggtitle(paste("rho =", rho, ", p =", pval))+
  theme(plot.title = element_text(hjust=0.5))
#   annotate("text", 
#            x = 0.2,
#            y = 0.95,
#            label = paste("rho =", rho, ", p =", pval),
#            size = 8, color = Canton::hue('NPG')[3], fontface = "bold")  # 红色加粗，一眼看到！
Canton::gs(p, name = i, outdir = outdir, format = "pdf", w =5, h = 5)
}
