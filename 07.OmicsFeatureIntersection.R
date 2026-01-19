Canton::using(tidyverse,data.table,clusterProfiler)
mRNA_CI_vs_HTN_DEG_Table = fread('结果/mRNA基因差异分析/Table_CI_vs_HTN_DEG_Table.csv')
mRNA_ICH_vs_HTN_DEG_Table = fread('结果/mRNA基因差异分析/Table_ICH_vs_HTN_DEG_Table.csv')

protein_CI_vs_HTN_DEG_Table = fread('结果/蛋白差异分析/Table_CI_vs_HTN_DEG_Table.csv')
protein_ICH_vs_HTN_DEG_Table = fread('结果/蛋白差异分析/Table_ICH_vs_HTN_DEG_Table.csv')

# 首先mrna 差异基因 ICH_vs_HTN - CI_vs_HTN
mRNA_CI_vs_HTN_Up = mRNA_CI_vs_HTN_DEG_Table %>% dplyr::filter(LFC >= 0.585, Pvalue<0.05) %>% pull(1)
mRNA_ICH_vs_HTN_Up = mRNA_ICH_vs_HTN_DEG_Table %>% dplyr::filter(LFC >= 0.585, Pvalue<0.05) %>% pull(1)
mrna_Up = setdiff(mRNA_ICH_vs_HTN_Up, mRNA_CI_vs_HTN_Up)

mRNA_CI_vs_HTN_Down = mRNA_CI_vs_HTN_DEG_Table %>% dplyr::filter(LFC <= -0.585, Pvalue<0.05) %>% pull(1)
mRNA_ICH_vs_HTN_Down = mRNA_ICH_vs_HTN_DEG_Table %>% dplyr::filter(LFC <= -0.585, Pvalue<0.05) %>% pull(1)
mrna_Down = setdiff(mRNA_ICH_vs_HTN_Down, mRNA_CI_vs_HTN_Down)

# 蛋白  ICH_vs_HTN - CI_vs_HTN
protein_CI_vs_HTN_Up = protein_CI_vs_HTN_DEG_Table %>% dplyr::filter(LFC >= 0.26, Pvalue<0.05) %>% pull(1)
protein_ICH_vs_HTN_Up = protein_ICH_vs_HTN_DEG_Table %>% dplyr::filter(LFC >= 0.26, Pvalue<0.05) %>% pull(1)
protein_Up = setdiff(protein_ICH_vs_HTN_Up, protein_CI_vs_HTN_Up)

protein_CI_vs_HTN_Down = protein_CI_vs_HTN_DEG_Table %>% dplyr::filter(LFC <= -0.26, Pvalue<0.05) %>% pull(1)
protein_ICH_vs_HTN_Down = protein_ICH_vs_HTN_DEG_Table %>% dplyr::filter(LFC <= -0.26, Pvalue<0.05) %>% pull(1)
protein_Down = setdiff(protein_ICH_vs_HTN_Down, protein_CI_vs_HTN_Down)
outdir='结果/marker筛选'
# 蛋白和基因表达一致的
intersect(mrna_Up,protein_Up) # "CSF2RA" "CXCR2"  "IL1RN"  "LYN"    "NAIP"   "PLAU" 
intersect(mrna_Down,protein_Down) # "ADAMTS13" "DKK3"     "ERBB2"    "TMPRSS5"  "VIPR2"

marker_list <- list(mrna_Up = mrna_Up, protein_Up = protein_Up)

p <- ggvenn::ggvenn(marker_list,
    show_percentage = F, text_size = 8,
    fill_color = c("#1B9E77", "#D95F02", "#7570B3"), fill_alpha = 0.7,
    stroke_size = 0.5, stroke_color = c("white")
)
p
ggsave(str_glue("{output_dir}/Up.pdf"), width = 6, height = 6)

marker_list <- list(mrna_Down = mrna_Down, protein_Down = protein_Down)

p <- ggvenn::ggvenn(marker_list,
    show_percentage = F, text_size = 8,
    fill_color = c("#1B9E77", "#D95F02", "#7570B3"), fill_alpha = 0.7,
    stroke_size = 0.5, stroke_color = c("white")
)
p
ggsave(str_glue("{output_dir}/Dwon.pdf"), width = 6, height = 6)

# "ADAM17"   "TMPRSS5"  "PLAU"     "ADAMTS13"

# 蛋白和基因表达一致的

# mRNA高蛋白低
intersect(mrna_Up,protein_Down) # "ADAM17"  IL18RAP



marker_list <- list(mrna_Up = mrna_Up, protein_Down = protein_Down)

p <- ggvenn::ggvenn(marker_list,
    show_percentage = F, text_size = 8,
    fill_color = c("#1B9E77", "#D95F02", "#7570B3"), fill_alpha = 0.7,
    stroke_size = 0.5, stroke_color = c("white")
)
p
ggsave(str_glue("{output_dir}/mrna_Up_protein_Down.pdf"), width = 6, height = 6)


data.table(Gene=c(intersect(mrna_Up,protein_Up), intersect(mrna_Down,protein_Down), "ADAM17")) %>% fwrite("结果/marker筛选/candidates.csv")


packageVersion('caret')