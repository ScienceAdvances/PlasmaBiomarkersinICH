source("/amax/home/wangcy/Alex/CODEHUB/config.R")
source("/amax/home/wangcy/Alex/CODEHUB/RCode/swr.R")
source("/amax/home/wangcy/Alex/CODEHUB/RCode/enrich.R")
source("/amax/home/wangcy/Alex/CODEHUB/RCode/deseq.R")
outdir <- "结果/ORA富集分析"
Canton::using(Canton, tidyverse, data.table, clusterProfiler)
mRNA_CI_vs_HTN_DEG_Table <- fread("结果/mRNA基因差异分析/Table_CI_vs_HTN_DEG_Table.csv")
mRNA_ICH_vs_HTN_DEG_Table <- fread("结果/mRNA基因差异分析/Table_ICH_vs_HTN_DEG_Table.csv")

protein_CI_vs_HTN_DEG_Table <- fread("结果/蛋白差异分析/Table_CI_vs_HTN_DEG_Table.csv")
protein_ICH_vs_HTN_DEG_Table <- fread("结果/蛋白差异分析/Table_ICH_vs_HTN_DEG_Table.csv")
# BiocManager::install('xlsx')
# 首先mrna 差异基因 ICH_vs_HTN-CI_vs_HTN--mrna_Up
mRNA_CI_vs_HTN_Up <- mRNA_CI_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC >= 0.585, Pvalue < 0.05) %>%
    pull(1)
mRNA_ICH_vs_HTN_Up <- mRNA_ICH_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC >= 0.585, Pvalue < 0.05) %>%
    pull(1)
mrna_Up <- setdiff(mRNA_ICH_vs_HTN_Up, mRNA_CI_vs_HTN_Up)
enrich(
    genelist = mrna_Up,
    suffix = "mrna_Up",
    outdir = outdir,
    fromType = "SYMBOL",
    simplify = FALSE,
    species = "hsa",
    orgdb = "org.Hs.eg.db",
    hue = c("#4CBAD4", "#E54B34")
)
conflicts_prefer(base::setdiff)
mRNA_CI_vs_HTN_Down <- mRNA_CI_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC <= -0.585, Pvalue < 0.05) %>%
    pull(1)
mRNA_ICH_vs_HTN_Down <- mRNA_ICH_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC <= -0.585, Pvalue < 0.05) %>%
    pull(1)
mrna_Down <- base::setdiff(mRNA_ICH_vs_HTN_Down, mRNA_CI_vs_HTN_Down)
fwrite(data.table(mrna_Down=mrna_Down),"mrna_Down.xls")

source("/amax/home/wangcy/Alex/CODEHUB/RCode/enrich.R")

enrich(
    genelist = mrna_Down,
    suffix = "mrna_Down",
    outdir = outdir,
    fromType = "SYMBOL",
    simplify = FALSE,
    species = "hsa",
    orgdb = "org.Hs.eg.db",
    hue = c("#4CBAD4", "#E54B34")
)

# 蛋白  ICH_vs_HTN - CI_vs_HTN
protein_CI_vs_HTN_Up <- protein_CI_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC >= 0.26, Pvalue < 0.05) %>%
    pull(1)
protein_ICH_vs_HTN_Up <- protein_ICH_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC >= 0.26, Pvalue < 0.05) %>%
    pull(1)
protein_Up <- setdiff(protein_ICH_vs_HTN_Up, protein_CI_vs_HTN_Up)
enrich(
    genelist = protein_Up,
    suffix = "protein_Up",
    outdir = outdir,
    fromType = "SYMBOL",
    simplify = FALSE,
    species = "hsa",
    orgdb = "org.Hs.eg.db",
    hue = c("#4CBAD4", "#E54B34")
)


protein_CI_vs_HTN_Down <- protein_CI_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC <= -0.26, Pvalue < 0.05) %>%
    pull(1)
protein_ICH_vs_HTN_Down <- protein_ICH_vs_HTN_DEG_Table %>%
    dplyr::filter(LFC <= -0.26, Pvalue < 0.05) %>%
    pull(1)
protein_Down <- setdiff(protein_ICH_vs_HTN_Down, protein_CI_vs_HTN_Down)
enrich(
    genelist = protein_Down,
    suffix = "protein_Down",
    outdir = outdir,
    fromType = "SYMBOL",
    simplify = FALSE,
    species = "hsa",
    orgdb = "org.Hs.eg.db",
    hue = c("#4CBAD4", "#E54B34")
)


c(base::intersect(mrna_Up, protein_Up), base::intersect(mrna_Down, protein_Down), "ADAM17", "IL18RAP")

cancidate <- fread("结果/marker筛选/candidates.csv") %>% pull(1)

enrich(cancidate,
    suffix = "筛选的13个特异基因",
    outdir = outdir,
    fromType = "SYMBOL",
    simplify = FALSE,
    species = "hsa",
    orgdb = "org.Hs.eg.db",
    hue = c("#4CBAD4", "#E54B34")
)
