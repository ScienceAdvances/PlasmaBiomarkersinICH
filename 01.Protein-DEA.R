# install required packages
install.packages("remotes")
remotes::install_github('ScienceAdvances/Canton',force = TRUE)
Canton::using(tidyverse,data.table,clusterProfiler)

source("config.R")

fdata=readxl::read_xlsx("Data/Set1/Set1_protein.xlsx")
d=clusterProfiler::bitr(geneID=fdata$entrezID, fromType="ENTREZID", toType='SYMBOL', OrgDb="org.Hs.eg.db", drop = TRUE)
d2=merge(d,fdata,by.y='entrezID',by.x='ENTREZID',all.y=T)
fwrite(d2[,c(1,2),3],"Set1_protein.xls",sep='\t')

d=fread('Data/Set1_protein.xls')
pdata=fread('Data/Set1/meta.tsv')

d=d %>% dplyr::select(-ENTREZID,-proteinID)
colnames(d)[1] = 'Feature'
d %<>% unique_exprs()

f2=d %>% column_to_rownames('Feature') %>% as.data.frame()
f3=log2(f2+1)

x=limma_deg(
    outdir = '结果/蛋白差异分析',
    fdata = f3,
    pdata = pdata %>% column_to_rownames('Sample') %>% dplyr::filter(Group %in%c('ICH','HTN')) %>% dplyr::select(Group) %>% as.data.frame(),
    control_label = 'HTN',
    case_label = 'ICH',
    lfc_threshold = 0.26,
    pvalue_threshold = 0.05,
    pvalue = c("Pvalue", "Padj")[1],
    text_gene = NULL,
    hue1 = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1], "#A6A6A6"),
    hue2 = c(Canton::hue("NPG")[2], "white", Canton::hue("NPG")[1]),
    top_hue = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1]),
    top_order = c("HTN", "ICH")
    )


source("limma_deg.R")
x=limma_deg(
    outdir = '结果/蛋白差异分析',
    fdata = f3,
    pdata = pdata %>% column_to_rownames('Sample') %>% dplyr::filter(Group %in%c('CI','HTN')) %>% dplyr::select(Group) %>% as.data.frame(),
    control_label = 'HTN',
    case_label = 'CI',
    lfc_threshold = 0.26,
    pvalue_threshold = 0.05,
    pvalue = c("Pvalue", "Padj")[1],
    text_gene = NULL,
    hue1 = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1], "#A6A6A6"),
    hue2 = c(Canton::hue("NPG")[2], "white", Canton::hue("NPG")[1]),
    top_hue = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1]),
    top_order = c("HTN", "CI")
    )
