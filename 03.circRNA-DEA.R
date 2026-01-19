source("/amax/home/wangcy/Alex/CODEHUB/RCode/limma_deg.R")
source('/amax/home/wangcy/Alex/CODEHUB/RCode/unique_exprs.R')
source('/amax/home/wangcy/Alex/CODEHUB/RCode/deseq.R')

Canton::using(tidyverse,data.table,clusterProfiler)
pdata=fread('Data/Set1/meta.tsv')
# circRNA差异分析
d_mran=readxl::read_xlsx('Data/cirRNAcount.xlsx') %>% column_to_rownames('circRNA_ID') %>% dplyr::select(-Location,-Source_Gene)
colnames(d_mran)
# fwrite(d_mran,"Data/Set1_mRNA.xls",sep='\t')
pdata=fread('Data/Set1/meta.tsv')
# common=intersect(pdata %>% dplyr::filter(Sample %in% colnames(d)) %>% pull(Sample2), colnames(d_mran))
# d_mran %<>% dplyr::select(dplyr::all_of(common))

d_mran=d_mran[apply(d_mran,1,sum) > 10,]
dim(d_mran)
x=deseq(
    outdir = '结果/circRNA差异分析',
    fdata = d_mran,
    pdata = pdata %>% column_to_rownames('Sample2') %>% dplyr::filter(Group %in%c('ICH','HTN')) %>% dplyr::select(Group) %>% as.data.frame(),
    control_label = 'HTN',
    case_label = 'ICH',
    lfc_threshold = 0.585,
    pvalue_threshold = 0.05,
    pvalue = c("Pvalue", "Padj")[1],
    text_gene = NULL,
    hue1 = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1], "#A6A6A6"),
    hue2 = c(Canton::hue("NPG")[2], "white", Canton::hue("NPG")[1]),
    top_hue = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1]),
    top_order = c("HTN", "ICH")
)

x=deseq(
    outdir = '结果/circRNA差异分析',
    fdata = d_mran,
    pdata = pdata %>% column_to_rownames('Sample2') %>% dplyr::filter(Group %in%c('CI','HTN')) %>% dplyr::select(Group) %>% as.data.frame(),
    control_label = 'HTN',
    case_label = 'CI',
    lfc_threshold = 0.585,
    pvalue_threshold = 0.05,
    pvalue = c("Pvalue", "Padj")[1],
    text_gene = NULL,
    hue1 = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1], "#A6A6A6"),
    hue2 = c(Canton::hue("NPG")[2], "white", Canton::hue("NPG")[1]),
    top_hue = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1]),
    top_order = c("HTN", "ICH")
)