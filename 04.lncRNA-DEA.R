source("/amax/home/wangcy/Alex/CODEHUB/RCode/limma_deg.R")
source('/amax/home/wangcy/Alex/CODEHUB/RCode/unique_exprs.R')
source('/amax/home/wangcy/Alex/CODEHUB/RCode/deseq.R')

Canton::using(tidyverse,data.table,clusterProfiler)
pdata=fread('Data/Set1/meta.tsv')
d1=fread('Data/ICH_HTN.lncRNA.anno.xls') %>% select(GeneName,ends_with('_count')) %>% rename(Feature=GeneName)

d1 %<>% unique_exprs %>% column_to_rownames('Feature')
d1=d1[apply(d1,1,sum) > 10,]

pdata = data.frame(row.names=colnames(d1),Group=ifelse(str_detect(colnames(d1),"^ICH"),'ICH','HTN'))
x=deseq(
    outdir = '结果/lncRNA差异分析',
    fdata = d1,
    pdata = pdata,
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


d2=fread('Data/CI_HTN.lncRNA.anno.xls') %>% select(GeneName,ends_with('_count')) %>% dplyr::rename(Feature=GeneName)
d2 %<>% unique_exprs %>% column_to_rownames('Feature')
d2=d2[apply(d2,1,sum) > 10,]

pdata = data.frame(row.names=colnames(d2),Group=ifelse(str_detect(colnames(d2),"^CI"),'CI','HTN'))
x=deseq(
    outdir = '结果/lncRNA差异分析',
    fdata = d2,
    pdata = pdata,
    control_label = 'HTN',
    case_label = 'CI',
    lfc_threshold = 0.585,
    pvalue_threshold = 0.05,
    pvalue = c("Pvalue", "Padj")[1],
    text_gene = NULL,
    hue1 = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1], "#A6A6A6"),
    hue2 = c(Canton::hue("NPG")[2], "white", Canton::hue("NPG")[1]),
    top_hue = c(Canton::hue("NPG")[2], Canton::hue("NPG")[1]),
    top_order = c("HTN", "CI")
)
