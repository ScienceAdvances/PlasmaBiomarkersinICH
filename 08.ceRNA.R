source("/amax/home/wangcy/Alex/CODEHUB/RCode/limma_deg.R")
source('/amax/home/wangcy/Alex/CODEHUB/RCode/unique_exprs.R')
source('/amax/home/wangcy/Alex/CODEHUB/RCode/deseq.R')

Canton::using(tidyverse,data.table,clusterProfiler)
Canton::using(multiMiR, igraph, conflicted, tidyverse, org.Hs.eg.db, org.Mm.eg.db, data.table)
library(data.table)
install.packages("multiMiR")
BiocManager::install("arrow")
library(arrow)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::between)
outdir <- "结果/ceRNA网络"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# 设置结果路径
conflicts_prefer(base::setdiff)
# mRNA高蛋白低 
mRNA <- c('ADAM17', 'IL18RAP')
# mRNA-miRNA
# mirecords mirtarbase tarbase
gene2mir <- get_multimir(
    org = "hsa",
    target = mRNA,
    table = "validated",
    summary = TRUE,
    predicted.cutoff = 1,
    predicted.cutoff.type = "p"
)
gene2mir@summary %>% as.data.frame() %>% dplyr::filter(target_symbol=='IL18RAP')
gene2mir@summary$target_symbol %>% unique()
miRNA <- gene2mir@summary %>%
    arrange(desc(all.sum)) %>%
    dplyr::filter(all.sum > 0) %>%
    dplyr::select(mature_mirna_id, target_symbol)
names(miRNA) <- c("source", "target")
miRNA$target %>% unique
fwrite(miRNA,file.path(outdir,'miRNA.csv'))
# dplyr::select(miRNA,source) %>% fwrite('/amax/home/wangcy/Alex/X261/结果/ceRNA网络/circRNA/miRNA.csv',col.names=F)
property <- tibble(
    nodes = c(
        miRNA$source,
        miRNA$target
    ),
    type = rep(c("miRNA", "mRNA"),
        times = c(
            nrow(miRNA),
            nrow(miRNA)
        )
    )
)
property <- distinct(property, nodes, .keep_all = T)
fwrite(property, str_glue("{outdir}/miRNA_property.csv"))

# miRNA-lncRNA
load("/amax/home/wangcy/Alex/DATAHUB/RNA/starbase/anno.Rdata")
lncRNA <- fread("/amax/home/wangcy/Alex/DATAHUB/RNA/starbase/starBaseV3_hg19_CLIP-seq_lncRNA_all.csv")
head(lncRNA, 2)
lncRNA %<>% dplyr::filter(
    # clipExpNum >= 4,
    miRNAname %in% miRNA$source,
    geneID %in% str_remove(lnc_anno$gene_id, "\\.\\d")
    # degraExpNum >= 1
) %>% dplyr::select(miRNAname, geneName) %>% distinct(.keep_all = TRUE)
nrow(lncRNA)
names(lncRNA) <- c("source", "target")
a1=fread('结果 copy/1.1差异分析/lncRNA差异分析/Table_ICH_vs_HTN_DEG_Table.csv') %>% dplyr::filter(Change == 'Up' ) %>% pull(1)
a2=fread('结果 copy/1.1差异分析/lncRNA差异分析/Table_CI_vs_HTN_DEG_Table.csv') %>% dplyr::filter(Change == 'Up' ) %>% pull(1)

lncRNA =lncRNA %>% dplyr::filter(target %in% base::setdiff(a1,a2)) %>% distinct()



property <- tibble(
    nodes = c(
        lncRNA$source,
        lncRNA$target
    ),
    type = rep(c("miRNA", "lncRNA"),
        times = c(
            nrow(lncRNA),
            nrow(lncRNA)
        )
    )
)
fwrite(lncRNA,file.path(outdir,'lncRNA.csv'))
property <- distinct(property, nodes, .keep_all = T)
fwrite(property, str_glue("{outdir}/lncRNA_property.csv"))


# circRNA
# 定义所需的预测程序
programs=("PITA" "RNA22" "miRmap" "DIANA-microT" "miRanda" "PicTar" "TargetScan")
prog='miRmap'
# 循环下载每个程序的数据
for prog in "${programs[@]}"
do
curl "https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=circRNA&miRNA=all&clipExpNum=5&degraExpNum=0&pancancerNum=0&programNum=1&program=${prog}&target=all&cellType=all" > "ENCORI_hg38_CLIP-seq_${prog}.circRNA.txt"
done


d1=fread('结果/ceRNA网络/ENCORI_hg38_CLIP-seq_miRmap.circRNA.txt',skip=3,sep='\t')
d1 %<>% dplyr::filter(miRNAname %in% miRNA$source) %>% tidyr::separate_longer_delim(cols=circID, delim=',')
head(d1,2)

a1=fread('结果 copy/1.1差异分析/circRNA差异分析/Table_ICH_vs_HTN_DEG_Table.csv') %>% dplyr::filter(Change == 'Up' ) %>% pull(1)
a2=fread('结果 copy/1.1差异分析/circRNA差异分析/Table_CI_vs_HTN_DEG_Table.csv') %>% dplyr::filter(Change == 'Up' ) %>% pull(1)

circRNA =d1 %>% dplyr::filter(circID %in% base::setdiff(a1,a2))  %>% group_by(miRNAname) %>% slice_max(TDMDScore,n=3,with_ties=F)%>% dplyr::select(miRNAname,circID) %>% distinct()


names(circRNA) <- c("source", "target")
property <- tibble(
    nodes = c(
        circRNA$source,
        circRNA$target
    ),
    type = rep(c("miRNA", "circRNA"),
        times = c(
            nrow(circRNA),
            nrow(circRNA)
        )
    )
)
fwrite(circRNA,file.path(outdir,'circRNA.csv'))
property <- distinct(property, nodes, .keep_all = T)
fwrite(property, str_glue("{outdir}/circRNA_property.csv"))

miRNA %<>% dplyr::filter(source %in% c(lncRNA$source,circRNA$source)) 

# to cytoscape
net <- rbind(miRNA, lncRNA, circRNA) %>% distinct
property <- tibble(
    nodes = c(
        miRNA$source,
        lncRNA$source,
        miRNA$target,
        lncRNA$target,
        circRNA$target
    ),
    type = rep(c("miRNA", "mRNA", "lncRNA",'circRNA'),
        times = c(
            nrow(miRNA) + nrow(lncRNA),
            nrow(miRNA),
            nrow(lncRNA),
            nrow(circRNA)
        )
    )
)
property <- distinct(property, nodes, .keep_all = T)
nrow(property)

fwrite(net, str_glue("{outdir}/net.csv"))
fwrite(property, str_glue("{outdir}/property.csv"))


network_data <- read.csv("data.csv", header = TRUE)
colnames(network_data) <- c("circRNA", "miRNA", "mRNA")
miRNA %>% head(2)
lncRNA %>% head(2)
colnames(miRNA) <- c('miRNA','mRNA')
colnames(lncRNA) <- c('miRNA','circRNA')

network_data <- merge(distinct(miRNA,.keep_all = T),distinct(lncRNA,.keep_all = T),by='miRNA')
# 创建空的网络对象
g <- graph.empty(n = length(c(unique(network_data$miRNA), unique(network_data$circRNA), unique(network_data$mRNA))), directed = TRUE)

# 添加节点
g <- g %>%
    set_vertex_attr("name", value = c(unique(network_data$circRNA), unique(network_data$miRNA), unique(network_data$mRNA))) %>%
    set_vertex_attr("type", value = c(
        rep("circRNA", length(unique(network_data$circRNA))),
        rep("miRNA", length(unique(network_data$miRNA))),
        rep("mRNA", length(unique(network_data$mRNA)))
    ))
g <- set_vertex_attr(g, "color", value = ifelse(V(g)$type == "circRNA", "#fb8072", ifelse(V(g)$type == "miRNA", "yellow3", "#80b1d3")))

# 添加边与边长
afedge <- c()
aflength <- c()
for (i in 1:nrow(network_data)) {
    circRNA_node <- which(V(g)$name == network_data[i, 1])
    miRNA_node <- which(V(g)$name == network_data[i, 2])
    mRNA_node <- which(V(g)$name == network_data[i, 3])
    aflength <- c(aflength, 20, 10)
    afedge <- c(afedge, circRNA_node, miRNA_node, miRNA_node, mRNA_node)
}
g <- g %>%
    add_edges(afedge) %>%
    set_edge_attr("edge.length", value = aflength)

# 添加节点大小
circRNA.size <- as.vector(scale(as.vector(table(network_data$circRNA)), center = F)) + 10
miRNA.size <- as.vector(scale(as.vector(table(network_data$miRNA)), center = F)) + 10
mRNA.size <- as.vector(scale(as.vector(table(network_data$mRNA)), center = F)) + 10
V(g)$size <- c(circRNA.size, miRNA.size, mRNA.size)


# 使用Graphopt算进行布局，保存为ceRNA.net.pdf文件
pdf(file = "ceRNA.net.pdf", height = 10, width = 10)
plot(g,
    layout = layout.graphopt(g),
    vertex.label = V(g)$name,
    vertex.label.family = "sans",
    vertex.label.cex = ifelse(V(g)$type == "circRNA", 0.8, ifelse(V(g)$type == "miRNA", 0.5, 0.5)),
    vertex.size = V(g)$size,
    vertex.color = V(g)$color,
    vertex.label.color = "black",
    edge.arrow.size = 0.5,
    edge.width = 1
)
dev.off()
