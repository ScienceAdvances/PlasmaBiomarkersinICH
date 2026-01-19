Canton::using(tidyverse, data.table, clusterProfiler)
mRNA_CI_vs_HTN_DEG_Table <- fread("结果 copy/1.1差异分析/mRNA基因差异分析/Table_CI_vs_HTN_DEG_Table.csv")
outdir <- "结果/GSEA富集分析/mRNA_CI_vs_HTN"

mRNA_ICH_vs_HTN_DEG_Table <- fread("结果 copy/1.1差异分析/mRNA基因差异分析/Table_ICH_vs_HTN_DEG_Table.csv")
outdir <- "结果/GSEA富集分析/mRNA_ICH_vs_HTN"

Canton::umkdir(outdir)
Canton::using(glue)

deg_table <- mRNA_ICH_vs_HTN_DEG_Table
deg_table <- deg_table$Feature %>%
    clusterProfiler::bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE) %>%
    merge(deg_table, by.y = "Feature", by.x = "SYMBOL") %>%
    dplyr::filter(ENTREZID != "", !is.na(ENTREZID), LFC != "", !is.na(LFC))

geneList <- deg_table$LFC
names(geneList) <- deg_table$ENTREZID
i <- "gsea"
gsea <- clusterProfiler::gseKEGG(geneList = sort(geneList, decreasing = TRUE), organism = "hsa", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1)
saveRDS(gsea, file = file.path(outdir, glue("gsea.rds")))

gsea %>%
    as.data.frame() %>%
    fwrite(file.path(outdir, glue("{i}.xls")), sep = "\t")

ids <- gsea %>%
    as.data.frame() %>%
    dplyr::filter(abs(NES) > 1, p.adjust < 0.05) %>%
    pull(ID)
i <- ids[1]
for (i in ids) {
    pos <- c(0.2, 0.2)
    p <- GseaVis::gseaNb(
        object = gsea,
        kegg = F,
        addPval = T,
        newGsea = F,
        curveCol = hue3[2],
        # rankCol=c(hue1,hue3),
        geneSetID = i %>% as.character(),
        legend.position = pos
    ) +
        ggplot2::annotate("text", x = 2000, y = 1, label = "ICH", size = 5) +
        ggplot2::annotate("text", x = 13000, y = -1, label = "HTN", size = 5) +
        theme(text = element_text(size = 18), legend.position.inside = pos, plot.margin = unit(c(0.5, 1, 1, 1), "cm"))
    Canton::gs(p, outdir = outdir, name = i, w = 7, h = 6)
}


gsea <- readRDS(file.path(outdir, "gsea.rds"))
gsea_df <- as.data.frame(gsea)

d <- gsea %>%
    as.data.frame() %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::mutate(pathway = Description, pval = pvalue, padj = p.adjust, log2err = NA, ES = enrichmentScore, size = setSize, leadingEdge = leading_edge) %>%
    dplyr::select(pathway, pval, padj, log2err, ES, NES, size, leadingEdge)

pathways <- gsea@geneSets[rownames(d)]
names(pathways) <- gsea_df[names(pathways), "Description"]

p <- fgsea::plotGseaTable(pathways = pathways, stats = geneList, fgseaRes = d, gseaParam = 0.5)
Canton::gs(p, outdir = outdir, name = "GseaTable", w = 9, h = 5)

# 山脊图
used_pathways <- fread("结果/GSEA富集分析/used_pathways.tsv") %>% pull(1)
p <- GseaVis::enrich_ridge_plot(
    object = gsea,
    terms_ID = used_pathways,
    gene_list = geneList
) +
    scale_fill_distiller(palette = "Spectral") +
    scale_color_distiller(palette = "RdBu", direction = -1)
Canton::gs(p, outdir = outdir, name = "enrich_ridge", w = 9, h = 5)
