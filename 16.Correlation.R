outdir="结果/基因与临床相关性"

m1 <- readxl::read_excel("Data/最终167例信息收集表.xlsx", sheet = "167例调查表") %>%
    dplyr::select(-`编号`, -Group1) %>%
    dplyr::mutate(
        Group =
            case_when(Group == "高血压组" ~ "HTN", Group == "脑梗组" ~ "CI", Group == "脑出血组" ~ "ICH", Group == "正常组" ~ "CTRL", TRUE ~ NA)
    )
m2 <- readxl::read_excel("Data/第二人群73例临床数据.xlsx", sheet = "Sheet2") %>%
    dplyr::select(-Group1)

m <- rbind(m1, m2)
colnames(m) %<>% str_replace_all(" ", "") %>% str_replace_all("-", "")

cancidate <- fread("结果/ICH_CI模型基因筛选/marker.csv") %>% pull(1)
d_mran <- readxl::read_xlsx("Data/Set1_mRNA.xlsx") %>% column_to_rownames("GeneName")
# d_protein=fread('Data/Set1_protein.xls')
m_mrna <- m %>%
    distinct(Sample, .keep_all = T) %>%
    column_to_rownames("Sample") %>%
    dplyr::select(-Group, -Sample1)
m_mrna <- m_mrna[colnames(d_mran), ]
x <- d_mran[cancidate, ]
y <- m_mrna

# 转换为矩阵格式
x <- as.matrix(x) %>% apply(c(1, 2), as.numeric)
y <- as.matrix(y) %>% apply(c(1, 2), as.numeric)

x[is.na(x)] <- 0
y[is.na(y)] <- 0

spec <- as.list(seq_along(cancidate))
names(spec) <- cancidate

mantel <- linkET::mantel_test(
    spec = t(x), env = y,
    spec_select = spec
) %>%
    dplyr::mutate(
        rd = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, 0.6, Inf), labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")),
        pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
    )

p <- linkET::qcorrplot(linkET::correlate(y, method = "spearman"), type = "upper", diag = TRUE) +
    linkET::geom_square() +
    linkET::geom_couple(aes(color = pd, size = rd, xend = .xend, yend = .yend), data = mantel, curvature = linkET::nice_curvature()) +
    geom_mark(sep = "\n", size = 2, sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, color = "black") +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
    scale_size_manual(values = c(0.2, 0.5, 1, 1.5, 2, 2.5)) +
    scale_colour_manual(values = color_pal(3)) +
    guides(
        size = guide_legend(title = "Mantel's r", override.aes = list(color = "grey35"), order = 2),
        color = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
        fill = guide_colorbar(title = "Spearman's rho", order = 3)
    )
Canton::gs(p,name="mRNA第一群", outdir = outdir, format = "pdf", w = 12, h = 12)
