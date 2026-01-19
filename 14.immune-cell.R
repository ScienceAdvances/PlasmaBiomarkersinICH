source("/amax/home/wangcy/Alex/CODEHUB/config.R")
source("/amax/home/wangcy/Alex/CODEHUB/RCode/unique_exprs.R")
source("/amax/home/wangcy/Alex/CODEHUB/RCode/immune_module.R")
source("/amax/home/wangcy/Alex/CODEHUB/RCode/immune_score.R")
source("/amax/home/wangcy/Alex/CODEHUB/RCode/catplot.R")

using(arrow)
outdir <- "结果/免疫浸润"
mkdir(outdir)
# 导入数据
fdata <- data.table::fread("Data/Set1_mRNA_tpm.xls", data.table = F) %>%
      dplyr::rename(Feature = Symbol) %>%
      unique_exprs() %>%
      tibble::column_to_rownames("Feature") %>%
      t() %>%
      as.data.frame()
pdata <- data.table::fread("Data/Set1/meta.tsv", data.table = F) %>%
      dplyr::select(Sample2, Group) %>%
      tibble::column_to_rownames("Sample2")
conflicts_prefer(dplyr::select)
conflicts_prefer(data.table::between)
res <- immune_module(t(fdata), tumor = FALSE, Group = select(pdata, Group), outdir = outdir, method = c("ssgsea", "cibersort"))

catplot(fdata=res$ssgsea %>% column_to_rownames('Sample'), pdata=pdata, outdir = outdir,filename='ssgsea')
catplot(fdata=res$cibersort %>% column_to_rownames('Sample'), pdata=pdata, outdir = outdir,filename='cibersort')

cancidate <- fread("结果/ICH_CI模型基因筛选/marker.csv") %>% pull(1)
res <- readRDS("结果/免疫浸润/immune_score_res.rds") %>% 
head(res$ssgsea,2)
fdata=res$ssgsea %>% column_to_rownames('Sample')
pdata=pdata

fdata <- t(fdata)
using(tidyverse, cols4all, stringr)
for (i in cancidate) {
      data <- fdata[ i,,drop=FALSE] %>%
            t() %>%
            as.data.table(keep.rownames = "Sample") %>%
            merge(res$ssgsea, by = "Sample")
      score <- pull(data, i)
      d2 <- map_dfr(setdiff(colnames(data), c("Sample", i)), function(x) {
            ct <- cor.test(pull(data, x), score)
            return(tibble(Feature = x, R = as.numeric(ct$estimate), P = ct$p.value))
      })
      data.table(N=c4a_palettes()) %>% fwrite('1.csv')
      dat <- d2 %>% arrange(R)
      dat$Feature <- factor(dat$Feature, levels = dat$Feature)
      p1 <- ggplot(dat, aes(x = Feature, y = R)) +
            coord_flip() +
            geom_col(aes(fill = R), width = 0.1) + # 添加条形图，收窄柱子为一条线
            geom_point(aes(size = -log10(P), color = R)) + # 添加散点/气泡
            scale_size_continuous(range = c(2, 6)) +
            cols4all::scale_color_continuous_c4a_div("kovesi.gn_wh_rd", mid = 0, limits = c(-1, 1), breaks = seq(-1, 1, 0.4), labels = seq(-1, 1, 0.4)) +
            cols4all::scale_fill_continuous_c4a_div("kovesi.gn_wh_rd", mid = 0, limits = c(-1, 1), breaks = seq(-1, 1, 0.4), labels = seq(-1, 1, 0.4)) +
            geom_hline(yintercept = c(-0, 0), color = "white", size = 0.5, lty = "dashed") +
            labs(
                  x = NULL,
                  y = "Spearman R",
                  # title = "KEGG Pathway Enrichment",
                  fill = "Spearman R",
                  color = "Spearman R",
                  size = bquote("-" ~ log[10] ~ "(P value)")
            ) +
            ## 负值坐标转换回来
            scale_y_continuous(
                  expand = expansion(add = c(0.1, 0.1)),
                  limits = c(-1, 1), breaks = seq(-1, 1, 0.4),
                  labels = seq(-1, 1, 0.4)
            )

      # 依次从下到上添加标签

      up_pathway <- length(which(dat$R > 0))
      down_pathway <- length(which(dat$R < 0))
      high <- nrow(dat)
      p2 <- p1 + theme_bw() + mytheme +
            geom_text(
                  data = dat[1:down_pathway, ], aes(x = Feature, y = 0.1, label = Feature),
                  hjust = 0, color = "black", size = 6
            ) +
            geom_text(
                  data = dat[(down_pathway + 1):high, ], aes(x = Feature, y = -0.1, label = Feature),
                  hjust = 1, color = "black", size = 6
            ) +
            scale_x_discrete(labels = NULL) +
            geom_text(x = 12, y = 10, label = "Up", size = 10, color = "#EE4E4A") +
            geom_text(x = 2, y = -10, label = "Down", size = 10, color = "#419A18") +
            guides(col = guide_colorbar(barwidth = 1.5, barheight = 10, ticks = TRUE))
      Canton::gs(p = p2, outdir = outdir, name = i, w = 12, h = 9)
}


mytheme <- theme(
      axis.text.x = element_text(hjust = 0.5, size = 20),
      axis.ticks.y = element_blank(), ## 删去y轴刻度线
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.line = element_line(size = 1),
      plot.margin = unit(c(1, 1, 1, 1), "cm"), # 画布边缘距离上(top)、右(right)、下(bottom)、左(left)
      plot.title = element_text(hjust = 0.5, size = 22),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20),
      legend.position = "right",
      legend.background = element_rect(fill = "transparent"),
      ## 删除网格线与纵坐标轴
      axis.line.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
)
