Canton::using(tidyverse,data.table,clusterProfiler,ggrepel)
mRNA_CI_vs_HTN_DEG_Table = fread('结果/mRNA基因差异分析/Table_CI_vs_HTN_DEG_Table.csv')
mRNA_ICH_vs_HTN_DEG_Table = fread('结果/mRNA基因差异分析/Table_ICH_vs_HTN_DEG_Table.csv')

protein_CI_vs_HTN_DEG_Table = fread('结果/蛋白差异分析/Table_CI_vs_HTN_DEG_Table.csv')
protein_ICH_vs_HTN_DEG_Table = fread('结果/蛋白差异分析/Table_ICH_vs_HTN_DEG_Table.csv')

d1=mRNA_ICH_vs_HTN_DEG_Table %>% dplyr::rename(log2FC_trans=LFC) %>% dplyr::select(Feature,log2FC_trans)
d2=protein_ICH_vs_HTN_DEG_Table %>% dplyr::rename(log2FC_prot=LFC) %>% dplyr::select(Feature,log2FC_prot)
df = merge(d1,d2,by='Feature')

d1=mRNA_CI_vs_HTN_DEG_Table %>% dplyr::rename(log2FC_trans=LFC) %>% dplyr::select(Feature,log2FC_trans)
d2=protein_CI_vs_HTN_DEG_Table %>% dplyr::rename(log2FC_prot=LFC) %>% dplyr::select(Feature,log2FC_prot)
df = merge(d1,d2,by='Feature')


name="CI_vs_HTN"
# === 分组逻辑（使用你的阈值） ===
df <- df %>%
  mutate(
    trans_group = cut(log2FC_trans, 
                     breaks = c(-Inf, -0.585, 0.585, Inf), 
                     labels = c("mRNA_Low", "mRNA_Medium", "mRNA_High")),
    prot_group = cut(log2FC_prot, 
                    breaks = c(-Inf, -0.26, 0.26, Inf), 
                    labels = c("Protein_Low", "Protein_Medium", "Protein_High")),
    quadrant = paste(trans_group, prot_group, sep = " & ")
  )

# === 1. 计算每个象限的基因数量 ===
quadrant_counts <- df %>%
  group_by(quadrant) %>%
  summarise(count = n()) %>%
  ungroup()

# === 2. 创建象限中心坐标（关键！） ===
# 定义9个象限的中心坐标（x,y）
quadrant_centers <- data.frame(
  quadrant = c(
    "mRNA_Low & Protein_Low", "mRNA_Medium & Protein_Low", "mRNA_High & Protein_Low",
    "mRNA_Low & Protein_Medium", "mRNA_Medium & Protein_Medium", "mRNA_High & Protein_Medium",
    "mRNA_Low & Protein_High", "mRNA_Medium & Protein_High", "mRNA_High & Protein_High"
  ),
  x = c(-1.5, 0, 1.5,    -1.5, 0, 1.5, -1.5, 0, 1.5),  # x轴中心坐标
  y = c(-1.5, -1.5, -1.5, 0, 0, 0, 1.5, 1.5, 1.5)   # y轴中心坐标
)

# === 3. 合并数量和坐标 ===
quadrant_counts <- quadrant_counts %>%
  left_join(quadrant_centers, by = "quadrant")

Canton::hue('NPG')
color_map=c("#00468B", "#EC0000", "#41B43F", "#0099B3", "#925E9F","#AC002A", "#9E3B4F", "#B2811E","#DA988C")

# === 4. 选10个最显著的基因（按log2FC绝对值和排序） ===
used_gene = c("ADAM17",'TMPRSS5','PLAU','ADAMTS13')
label_df <- df %>%
  dplyr::filter(Feature %in% used_gene)

# === 5. 修正颜色映射 ===

names(color_map) <- c(
  "mRNA_Low & Protein_Low",
  "mRNA_Low & Protein_Medium",
  "mRNA_Low & Protein_High",
  "mRNA_Medium & Protein_Low",
  "mRNA_Medium & Protein_Medium",
  "mRNA_Medium & Protein_High",
  "mRNA_High & Protein_Low",
  "mRNA_High & Protein_Medium",
  "mRNA_High & Protein_High"
)

# === 6. 画图（含象限计数标注 + 基因标注） ===
p=ggplot(df, aes(x = log2FC_trans, y = log2FC_prot, color = quadrant)) +
  geom_point(size = 2.5, alpha = 0.7) +  # 所有点
  # === 关键：添加象限计数标注（每个象限中心） ===
  geom_text(
    data = quadrant_counts,
    aes(x = x, y = y, label = count),
    size = 6, 
    fontface = "bold",
    color = "black",
    check_overlap = TRUE  # 避免重叠
  ) +
  # === 保留基因标注 ===
 geom_text_repel(
    data = label_df, 
    aes(label = Feature),
    size = 7, 
    box.padding = 0.5, 
    point.padding = 0.5,
    # color = "black",
    segment.color = "black",  # 箭头颜色
    segment.size = 1,        # 箭头线粗细
    segment.alpha = 1,       # 箭头透明度
    force = 2.0,               # 避免重叠力度
    arrow = arrow(
      length = unit(0.1, "npc"),  # 箭头长度（小单位）
      type = "closed"              # "closed"表示封闭箭头，"open"表示开口箭头
    ),
    seed = 123                 # 固定随机种子保证位置稳定
  ) +
  labs(
    title = "",
    x = "Log2FC (Transcriptome)",
    y = "Log2FC (Proteome)",
    color = "Quadrant"
  ) +
  scale_color_manual(values = color_map) +
  # 参考线更新为你的阈值
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = c(-0.26, 0.26), linetype = "dashed", color = "gray30") +theme_bw(base_size = 22) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(linetype = "dotted")
  )

# 保存高清图
Canton::gs(p, name = name, outdir = base::getwd(), format = "pdf", w = 14, h = 9)
