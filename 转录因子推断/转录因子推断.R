#https://saezlab.github.io/decoupleR/articles/pw_sc.html
#BiocManager::install("dorothea")
library(Seurat)
library(dorothea)
library(decoupleR)
library(dplyr)
library(tibble)
library(pheatmap)
library(tidyverse)
library(qs)


# DoRothEA 需要标准化的基因表达矩阵
#模块一：环境准备与转录因子活性推断 (Core Inference)
#这是最基础的步骤，目的是将基因表达矩阵（Gene Expression）转换为转录因子活性矩阵（TF Activity）。
# ==============================================================================
# 模块一：环境准备与转录因子活性推断
# 目的：利用 DoRothEA 数据库和 VIPER 算法，从基因表达数据推断转录因子(TF)的活性。
# ==============================================================================

# 1. 加载必要的 R 包
library(Seurat)
library(dorothea)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr) # 用于 pivot_wider

# 2. 准备输入数据
Epi = qread("sc.obj.qs")

# DoRothEA 需要标准化的基因表达矩阵作为输入。
# 这里提取 Seurat 对象中的 "RNA" assay 的 "data" slot (即 LogNormalize 后的数据)。
# 请确保 'Epi' 是你已经加载好的 Seurat 对象。
mat <- GetAssayData(Epi, assay = "RNA", slot = "data")

# 3. 加载 DoRothEA 调控网络数据库
# dorothea_hs 对应人类 (Human), dorothea_mm 对应小鼠 (Mouse)。
# 这里以人类为例。
data(dorothea_hs, package = "dorothea")

# 4. 筛选高置信度的调控关系
# DoRothEA 对 TF-Target 关系有 A/B/C/D/E 分级。
# A 为最高置信度。通常建议保留 A, B, C 级以保证结果的准确性。
regulons <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# 5. 运行 VIPER 算法推断 TF 活性
# run_viper 是推断 TF 活性的经典算法。
# .source: 调控因子列 (TF)
# .target: 靶基因列 (Target)
# .mor: 作用模式列 (Mode of Regulation, 激活或抑制)
# 注意：cores 参数控制并行计算的核心数，数据量大时建议调大 (如 cores = 4)。
tf_activities <- run_viper(mat, regulons, 
                           .source = 'tf', 
                           .target = 'target', 
                           .mor = 'mor', 
                           cores = 1) 

# 6. 数据格式转换
# VIPER 输出的是长格式 (Long format)，需要转换为 Seurat 可用的宽格式矩阵 (Wide matrix)。
# 转换后格式：行 = 转录因子 (TF), 列 = 细胞 (Cell)
tf_activities_mat <- tf_activities %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  as.matrix()

# 确保活性矩阵的细胞顺序与原始 Seurat 对象完全一致，防止元数据错位。
tf_activities_mat <- tf_activities_mat[, colnames(Epi)]

# 7. 将结果存回 Seurat 对象
# 创建一个新的 Assay 专门存放 TF 活性数据，命名为 'dorothea'。
# 这样后续可以使用 Seurat 的标准函数 (如 FindMarkers, VlnPlot) 对 TF 活性进行分析。
Epi[['dorothea']] <- CreateAssayObject(data = tf_activities_mat)
#警告: Layer counts isn't present in the assay object; returning NULL
# 将默认 Assay 切换为 'dorothea'，以便进行后续分析。
DefaultAssay(Epi) <- 'dorothea'

# 对活性数据进行标准化 (Scale)，这是绘制热图和降维前的必要步骤。
Epi <- ScaleData(Epi)

# (可选) 保存中间结果，防止 R Session 崩溃导致重新计算
# qs::qsave(Epi, "Epi_with_dorothea.qs")

###################
#模块二：Top 转录因子全局热图 (Global Visualization)
#使用 ComplexHeatmap 绘制发表级的热图，展示不同细胞亚群中特异性激活的转录因子。
# ==============================================================================
# 模块二：Top 转录因子全局热图
# 目的：寻找不同细胞群的 Marker TF，并绘制热图展示其活性强度。
# ==============================================================================

# 1. 加载绘图包
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# --- 用户配置区 ---
my_seurat_object <- Epi           # 你的 Seurat 对象
my_grouping_variable <- "celltype.1" # 你想要展示的分组列名 (如细胞类型)
my_top_n_tfs <- 10                # 每个分组展示前多少个显著 TF

# 2. 寻找差异激活的 TF (Marker TFs)
# 确保分组变量存在
if(!my_grouping_variable %in% colnames(my_seurat_object@meta.data)) { 
  stop("错误: 在 metadata 中找不到指定的分组变量。") 
}
Idents(my_seurat_object) <- my_grouping_variable
DefaultAssay(my_seurat_object) <- 'dorothea'

# 计算差异 TF
# only.pos = TRUE: 只关注活性升高的 TF
# logfc.threshold: 过滤掉变化微小的 TF
marker_tfs <- FindAllMarkers(my_seurat_object, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.15, verbose = FALSE)

# 提取每个分组 Top N 的 TF 名称
top_tfs_list <- marker_tfs %>% 
  group_by(cluster) %>% 
  slice_max(n = my_top_n_tfs, order_by = avg_log2FC) %>% 
  pull(gene) %>% 
  unique()

# 3. 准备绘图矩阵
# 计算选定 TF 在各分组的平均活性
avg_activity <- AverageExpression(my_seurat_object, features = top_tfs_list, assay = 'dorothea', slot = 'data', verbose = FALSE)$dorothea

# 对矩阵进行 Z-score 标准化 (Scale)，使热图颜色对比更鲜明
# t(scale(t(...))) 是对行 (TF) 进行标准化
matrix_activity_scaled <- t(scale(t(avg_activity)))

# 4. 构建热图注释 (Annotation)
final_cluster_names <- colnames(matrix_activity_scaled)
unique_cluster_levels <- sort(unique(final_cluster_names))
num_colors <- length(unique_cluster_levels)

# 生成颜色板
color_palette <- colorRampPalette(brewer.pal(min(num_colors, 8), "Set2"))(num_colors)
names(color_palette) <- unique_cluster_levels

# 创建顶部注释条
ha_top <- HeatmapAnnotation(
  df = data.frame(cluster = final_cluster_names),
  col = list(cluster = color_palette),
  show_annotation_name = FALSE,
  annotation_legend_param = list(title = "Cluster", direction = "horizontal", nrow = 1)
)

# 5. 定义热图颜色映射
# 使用 98% 分位数来处理极端值，避免颜色过饱和
activity_q <- quantile(abs(matrix_activity_scaled), 0.98)
col_activity <- colorRamp2(c(-activity_q, 0, activity_q), c("#0f2c4b", "white", "#a70b0b"))

# 6. 绘制热图
ht1 <- Heatmap(as.matrix(matrix_activity_scaled),
               name = "TF Activity", 
               col = col_activity,
               top_annotation = ha_top,
               column_title = "Top TF Regulon Activity", 
               column_title_gp = gpar(fontsize = 12, fontface = "bold"),
               cluster_columns = FALSE,  # 不对列聚类，保持分组顺序
               cluster_rows = TRUE,      # 对行 (TF) 进行聚类
               show_row_dend = TRUE,
               row_names_gp = gpar(fontsize = 9),
               column_names_rot = 45)

# 输出图像
draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "bottom")

# (可选) 保存为 PDF
# pdf("TF_Activity_Heatmap.pdf", width = 8, height = 10)
# draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "bottom")
# dev.off()

###################

#差异活性分析 (Differential Activity Analysis)
DefaultAssay(Epi) <- "dorothea"

# 假设我们要对比 肿瘤(Tumor) vs 正常(Normal)
# 先设置 idents
Idents(Epi) <- "sample.type" 

# 寻找差异 TF
# test.use = "wilcox" 是常用的非参数检验
da_tfs <- FindMarkers(Epi, 
                      ident.1 = "Tumor", 
                      ident.2 = "Normal", 
                      only.pos = FALSE, # 这里要设为 FALSE，因为我们要看下调的 TF
                      min.pct = 0.1, 
                      logfc.threshold = 0.1)

# 整理结果，添加显著性标签
da_tfs$TF <- rownames(da_tfs)
head(da_tfs %>% arrange(desc(avg_log2FC))) # 查看激活最显著的 TF
head(da_tfs %>% arrange(avg_log2FC))       # 查看抑制最显著的 TF


# --- 对比视角：RNA 表达 vs. 实际活性 ---

# 1. 切换回 RNA 看基因表达量
DefaultAssay(Epi) <- "RNA"
p1 <- FeaturePlot(Epi, features = "HIF1A", order = TRUE) + 
      ggtitle("HIF1A mRNA Expression")

# 2. 切换到 dorothea 看蛋白活性
DefaultAssay(Epi) <- "dorothea"
# 注意：VIPER 算出来的活性值有正负，建议用比较鲜明的红蓝配色
p2 <- FeaturePlot(Epi, features = "HIF1A", order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90') + 
      ggtitle("HIF1A Protein Activity") +
      scale_color_gradient2(low = "blue", mid = "white", high = "red")

# 3. 拼图对比
p1 | p2