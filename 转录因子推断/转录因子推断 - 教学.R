#https://saezlab.github.io/decoupleR/articles/pw_sc.html
# ==============================================================================
# 教程名称：基于 DoRothEA 的单细胞转录因子活性推断与可视化
# 适用人群：R 语言与单细胞分析新手
# ==============================================================================
# 如果尚未安装 dorothea，请取消下一行的注释并运行（仅需运行一次）
# BiocManager::install("dorothea")

# 1. Seurat: 单细胞分析的标准工具包，负责管理数据、降维聚类等。
library(Seurat)
# 2. dorothea: 包含人类和小鼠的“转录因子-靶基因”调控网络数据库。
library(dorothea)
# 3. decoupleR: 提供算法（如 run_viper），用于从基因表达推断生物学活动。
library(decoupleR)
# 4. dplyr, tibble, tidyr, tidyverse: 这一组叫 "Tidyverse" 系列，是 R 语言中处理表格（筛选、排序、变形）的神器。
library(dplyr)
library(tibble)
library(tidyr) # 专门用于长宽数据格式转换
library(tidyverse)
# 5. pheatmap / ComplexHeatmap: 专门用于画热图的包。
library(pheatmap)
# 6. qs: 用于快速读取和保存大型数据文件（比传统的 readRDS 快很多）。
library(qs)


# ==============================================================================
# 模块一：环境准备与转录因子活性推断 (Core Inference)
# 核心逻辑：
# 基因表达量 (Gene Expression) != 蛋白活性 (Protein Activity)。
# 我们利用“靶基因的表达变化”来反推“转录因子的活性”。
# ==============================================================================

# 1. 准备输入数据
# qread 是读取 .qs 格式文件的函数。
# "sc.obj.qs" 是你的单细胞数据文件，读取进来后赋值给变量 Epi。
# Epi 是一个 Seurat 对象，它是存储单细胞所有信息（矩阵、分组、降维结果）的容器。
Epi = qread("sc.obj.qs")

# 2. 提取基因表达矩阵
# DoRothEA 算法需要“标准化后”的数据，而不是原始的测序计数（Counts）。
# GetAssayData: 从 Seurat 对象中取数据的函数。
# assay = "RNA": 指定取 RNA 测序数据。
# slot = "data": 这里存放的是经过 LogNormalize（对数标准化）的数据，消除了测序深度的影响。
mat <- GetAssayData(Epi, assay = "RNA", slot = "data")

# 3. 加载 DoRothEA 调控网络数据库
# dorothea_hs 是该包内置的一个数据框，存储了人类（Homo sapiens）的调控关系。
# 如果是小鼠数据，这里应该用 data(dorothea_mm, package = "dorothea")。
data(dorothea_hs, package = "dorothea")

# 4. 筛选高置信度的调控关系
# DoRothEA 数据库里的关系分为 A, B, C, D, E 五个等级。
# A 级最可信（有文献支持），E 级最不可信（仅基于预测）。
# 为了计算结果准确，我们只保留 A, B, C 级。
# %>% 是管道符，意思是把左边的数据传给右边的函数。
# filter: 筛选行。
regulons <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# 5. 运行 VIPER 算法推断 TF 活性
# run_viper 是核心计算步骤。它会查看每个 TF 对应的靶基因在你的数据(mat)里表达高不高。
# 输入参数解释：
# - mat: 你的基因表达矩阵。
# - regulons: 刚才筛选好的调控规则表。
# - .source: 指明 regulons 表里哪一列是转录因子名字 ('tf')。
# - .target: 指明哪一列是靶基因名字 ('target')。
# - .mor: Mode of Regulation (调控模式)，指明是激活(+1)还是抑制(-1)。
# - cores: 电脑 CPU 核心数。设为 1 最安全；如果数据很大且电脑配置好，可设为 4 加速。
tf_activities <- run_viper(mat, regulons, 
                           .source = 'tf', 
                           .target = 'target', 
                           .mor = 'mor', 
                           cores = 1) 

# 此时 tf_activities 是一个“长格式”表格（Long Format）：
# 每一行记录：某一个 TF 在某一个 细胞 中的活性得分。

# 6. 数据格式转换
# Seurat 对象要求数据必须是“宽格式矩阵”：行是特征(TF)，列是细胞(Cell)。
# pivot_wider: 将长表变宽。
# - id_cols: 保持 TF 名字作为行标识。
# - names_from: 将细胞 ID ('condition') 变成列名。
# - values_from: 将活性得分 ('score') 填入表格中间。
tf_activities_mat <- tf_activities %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>% # 把 TF 名字这一列变成真正的“行名”
  as.matrix() # 转换为纯数字矩阵，以便后续计算

# 重要检查：确保活性矩阵的列（细胞）顺序，和原始 Seurat 对象完全一致。
# 否则后续分析会把细胞 A 的活性安到细胞 B 头上，导致严重错误。
tf_activities_mat <- tf_activities_mat[, colnames(Epi)]

# 7. 将结果存回 Seurat 对象
# Seurat 对象可以存放多种类型的数据（称为 Assay）。
# 默认只有 'RNA' Assay。现在我们新建一个名为 'dorothea' 的 Assay 来存 TF 活性。
# 这样 Seurat 就会把 TF 活性当成一种“类似于基因表达”的数据来处理。
Epi[['dorothea']] <- CreateAssayObject(data = tf_activities_mat)

# 警告解释：CreateAssayObject 可能会提示 "Layer counts isn't present"。
# 这是正常的，因为 TF 活性是计算出来的分数，不是测序得到的整数 Counts，忽略即可。

# 切换默认分析对象
# 这一步之后，你运行 FindMarkers 或 VlnPlot 时，默认都是针对 TF 活性进行分析，而不是基因表达。
DefaultAssay(Epi) <- 'dorothea'

# 标准化 (Scale)
# ScaleData 的作用是将数据转换为 Z-score（均值为 0，标准差为 1）。
# 目的：消除不同 TF 基础活性值的差异，让它们可以在同一张热图上进行比较。
Epi <- ScaleData(Epi)

# (可选) 保存中间结果
# 这是一个好习惯，防止后面画图出错导致 R 崩溃，需要重新跑上面的计算（很耗时）。
# qs::qsave(Epi, "Epi_with_dorothea.qs")

# ==============================================================================
# 模块二：Top 转录因子全局热图 (Global Visualization)
# 核心逻辑：
# 1. 既然算出活性了，我们要看哪些 TF 在哪种细胞里最活跃。
# 2. 找出每种细胞特有的 TF (Marker TFs)。
# 3. 画热图展示这些 TF 在所有细胞群中的表现。
# ==============================================================================

# 1. 加载高级绘图包
library(ComplexHeatmap) # 绘制复杂热图的标准工具
library(RColorBrewer)   # 提供好看的配色板
library(circlize)       # 用于更精细的颜色映射控制

# --- 用户配置区 (新手请关注这里) ---
# 这里的变量可以根据你的实际需求修改
my_seurat_object <- Epi           # 你的 Seurat 对象变量名
my_grouping_variable <- "celltype.1" # *重要*：这是你的细胞分类列名，必须在 meta.data 里存在
my_top_n_tfs <- 10                # 每个细胞群展示前多少个最显著的 TF

# 2. 寻找差异激活的 TF (Marker TFs)
# 安全检查：确认你写的分类列名真的存在
if(!my_grouping_variable %in% colnames(my_seurat_object@meta.data)) { 
  stop("错误: 在 metadata 中找不到指定的分组变量，请检查 'my_grouping_variable' 的拼写。") 
}

# 设置细胞身份 (Idents)
# 告诉 Seurat，我们要按照 "celltype.1" 来把细胞分组。
Idents(my_seurat_object) <- my_grouping_variable
# 再次确认为 'dorothea' assay
DefaultAssay(my_seurat_object) <- 'dorothea'

# 计算差异 TF (FindAllMarkers)
# 这是一个统计检验过程，对比“这一群细胞”和“其他所有细胞”的 TF 活性。
# - only.pos = TRUE: 我们只关心活性“升高”的 TF，不关心降低的。
# - min.pct = 0.1: 加速计算，忽略那些在极少数细胞里有活性的 TF。
# - logfc.threshold = 0.15: 过滤掉变化幅度太小的 TF。
marker_tfs <- FindAllMarkers(my_seurat_object, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.15, verbose = FALSE)

# 提取 Top N 列表
# 这里的逻辑是：按分组(cluster)打包 -> 按差异倍数(avg_log2FC)排序 -> 取前 10 个 -> 提取名字 -> 去重
top_tfs_list <- marker_tfs %>% 
  group_by(cluster) %>% 
  slice_max(n = my_top_n_tfs, order_by = avg_log2FC) %>% 
  pull(gene) %>% 
  unique()

# 3. 准备绘图矩阵
# 热图不能画几万个细胞，那样看不清。
# AverageExpression: 计算每种细胞类型中 TF 活性的“平均值”。
# 结果是一个矩阵：行是 TF，列是细胞类型。
avg_activity <- AverageExpression(my_seurat_object, features = top_tfs_list, assay = 'dorothea', slot = 'data', verbose = FALSE)$dorothea

# 再次 Z-score 标准化
# 为什么前面 ScaleData 了一次，这里还要 Scale？
# 前面的 ScaleData 是针对“单细胞”层面的。这里是针对“平均值矩阵”层面的。
# t() 是转置（行列互换）。因为 scale() 函数默认只对列进行操作，我们要对行(TF)标准化，所以要转置两次。
# 结果：让每个 TF 在不同细胞群之间的对比更强烈（高为红，低为蓝）。
matrix_activity_scaled <- t(scale(t(avg_activity)))

# 4. 构建热图注释 (Annotation)
# 这一步是为了给热图上方加一条彩色的条，用来标记每一列属于哪个细胞群。
final_cluster_names <- colnames(matrix_activity_scaled)
unique_cluster_levels <- sort(unique(final_cluster_names))
num_colors <- length(unique_cluster_levels)

# 生成颜色板
# colorRampPalette: 创建渐变色生成器。
# brewer.pal: 从 Set2 色板取色。
# 结果：为每个细胞群分配一个唯一的颜色。
color_palette <- colorRampPalette(brewer.pal(min(num_colors, 8), "Set2"))(num_colors)
names(color_palette) <- unique_cluster_levels

# 创建顶部注释对象
ha_top <- HeatmapAnnotation(
  df = data.frame(cluster = final_cluster_names),
  col = list(cluster = color_palette),
  show_annotation_name = FALSE, # 不显示注释条的名字，为了美观
  annotation_legend_param = list(title = "Cluster", direction = "horizontal", nrow = 1) # 图例设置
)

# 5. 定义热图颜色映射
# 这是一个处理“异常值”的技巧。
# 如果有一个极大值（比如 10），其他值都是 1 或 2，那么热图大部分会是白色的，看不清差异。
# quantile(..., 0.98): 找出绝对值排在 98% 位置的那个数。
activity_q <- quantile(abs(matrix_activity_scaled), 0.98)

# colorRamp2: 定义颜色映射规则。
# 负值(蓝) -> 零(白) -> 正值(红)
# 大于 98% 分位数的极值，颜色将被“截断”在这个最深色，保证整体对比度。
col_activity <- colorRamp2(c(-activity_q, 0, activity_q), c("#0f2c4b", "white", "#a70b0b"))

# 6. 绘制热图
# Heatmap: ComplexHeatmap 包的核心函数。
ht1 <- Heatmap(as.matrix(matrix_activity_scaled),
               name = "TF Activity",          # 图例的标题
               col = col_activity,            # 使用刚才定义的红白蓝颜色
               top_annotation = ha_top,       # 加上顶部的彩色注释条
               column_title = "Top TF Regulon Activity", 
               column_title_gp = gpar(fontsize = 12, fontface = "bold"),
               cluster_columns = FALSE,       # FALSE: 不重新排序细胞群，按默认顺序展示
               cluster_rows = TRUE,           # TRUE: 将表达模式相似的 TF 聚在一起
               show_row_dend = TRUE,          # 显示左侧的聚类树
               row_names_gp = gpar(fontsize = 9), # TF 名字的字体大小
               column_names_rot = 45)         # 列名旋转 45 度，防止重叠

# 输出图像到 R 的绘图窗口
draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "bottom")

# (可选) 保存为 PDF 文件
# 取消注释以下三行可以将图片保存到电脑
# pdf("TF_Activity_Heatmap.pdf", width = 8, height = 10)
# draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "bottom")
# dev.off() # 关闭绘图设备，完成保存