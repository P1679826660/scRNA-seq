library(Seurat)
library(qs)
library(harmony)
library(ggplot2)


sc.obj = qread("GSE149614.Tumor.qs")
#qsave(sc.obj,"GSE149614.Tumor.qs")

#head(sc.obj@meta.data)
#table(sc.obj$sample)

sc.obj[["percent.mt"]] <- PercentageFeatureSet(object = sc.obj, pattern = "^MT-")
VlnPlot(object = sc.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by ="sample" , ncol = 3)
sc.obj <- subset(x = sc.obj, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15) 
VlnPlot(object = sc.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by ="sample" , ncol = 3)

sc.obj <- NormalizeData(object = sc.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc.obj <- FindVariableFeatures(object = sc.obj, selection.method = "vst", nfeatures = 2000) 
sc.obj<- ScaleData(sc.obj)   
gc()



sc.obj <- RunPCA(sc.obj, pc.genes = VariableFeatures(object = sc.obj))
sc.obj <- RunHarmony(sc.obj, group.by.vars = "sample")
ElbowPlot(sc.obj)

pcSelect=30    #这个地方根据上个PC纳入
sc.obj <- FindNeighbors(sc.obj, reduction = "harmony", dims = 1:pcSelect)             #计算邻接距离
sc.obj <- FindClusters(object = sc.obj, resolution = 1)                  
sc.obj <- RunUMAP(sc.obj, reduction = "harmony", dims = 1:pcSelect)
DimPlot(sc.obj, reduction = "umap", label = TRUE, label.size = 2)       
#######################################################################
#去双细胞 001
library(scDblFinder)
name = sc.obj@meta.data$sample
sce <- as.SingleCellExperiment(sc.obj)
sce <- scDblFinder(sce,samples = name) #, dbr = 0.08 设置8%

sc.obj$scDblFinder.score <- sce$scDblFinder.score
sc.obj$scDblFinder.class <- sce$scDblFinder.class

rm(sce)
table(sc.obj@meta.data$scDblFinder.class)

DimPlot(sc.obj, group.by = "scDblFinder.class", cols = c("doublet" = "red", "singlet" = "grey"))

sc.obj <- subset(sc.obj, subset = scDblFinder.class == "singlet")

#去RNA污染 002
library(Seurat)
library(decontX)
library(ggplot2)
library(viridis) # 用于 scale_color_viridis_c

# 加载您的 Seurat 对象
# 从 Seurat 对象的 RNA assay 中提取原始计数矩阵
# 注意：decontX 需要的是未经标准化的原始 count 数据
counts <-  GetAssayData(sc.obj, assay = "RNA", layer = "count")

# 检查矩阵维度和行列名，确保是 基因 x 细胞 格式
# dim(counts)
# head(colnames(counts))
# head(rownames(counts))

# 运行 decontX，结果将储存在一个新的变量中
# 这一步可能会消耗一些时间和内存，具体取决于您的细胞数量
decontX_results <- decontX(counts)

# 将 decontX 计算出的污染分数 ('contamination') 添加到 metadata
sc.obj$Contamination <- decontX_results$contamination

# 查看 metadata，确认新列已添加成功
head(sc.obj@meta.data)

# 在 UMAP 图上可视化 Contamination 分数
FeaturePlot(sc.obj, 
            features = 'Contamination', 
            raster = FALSE) + # 如果细胞数非常多，建议设置为 TRUE
  scale_color_viridis_c() +   # 使用 viridis 色彩方案，更美观
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle("Cell Contamination Score by decontX") # 添加标题
############
# 加载必要的包
library(dplyr)
library(ggplot2)

# ---
# 步骤 1: 提取并排序污染分数
# ---

# 从 Seurat 对象的 metadata 中提取 Contamination 分数
contamination_scores <- sc.obj@meta.data$Contamination

# 从高到低排序
sorted_scores <- sort(contamination_scores, decreasing = TRUE)

# 创建一个用于绘图的数据框
# X轴是细胞的排名 (1, 2, 3, ...)
# Y轴是排序后的污染分数
elbow_df <- data.frame(
  cell_rank = 1:length(sorted_scores),
  contamination = sorted_scores
)

# ---
# 步骤 2: 绘制拐点图
# ---

# 使用 ggplot2 绘制精美的拐点图
elbow_plot <- ggplot(elbow_df, aes(x = cell_rank, y = contamination)) +
  geom_point(size = 1.5, alpha = 0.5, color = "dodgerblue4") + # 用点图更清晰
  geom_line(color = "dodgerblue4", alpha = 0.8) +             # 用线连接趋势
  theme_bw(base_size = 14) +
  labs(
    title = "Contamination Score Elbow Plot",
    x = "Cell Rank (sorted by contamination)",
    y = "Contamination Score"
  ) +
  theme(panel.grid.minor = element_blank()) # 精简背景

# 显示图表
print(elbow_plot)


# ---
# 步骤 3: 使用“距离法”自动寻找拐点
# ---
library(dplyr) # 确保 dplyr 已加载

# 1. 数据归一化 (将 X 和 Y 缩放到 [0, 1] 区间)
# 这是必须的，因为 X 和 Y 轴的单位和范围完全不同
df_norm <- elbow_df %>%
  mutate(
    rank_norm = (cell_rank - min(cell_rank)) / (max(cell_rank) - min(cell_rank)),
    cont_norm = (contamination - min(contamination)) / (max(contamination) - min(contamination))
  )

# 2. 计算每个点到“首-尾”连线的距离
# 归一化后，首点是 (0, 1)，尾点是 (1, 0)
# 这条线的方程是 y = -x + 1，或者 x + y - 1 = 0
# 点 (x0, y0) 到线 Ax + By + C = 0 的距离公式为：
# d = |A*x0 + B*y0 + C| / sqrt(A^2 + B^2)
# 在这里, A=1, B=1, C=-1, x0=rank_norm, y0=cont_norm
df_norm <- df_norm %>%
  mutate(
    distance = abs(rank_norm + cont_norm - 1) / sqrt(1^2 + 1^2)
  )

# 3. 找到距离最远的点（即拐点）
# which.max() 会返回第一个达到最大值的索引
elbow_index <- which.max(df_norm$distance)

# 4. 获取原始数据中的拐点信息
auto_rank <- elbow_df$cell_rank[elbow_index]
auto_threshold <- elbow_df$contamination[elbow_index]

# 5. 自动输出结果
cat("--- 自动拐点检测结果 ---\n")
cat("自动检测到的拐点位于细胞排名:", auto_rank, "\n")
cat("对应的建议阈值 (Contamination Score) 大约是:", round(auto_threshold, 3), "\n")


# ---
# 步骤 4: 可视化自动确定的阈值
# ---

print(
  elbow_plot + 
    # 使用自动计算的阈值
    geom_hline(yintercept = auto_threshold, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = auto_rank, linetype = "dashed", color = "red", size = 1) +
    
    # 在标题中显示自动计算的值
    ggtitle("Elbow Plot with Automated Threshold",
            subtitle = paste("Automated Rank:", auto_rank, 
                             "| Automated Threshold:", round(auto_threshold, 3)))
)

#############
# 设置过滤阈值
# 注意：这个阈值需要根据您自己的数据进行评估和调整。
# 官方没有给出固定标准，这里参考 SoupX 和经验选择了 0.2，请务必结合您的可视化结果进行判断。
contamination_threshold <- round(auto_threshold, 3)

# 使用 subset 函数根据阈值筛选细胞
sc.obj <- subset(sc.obj, subset = Contamination < contamination_threshold)

# 查看过滤前后的细胞数量
cat("过滤后细胞数:", ncol(sc.obj), "\n")

##########################################

genes_to_check <- c("PTPRC","IGKC",                                      # 免疫
                    "EPCAM", "KRT8", "KRT18","ALB" ,"SERPINA1","HNF4A",                    # 上皮
                    "PECAM1","PLVAP","VWF", "COL1A1","DCN"
)                   # 基质细胞 (前三个内皮细胞  COL1A1-成纤维细胞)


DotPlot(sc.obj, features = genes_to_check,cols = c("RdYlBu")) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank())


##################################################
#提取肿瘤+一正常
sc.obj.gl=sc.obj
sc.obj= subset(sc.obj,subset=seurat_clusters %in% c(1,2,3,7,12,13,14,15,21,23,6))

rm(list = ls()[!ls() %in% "sc.obj"])

Epi= subset(sc.obj,subset=seurat_clusters %in% c(7,16,20))



qsave(Epi,"正常上皮.qs")















#############################
genes_to_check <- c(
  # --- 细胞周期阻滞核心通路 ---
  "CDKN2A",  # p16(INK4a)，最经典的衰老marker之一，作用于RB通路
  "CDKN1A",  # p21(Cip1)，另一个核心marker，p53通路的关键下
  # 其他分泌蛋白
  "GDF15",
  # --- 细胞内生物标志物 ---
  "GLB1",      # 编码衰老相关β-半乳糖苷酶(SA-β-gal)的基因，经典染色实验的靶点
  "LMNB1",      # Lamin B1，一个经典的【下调】marker，其表达丢失是衰老的标志
  "MKI67", "PCNA"
)

colnames(sc.obj@meta.data)

DotPlot(sc.obj, features = genes_to_check,cols = c("RdYlBu")) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank())




#######





# 步骤 1: 从Seurat的内置数据 cc.genes 中加载S期和G2/M期的基因列表
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# s.genes 和 g2m.genes 是Seurat自带的基因列表
sc.obj <- CellCycleScoring(sc.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(sc.obj, group.by = "Phase")
# 然后看看您的“功能失调群”里，S和G2M期的细胞比例是不是很高
table(sc.obj$senescence_stage, sc.obj$Phase)

# 确保库已加载
library(dplyr)
library(ggplot2)

# 1. 一步完成数据准备和增殖指数计算
proliferative_index_df <- sc.obj@meta.data %>%
  select(senescence_stage, Phase) %>%
  na.omit() %>%
  group_by(senescence_stage) %>%
  summarise(
    proliferating_cells = sum(Phase %in% c("S", "G2M")),
    total_cells = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    proliferative_pct = (proliferating_cells / total_cells) * 100,
    # 确保绘图时的X轴顺序正确
    senescence_stage = factor(senescence_stage, levels = c("Normal", "Transition", "Aged"))
  )

# 打印计算结果的表格
print(proliferative_index_df)

# 2. 绘制紧凑版的柱状图
prolif_plot <- ggplot(proliferative_index_df, aes(x = senescence_stage, y = proliferative_pct, fill = senescence_stage)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(aes(label = paste0(round(proliferative_pct, 1), "%")), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Normal" = "grey", "Transition" = "#FDBF6F", "Aged" = "red")) +
  labs(
    title = "增殖指数揭示'Aged'群体的功能失调状态",
    subtitle = "尽管衰老分数高，'Aged'群体的增殖活性并未被有效抑制",
    x = "细胞衰老阶段 (GMM分类)",
    y = "增殖指数 (%) [S + G2M 期细胞比例]"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  ) +
  ylim(0, max(proliferative_index_df$proliferative_pct, na.rm = TRUE) * 1.2)

# 打印图表
print(prolif_plot)
#####################################
# ==============================================================================
# 绘制堆叠百分比柱状图，展示细胞状态比例随肿瘤分期的动态变化
# ==============================================================================

# 确保库已加载
library(ggplot2)
library(dplyr)

# ---
# 步骤 1: 创建一个【正确的】长格式数据框
# ---

# a. 使用 table() 计算每个分期和每个状态组合的细胞数
stage_table <- table(sc.obj@meta.data$stage, sc.obj@meta.data$senescence_stage)

# b. 将 table 的输出转换为一个数据框
#    table()生成的数据框默认列名为 Var1, Var2, Freq
stage_data <- data.frame(stage_table)

# c. 【关键修正】为列重命名，使其更具可读性
colnames(stage_data) <- c("Stage", "Status", "Count")

# 检查一下数据框的结构，您会看到它已经是长格式了
print(head(stage_data))

# ---
# 步骤 2: 【关键修正】确保分期和状态按正确的生物学顺序排列
# ---
# 这能保证您X轴的顺序是 I -> II -> IIIA -> ...
stage_data$Stage <- factor(stage_data$Stage, levels = c("I", "II", "IIIA", "IIIB", "IV"))

# 这能保证图例和堆叠的顺序是 Normal -> Transition -> Aged
stage_data$Status <- factor(stage_data$Status, levels = c("Normal", "Transition", "Aged"))


# ---
# 步骤 3: 绘制堆叠百分比柱状图
# ---
# pivot_longer() 这一步被完全删除了，因为不再需要

proportions_plot <- ggplot(stage_data, aes(x = Stage, y = Count, fill = Status)) +
  
  # 【关键】使用 position = "fill" 来自动计算并绘制百分比堆叠
  geom_bar(stat = "identity", position = "fill", color = "black") +
  
  # 将Y轴的标签格式化为百分比
  scale_y_continuous(labels = scales::percent) +
  
  # 自定义颜色，确保与您的UMAP图一致
  scale_fill_manual(values = c("Aged" = "red", "Transition" = "yellow", "Normal" = "grey")) +
  
  # 美化图表
  labs(
    title = "细胞状态比例随肿瘤分期的动态变化",
    subtitle = "在II期观察到'功能失调/增殖'细胞比例的显著峰值",
    x = "肿瘤分期 (Stage)",
    y = "细胞比例 (%)",
    fill = "细胞状态"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) # 让X轴标签倾斜，避免重叠
  )

# 打印最终的图表
print(proportions_plot)


###############################################
aged_cells_obj <- subset(sc.obj, subset = senescence_stage == "Aged")

# --- 1. Define Parameters ---
# Define the number of principal components based on the ElbowPlot
pcSelect <- 20

# --- 2. Main Processing Pipe ---
# This single chain of commands performs all data processing steps
aged_cells_obj <- aged_cells_obj %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(.)) %>%
  harmony::RunHarmony(group.by.vars = "sample", reduction.save = "harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:pcSelect) %>%
  FindClusters(resolution = 1.5) %>%
  RunUMAP(reduction = "harmony", dims = 1:pcSelect)

gc()

ElbowPlot(aged_cells_obj, ndims = 50, reduction = "harmony")
DimPlot(aged_cells_obj, reduction = "umap", label = TRUE)

genes_to_check <- c(
  # --- 细胞周期阻滞核心通路 ---
  "CDKN2A",  # p16(INK4a)，最经典的衰老marker之一，作用于RB通路
  "CDKN1A",  # p21(Cip1)，另一个核心marker，p53通路的关键下
  # 其他分泌蛋白
  "GDF15",
  # --- 细胞内生物标志物 ---
  "GLB1",      # 编码衰老相关β-半乳糖苷酶(SA-β-gal)的基因，经典染色实验的靶点
  "LMNB1",      # Lamin B1，一个经典的【下调】marker，其表达丢失是衰老的标志
  "MKI67", "TOP2A","PCNA"
)


DotPlot(aged_cells_obj, features = genes_to_check,cols = c("RdYlBu"),group.by = "seurat_clusters") + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank())

table(aged_cells_obj$seurat_clusters,aged_cells_obj$sample)

#qsave(aged_cells_obj,"aged_cells_obj.qs")

# ==============================================================================
# 完整流程: 计算并可视化每个 Seurat Cluster 的细胞周期比例
# ==============================================================================

# ---
# 步骤 0: 加载必要的库
# ---
library(Seurat)
library(dplyr)
library(ggplot2)

# ---
# 步骤 1: 对 Seurat 对象进行细胞周期评分
# ---
cat("--- 步骤 1: 正在进行细胞周期评分 ---\n")

# a. 从Seurat的内置数据中加载S期和G2/M期的基因列表
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# b. 运行 CellCycleScoring 函数
#    这会在 sc.obj@meta.data 中添加三列: S.Score, G2M.Score, 和 Phase
aged_cells_obj <- CellCycleScoring(
  aged_cells_obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE # 保持我们现有的Idents不变
)

# c. 验证：查看元数据的前几行，确认新列已添加
cat("细胞周期评分完成！元数据已更新，预览如下:\n")
print(head(aged_cells_obj@meta.data[, c("seurat_clusters", "Phase", "S.Score", "G2M.Score")]))


# ---
# 步骤 2: 定量分析 - 计算每个 seurat_cluster 的细胞周期比例
# ---
cat("\n--- 步骤 2: 正在计算每个Cluster的细胞周期比例 ---\n")

# 使用dplyr进行分组和汇总计算
cycle_proportions_by_cluster <- aged_cells_obj@meta.data %>%
  # 按 seurat_clusters 分组
  group_by(seurat_clusters) %>%
  # 统计每个cluster的总细胞数，以及G1, S, G2M期的细胞数
  summarise(
    Total = n(),
    G1_count = sum(Phase == "G1"),
    S_count = sum(Phase == "S"),
    G2M_count = sum(Phase == "G2M"),
    .groups = 'drop'
  ) %>%
  # 计算增殖指数 (S+G2M) 的百分比
  mutate(
    Proliferative_Index_Pct = ((S_count + G2M_count) / Total) * 100
  )

# 打印计算结果的表格
cat("每个Seurat Cluster的细胞周期定量分析结果:\n")
print(cycle_proportions_by_cluster)


# ---
# 步骤 3: 可视化 - 绘制堆叠百分比柱状图
# ---
cat("\n--- 步骤 3: 正在绘制细胞周期比例图 ---\n")

# a. 准备用于ggplot2的“长格式”数据
plot_data <- aged_cells_obj@meta.data %>%
  # 按cluster和Phase分组，并计数
  group_by(seurat_clusters, Phase) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  # 确保Phase的顺序是 G1 -> S -> G2M
  mutate(Phase = factor(Phase, levels = c("G1", "S", "G2M")))

# b. 绘制图表
cycle_plot_by_cluster <- ggplot(plot_data, aes(x = seurat_clusters, y = Count, fill = Phase)) +
  # 使用 position = "fill" 来自动计算并绘制百分比堆叠
  geom_bar(stat = "identity", position = "fill", color = "black") +
  
  # 将Y轴的标签格式化为百分比
  scale_y_continuous(labels = scales::percent) +
  
  # 自定义颜色
  scale_fill_manual(values = c("G1" = "grey", "S" = "blue", "G2M" = "red")) +
  
  # 美化图表
  labs(
    title = "各 Seurat Cluster 的细胞周期组成",
    subtitle = "展示了不同细胞亚群的增殖状态异质性",
    x = "Seurat Cluster",
    y = "细胞比例 (%)",
    fill = "细胞周期阶段"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) # 让X轴标签倾斜
  )

# 打印最终的图表
print(cycle_plot_by_cluster)



# 使用 mutate 和 case_when 添加新列
aged_cells_obj@meta.data <- aged_cells_obj@meta.data %>%
  mutate(proliferation = case_when(
    Phase %in% c('S', 'G2M') ~ 'prolif',
    Phase == 'G1' ~ 'no.prolif',
    TRUE ~ NA_character_  # 对于其他可能的值（如果有的话），标记为NA
  ))

# --- 验证一下结果 ---
# 查看新列是否添加成功
head(aged_cells_obj@meta.data[, c("Phase", "proliferation")])

# 统计一下各个状态的细胞数
table(aged_cells_obj$proliferation)
##########
# 加载需要的库 (如果还未加载)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. 提取绘图所需的数据
plot_data <- aged_cells_obj@meta.data

# --- 图1: Proliferation 状态比例柱状图 ---
p1_bar <- ggplot(plot_data, aes(x = "Status", fill = proliferation)) +
  geom_bar(position = "fill", stat = "count") + # position="fill" 自动计算百分比
  scale_y_continuous(labels = scales::percent) + # Y轴显示为百分比
  labs(
    title = "Proportion of Proliferation Status",
    x = "",
    y = "Proportion",
    fill = "Proliferation"
  ) +
  theme_minimal() +
  # 添加百分比标签
  geom_text(aes(label = scales::percent(..count.. / sum(..count..), 1)),
            stat = "count",
            position = position_fill(vjust = 0.5))

# --- 图2: Phase 细胞周期比例柱状图 ---
p2_bar <- ggplot(plot_data, aes(x = "Phase", fill = Phase)) +
  geom_bar(position = "fill", stat = "count") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportion of Cell Cycle Phase",
    x = "",
    y = "Proportion",
    fill = "Cell Cycle Phase"
  ) +
  theme_minimal() +
  # 添加百分比标签
  geom_text(aes(label = scales::percent(..count.. / sum(..count..), 1)),
            stat = "count",
            position = position_fill(vjust = 0.5))

# 使用 patchwork 将两个图并排显示
p1_bar + p2_bar
############################
VlnPlot(aged_cells_obj, features = "SenePy", group.by = "proliferation")

library(openxlsx)
library(tidyverse)

# 2. 设置身份并查找差异基因
deg_results <- FindMarkers(
  aged_cells_obj,
  ident.1 = "prolif",
  ident.2 = "no.prolif",
  logfc.threshold = 0,
  min.pct = 0.1,
  verbose = FALSE,group.by = "proliferation"
)

# 3. 筛选并保存到Excel
filtered_degs <- deg_results %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC))

write.xlsx(filtered_degs, file = "DEG_Proliferation_vs_NoProliferation_compact.xlsx")

#########
library(Seurat)
library(dplyr)
library(tibble)
library(openxlsx)

# 1. 准备实验组：从 aged_cells_obj 中提取 prolif 细胞
prolif_cells <- subset(aged_cells_obj, subset = proliferation == 'no.prolif')
# 添加分组标签
prolif_cells$deg_group <- "Aged_Prolif"

# 2. 准备对照组：创建 normal_cells_obj
normal_cells_obj <- subset(sc.obj, subset = senescence_stage == "Normal")
# 添加分组标签
normal_cells_obj$deg_group <- "Normal"

# 3. 合并两个对象
# 为了避免cell name冲突，添加前缀
merged_obj <- merge(prolif_cells, y = normal_cells_obj, add.cell.ids = c("Aged", "Normal"))

merged_obj=JoinLayers(merged_obj)

# 4. 在合并后的对象上进行差异分析
Idents(merged_obj) <- "deg_group"
deg_results <- FindMarkers(
  merged_obj,
  ident.1 = "Aged_Prolif",
  ident.2 = "Normal",
  logfc.threshold = 0, # 先不过滤
  min.pct = 0.1,
  verbose = FALSE
)

# 5. 筛选结果并保存到Excel
filtered_degs <- deg_results %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  arrange(desc(avg_log2FC))

write.xlsx(filtered_degs, file = "DEG_AgedProlif_vs_Normal.xlsx")













genes_to_check <- c("CTSC","GABRA3","SLC12A8","SLC12A2","CFTR")
group.by = "proliferation"

DotPlot(aged_cells_obj, features = genes_to_check,cols = c("RdYlBu"),group.by = group.by) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank())

VlnPlot(aged_cells_obj, features = genes_to_check, group.by = group.by)


