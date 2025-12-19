# ==============================================================================
# Seurat 单细胞转录组标准分析流程 (优化版)
# 包括：QC -> Harmony去批次 -> 聚类 -> 去双细胞 -> 去背景噪音 -> 注释 -> 找Marker
# ==============================================================================

# ---
# 1. 环境准备与包加载
# ---
# 建议一次性加载所有需要的包，避免中途报错
library(Seurat)
library(qs)            # 极速读写大数据
library(harmony)       # 去除批次效应
library(ggplot2)       # 绘图基础
library(scDblFinder)   # 去除双细胞
library(SingleCellExperiment) # scDblFinder需要的数据格式
library(decontX)       # 去除游离RNA污染 (背景噪音)
library(viridis)       # 美观的配色方案
library(tidyverse)         # 数据处理神器 (mutate, filter等)

# ---
# 2. 数据加载与初步质控 (QC)
# ---

#第一种格式，经典10X
Cell Ranger输出的标准结果通常在一个名为 filtered_feature_bc_matrix 的文件夹里。 
重要：你不需要解压内部的文件，只需要提供这个文件夹的路径。
该文件夹内必须包含以下三个文件（文件名必须完全一致，可以是压缩包 .gz）：
barcodes.tsv.gz
features.tsv.gz (旧版本可能是 genes.tsv)
matrix.mtx.gz
# 读取数据
# Read10X 会自动寻找上述三个文件并构建稀疏矩阵
counts_matrix <- Read10X("filtered_feature_bc_matrix")
# 创建 Seurat 对象
sc.obj <- CreateSeuratObject(
  counts = counts_matrix,
  project = "Sample1",       # 给你的项目或样本起个名字
  min.cells = 3,             # 过滤掉：只在少于3个细胞中表达的基因（极低表达基因）
  min.features = 200         # 过滤掉：检测到基因数少于200的细胞（通常是空液滴或死细胞）
)

#第二种格式，如果是filtered_feature_bc_matrix.h5，h5格式
counts_matrix <- Read10X_h5(filename = h5_file)

#第三种格式，TXT 表达矩阵读取
library(data.table) # 用于 fread 快速读取大数据
library(limma)      # 用于 avereps 处理重复基因
library(Seurat)     # 单细胞分析核心包
library(Matrix)     # 用于处理矩阵

# ---
# 步骤 2: 读取数据 (使用 fread 提速)
# ---
# data.table::fread 读取速度远快于 read.table，特别适合几百兆以上的大文件
# check.names=F 保证细胞名中的特殊字符（如 "-"）不被自动修改为 "."
file_path <- "sc.obj.txt"

rt <- fread(file_path, sep = "\t", header = TRUE, check.names = FALSE)
# 转换格式：fread读进来是data.table，我们需要分离基因列和表达量矩阵
# 假设第一列是基因名 (Gene Symbol)
gene_names <- rt[[1]]   # 提取第一列作为基因名向量
expr_data  <- as.matrix(rt[, -1, with = FALSE]) # 提取除第一列外的所有列，转为矩阵

# 立即删除原始读取的大对象，释放内存
rm(rt)
gc() 

# 步骤 3: 处理重复基因 (去重)
# 在单细胞或转录组数据中，经常会出现同名基因（由于注释版本或多转录本原因）。
# 如果不处理，Seurat 会报错 "Duplicate feature names allowed"。
# 使用 limma::avereps 对重复的基因名取平均值 (Average)
# ID 参数指定了每一行对应的基因名
# 这一步会合并重复行，并返回一个行名唯一的矩阵
final_matrix <- limma::avereps(expr_data, ID = gene_names)

# 再次清理内存
rm(expr_data, gene_names)
gc()

# 步骤 4: 构建 Seurat 对象
# min.cells = 5:    过滤掉在少于5个细胞中表达的基因（低表达基因）
# min.features = 300: 过滤掉检测到少于300个基因的细胞（低质量细胞/空液滴）
# names.delim = "_": 如果细胞名是 "Sample1_Barcode" 格式，这会告诉Seurat前缀是样本名
pbmcT <- CreateSeuratObject(
  counts = final_matrix,
  project = "seurat",
  min.cells = 5,
  min.features = 300,
  names.delim = "_"
)



#读取qs格式
# 读取数据 (假设是qs格式已经转换完毕的seurat)
# sc.obj <- qread("path/to/GSE149614.Tumor.qs") 

###########################################################
library(Seurat)
library(dplyr) # 用于处理数据

# ==============================================================================
# 步骤 1: 读取外部 Metadata 文件
# ==============================================================================
# 假设 metadata 是 csv 或 txt 文件
# header=T: 第一行是列名
# row.names=1: 【关键】直接尝试把第一列作为行名读取
# check.names=F: 防止 R 把细胞名里的 '-' 变成 '.'
meta_file <- "GSEXXXX_metadata.csv" 

# 如果是 csv
metadata <- read.csv(meta_file, header = TRUE, row.names = 1, check.names = FALSE)

# 如果是 txt (制表符分隔)
# metadata <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

head(metadata)

# ==============================================================================
# 步骤 2: 【至关重要】核对细胞名格式
# ==============================================================================
# 很多时候，外部文件的细胞名和 Seurat 里的细胞名长得不一样
# 比如：Seurat 里是 "Sample1_AAACCC...", 而文件里只有 "AAACCC..."

head(colnames(sc.obj))
head(rownames(metadata))

# --- 检查交集 ---
# 计算有多少细胞是名字完全一样的
common_cells <- intersect(colnames(sc.obj), rownames(metadata))
cat("Seurat 细胞总数:", ncol(sc.obj), "\n")
cat("Metadata 记录数:", nrow(metadata), "\n")
cat("名字能匹配上的数量:", length(common_cells), "\n")

if (length(common_cells) == 0) {
  stop("错误：没有一个细胞名能匹配上！请检查两个数据的命名格式差异。")
} else if (length(common_cells) < ncol(sc.obj)) {
  warning("警告：只有部分细胞匹配成功，未匹配的细胞其新Metadata将为 NA。")
}

# ==============================================================================
# 步骤 3: 写入 Metadata
# ==============================================================================
# 只要行名对得上，Seurat 会自动匹配，不需要顺序一致
sc.obj <- AddMetaData(object = sc.obj, metadata = metadata)
# ==============================================================================
# 步骤 4: 验证结果
# ==============================================================================
print("写入后的 Seurat Metadata:")
head(sc.obj@meta.data)
################################################



# 2.1 计算线粒体基因比例
# 线粒体比例过高通常意味着细胞破裂或死细胞
sc.obj[["percent.mt"]] <- PercentageFeatureSet(object = sc.obj, pattern = "^MT-")

# 2.2 绘制质控小提琴图 (过滤前)
# nFeature_RNA: 基因数量 (代表细胞复杂度)
# nCount_RNA:   UMI总数 (代表测序深度)
# percent.mt:   线粒体比例
VlnPlot(object = sc.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "sample", ncol = 3)

# 2.3 执行过滤
# 逻辑：保留基因数在300-6000之间且线粒体小于15%的细胞
# 注意：这些阈值需根据实际VlnPlot图调整
sc.obj <- subset(x = sc.obj, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 15) 

# 绘制过滤后的质控图以供对比
VlnPlot(object = sc.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by = "sample", ncol = 3)

# ---
# 3. 标准化、降维与去批次 (Harmony)
# ---

# 3.1 标准化与高变基因寻找
# LogNormalize: 将表达量对数化，消除测序深度差异
sc.obj <- NormalizeData(object = sc.obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 寻找高变基因 (Top 2000)
# 这些基因在细胞间差异最大，用于后续的主成分分析
sc.obj <- FindVariableFeatures(object = sc.obj, selection.method = "vst", nfeatures = 2000) 

# 归一化 (Scale)
# 将数据调整为均值为0，方差为1，防止高表达基因主导PCA
sc.obj <- ScaleData(sc.obj)   
gc() # 释放内存

# 3.2 PCA 主成分分析
sc.obj <- RunPCA(sc.obj, features = VariableFeatures(object = sc.obj))

# 3.3 Harmony 去批次效应
# group.by.vars = "sample": 告诉算法根据 sample 列来校正批次
sc.obj <- RunHarmony(sc.obj, group.by.vars = "sample")

# 3.4 确定后续分析使用的维度
# 拐点图：帮助决定用多少个PC (Principal Components)
ElbowPlot(sc.obj)

# ---
# 4. 聚类与可视化 (初次聚类)
# ---

pcSelect <- 30  # 根据ElbowPlot选择，通常选20-30
# 注意：这里 reduction 必须选 "harmony"，否则没去批次
sc.obj <- FindNeighbors(sc.obj, reduction = "harmony", dims = 1:pcSelect)
sc.obj <- FindClusters(object = sc.obj, resolution = 1) # 分辨率越高，亚群分得越细

# UMAP 可视化降维
sc.obj <- RunUMAP(sc.obj, reduction = "harmony", dims = 1:pcSelect)
DimPlot(sc.obj, reduction = "umap", label = TRUE, label.size = 2) + ggtitle("Initial Clustering")

# ==============================================================================
# 5. 高级清洗：去除双细胞 (Doublets)
# ==============================================================================
# 原理：模拟双细胞的表达谱，看现有细胞是否与模拟数据高度相似

print("正在运行 scDblFinder 去除双细胞...")
# 提取样本信息，因为双细胞主要在同一样本内产生
name <- sc.obj@meta.data$sample

# 转换为 SingleCellExperiment 格式
sce <- as.SingleCellExperiment(sc.obj)

# 运行检测程序
# dbr (Doublet Rate) 默认为自动估算，如果知道10x的上样量，可手动指定 (如 0.08)
sce <- scDblFinder(sce, samples = name)

# 将结果导回 Seurat 对象
sc.obj$scDblFinder.score <- sce$scDblFinder.score # 双细胞得分
sc.obj$scDblFinder.class <- sce$scDblFinder.class # 分类结果 (singlet/doublet)
rm(sce); gc() # 清理内存

# 查看双细胞比例
print(table(sc.obj$scDblFinder.class))

# 可视化
DimPlot(sc.obj, group.by = "scDblFinder.class", cols = c("doublet" = "red", "singlet" = "grey"))

# 执行过滤：只保留单细胞 (singlet)
sc.obj <- subset(sc.obj, subset = scDblFinder.class == "singlet")

# ==============================================================================
# 6. 高级清洗：去除环境游离RNA (decontX)
# ==============================================================================
# 原理：悬液中漂浮的RNA (汤) 会被所有细胞吸附，导致背景噪音。

print("正在运行 decontX 去除环境RNA污染...")
# 提取原始计数矩阵 (decontX 要求 raw counts)
counts_matrix <- GetAssayData(sc.obj, assay = "RNA", layer = "count")

# 运行 decontX
decontX_results <- decontX(counts_matrix)

# 将污染分数添加回 metadata
sc.obj$Contamination <- decontX_results$contamination

# 可视化污染情况
FeaturePlot(sc.obj, features = 'Contamination', raster = FALSE) + 
  scale_color_viridis_c() + 
  theme_bw() + 
  ggtitle("Ambient RNA Contamination Score")

# --- 6.1 自动计算截断阈值 (几何距离法) ---
# 这是一个非常经典的算法，用于寻找曲线的"拐点" (Elbow point)

# 准备数据
contamination_scores <- sc.obj@meta.data$Contamination
sorted_scores <- sort(contamination_scores, decreasing = TRUE)
elbow_df <- data.frame(rank = 1:length(sorted_scores), score = sorted_scores)

# 数据归一化 (Normalization)
# 必须做！因为 Rank (X轴，几万) 和 Score (Y轴，0-1) 量级不同，无法直接算距离
df_norm <- elbow_df %>%
  mutate(
    rank_norm = (rank - min(rank)) / (max(rank) - min(rank)),
    score_norm = (score - min(score)) / (max(score) - min(score))
  )

# 计算每个点到“对角线”的距离
# 对角线连接了 (0,1) 和 (1,0)。距离最大的点，就是曲线最弯的那个点（拐点）
df_norm <- df_norm %>%
  mutate(distance = abs(rank_norm + score_norm - 1) / sqrt(2))

# 找到最大距离对应的索引
elbow_index <- which.max(df_norm$distance)
auto_threshold <- elbow_df$score[elbow_index]
auto_rank <- elbow_df$rank[elbow_index]

cat("自动检测阈值结果:\n")
cat("建议过滤阈值 (Contamination):", round(auto_threshold, 3), "\n")

# 绘制带有阈值线的拐点图
p_elbow <- ggplot(elbow_df, aes(x = rank, y = score)) +
  geom_point(size = 1, alpha = 0.5, color = "dodgerblue4") + 
  geom_vline(xintercept = auto_rank, linetype = "dashed", color = "red") +
  geom_hline(yintercept = auto_threshold, linetype = "dashed", color = "red") +
  labs(title = "Contamination Elbow Plot", 
       subtitle = paste("Threshold:", round(auto_threshold, 3))) +
  theme_bw()
print(p_elbow)

# --- 6.2 执行过滤 ---
# 使用自动阈值进行过滤 (您也可以手动指定一个值，如 0.2)
sc.obj <- subset(sc.obj, subset = Contamination < auto_threshold)
cat("过滤后剩余细胞数:", ncol(sc.obj), "\n")

# ==============================================================================
# 7. 细胞类型注释 (Annotation)
# ==============================================================================

# 7.1 检查 Marker 基因气泡图
# PTPRC(CD45)=免疫, EPCAM/KRT=上皮, PECAM1/VWF=内皮, COL1A1=成纤维
genes_to_check <- c("PTPRC", "IGKC",                      # 免疫细胞
                    "EPCAM", "KRT8", "KRT18", "ALB", "SERPINA1", "HNF4A", # 上皮/肝细胞
                    "PECAM1", "PLVAP", "VWF",             # 内皮细胞
                    "COL1A1", "DCN")                      # 成纤维细胞

DotPlot(sc.obj, features = genes_to_check, cols = c("RdYlBu")) + 
  RotatedAxis() +
  theme(axis.title = element_blank())

# 7.2 手动注释 (根据上面的DotPlot结果)
# 注意：这里的 new.cluster.ids 必须严格对应 cluster 0, 1, 2... 的顺序
# 请根据您最新的 FindClusters 结果核对数量
new.cluster.ids <- c(
  "PTPRC", "Epi", "Epi", "PTPRC", "PTPRC", 
  "PTPRC", "Endo", "Epi", "Epi", "Fib", 
  "PTPRC", "PTPRC", "PTPRC", "Epi", "PTPRC", 
  "Epi", "PTPRC"
)

# 获取当前 cluster ID
clusters <- as.character(sc.obj$seurat_clusters)
max_cluster_num <- max(as.numeric(clusters))

# 安全检查：确保ID数量匹配
if (length(new.cluster.ids) != (max_cluster_num + 1)) {
  stop(paste0("错误：您提供了 ", length(new.cluster.ids), " 个名字，但数据中有 ", 
              (max_cluster_num + 1), " 个 Cluster。请重新检查 new.cluster.ids"))
}

# 创建新的注释列 celltype.1
# 使用 factor 并指定 levels，保证顺序不乱
sc.obj$celltype.1 <- factor(x = clusters,
                            levels = as.character(0:max_cluster_num),
                            labels = new.cluster.ids)

# 可视化最终注释结果
DimPlot(sc.obj, reduction = "umap", label = TRUE, group.by = "celltype.1") + 
  ggtitle("Annotated Clusters")

# 保存阶段性成果
# qsave(sc.obj, "GSE149614.Tumor.Annotated.qs")

# ==============================================================================
# 8. 差异基因分析 (FindMarkers)
# ==============================================================================

print("正在计算各细胞类型的 Marker 基因...")

# FindAllMarkers: 计算每个群相对于其他所有群的差异基因
markers <- FindAllMarkers(
  object = sc.obj,
  group.by = "celltype.1", # 按刚才注释的细胞类型分组
  only.pos = TRUE,         # 只看高表达基因
  min.pct = 0.25,          # 至少在25%的细胞中表达 (由0.1提高到0.25可减少噪音)
  logfc.threshold = 1      # 差异倍数至少2倍 (loge(2) ≈ 0.69, 这里设1更严格)
) %>% 
  filter(p_val_adj < 0.05) # 过滤掉不显著的

# 保存结果
write.csv(markers, file = "celltype_markers.csv", row.names = FALSE)

print("全部分析流程结束！")


