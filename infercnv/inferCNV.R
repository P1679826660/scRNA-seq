# ==============================================================================
# inferCNV 全流程标准化分析代码 
# ==============================================================================

# --- 步骤 0: 环境准备 ---
library(infercnv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)
set.seed(1234)

# ==============================================================================
# 步骤 1: 参数配置
# ==============================================================================

# 1.1 设置路径
# sc.obj <- qread("你的数据路径.qs") 
# ！！！关键修改：改回相对路径，不要使用 normalizePath，这能解决文件乱跑的问题
output_dir <- "./inferCNV_output" 
gene_order_file <- "hg38_gencode_v27.txt" 

# 1.2 分组定义
annotation_column  <- "seurat_clusters" 
observation_groups <- c("1", "2", "3", "7", "12", "13", "14", "15", "21", "23")
reference_groups   <- c("6") # T细胞, NK等

# 1.3 目录初始化
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# 注意：这里删除了 abs_path <- normalizePath(...) 这一行

# ==============================================================================
# 步骤 2: 数据清洗与输入准备
# ==============================================================================
print("正在清洗数据：剔除未定义的细胞并提取Counts...")
sc.obj@meta.data$seurat_clusters <- as.character(sc.obj@meta.data$seurat_clusters)
# 2.1 筛选细胞：只保留观测组和参考组
keep_cells_types <- c(observation_groups, reference_groups)
sc.infer <- subset(sc.obj, subset = !!sym(annotation_column) %in% keep_cells_types)

# 2.2 提取原始计数矩阵
raw_counts_matrix <- GetAssayData(sc.infer, slot = "counts")


# 2.3 构建并保存注释文件
annotations_df <- data.frame(
  cell_id = colnames(sc.infer),
  cell_type = sc.infer@meta.data[[annotation_column]],
  stringsAsFactors = FALSE
)

# 统一参考组名称
annotations_df$cell_type[annotations_df$cell_type %in% reference_groups] <- "Reference"

# ！！！关键修改：直接使用 output_dir
annotations_path <- file.path(output_dir, "cell_annotations.txt")

write.table(annotations_df, file = annotations_path, sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

print(paste0("数据准备完成。输出目录：", output_dir))
print(paste0("分析细胞数：", ncol(sc.infer)))

# ==============================================================================
# 步骤 3: 运行 inferCNV
# ==============================================================================

# 3.1 创建对象
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file = annotations_path,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = "Reference"
)

# 3.2 运行主程序
infercnv_obj_run <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,             # 10x数据推荐0.1
  out_dir = output_dir,     # ！！！关键修改：这里直接用相对路径
  cluster_by_groups = TRUE, 
  denoise = TRUE,           
  HMM = F,               
  num_threads = parallel::detectCores() - 2, 
  write_expr_matrix = TRUE,
  output_format = "pdf"
)

qsave(infercnv_obj_run, file.path(output_dir, "infercnv_obj_run.qs"))





# 计算每个细胞的CNV总分 (简单逻辑：偏离1的绝对值之和)
cnv_matrix <- infercnv_obj_run@expr.data
cnv_score <- colSums((cnv_matrix - 1)^2)

# 将分值加回到你的Seurat对象中查看分布
sc.obj$cnv_score <- cnv_score[colnames(sc.obj)]

# 按照你的 'infer' 分组查看 CNV 指数分布
VlnPlot(sc.obj, features = "cnv_score", group.by = "infer", pt.size = 0) +
  geom_hline(yintercept = mean(sc.obj$cnv_score[sc.obj$infer == "Endo"]), 
             linetype = "dashed", color = "red") +
  ggtitle("CNV Score by Cluster")






下面为旧代码
# ==============================================================================
# 步骤 4: CNV 评分计算与性质定义
# ==============================================================================

# 4.1 读取矩阵
#obs_data <- read.table(file.path(output_dir, "infercnv.observations.txt"), header = TRUE, check.names = FALSE)
#ref_data <- read.table(file.path(output_dir, "infercnv.references.txt"), header = TRUE, check.names = FALSE)
#all_cnv_data <- cbind(obs_data, ref_data)

# 4.1 直接从对象提取矩阵，避免读取文件的各种符号错误
all_data_matrix <- infercnv_obj_run@expr.data
obs_cells <- colnames(raw_counts_matrix)[!(colnames(raw_counts_matrix) %in% infercnv_obj_run@reference_grouped_cell_indices)]

# 4.2 计算分数 (RSS方法)
# 计算每个细胞所有基因CNV值的平方误差均值
cnv_score <- colMeans((all_cnv_data - 1)^2)
cnv_score_df <- data.frame(cell_id = gsub("\\.", "-", names(cnv_score)), cnv_score = cnv_score)

# 4.3 映射回原始 Seurat 对象 (sc.obj)
# 注意：这里我们要把分数映射回最原始的大对象 sc.obj，而不是 sc.infer，方便你后续看全景
sc.obj$cnv_score <- NA
common_cells <- intersect(rownames(sc.obj@meta.data), cnv_score_df$cell_id)
sc.obj@meta.data[common_cells, "cnv_score"] <- cnv_score_df[match(common_cells, cnv_score_df$cell_id), "cnv_score"]

# 4.4 阈值判读 (Mean + 3*SD)
ref_scores <- sc.obj@meta.data %>%
  filter(!!sym(annotation_column) %in% reference_groups) %>%
  pull(cnv_score) %>% na.omit()

# 如果参考细胞太少，可能会报错，建议检查 length(ref_scores)
if(length(ref_scores) < 10) warning("警告：参考细胞数量过少，阈值计算可能不准确")

threshold <- mean(ref_scores) + 3 * sd(ref_scores)
print(paste0("恶性判定阈值设定为: ", round(threshold, 4)))

# 4.5 标记状态
sc.obj$Malignant_Status <- "Undefined"
# 先标记参考细胞
sc.obj$Malignant_Status[sc.obj@meta.data[[annotation_column]] %in% reference_groups] <- "Reference"

# 标记观测细胞
target_cells <- !is.na(sc.obj$cnv_score) & sc.obj$Malignant_Status != "Reference"
sc.obj$Malignant_Status[target_cells & sc.obj$cnv_score >= threshold] <- "Malignant"
sc.obj$Malignant_Status[target_cells & sc.obj$cnv_score < threshold]  <- "Non-Malignant"

# ==============================================================================
# 步骤 5: 可视化保存与提取
# ==============================================================================

# 小提琴图
p <- VlnPlot(sc.obj, features = "cnv_score", group.by = annotation_column, pt.size = 0) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  ggtitle(paste0("CNV Score Distribution (Threshold: ", round(threshold, 3), ")"))

ggsave(file.path(output_dir, "Malignant_VlnPlot.pdf"), p, width = 10, height = 6)

# 保存最终对象
qsave(sc.obj, file.path(output_dir, "sc_obj_final.qs"))

# 提取恶性细胞子集并保存
print(table(sc.obj$Malignant_Status))
malignant_obj <- subset(sc.obj, subset = Malignant_Status == "Malignant")
qsave(malignant_obj, file.path(output_dir, "malignant_obj.qs"))

print("分析全部结束。所有文件已保存至 inferCNV_output 文件夹。")




#############
#############
#############
#针对肝癌（Liver Cancer），TCGA 的标准缩写是 LIHC (Liver Hepatocellular Carcinoma)。
#以下是针对你的对象 vvv 的完整分析代码：
#第一步：从 vvv 中提取数据并创建 SpaCET 对象
#SpaCET 需要两个核心输入：表达矩阵和空间坐标。
library(SpaCET)
library(Seurat)
# 1. 从你的对象 vvv 中提取 Count 计数矩阵
# 注意：使用原始 counts，不要用 scale.data
counts_matrix <- GetAssayData(vvv, slot = "counts", assay = "Spatial")
# 2. 提取空间坐标 (适配 Seurat V4/V5)
# 这步是为了让 SpaCET 知道每个点在图上的位置
spatial_coords <- GetTissueCoordinates(vvv)
# 3. 创建 SpaCET 对象
# 既然是从对象转换，我们用 create.SpaCET.object 而不是 create.SpaCET.object.10X
spaCET_obj <- create.SpaCET.object(
  counts = counts_matrix,
  spotCoordinates = spatial_coords
)
# 4. 基础质量控制 (过滤掉基因数极少的空点)
spaCET_obj <- SpaCET.quality.control(spaCET_obj, min.gene = 100)

#第二步：识别肝癌细胞 (LIHC)
#这是核心步骤。我们将 cancerType 设置为 "LIHC"。
# 5. 运行去卷积分析
# coreNo = 4 表示用4个核并行计算，速度更快
spaCET_obj <- SpaCET.deconvolution(
  spaCET_obj, 
  cancerType = "LIHC", 
  coreNo = 4
)
# 6. (可选) 可视化 SpaCET 自带的热图
SpaCET.visualize.spatialFeature(
  spaCET_obj, 
  spatialType = "cellFraction", 
  spatialFeatures = "Malignant"
)


#第三步：将鉴定结果“装回” vvv 对象 (最重要的一步)
#为了方便你继续使用 Seurat 的功能（比如画更好看的图，或者做亚群分析），我们把算出来的“恶性比例”存回到 vvv 的 meta.data 里。
# 7. 提取恶性细胞比例 (Malignant Fraction)
# SpaCET 的结果存储在 @results$deconvolution$prop 矩阵中
prop_mat <- spaCET_obj@results$deconvolution$prop
malignant_frac <- prop_mat["Malignant", ]

# 8. 映射回 vvv 对象
# 确保细胞名称一致
common_spots <- intersect(colnames(vvv), names(malignant_frac))
vvv$Malignant_Score <- 0 # 初始化
vvv@meta.data[common_spots, "Malignant_Score"] <- malignant_frac[common_spots]

# 9. 使用 Seurat 画图
# 这样你就能用大家熟悉的 Seurat 语法画图了
SpatialFeaturePlot(vvv, features = "Malignant_Score") + 
  ggtitle("Liver Cancer Malignant Cells Distribution")





#################
#################
#################
#Cottrazm包，自动识别肿瘤边界、自动识别肿瘤，空间组学分析
#精准地画出“肿瘤边界”


