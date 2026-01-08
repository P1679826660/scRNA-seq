# ==============================================================================
# inferCNV 全流程标准化分析代码 (修正版)
# ==============================================================================

# --- 步骤 0: 环境准备 ---
library(infercnv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)
set.seed(1234)

# ==============================================================================
# 步骤 1: 参数配置 (在此修改)
# ==============================================================================

# 1.1 设置路径 (建议使用绝对路径)
# sc.obj <- qread("你的数据路径.qs") 
output_dir <- "./inferCNV_output"
gene_order_file <- "hg38_gencode_v27.txt" # 确保此文件在工作目录下

# 1.2 分组定义
annotation_column  <- "seurat_clusters" 
observation_groups <- c("1", "2", "3", "7", "12", "13", "14", "15", "21", "23")
reference_groups   <- c("6") # T细胞, NK等

# 1.3 目录初始化
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
abs_path <- normalizePath(output_dir)

# ==============================================================================
# 步骤 2: 数据清洗与输入准备 (核心修正区)
# ==============================================================================
print("正在清洗数据：剔除未定义的细胞并提取Counts...")

# 2.1 筛选细胞：只保留观测组和参考组，剔除B细胞等无关干扰
keep_cells_types <- c(observation_groups, reference_groups)
sc.infer <- subset(sc.obj, subset = !!sym(annotation_column) %in% keep_cells_types)

# 2.2 提取原始计数矩阵 (必须是 counts)
if (packageVersion("Seurat") >= "5.0.0") {
  raw_counts_matrix <- Seurat::LayerData(sc.infer, assay = "RNA", layer = "counts")
} else {
  raw_counts_matrix <- GetAssayData(sc.infer, slot = "counts")
}

# 2.3 构建并保存注释文件
annotations_df <- data.frame(
  cell_id = colnames(sc.infer),
  cell_type = sc.infer@meta.data[[annotation_column]],
  stringsAsFactors = FALSE
)

# 统一参考组名称
annotations_df$cell_type[annotations_df$cell_type %in% reference_groups] <- "Reference"

annotations_path <- file.path(abs_path, "cell_annotations.txt")
write.table(annotations_df, file = annotations_path, sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

print(paste0("数据准备完成。输出目录：", abs_path))
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

# 3.2 运行主程序 (耗时操作)
infercnv_obj_run <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,             # 10x数据用0.1
  out_dir = abs_path,
  cluster_by_groups = TRUE, 
  denoise = TRUE,           
  HMM = TRUE,               
  num_threads = parallel::detectCores() - 2, 
  write_expr_matrix = TRUE,
  output_format = "pdf"
)

qsave(infercnv_obj_run, file.path(abs_path, "infercnv_obj_run.qs"))

# ==============================================================================
# 步骤 4: CNV 评分计算与性质定义
# ==============================================================================

# 4.1 读取矩阵
obs_data <- read.table(file.path(abs_path, "infercnv.observations.txt"), header = TRUE, check.names = FALSE)
ref_data <- read.table(file.path(abs_path, "infercnv.references.txt"), header = TRUE, check.names = FALSE)
all_cnv_data <- cbind(obs_data, ref_data)

# 4.2 计算分数 (RSS方法)
cnv_score <- colMeans((all_cnv_data - 1)^2)
cnv_score_df <- data.frame(cell_id = gsub("\\.", "-", names(cnv_score)), cnv_score = cnv_score)

# 4.3 映射回原始 Seurat 对象 (sc.obj)
sc.obj$cnv_score <- NA
common_cells <- intersect(rownames(sc.obj@meta.data), cnv_score_df$cell_id)
sc.obj@meta.data[common_cells, "cnv_score"] <- cnv_score_df[match(common_cells, cnv_score_df$cell_id), "cnv_score"]

# 4.4 阈值判读
ref_scores <- sc.obj@meta.data %>%
  filter(!!sym(annotation_column) %in% reference_groups) %>%
  pull(cnv_score) %>% na.omit()

threshold <- mean(ref_scores) + 3 * sd(ref_scores)

# 4.5 标记状态
sc.obj$Malignant_Status <- "Undefined"
sc.obj$Malignant_Status[sc.obj@meta.data[[annotation_column]] %in% reference_groups] <- "Reference"
target_cells <- !is.na(sc.obj$cnv_score) & sc.obj$Malignant_Status != "Reference"
sc.obj$Malignant_Status[target_cells & sc.obj$cnv_score >= threshold] <- "Malignant"
sc.obj$Malignant_Status[target_cells & sc.obj$cnv_score < threshold]  <- "Non-Malignant"

# ==============================================================================
# 步骤 5: 可视化保存
# ==============================================================================

p <- VlnPlot(sc.obj, features = "cnv_score", group.by = annotation_column, pt.size = 0) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red")
ggsave(file.path(abs_path, "Malignant_VlnPlot.pdf"), p, width = 10, height = 6)

qsave(sc.obj, file.path(abs_path, "sc_obj_final.qs"))

print("分析全部结束。")
