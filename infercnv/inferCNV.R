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
  HMM = TRUE,               
  num_threads = parallel::detectCores() - 2, 
  write_expr_matrix = TRUE,
  output_format = "pdf"
)

qsave(infercnv_obj_run, file.path(output_dir, "infercnv_obj_run.qs"))

# ==============================================================================
# 步骤 4: CNV 评分计算与性质定义
# ==============================================================================

# 4.1 读取矩阵
obs_data <- read.table(file.path(output_dir, "infercnv.observations.txt"), header = TRUE, check.names = FALSE)
ref_data <- read.table(file.path(output_dir, "infercnv.references.txt"), header = TRUE, check.names = FALSE)
all_cnv_data <- cbind(obs_data, ref_data)

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

