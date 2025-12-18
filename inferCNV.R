# ==============================================================================
# inferCNV 全流程标准化分析代码 (新手友好版)
# ==============================================================================

# ---
# 步骤 0: 加载必要的R包
# ---
library(infercnv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)       # 用于快速读写大数据
library(mixtools) # 如果需要做高斯混合模型(GMM)会用到

# 设置随机种子，保证结果可重复
set.seed(1234)

# ==============================================================================
# 步骤 1: 参数配置与数据准备 (!! 仅需修改此部分 !!)
# ==============================================================================

# 1.1 输入与输出设置
# 假设您的Seurat对象已经在环境中，名字叫 sc.obj
# 如果不在，请取消下面注释加载：
# sc.obj <- qread("path/to/your/seurat_object.qs") 

# 基因位置文件 (人类用hg38，小鼠用mm10，确保格式正确)
gene_order_file <- "hg38_gencode_v27.txt" 
output_dir <- "./inferCNV_output"

# 1.2 分组定义
# 定义哪些 cluster 是观测组（可能是肿瘤），哪些是参考组（正常细胞，如T细胞、内皮细胞）
# 注意：这里必须是字符(character)格式
observation_groups <- c("1", "2", "3", "7", "12", "13", "14", "15", "21", "23")
reference_groups   <- c("6") 
annotation_column  <- "seurat_clusters" # Seurat元数据中存储分类的列名

# 1.3 自动检查目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---
# 步骤 1.4: 生成 inferCNV 所需的注释文件
# ---
print("正在准备输入矩阵和注释文件...")

# 提取原始计数矩阵 (Raw Counts) - inferCNV 要求使用原始count
raw_counts_matrix <- GetAssayData(sc.obj, slot = "counts")

# 构建注释数据框
annotations_df <- data.frame(
  cell_id = rownames(sc.obj@meta.data),
  cell_type = sc.obj@meta.data[[annotation_column]],
  stringsAsFactors = FALSE
)

# 【重要技巧】将所有参考组的名称统一修改为 "Reference"
# 这样inferCNV在画热图时会把它们放在一起，且作为基准
annotations_df$cell_type[annotations_df$cell_type %in% reference_groups] <- "Reference"

# 检查是否有未定义的细胞
# 如果想只跑部分细胞，可以在这里做subset，但通常建议跑全集
annotations_path <- file.path(output_dir, "cell_annotations.txt")
write.table(annotations_df, file = annotations_path, sep = "\t", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

print("数据准备完成。")
gc() # 释放内存

# ==============================================================================
# 步骤 2: 创建对象并运行 inferCNV
# ==============================================================================

# 2.1 创建 inferCNV 对象
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file = annotations_path,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = "Reference" # 这里对应上面我们统一修改的名字
)

# 2.2 运行 inferCNV 主程序
# 注意：这一步非常耗时，可能需要几小时
infercnv_obj_run <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,      # 10x数据用0.1, Smart-seq2数据用1
  out_dir = output_dir,
  cluster_by_groups = TRUE, # 是否按组聚类
  denoise = TRUE,    # 降噪，推荐开启
  HMM = TRUE,        # 【推荐开启】用于预测CNV区域，发文章通常需要
  num_threads = parallel::detectCores() - 2, # 使用多核加速
  write_expr_matrix = TRUE, # 输出处理后的表达矩阵，后续算分必须用
  output_format = "pdf"
)

# 保存运行结果对象，防止R崩溃白跑
qsave(infercnv_obj_run, file.path(output_dir, "infercnv_obj_run.qs"))
print("inferCNV 运行结束。")

# ==============================================================================
# 步骤 3: 结果解析与恶性细胞识别 (关键步骤)
# ==============================================================================

print("开始计算CNV评分并定义细胞性质...")

# --- 3.1 读取 inferCNV 处理后的表达矩阵 ---
# 这些矩阵代表了相对于正常参考细胞的表达倍数变化

# 读取观测组数据
obs_data <- read.table(file.path(output_dir, "infercnv.observations.txt"), header = TRUE, check.names = FALSE)
# 读取参考组数据
ref_data <- read.table(file.path(output_dir, "infercnv.references.txt"), header = TRUE, check.names = FALSE)

# 【核心修正】将两者合并，确保所有细胞都在矩阵里
all_cnv_data <- cbind(obs_data, ref_data)

# --- 3.2 计算 CNV Score (残差平方和方法) ---
# 原理：正常细胞的相对表达量应接近1。
# (表达量 - 1)^2 的总和越大，说明基因组变异越剧烈。
# 这种方法比解析HMM文件更简单且极其稳健。

cnv_score <- colMeans((all_cnv_data - 1)^2)

# 构建分数数据框
cnv_score_df <- data.frame(cell_id = names(cnv_score), cnv_score = cnv_score)
# 修正细胞名：R读取文件时会将 "-" 变为 "."，需要变回来以匹配 Seurat
cnv_score_df$cell_id <- gsub("\\.", "-", cnv_score_df$cell_id)

# --- 3.3 将分数映射回 Seurat 对象 ---
sc.obj$cnv_score <- NA 
# 找到交集细胞
common_cells <- intersect(rownames(sc.obj@meta.data), cnv_score_df$cell_id)

if(length(common_cells) == 0) stop("错误：细胞ID无法匹配，请检查inferCNV输出文件与Seurat对象命名格式。")

sc.obj@meta.data[common_cells, "cnv_score"] <- cnv_score_df[common_cells, "cnv_score"]

# --- 3.4 定义恶性/非恶性阈值 ---
# 方法：利用参考组（正常细胞）分数的统计特征。
# 逻辑：如果一个细胞的分数显著高于参考组的平均水平，它就是恶性的。
# 阈值 = 参考组均值 + 3倍标准差 (99.7%置信区间)

# 提取参考组的分数
ref_scores <- sc.obj@meta.data %>%
  filter(!!sym(annotation_column) %in% reference_groups) %>%
  pull(cnv_score) %>%
  na.omit()

if(length(ref_scores) < 10) warning("警告：参考细胞数量过少，计算出的阈值可能不准确！")

# 计算阈值
threshold <- mean(ref_scores) + 3 * sd(ref_scores)
print(paste0("参考组均值: ", round(mean(ref_scores), 4)))
print(paste0("判读阈值 (Mean + 3SD): ", round(threshold, 4)))

# --- 3.5 分类细胞 ---
sc.obj$Malignant_Status <- "Undefined"

# 标记参考细胞
sc.obj$Malignant_Status[sc.obj@meta.data[[annotation_column]] %in% reference_groups] <- "Reference"

# 标记非参考细胞 (基于阈值)
# 只有当 cnv_score 存在 且 不是参考组时才进行判断
target_cells <- !is.na(sc.obj$cnv_score) & sc.obj$Malignant_Status != "Reference"

sc.obj$Malignant_Status[target_cells & sc.obj$cnv_score >= threshold] <- "Malignant"
sc.obj$Malignant_Status[target_cells & sc.obj$cnv_score < threshold]  <- "Non-Malignant"

# 输出统计结果
print("最终细胞分类统计：")
print(table(sc.obj$Malignant_Status))

# ==============================================================================
# 步骤 4: 可视化与保存
# ==============================================================================

print("正在生成可视化图表...")

# 4.1 CNV分数小提琴图 (查看分数分布)
p1 <- VlnPlot(sc.obj, features = "cnv_score", group.by = annotation_column, pt.size = 0) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  ggtitle(paste0("CNV Scores (Red Line = Threshold: ", round(threshold, 3), ")"))
ggsave(file.path(output_dir, "1_VlnPlot_CNV_Score.pdf"), p1, width = 12, height = 6)

# 4.2 UMAP图 - 展示 CNV 分数高低
p2 <- FeaturePlot(sc.obj, features = "cnv_score", order = TRUE) +
  scale_color_gradientn(colors = c("navy", "white", "firebrick"), name = "CNV Score") +
  ggtitle("CNV Score FeaturePlot")
ggsave(file.path(output_dir, "2_FeaturePlot_CNV_Score.pdf"), p2, width = 8, height = 7)

# 4.3 UMAP图 - 展示最终分类结果 (恶性 vs 非恶性)
p3 <- DimPlot(sc.obj, group.by = "Malignant_Status", 
              cols = c("Malignant" = "#E41A1C",       # 红色
                       "Non-Malignant" = "#377EB8",   # 蓝色
                       "Reference" = "#4DAF4A",       # 绿色
                       "Undefined" = "grey80")) +     # 灰色
  ggtitle("Final Malignant Classification")
ggsave(file.path(output_dir, "3_DimPlot_Malignancy.pdf"), p3, width = 8, height = 7)

# 4.4 密度分布图 (验证阈值合理性)
# 理想情况下，Reference在左侧，Malignant在右侧，阈值在中间
p4 <- ggplot(sc.obj@meta.data, aes(x = cnv_score, fill = Malignant_Status)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
  theme_classic() +
  ggtitle("Distribution of CNV Scores")
ggsave(file.path(output_dir, "4_Density_CNV_Distribution.pdf"), p4, width = 8, height = 6)

# --- 保存结果 ---
# 保存含有CNV信息和分类信息的完整Seurat对象
# qsave(sc.obj, file.path(output_dir, "sc_obj_with_cnv.qs"))

# 如果只想提取恶性细胞用于后续分析：
malignant_subset <- subset(sc.obj, subset = Malignant_Status == "Malignant")
qsave(malignant_subset, file.path(output_dir, "malignant_subset.qs"))

print("全部分析完成！")