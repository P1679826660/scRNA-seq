# 步骤 1: 安装和加载 DESeq2 包
# if (!requireNamespace("DESeq2", quietly = TRUE))
#   BiocManager::install("DESeq2")
###############################################
#第一步运行这个
###############################################
# 加载包
library(DESeq2)
library(dplyr) 

# 步骤 2: 读取基因计数文件
# 确保 GSE247245_raw_counts.csv 文件在您的 R 工作目录中

###############################################
#这里有两个要注意1.名称   2. sep = "," 分隔符
###############################################
countData_raw <- read.table("GSE215011_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t")
head(countData_raw)

# 使用 make.unique 函数强制使第一列的基因名唯一
unique_gene_names <- make.unique(as.character(countData_raw[, 1]))

# 将唯一的基因名设为行名
rownames(countData_raw) <- unique_gene_names

# 移除原始的基因名列，得到最终的计数矩阵
countData <- countData_raw[, -1]

# 查看处理后的数据的前几行
head(countData)

#############################################################
#改下面的两个txt 设置对照和实验组
##############################################################
# 步骤 3: 从外部文件动态创建样本信息表 (colData)
# 读取对照组和实验组的样本名
# 确保 control_groups.txt 和 experimental_groups.txt 文件在您的工作目录中
control_samples <- readLines("对照组.txt")
experimental_samples <- readLines("实验组.txt")

# 合并所有需要分析的样本名
all_samples <- c(control_samples, experimental_samples)

# 根据读取到的样本名，从原始计数矩阵中提取对应的列
countData.1 <- countData[, all_samples]

# 创建条件因子 (conditions)
# "Control" 对应对照组样本, "Dox" 对应实验组样本
conditions <- factor(c(rep("Control", length(control_samples)), rep("Dox", length(experimental_samples))))

# 创建 colData 数据框
colData <- data.frame(row.names = all_samples, condition = conditions)

# 检查 colData 是否正确
print(colData)

# 确保 colData 的行名和 countData.1 的列名完全一致
all(rownames(colData) == colnames(countData.1))

# 步骤 4: 创建 DESeqDataSet 对象
# design 参数告诉 DESeq2 如何根据 condition 列进行比较
dds <- DESeqDataSetFromMatrix(countData = countData.1,
                              colData = colData,
                              design = ~ condition)

# 预过滤：去除那些在所有样本中计数都极低的基因（例如，总计数小于10）
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 步骤 5: 运行 DESeq2 差异表达分析
# 这个函数包含了归一化、离散度估计和统计检验
dds <- DESeq(dds)

# 步骤 6: 提取差异表达结果
# contrast 参数指定了比较的分组, 这里比较 Dox 相对于 Control 组的变化
res <- results(dds, contrast=c("condition", "Dox", "Control"))

# 将结果转换为数据框，方便后续处理
res_df <- as.data.frame(res)

# 查看原始结果的前几行
head(res_df)

# 步骤 7: 根据您提供的阈值筛选差异基因
# pvalue < 0.05
# log2FoldChange > 0.5 (上调) 或 < -0.5 (下调)

significant_genes <- res_df %>%
  filter(pvalue < 0.05 & abs(log2FoldChange) > 0.5) %>%
  # 按 padj 从小到大排序
  arrange(padj)

# 查看筛选出的差异基因数量
cat(paste("\n找到", nrow(significant_genes), "个显著差异基因。\n"))

# 查看筛选结果的前几行
head(significant_genes)

# 步骤 8: 保存结果到文件
# 保存显著差异的基因
write.csv(significant_genes, file = "significant_differential_genes.csv")

cat("\n分析完成！差异基因结果已保存到 significant_differential_genes.csv 文件中。\n")

