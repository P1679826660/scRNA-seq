library(limma)
library(pheatmap)

# --- 参数设置 ---
inputFile  <- "geneMatrix.txt"
logFCfilter <- 0.5
adjPfilter <- 0.05

# --- 数据读取与处理 ---
rt <- read.table(inputFile, header=T, sep="\t", row.names=1, check.names=F)
data <- as.matrix(rt)

# 处理重复基因名 (取平均值)，并过滤全为0的低表达数据
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]

# 【关键修改】如果你的输入数据是原始的 FPKM 或 TPM，请取消下面这行的注释
# data <- log2(data + 1)

# --- 获取分组信息 ---
getSampleNames <- function(pattern) {
  files <- list.files(pattern = pattern)
  samples <- c()
  for (f in files) {
    tmp <- read.table(f, header=F, sep="\t", check.names=F)
    samples <- c(samples, as.vector(tmp[,1]))
  }
  return(unique(samples))
}

sampleName1 <- getSampleNames("s1.txt$") # 对照组
sampleName2 <- getSampleNames("s2.txt$") # 处理组

# --- 截取有效数据 ---
valid_s1 <- intersect(sampleName1, colnames(data))
valid_s2 <- intersect(sampleName2, colnames(data))
data <- data[, c(valid_s1, valid_s2)]

# --- 差异分析 ---
Type <- c(rep("con", length(valid_s1)), rep("treat", length(valid_s2)))
design <- model.matrix(~0 + factor(Type))
colnames(design) <- c("con", "treat")

fit <- lmFit(data, design)
cont.matrix <- makeContrasts(treat - con, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# --- 结果输出 (CSV) ---
allDiff <- topTable(fit2, adjust='fdr', number=Inf)
write.csv(allDiff, file="all.csv", quote=F)

diffSig <- allDiff[abs(allDiff$logFC) > logFCfilter & allDiff$adj.P.Val < adjPfilter, ]
write.csv(diffSig, file="diff.csv", quote=F)