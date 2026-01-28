# ==============================================================================
# KEGG 分析、筛选与可视化一体化流程 (优化版)
# ==============================================================================

# --- 0. 载入包 ---
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(showtext)
})

# --- 1. 全局配置区域 (请在这里修改参数) ---

# 1.1 输入输出文件
input_file <- "要区分上下调.txt"
output_file <- "KEGG_Final_Result.txt" # 最终保存的文件名

# 1.2 统计阈值
p_cutoff <- 0.05
q_cutoff <- 0.05

# 1.3 筛选关键词设置 (解决筛选不方便的问题)
# 模式A：剔除模式 (最常用)。保留大部分，剔除包含以下词的通路。
# 例如：剔除所有包含 "Disease", "Cancer", "Infection" 的通路
exclude_keywords <- c("Disease", "Cancer", "Infection", "Leukemia", "Carcinoma")

# 模式B：仅保留模式。只保留包含以下词的通路。
# 如果想用这个模式，请将下方变量设为 TRUE，并在 include_keywords 填入内容
use_include_mode <- FALSE 
include_keywords <- c("Signaling pathway", "Metabolism") 

# 1.4 绘图设置
plot_color <- "#F8A2A4"
top_n_plot <- 30 # 画前多少个

# 字体设置
font_add("arial", "arial.ttf") 
showtext_auto()


# --- 2. 数据读取与ID转换 ---

cat(">>> 正在读取数据并转换 ID...\n")
rt <- read.table(input_file, sep="\t", check.names=F, header=TRUE) # 建议加上header=TRUE，除非你确定没有表头

# 提取基因列 (假设第一列是Symbol)
genes <- unique(as.vector(rt[,1]))

# ID转换 (Symbol -> Entrez ID)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)

# 创建映射表供后续使用
rt_map <- data.frame(genes = genes, entrezID = entrezIDs, stringsAsFactors = FALSE)

# 过滤掉无法转换ID的基因
gene_input <- entrezIDs[entrezIDs != "NA" & !is.na(entrezIDs)]

cat(sprintf(">>> 输入基因 %d 个，成功转换 ID %d 个。\n", length(genes), length(gene_input)))


# --- 3. 运行 KEGG 富集分析 ---

cat(">>> 正在运行 KEGG 富集分析 (请耐心等待)...\n")
options(timeout = 600)

kk <- enrichKEGG(gene = gene_input, 
                 organism = "hsa", 
                 pvalueCutoff = 1, # 先放宽阈值，后续在本地筛选
                 qvalueCutoff = 1)

if (is.null(kk) || nrow(kk) == 0) {
  stop("KEGG 分析未找到任何富集通路，请检查基因 ID 是否正确。")
}

# 将结果转换为 Dataframe 并进行初步处理
kegg_df <- as.data.frame(kk)

# 将 GeneID 转换回 Symbol (利用之前的映射表)
# 这里优化了之前的 apply 写法，提高可读性
kegg_df$geneID <- sapply(kegg_df$geneID, function(x) {
  ids <- strsplit(x, "/")[[1]]
  symbols <- rt_map$genes[match(ids, rt_map$entrezID)]
  paste(symbols, collapse = "/")
})

# 初步阈值过滤
kegg_df <- kegg_df %>% 
  filter(pvalue < p_cutoff, qvalue < q_cutoff)

cat(sprintf(">>> 初步富集到 %d 条通路。\n", nrow(kegg_df)))


# --- 4. 自动化筛选 (核心优化) ---

cat(">>> 正在根据关键词进行筛选...\n")

if (use_include_mode) {
  # 模式B：只保留匹配白名单的
  pattern <- paste(include_keywords, collapse = "|")
  kegg_filtered <- kegg_df %>% filter(str_detect(Description, regex(pattern, ignore_case = TRUE)))
  cat(">>> 已执行【仅保留】模式，关键词：", pattern, "\n")
  
} else {
  # 模式A：剔除匹配黑名单的 (默认)
  if (length(exclude_keywords) > 0) {
    pattern <- paste(exclude_keywords, collapse = "|")
    kegg_filtered <- kegg_df %>% filter(!str_detect(Description, regex(pattern, ignore_case = TRUE)))
    cat(">>> 已执行【剔除】模式，去除了包含以下词的通路：", pattern, "\n")
  } else {
    kegg_filtered <- kegg_df
  }
}

cat(sprintf(">>> 筛选后剩余 %d 条通路。\n", nrow(kegg_filtered)))

# 保存最终筛选结果
write.table(kegg_filtered, file = output_file, sep = "\t", quote = F, row.names = F)


# --- 5. 绘图 ---

if (nrow(kegg_filtered) > 0) {
  cat(">>> 正在绘图...\n")
  
  plot_data <- kegg_filtered %>%
    mutate(logP = -log10(p.adjust)) %>%
    arrange(p.adjust) %>%
    head(top_n_plot) %>%
    mutate(Description = fct_reorder(Description, logP))
  
  p_final <- ggplot(plot_data, aes(x = logP, y = Description)) +
    geom_segment(aes(x = 0, y = Description, xend = logP, yend = Description),
                 color = plot_color,
                 linewidth = 0.6) +
    geom_point(color = plot_color,
               size = 7) +
    geom_text(aes(label = Count),
              color = "white",
              size = 3.5,
              fontface = "bold") +
    labs(x = bquote(-log[10]("Adjusted P-value")),
         y = NULL,
         title = "KEGG Pathways (Filtered)") +
    theme_classic(base_size = 16) +
    theme(
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(family = "arial", face = "bold", size = 16, margin = margin(t = 10)),
      axis.text = element_text(family = "arial", color = "black", size = 14),
      plot.title = element_text(family = "arial", hjust = 0.5, face = "bold", size = 20, margin = margin(b = 15)),
      legend.position = "none"
    ) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  print(p_final)
  
  # 如果需要自动保存图片，取消下面注释
  # ggsave("KEGG_Plot.pdf", width = 10, height = 7)
  
} else {
  warning("筛选后没有剩余通路，无法绘图。请检查筛选关键词是否过于严格。")
}

cat(">>> 流程结束。\n")