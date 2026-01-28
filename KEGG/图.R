# ==============================================================================
# KEGG 自动化分析流程：一键双图 (棒棒糖图 + 网络图)
# 功能：自动分析、自动去除疾病通路、自动保存表格和两张图片
# ==============================================================================

# --- 0. 加载必要的包 ---
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot) # 用于画网络图
library(forcats)
library(showtext)

# 配置字体 (如果报错找不到字体，可以注释掉这部分)
font_add("arial", "arial.ttf") 
showtext_auto()

# ==============================================================================
# 步骤 1: 数据读取与 KEGG 分析
# ==============================================================================
cat(">>> 1. 正在读取数据并转换 ID...\n")

# 1.1 读取数据
rt <- read.table("要区分上下调.txt", sep="\t", check.names=F)
genes <- unique(as.vector(rt[,1]))

# 1.2 ID 转换 (Symbol -> Entrez ID)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
rt_map <- data.frame(genes, entrezID=entrezIDs)
gene_list <- entrezIDs[entrezIDs != "NA"]

# 1.3 运行 KEGG (放宽阈值，先全抓下来，后面再筛)
options(timeout = 600)
cat(">>> 2. 正在连接 KEGG 数据库进行富集分析...\n")
kk <- enrichKEGG(gene = gene_list, 
                 organism = "hsa", 
                 pvalueCutoff = 1, 
                 qvalueCutoff = 1)

if (is.null(kk) || nrow(kk) == 0) {
  stop("错误：没有富集到任何通路，请检查基因名是否正确。")
}

# 1.4 ID 转回 Symbol (为了生成可读的表格)
# 注意：这里我们保留 kk 对象用于网络图，另外生成一个 dataframe 用于表格和棒棒糖图
kk_readable <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
KEGG_df <- as.data.frame(kk_readable)

# ==============================================================================
# 步骤 2: 自动化筛选 (剔除疾病通路)
# ==============================================================================
cat(">>> 3. 正在自动筛选通路 (剔除 Human Diseases)...\n")

# 设定显著性阈值
pvalueFilter <- 0.05
qvalueFilter <- 0.05

# 筛选逻辑：
# 1. pvalue 和 qvalue 小于 0.05
# 2. Count (基因数) 大于等于 3
# 3. 剔除 ID 包含 "hsa05" 的通路 (05xxx 代表疾病)
# 4. 剔除 ID 包含 "hsa06" 的通路 (06xxx 代表药物开发，通常也不要)

clean_result <- KEGG_df %>%
  filter(pvalue < pvalueFilter & qvalue < qvalueFilter) %>%
  filter(Count >= 3) %>%
  # 使用 ID 进行精确剔除，这比用文字匹配 "Human Diseases" 更准
  filter(!grepl("hsa05\\d{3}", ID)) %>% 
  filter(!grepl("hsa06\\d{3}", ID)) %>%
  arrange(p.adjust)

# 保存筛选后的结果
write.table(clean_result, file="KEGG_Filtered_Result.txt", sep="\t", quote=F, row.names=F)
cat(paste0(">>> 筛选完成！剩余通路数: ", nrow(clean_result), "。结果已保存。\n"))

# ==============================================================================
# 步骤 3: 绘制图一 —— 棒棒糖图 (Lollipop Plot)
# ==============================================================================
cat(">>> 4. 正在绘制图一：棒棒糖图...\n")

if (nrow(clean_result) > 0) {
  
  top_n_lollipop <- 30 # 棒棒糖图显示前30个
  
  plot_data <- clean_result %>%
    mutate(logP = -log10(p.adjust)) %>%
    arrange(p.adjust) %>%
    head(top_n_lollipop) %>%
    mutate(Description = fct_reorder(Description, logP))
  
  plot_color <- "#F8A2A4"
  
  p1 <- ggplot(plot_data, aes(x = logP, y = Description)) +
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
         title = "KEGG Pathways") +
    theme_classic(base_size = 16) +
    theme(
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(family = "arial", face = "bold", size = 16, margin = margin(t = 10)),
      axis.text = element_text(family = "arial", color = "black", size = 12), # 字体稍微调小一点防重叠
      plot.title = element_text(family = "arial", hjust = 0.5, face = "bold", size = 20, margin = margin(b = 15)),
      legend.position = "none"
    ) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  # 保存图一
  ggsave("KEGG_Plot1_Lollipop.pdf", plot = p1, width = 10, height = 8)
  print(p1)
  
} else {
  cat("没有符合条件的通路，无法绘制棒棒糖图。\n")
}

# ==============================================================================
# 步骤 4: 绘制图二 —— 网络图 (Network Plot)
# ==============================================================================
cat(">>> 5. 正在绘制图二：网络图...\n")

if (nrow(clean_result) > 0) {
  
  # 技巧：我们将筛选后的 dataframe 塞回 kk 对象里
  # 这样 cnetplot 就会只画我们筛选剩下的通路
  kk_readable@result <- clean_result
  
  # 提取前 N 个通路的名称用于画图
  top_n_network <- 15 # 网络图不能太多，太挤了看不清，建议 10-15 个
  top_pathways_names <- head(clean_result$Description, top_n_network)
  
  tryCatch({
    p2 <- cnetplot(kk_readable,
                   showCategory = top_pathways_names, 
                   color_gene = "#CCCCCC",    # 基因点：浅灰
                   color_category = "#F8A2A4",# 通路点：粉红
                   node_label = "category",   # 只显示通路名，不显示基因名
                   circular = FALSE,
                   layout = "kk")
    
    # 保存图二
    ggsave("KEGG_Plot2_Network.pdf", plot = p2, width = 12, height = 10)
    print(p2)
    p3 <- cnetplot(kk_readable,
                   showCategory = top_pathways_names, 
                   color_gene = "#CCCCCC",    # 基因点：浅灰
                   color_category = "#F8A2A4",# 通路点：粉红
                   #node_label = "category",   # 只显示通路名，不显示基因名
                   circular = FALSE,
                   layout = "kk")
    
    # 保存图二
    ggsave("KEGG_Plot3_Network.pdf", plot = p2, width = 12, height = 10)
    print(p3)
    cat(">>> 所有绘图完成！PDF文件已生成。\n")
    
  }, error = function(e) {
    cat("网络图绘制失败，可能是通路名过长或网络结构问题。\n")
    print(e)
  })
  
} else {
  cat("没有符合条件的通路，无法绘制网络图。\n")
}
