# ==============================================================================
# 模块 0：环境加载与工作目录设置
# ==============================================================================

library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(qs) # 用于快速读写数据

# 1. 创建结果输出目录 (所有图片和文件都会保存在这个文件夹里，保持根目录整洁)
out_dir <- "CellChat_Results"
if(!dir.exists(out_dir)) dir.create(out_dir)

# 2. 设置并行运算 (可选，如果报错可注释掉)
# future::plan("multisession", workers = 4) 
options(stringsAsFactors = FALSE)

# ==============================================================================
# 模块 1：创建 CellChat 对象与数据库设置
# ==============================================================================

# [用户需修改]: 确认你的 Seurat 对象名，以及细胞类型所在的列名 (group.by)
# 假设你的 Seurat 对象叫 seurat_obj，细胞类型列叫 "cell_type"
#cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")

# 1. 提取数据矩阵（Seurat 5 用 layer）
data.input <- GetAssayData(sc.obj, assay = "RNA", layer = "data") 
# 2. 提取元数据
meta <- sc.obj@meta.data 
# 3. 通过矩阵创建 CellChat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")


# 检查细胞类型对应的细胞数量 (如果不平衡，太少的细胞群可能会报错)
groupSize <- as.numeric(table(cellchat@idents))
print(table(cellchat@idents)) 

# 设置配体-受体数据库
# [用户需修改]: 如果是小鼠，请改为 CellChatDB.mouse
CellChatDB.use <- CellChatDB.human 

# 这里的 search参数 控制分析范围：
# "Secreted Signaling": 分泌信号 (最常用)
# "ECM-Receptor": 细胞外基质受体
# "Cell-Cell Contact": 细胞接触
CellChatDB.use <- subsetDB(CellChatDB.use, search = "Secreted Signaling")

# 将筛选后的数据库赋值给对象
cellchat@DB <- CellChatDB.use

# ==============================================================================
# 模块 2：预处理与核心计算 (标准流程，通常无需修改)
# ==============================================================================

message("正在进行预处理...")

# 1. 提取表达数据 (subsetData) 并处理
# 这一步是为了提取信号通路相关的基因数据，减少计算量
cellchat <- subsetData(cellchat) 

# 2. 识别高表达基因和互作
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 3. 推断通讯概率 (核心步骤)
# type = "triMean" 是较稳健的计算方法，减少极值影响
# 这一步根据数据大小可能需要几分钟
cellchat <- computeCommunProb(cellchat, type = "triMean")

# 4. 过滤低质量通讯
# min.cells = 10 表示：如果一个细胞群参与通讯的细胞少于10个，则忽略该通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 5. 推断信号通路层级的通讯概率
cellchat <- computeCommunProbPathway(cellchat)

# 6. 整合网络 (将所有配体-受体对聚合为通路)
cellchat <- aggregateNet(cellchat)

message("计算完成，开始生成并保存图表...")

# ==============================================================================
# 模块 3：可视化与自动保存 (结果输出)
# ==============================================================================

# --------------------------------------------------------
# 3.1 全局概览图 (保存为 PDF)
# --------------------------------------------------------
pdf(file.path(out_dir, "01_Global_Interactions.pdf"), width = 10, height = 8)

# 这里的 par 设置是为了在一页画两张图
par(mfrow = c(1, 2), xpd = TRUE)

# 图1：互作数量 (Number of interactions) - 连线越多，通讯关系越复杂
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of interactions")

# 图2：互作强度 (Interaction weights) - 连线越粗，通讯概率/强度越高
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction strength")

dev.off() # 关闭画板，保存文件

# --------------------------------------------------------
# 3.2 单个细胞发出的信号概览 (保存为 PDF)
# --------------------------------------------------------
# 这张图展示每种细胞作为“发送者”时，主要和谁通讯
mat <- cellchat@net$weight
pdf(file.path(out_dir, "02_Per_Cell_Type_Sender.pdf"), width = 12, height = 12)
par(mfrow = c(3, 4), xpd = TRUE) # 根据你的细胞类型数量调整行列，这里预设3行4列

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# --------------------------------------------------------
# 3.3 具体信号通路的详细分析 (气泡图与弦图)
# --------------------------------------------------------

# 第一步：查看有哪些显著的信号通路
# 运行这行代码在控制台查看结果，然后选一个你感兴趣的
print(cellchat@netP$pathways)

# [用户需修改]: 将下面的 "TGFb" 替换为你数据中存在的、感兴趣的通路名称
pathway.show <- "TGFb" 

# 检查该通路是否存在于结果中，避免报错
if (pathway.show %in% cellchat@netP$pathways) {
  
  # A. 层级图 (Hierarchy plot) 和 弦图 (Chord diagram)
  pdf(file.path(out_dir, paste0("03_Pathway_", pathway.show, "_Details.pdf")), width = 10, height = 8)
  
  # 画弦图
  p_chord <- netVisual_aggregate(cellchat, signaling = pathway.show, layout = "chord")
  print(p_chord) # PDF中需要显式print
  
  # 画热图 (Heatmap) - 展示谁发给谁最强
  p_heatmap <- netVisual_heatmap(cellchat, signaling = pathway.show, color.heatmap = "Reds")
  print(p_heatmap)
  
  dev.off()
  
  # B. 基因表达贡献图 (查看是哪对 Ligand-Receptor 起主要作用)
  pdf(file.path(out_dir, paste0("04_Pathway_", pathway.show, "_Contribution.pdf")), width = 8, height = 6)
  p_contrib <- netAnalysis_contribution(cellchat, signaling = pathway.show)
  print(p_contrib)
  dev.off()
  
} else {
  message(paste("警告:", pathway.show, "通路在数据中未检测到，请更换其他通路名称。"))
}

# --------------------------------------------------------
# 3.4 气泡图 (Bubble Plot) - 最直观的展示方式
# --------------------------------------------------------
# 展示所有细胞群之间的配体-受体对

pdf(file.path(out_dir, "05_All_Bubble_Plot.pdf"), width = 12, height = 18)

# remove.isolate = FALSE 表示即使某些细胞没有通讯也显示在图上
# [用户提示]: 如果图太大，可以通过 sources.use 和 targets.use 指定特定细胞群
p_bubble <- netVisual_bubble(cellchat, remove.isolate = FALSE)
print(p_bubble)

dev.off()

# ==============================================================================
# 模块 4：保存结果对象
# ==============================================================================

# 保存处理好的对象，方便下次直接读取，不用重新跑 computeCommunProb
qsave(cellchat, file.path(out_dir, "cellchat_final.qs"))

message("分析结束！所有结果已保存在 CellChat_Results 文件夹中。")

############################################################################
############################################################################
#Compare_Function.R
library(tidyverse)
library(CellChat)

# ==========================================================
# 1. 设置你的比较场景 (在这里修改名称即可复用于任何任务)
# ==========================================================
# 必须确保这些名字在 levels(cellchat@idents) 中完全一致

source_cell <- "LiverCancer"  # 谁在发信号？(肝癌细胞)
target_1    <- "CD8_T"        # 接收者1 (例如杀伤性T细胞)
target_2    <- "Treg"         # 接收者2 (例如调节性T细胞)

# ==========================================================
# 2. 自动化计算函数 (无需修改，直接运行)
# ==========================================================

compare_communication <- function(cellchat_obj, source, t1, t2) {
  
  # 获取所有计算过的通路名称
  all_pathways <- cellchat_obj@netP$pathways
  
  # 创建一个空列表来存结果
  results_list <- list()
  
  # 遍历每一个通路
  for (path in all_pathways) {
    # 提取该通路的三维概率矩阵 [发送者, 接收者, 通路]
    # 注意：cellchat@netP$prob 是核心数据所在
    # 有些旧版本可能是 cellchat@net$prob，如果报错请检查版本，通常 netP 用于通路聚合数据
    
    # 提取特定通路的矩阵
    prob_matrix <- cellchat_obj@netP$prob[,, path] 
    
    # 获取强度数值 (如果该通路不活跃，数值为0)
    # 逻辑：矩阵的行是发送者，列是接收者
    val_t1 <- prob_matrix[source, t1]
    val_t2 <- prob_matrix[source, t2]
    
    # 如果两个值都是0，说明这个通路在这三者间完全不工作，跳过
    if (sum(val_t1, val_t2) > 0) {
      results_list[[path]] <- data.frame(
        Pathway = path,
        Source = source,
        To_Target_1 = val_t1, # 发给 T细胞1 的强度
        To_Target_2 = val_t2, # 发给 T细胞2 的强度
        Difference = val_t1 - val_t2, # 差异
        Dominant = ifelse(val_t1 > val_t2, t1, t2) # 谁接收得更多
      )
    }
  }
  
  # 合并结果
  if(length(results_list) == 0) {
    return(NULL)
  }
  final_df <- do.call(rbind, results_list)
  return(final_df)
}

# ==========================================================
# 3. 执行分析与保存结果
# ==========================================================

# 运行计算
comparison_df <- compare_communication(cellchat, source_cell, target_1, target_2)

# 检查是否有结果
if (!is.null(comparison_df)) {
  
  # 按差异绝对值排序，把差异最大的通路排前面
  comparison_df <- comparison_df %>% arrange(desc(abs(Difference)))
  
  # 打印前几行看看
  print(head(comparison_df))
  
  # 保存为 CSV (这是最实用的东西，直接用 Excel 打开看)
  write.csv(comparison_df, paste0("Compare_", source_cell, "_to_", target_1, "_vs_", target_2, ".csv"), row.names = FALSE)
  
  message("计算完成！CSV表格已保存。")
  
} else {
  message("未检测到这两个细胞群之间的有效通讯差异。")
}

# ==========================================================
# 4. 绘图：直观展示差异 (双向柱状图)
# ==========================================================

if (!is.null(comparison_df)) {
  
  # 为了画图整洁，只取差异最大的前 20 个通路
  plot_data <- head(comparison_df, 20)
  
  p <- ggplot(plot_data, aes(x = reorder(Pathway, Difference), y = Difference, fill = Dominant)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() + # 翻转坐标轴，横着看方便
    theme_bw() +
    scale_fill_manual(values = c("red", "blue")) + # 你可以自定义颜色
    labs(
      title = paste0("Signal Difference: ", source_cell, " sending to..."),
      subtitle = paste0("Positive (Right) = Stronger in ", target_1, "\nNegative (Left) = Stronger in ", target_2),
      x = "Signaling Pathway",
      y = "Interaction Strength Difference"
    ) +
    theme(axis.text = element_text(size = 10, color = "black"))
  
  # 自动保存图片
  ggsave(paste0("Diff_Plot_", source_cell, ".pdf"), p, width = 8, height = 6)
  print(p)
}
