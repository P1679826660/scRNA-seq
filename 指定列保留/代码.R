# Check if packages are installed, install if not
if (!require("readxl")) install.packages("readxl")
if (!require("writexl")) install.packages("writexl")

# Load libraries
library(readxl)
library(writexl)

# ================= 配置区域 =================
# 请在这里修改你的文件名和列名
excel_file  <- "OEvsWT_deg_all.xls"        # 输入的Excel文件名
txt_file    <- "save.txt"   # 输入的txt文件名
target_col  <- "gene_name"        # Excel中用来比对的那一列的列名
output_file <- "result.xlsx"     # 输出结果的文件名
# ===========================================

# 1. 读取 Excel 文件
# suppressMessages 用于减少不必要的控制台输出
raw_data <- read.delim(excel_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# 2. 读取 TXT 文件
# readLines 会把每一行读取为一个字符串
match_list <- readLines(txt_file, warn = FALSE)

# 3. 数据清洗
# 去除 TXT 名单中可能存在的首尾空格
match_list <- trimws(match_list)
# 去除空行（如果有）
match_list <- match_list[match_list != ""]

# 4. 执行筛选
# 必须确保 Excel 中那一列也是字符类型，否则数字和字符串无法匹配
# raw_data[[target_col]] 用于提取指定列的数据
col_values <- as.character(raw_data[[target_col]])

# 判断 col_values 中的每个元素是否在 match_list 中
# %in% 操作符返回 TRUE 或 FALSE
is_matched <- col_values %in% match_list

# 根据判断结果筛选行
final_result <- raw_data[is_matched, ]

# 5. 保存结果
write_xlsx(final_result, output_file)

# 6. 输出简报
cat("处理完成。\n")
cat(sprintf("原始数据行数: %d\n", nrow(raw_data)))
cat(sprintf("txt名单数量: %d\n", length(match_list)))
cat(sprintf("匹配成功行数: %d\n", nrow(final_result)))
cat(sprintf("结果已保存至: %s\n", output_file))

