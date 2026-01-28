# @Author: WenCong Wu (Optimized for Desktop)
# @Timestamp: 2026-01-16

# 1. 环境与包加载
options(echo = TRUE)
if (!"meta" %in% installed.packages()[, "Package"]) install.packages("meta")
library(meta)

# 2. 路径设置
desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
input_dir    <- file.path(desktop_path, "数据")
output_dir   <- file.path(desktop_path, "002")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 智能文件查找函数
find_file <- function(pattern) {
  files <- list.files(input_dir, full.names = TRUE)
  matched <- files[grep(paste0("^", pattern, "$"), basename(files), ignore.case = TRUE)]
  if (length(matched) == 0) stop(paste("找不到文件:", pattern))
  return(matched[1])
}

# 3. 核心分析函数 (补全了所有漏掉的绘图逻辑)
run_full_analysis <- function(data_list, model_label, cols) {
  target_dir <- file.path(output_dir, model_label)
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
  
  # --- A. 基础 Meta 分析模块 ---
  df <- data_list$Main
  m <- metabin(df[[cols$e]], df[[cols$ne]], df[[cols$c]], df[[cols$nc]], 
               data = df, method = "I", sm = "OR", common = TRUE, random = TRUE,
               studlab = paste(df$Author, df$Year))
  
  # 1. 基础森林图
  pdf(file.path(target_dir, paste0(model_label, " forest.PDF")), height=5, width=12)
  forest(m, lab.e="Exposure", lab.c="Control")
  dev.off()
  
  # 2. 敏感性分析 (补回)
  pdf(file.path(target_dir, paste0(model_label, " sensetivity_analysis_random.PDF")), height=4, width=10)
  forest(metainf(m, pooled = "random"))
  dev.off()
  
  # 3. 漏斗图
  pdf(file.path(target_dir, paste0(model_label, " funnel.PDF")), height=7, width=7)
  funnel(m)
  dev.off()
  
  # 4. 甘氏图 (补回)
  pdf(file.path(target_dir, paste0(model_label, " galbraith.PDF")), height=7, width=7)
  radial(m, text=df$Author, level=0.95)
  dev.off()
  
  # --- B. 亚组分析模块 ---
  # 国家亚组
  pdf(file.path(target_dir, paste0(model_label, " subgroup_analysis_Country.PDF")), height=9, width=14)
  forest(update(m, byvar = df$Country))
  dev.off()
  
  # 年龄亚组 (基于 Sub 列表)
  sub_df <- data_list$Sub
  m_age <- metabin(sub_df[[cols$e]], sub_df[[cols$ne]], sub_df[[cols$c]], sub_df[[cols$nc]], 
                   data = sub_df, method = "I", sm = "OR", byvar = sub_df$Sample_age,
                   studlab = paste(sub_df$Author, sub_df$Year))
  pdf(file.path(target_dir, paste0(model_label, " subgroup_analysis_Sample_age.PDF")), height=9, width=14)
  forest(m_age)
  dev.off()
  
  # --- C. 严重程度分析模块 (补全所有相关图表) ---
  sev_df <- data_list$Severity
  m_sev <- metabin(sev_df[[cols$se]], sev_df[[cols$sne]], sev_df[[cols$sc]], sev_df[[cols$snc]], 
                   data = sev_df, method = "I", sm = "OR", 
                   studlab = paste(sev_df$Author, sev_df$Year))
  
  # 1. 严重程度森林图
  pdf(file.path(target_dir, paste0(model_label, " forest_severity.PDF")), height=5, width=12)
  forest(m_sev, lab.e="Not_Severe", lab.c="Severe")
  dev.off()
  
  # 2. 严重程度敏感性分析 (补回)
  pdf(file.path(target_dir, paste0(model_label, " sensetivity_analysis_random_severity.PDF")), height=4, width=10)
  forest(metainf(m_sev, pooled = "random"))
  dev.off()
  
  # 3. 严重程度漏斗图 (补回)
  pdf(file.path(target_dir, paste0(model_label, " funnel_severity.PDF")), height=7, width=7)
  funnel(m_sev)
  dev.off()
  
  # 4. 严重程度甘氏图 (补回)
  pdf(file.path(target_dir, paste0(model_label, " galbraith_severity.PDF")), height=7, width=7)
  radial(m_sev, text=sev_df$Author, level=0.95)
  dev.off()
}

# 4. 读取数据（适配您奇怪的文件名）
data_genotype <- list(
  Main     = read.csv(find_file("mian.CSV")),
  Sub      = read.csv(find_file("subgroup_children_adults.CSV")),
  Severity = read.csv(find_file("severity.CSV"))
)

# 常规模型配置
configs <- list(
  "Dominant_Model"     = list(e="AG_GG_case", ne="AG_GG_case", c="AG_GG_control", nc="AA_control", 
                              se="AG_GG_not_severe", sne="AG_GG_not_severe", sc="AG_GG_severe", snc="AA_severe"),
  "Heterozygote_Model" = list(e="AG_case", ne="AA_AG_case", c="AG_control", nc="AA_AG_control",
                              se="AG_not_severe", sne="AA_AG_not_severe", sc="AG_severe", snc="AA_AG_severe"),
  "Homozygote_Model"   = list(e="GG_case", ne="AA_case", c="GG_control", nc="AA_control",
                              se="GG_not_severe", sne="AA_not_severe", sc="GG_severe", snc="AA_severe"),
  "Recessive_Model"    = list(e="GG_case", ne="AA_AG_case", c="GG_control", nc="AA_AG_control",
                              se="GG_not_severe", sne="AA_AG_not_severe", sc="GG_severe", snc="AA_AG_severe")
)

# 运行前4个模型
for (label in names(configs)) run_full_analysis(data_genotype, label, configs[[label]])

# 5. Allele 模型特殊读取
data_allele <- list(
  Main     = read.csv(find_file("mian_allele.csv.CSV")),
  Sub      = read.csv(find_file("subgroup_children_adults_allele.csv.CSV")),
  Severity = read.csv(find_file("severity_allele.csv.CSV"))
)

allele_cols <- list(e="G_case", ne="Sample_case", c="A_control", nc="Sample_control",
                    se="G_not_severe", sne="Sample_not_severe", sc="A_severe", snc="Sample_severe")

run_full_analysis(data_allele, "Allele_Model", allele_cols)

cat("\n>>> [全部完成] 每个模型已输出约 10 个 PDF 文件至桌面 002 文件夹。\n")