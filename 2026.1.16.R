# @Author: Zhaikai Liu (Optimized for Desktop)
# @Timestamp: 2026-01-26

# 1. 路径自动配置 ----------------------------------------------------------
# 自动获取当前用户的桌面路径（兼容Windows/Mac/Linux）
desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
if (desktop_path == "Desktop") { # 兼容部分系统环境
  desktop_path <- "~/Desktop"
}

input_dir <- file.path(desktop_path, "数据")
output_root <- file.path(desktop_path, "001")

# 定义模型子文件夹
models <- c("Dominant_Model", "Heterozygote_Model", "Homozygote_Model", "Recessive_Model", "Allele_Model")

# 创建输出目录结构
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
for (m in models) {
  subdir <- file.path(output_root, m)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
}

# 2. 环境准备 --------------------------------------------------------------
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

load_package <- function(p) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

load_package("meta")

# 3. 数据读取 --------------------------------------------------------------
# 假设文件名为固定名称，您可以根据实际文件名修改
main_file <- file.path(input_dir, "mian.csv")
sub_file <- file.path(input_dir, "subgroup_children_adults.csv")
sev_file <- file.path(input_dir, "severity.csv")

# Allele模型专用的衍生文件名
main_allele_file <- file.path(input_dir, "mian.csv_allele.csv")
sub_allele_file <- file.path(input_dir, "subgroup_children_adults.csv_allele.csv")
sev_allele_file <- file.path(input_dir, "severity.csv_allele.csv")

Main <- read.csv(main_file, header = TRUE, dec = ".")
Sub <- read.csv(sub_file, header = TRUE, dec = ".")
Severity <- read.csv(sev_file, header = TRUE, dec = ".")

# 4. 自动化分析函数 ----------------------------------------------------------
# 定义一个函数来减少重复代码，提高可维护性
run_meta_analysis <- function(data_obj, ev_e, n_e, ev_c, n_c, model_name, sub_folder, lab_e="Exposure", lab_c="Control", is_subgroup=FALSE, by_var=NULL) {
  
  # 执行Meta分析
  m_bin <- metabin(ev_e, n_e, ev_c, n_c, data = data_obj,
                   method = "I", sm = "OR", 
                   common = TRUE, random = TRUE,
                   byvar = by_var,
                   studlab = paste(data_obj$Author, data_obj$Year))
  
  base_path <- file.path(output_root, sub_folder)
  
  # 1. Forest Plot
  pdf(file.path(base_path, paste0(model_name, " forest.pdf")), height=6, width=12)
  forest(m_bin, lab.e=lab_e, lab.c=lab_c)
  dev.off()
  
  # 2. Sensitivity (Random)
  pdf(file.path(base_path, paste0(model_name, " sensitivity_random.pdf")), height=5, width=10)
  forest(metainf(m_bin, pooled="random"), lab.e=lab_e, lab.c=lab_c)
  dev.off()
  
  # 3. Funnel Plot
  pdf(file.path(base_path, paste0(model_name, " funnel.pdf")), height=8, width=8)
  funnel(m_bin)
  dev.off()
  
  # 4. Galbraith Plot
  pdf(file.path(base_path, paste0(model_name, " galbraith.pdf")), height=8, width=8)
  radial(m_bin, text=data_obj$Author, level=0.95)
  dev.off()
}

# 5. 执行各模型分析 ----------------------------------------------------------

## --- Dominant Model ---
run_meta_analysis(Main, Main$AG_GG_case, Main$AG_GG_case+Main$AG_case, 
                  Main$AG_GG_control, Main$AG_GG_control+Main$AA_control, 
                  "Dominant", "Dominant_Model")

## --- Heterozygote Model ---
run_meta_analysis(Main, Main$AG_case, Main$AA_AG_case, 
                  Main$AG_control, Main$AA_AG_control, 
                  "Heterozygote", "Heterozygote_Model")

## --- Homozygote Model ---
run_meta_analysis(Main, Main$GG_case, Main$GG_case+Main$AA_case, 
                  Main$GG_control, Main$GG_control+Main$AA_control, 
                  "Homozygote", "Homozygote_Model")

## --- Recessive Model ---
run_meta_analysis(Main, Main$GG_case, Main$GG_case+Main$AA_AG_case, 
                  Main$GG_control, Main$GG_control+Main$AA_AG_control, 
                  "Recessive", "Recessive_Model")

## --- Allele Model (需读取专用文件) ---
if(file.exists(main_allele_file)){
  MainA <- read.csv(main_allele_file, header=TRUE)
  run_meta_analysis(MainA, MainA$G_case, MainA$Sample_case, 
                    MainA$A_control, MainA$Sample_control, 
                    "Allele", "Allele_Model")
}

cat("所有分析已完成，请检查桌面 '001' 文件夹。\n")