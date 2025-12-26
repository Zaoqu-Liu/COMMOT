# 安装修改后的 fastCNV 包
# 确保使用 Seurat 4.4.0 和 SeuratObject 4.1.4

# 1. 首先确保安装 remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# 2. 安装 Seurat 4.4.0 和 SeuratObject 4.1.4（如果需要）
if (!requireNamespace("SeuratObject", quietly = TRUE) || 
    packageVersion("SeuratObject") != "4.1.4") {
  message("Installing SeuratObject 4.1.4...")
  remotes::install_version("SeuratObject",
    version = "4.1.4",
    repos = c(
      "https://satijalab.r-universe.dev",
      "https://cloud.r-project.org"
    ),
    upgrade = "never"
  )
}

if (!requireNamespace("Seurat", quietly = TRUE) || 
    packageVersion("Seurat") != "4.4.0") {
  message("Installing Seurat 4.4.0...")
  remotes::install_version("Seurat",
    version = "4.4.0",
    repos = c(
      "https://satijalab.r-universe.dev",
      "https://cloud.r-project.org"
    ),
    upgrade = "never"
  )
}

# 3. 安装修改后的 fastCNV
message("Installing modified fastCNV from local path...")
remotes::install_local(
  path = "/Users/liuzaoqu/Downloads/fastCNV",
  upgrade = "never",
  force = TRUE
)

# 4. 验证安装
if (requireNamespace("fastCNV", quietly = TRUE)) {
  message("\n=== fastCNV installed successfully! ===")
  message("fastCNV version: ", packageVersion("fastCNV"))
  message("Seurat version: ", packageVersion("Seurat"))
  message("SeuratObject version: ", packageVersion("SeuratObject"))
  
  # 测试能否加载数据
  tryCatch({
    data <- fastCNV::geneMetadata
    message("✓ geneMetadata accessible: ", nrow(data), " genes")
  }, error = function(e) {
    message("✗ Warning: geneMetadata not accessible")
  })
} else {
  message("✗ Installation failed!")
}

