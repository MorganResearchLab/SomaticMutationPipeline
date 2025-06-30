#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install CRAN packages first
cran_packages <- c("optparse", "stringr", "reshape2", "dplyr", "devtools")

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "from CRAN\n")
        install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    } else {
        cat(pkg, "already installed\n")
    }
}

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager\n")
    install.packages("BiocManager", quiet = TRUE)
}

# Set Bioconductor version to 3.21 (for R 4.5.1)
cat("Setting Bioconductor version to 3.21\n")
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)

# Install Bioconductor packages
bioc_packages <- c(
    "VariantAnnotation",
    "GenomicFeatures", 
    "GenomicRanges",
    "BSgenome.Hsapiens.UCSC.hg38",
    "biomaRt",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "Biostrings",
    "org.Hs.eg.db",
    "Homo.sapiens"
)

cat("Installing Bioconductor packages...\n")
for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "from Bioconductor\n")
        BiocManager::install(pkg, ask = FALSE, update = FALSE, dependencies = TRUE)
    } else {
        cat(pkg, "already installed\n")
    }
}

cat("All packages installed successfully!\n")

# Test loading key packages
cat("Testing package loading...\n")
test_packages <- c("optparse", "stringr", "dplyr", "VariantAnnotation", "GenomicRanges")
for (pkg in test_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("✓", pkg, "loads successfully\n")
    } else {
        cat("✗", pkg, "failed to load\n")
    }
}

message("All done")