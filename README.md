# SomaticMutationPipeline
A comprehensive Snakemake pipeline for calling somatic mutations from single-cell RNA sequencing data using SComatic.

## Overview
This pipeline processes single-cell RNA-seq BAM files to identify somatic mutations by:

1. Splitting BAM files by cell type and chromosome
2. Counting bases at each genomic position
3. Calling variants using statistical methods
4. Filtering out germline variants and technical artifacts
5. Genotyping individual cells for identified variants

## Requirements

* Conda/Mamba for environment management
* Snakemake (≥7.0) with SLURM executor plugin
* Python 3.7+
* R 4.2+ with Bioconductor packages
* Access to reference genome files (hg38)
* SComatic repository (see below)

## Installation
### 1. Clone the Repository
```
git clone https://github.com/MorganResearchLab/SomaticMutationPipeline.git
cd SomaticMutationPipeline
```

### 2. Create the Main Environment
```
# Create the main environment
conda env create -f environment.yaml

# Activate the environment
conda activate somatic-mutation-pipeline
```

### 3. Install Snakemake SLURM Plugin
__NB:__ Select the executor plugin appropriate to your HPC. The example below is for SLURM.
A list of all snakemake executor plugins can be found online: https://github.com/snakemake/snakemake-interface-executor-plugins

```bash
# Install the SLURM executor plugin for Snakemake
pip install snakemake-executor-plugin-slurm
```

### 4. Clone SComatic Repository

If you don't have access to a shared SComatic installation, clone it:

```bash
git clone https://github.com/cortes-ciriano-lab/SComatic.git
```

## Configuration Setup

### Main Configuration File: `snakemake.dir/config.yaml`

This is the primary configuration file you need to modify for your analysis:

#### **Required Paths to Update:**

```yaml
# Update these paths for your system:
conda_prefix: "/path/to/your/miniforge3/envs"  # Path to your conda environments
src_dir: "/path/to/SomaticMutationPipeline/scripts"  # Path to this pipeline's scripts
scomatic_dir: "/path/to/SComatic"  # Path to SComatic repository
ref_dir: "/path/to/reference/chromosomes"  # Directory containing genome (e.g. hg38) chromosome files (chr1.fa.gz, chr2.fa.gz, etc.)
bam_dir: "/path/to/your/bam/files"  # Directory containing your input BAM files
```

#### **Required Input Files:**

```yaml
# Cell type annotation file (TSV format)
cell_type_file: "/path/to/cell_metadata.tsv"
# Required columns: Index (cell barcode), Cell_type (cell type annotation)

# Sample metadata file (CSV format)
sra_metadata: "/path/to/sample_metadata.csv"
# Should contain sample IDs and batch information
```

#### **Analysis Parameters:**

```yaml
# Batch handling
batch_position: "None"  # Options: "prefix", "suffix", "None"
batch_column: "batch"   # Column name in sra_metadata containing batch info

# Cell types to analyze (comma-separated, no spaces)
celltypes: "B,CD4,CD8,Monocytes,NK"  # Customize for your data

# Parallel processing
window: 5000000  # Genomic chunk size for parallel BAM splitting (5Mb recommended)

# SComatic filtering parameters (usually don't need to change)
min_cells: 5
min_ac_cells: 2
min_ac_reads: 3
max_cell_types: 6
min_cell_types: 1
fisher_cutoff: 0.001
```

#### **Reference Files (Pre-configured):**

These paths point to SComatic's reference files and typically don't need modification:

```yaml
panel_of_normals: "/path/to/SComatic/PoNs/PoN.scRNAseq.hg38.tsv"
rna_editing_sites: "/path/to/SComatic/RNAediting/AllEditingSites.hg38.txt"
```

### SLURM Profile Configuration: `snakeprofile/config.yaml`

Configure for your HPC cluster:

```yaml
executor: slurm
jobs: 100                    # Maximum concurrent jobs
latency-wait: 250           # Wait time for file system

default-resources:
	slurm_partition: "your-partition"  # Your SLURM partition
	slurm_account: "your-account"      # Your SLURM account
	runtime: "12d"                     # Default job runtime
	mem_mb: 4000                      # Default memory
	cpus_per_task: 1                  # Default CPUs
	nodes: 1                          # Default nodes
```

### Job-Specific Resources: `snakemake.dir/Maxwell.json`

This file contains resource specifications for each pipeline step. You may need to adjust based on your cluster's limits:

```json
{
    "split_by_celltype": {
	"memory": 12000,      # Memory in MB
	"time": "95:00:00",   # Time limit
	"nCPUs": 8            # CPU cores
},
    "base_counter": {
        "memory": 12000,
		"time": "47:00:00"
},
    // ... other rules
}
```

## Input Data Requirements

### 1. BAM Files
- **Location**: Place in the directory specified by `bam_dir`
- **Format**: Standard BAM files with cell barcode (CB) tags
- **Naming**: Files should be named consistently (e.g., `Sample1.bam`, `Sample2.bam`)

### 2. Cell Type Metadata File
- **Format**: Tab-separated values (.tsv)
- **Required columns**:
  - `Index`: Cell barcodes (must match CB tags in BAM files)
    - `Cell_type`: Cell type annotations
    - **Example for 10X Genomics data**:
    ```
    Index	Cell_type
    AAACCTGAGAAACCAT-1	CD4
    AAACCTGAGAAACGAG-1	B
    AAACCTGAGAAATCCC-1	Monocytes
    ```

### 3. Sample Metadata File (Optional)
- **Format**: Comma-separated values (.csv)
- **Purpose**: Maps sample names to batch IDs
- **Example**:
```
Sample,batch
Sample1,Batch_A
Sample2,Batch_B
```

### 4. Reference Genome Files
- **Required**: genome build chromosome FASTA files (compressed)
- **Format**: `chr1.fa.gz`, `chr2.fa.gz`, ..., `chr22.fa.gz`
- **Location**: Directory specified by `ref_dir`

## Running the Pipeline

### 1. Test Configuration (Dry Run)

```bash
# Test that everything is configured correctly
snakemake -s snakemake.dir/snakefile \
    --configfile snakemake.dir/config.yaml \
        -n --rerun-triggers mtime
	```

### 2. Create Conda Environments (Optional)

Pre-create all required environments:

```bash
snakemake --sdm conda --conda-create-envs-only \
    -s snakemake.dir/snakefile \
        --configfile snakemake.dir/config.yaml
	```

### 3. Execute the Pipeline

```bash
nohup nice -19 snakemake \
    -s snakemake.dir/snakefile \
        --configfile snakemake.dir/config.yaml \
	    --sdm conda \
	        --profile snakeprofile \
		    --latency-wait 650 \
		        --slurm-init-seconds-before-status-checks 500 \
			    -j 25 \
			        --rerun-triggers mtime \
				    --rerun-incomplete &
				    ```

### 4. Monitor Progress

```bash
# Watch the pipeline progress
tail -f -n20 nohup.out

# Check SLURM job status
squeue -u $USER
```

## Output Structure

The pipeline creates the following output directories:

```
results/
├── split_bams/           # Cell-type and chromosome-specific BAM files
├── base_counts/          # Base counts for each position
├── merged_counts/        # Merged counts across cell types
├── variant_calls/        # Initial variant calls
├── filtered_calls/       # Filtered variant calls
├── qc_filtered/          # Quality-filtered variants
├── genotyped/           # Single-cell genotypes
└── germline/            # Final somatic variants (germline filtered)
```

### Key Output Files

- **`results/germline/{sample}_{celltype}_not_germline.tsv`**: Final somatic mutations for each sample and cell type
- **`results/genotyped/{sample}_{celltype}_{chrom}.single_cell_genotype.tsv`**: Single-cell genotype data
- **Log files**: Located in `logs/` directory for troubleshooting

## Troubleshooting

### Common Issues

1. **Environment Creation Fails**
   - Ensure you have conda/mamba installed
      - Try using mamba instead of conda for faster environment resolution

2. **BAM Files Not Found**
   - Verify `bam_dir` path is correct
      - Check file permissions
         - Ensure BAM files have been index, i.e. there is a *bam.bai for each BAM file	 

3. **Reference Files Missing**
   - Download genome build chromosome files
      - Ensure files are compressed (.gz) and properly indexed

4. **Memory/Time Limits Exceeded**
   - Adjust resources in `snakemake.dir/hpc.json`
      - Consider reducing `window` size for smaller memory usage

5. **Cell Barcodes Don't Match**
   - Check cell barcode format in metadata vs BAM files
      - Verify batch handling settings (`batch_position`, `batch_column`)

### Log Files

Check specific log files for detailed error messages:
- `logs/{rule_name}/{sample}_{chrom}.log`
- Main pipeline log: `nohup.out`

## Performance Optimization

### For Large Datasets:
- Increase `window` size for fewer, larger parallel jobs
- Adjust memory allocations based on your data size
- Consider using faster storage for temporary files

### For Many Samples:
- Increase `jobs` parameter in SLURM profile
- Use cluster's high-memory nodes for memory-intensive steps

## Citation

If you use this pipeline, please cite:
- **SComatic**: Muyas et al. (2023) Nature Genetics, https://doi.org/10.1038/s41587-023-01863-z

## Support

For questions and issues:
- **GitHub Issues**: https://github.com/MorganResearchLab/SomaticMutationPipeline/issues

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.~
