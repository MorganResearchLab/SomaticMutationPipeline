# ==============================================================================
# Rules for Variant Calling
# ==============================================================================

rule base_counter:
    input:
        bam = "results/split_bams/{sample}/{sample}_{celltype}_{chrom}.bam",
        bai = "results/split_bams/{sample}/{sample}_{celltype}_{chrom}.bam.bai"
    output:
        counts = "results/base_counts/{sample}/{chrom}/{sample}_{celltype}_{chrom}.tsv",
        temp = directory("results/temp/{sample}_{chrom}_{celltype}")
    params:
        fasta = lambda wildcards: f"{config['paths']['ref_dir']}/{wildcards.chrom}.fa.gz",
        out_folder = lambda wildcards: f"results/base_counts/{wildcards.sample}/{wildcards.chrom}/",
        script = os.path.join(config["paths"]["scomatic_dir"], "scripts/BaseCellCounter/BaseCellCounter.py")
    threads: 1
    resources:
        mem_mb = 12000,
        runtime = 2820,
        partition = "uoa-compute"
    conda:
        "../envs/scomatic.yaml"
    log:
        "logs/base_counter/{sample}_{chrom}_{celltype}.log"
    shell:
        """
        python {params.script} \
            --bam {input.bam} \
            --ref {params.fasta} \
            --chrom {wildcards.chrom} \
            --out_folder {params.out_folder} \
            --min_bq 30 \
            --tmp_dir {output.temp} \
            --nprocs {threads} \
            > {log} 2>&1
        """

rule merge_base_counts:
    input:
        expand("results/base_counts/{{sample}}/{{chrom}}/{{sample}}_{celltype}_{{chrom}}.tsv",
               celltype=CELL_TYPES)
    output:
        merged = "results/merged_counts/{sample}_{chrom}_AllCellTypes.tsv"
    params:
        tsv_dir = "results/base_counts/{sample}/{chrom}/",
        script = os.path.join(config["paths"]["scomatic_dir"], "scripts/MergeCounts/MergeBaseCellCounts.py")
    resources:
        mem_mb = 12000,
        runtime = 2820,
        partition = "uoa-compute"
    conda:
        "../envs/scomatic.yaml"
    log:
        "logs/merge_counts/{sample}_{chrom}.log"
    shell:
        """
        python {params.script} \
            --tsv_folder {params.tsv_dir} \
            --outfile {output.merged} \
            > {log} 2>&1
        """

rule variant_calls:
    input:
        "results/merged_counts/{sample}_{chrom}_AllCellTypes.tsv"
    output:
        "results/variant_calls/{sample}_{chrom}.calling.step1.tsv"
    params:
        ref = lambda wildcards: f"{config['paths']['ref_dir']}/{wildcards.chrom}.fa.gz",
        outfile_prefix = "results/variant_calls/{sample}_{chrom}",
        script = os.path.join(config["paths"]["scomatic_dir"], "scripts/BaseCellCalling/BaseCellCalling.step1.py"),
        # Assign each parameter directly to avoid the dictionary access error
        min_cells = config["params"]["variant_calling"]["min_cells"],
        min_ac_cells = config["params"]["variant_calling"]["min_ac_cells"],
        min_ac_reads = config["params"]["variant_calling"]["min_ac_reads"],
        max_cell_types = config["params"]["variant_calling"]["max_cell_types"],
        min_cell_types = config["params"]["variant_calling"]["min_cell_types"],
        fisher_cutoff = config["params"]["variant_calling"]["fisher_cutoff"]
    resources:
        mem_mb = 12000,
        runtime = 2820,
        partition = "uoa-compute"
    conda:
        "../envs/scomatic.yaml"
    log:
        "logs/variant_calls/{sample}_{chrom}.log"
    shell:
        """
        python {params.script} \
            --infile {input} \
            --outfile {params.outfile_prefix} \
            --ref {params.ref} \
            --min_cells {params.min_cells} \
            --min_ac_cells {params.min_ac_cells} \
            --min_ac_reads {params.min_ac_reads} \
            --max_cell_types {params.max_cell_types} \
            --min_cell_types {params.min_cell_types} \
            --fisher_cutoff {params.fisher_cutoff} \
            > {log} 2>&1
        """

rule filter_calls:
    input:
        "results/variant_calls/{sample}_{chrom}.calling.step1.tsv"
    output:
        "results/filtered_calls/{sample}_{chrom}.calling.step2.tsv"
    params:
        pon = config["inputs"]["panel_of_normals"],
        editing = config["inputs"]["rna_editing_sites"],
        outfile_prefix = "results/filtered_calls/{sample}_{chrom}",
        script = os.path.join(config["paths"]["scomatic_dir"], "scripts/BaseCellCalling/BaseCellCalling.step2.py")
    resources:
        mem_mb = 8000,
        runtime = 2040,
        partition = "spot-compute"
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/filter_calls/{sample}_{chrom}.log"
    shell:
        """
        python {params.script} \
            --infile {input} \
            --outfile {params.outfile_prefix} \
            --pon {params.pon} \
            --editing {params.editing} \
            > {log} 2>&1
        """

rule qc_filter:
    input:
        "results/filtered_calls/{sample}_{chrom}.calling.step2.tsv"
    output:
        "results/qc_filtered/{sample}_{chrom}.filtered.vcf"
    params:
        bed_file = config["inputs"]["umap_bed_file"]
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/qc_filter/{sample}_{chrom}.log"
    resources:
        mem_mb = 8000,
        runtime = 2040,
        partition = "spot-compute"
    shell:
        """
        bedtools intersect -header \
            -a {input} \
            -b {params.bed_file} | \
            awk -F "\\t" '$1 ~ /^#/ || $6 == ""' > {output} 2> {log}
        """
