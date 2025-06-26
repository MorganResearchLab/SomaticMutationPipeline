# ==============================================================================
# Rules for Splitting and Processing BAMs
# ==============================================================================

rule split_by_celltype:
    input:
       bam = lambda wildcards: os.path.join(config["paths"]["bam_dir"], f"{wildcards.sample.split('_')[0]}.bam")
    output:
        bam = "results/split_bams/{sample}/{sample}_{celltype}_{chrom}.bam"
    params:
        outdir = "results/split_bams/{sample}",
        script = os.path.join(config["paths"]["scripts_dir"], "split_by_celltype-parallel.py"),
        meta = config["inputs"]["cell_type_annotations"],
        sra_meta = config["inputs"]["sra_metadata"],
        window = config["params"]["split_bams"]["window"]
    threads: 8
    resources:
        mem_mb = 12000,
        runtime = "12d",
        partition = "uoa-compute"
    conda:
        "../envs/scomatic.yaml"
    log:
        "logs/split_bams/{sample}_{celltype}_{chrom}.log"
    shell:
        """
        python {params.script} \
            --bamfile {input.bam} \
            --outdir {params.outdir} \
            --chrom {wildcards.chrom} \
            --celltype {wildcards.celltype} \
            --threads {threads} \
            --metafile {params.meta} \
            --chunk_size {params.window} > {log} 2>&1
        """

rule index_split_bams:
    input:
        "results/split_bams/{sample}/{sample}_{celltype}_{chrom}.bam"
    output:
        "results/split_bams/{sample}/{sample}_{celltype}_{chrom}.bam.bai"
    resources:
        mem_mb = 4000,
        runtime = 180,
        partition = "uoa-compute"
    conda:
        "../envs/scomatic.yaml"
    log:
        "logs/index_bam/{sample}_{celltype}_{chrom}.log"
    shell:
        "samtools index {input} > {log} 2>&1"