# ==============================================================================
# Rule for Final Genotyping
# ==============================================================================

rule genotype_cells:
    input:
        vcf = "results/qc_filtered/{sample}_{chrom}.filtered.vcf",
        bam = "results/split_bams/{sample}/{sample}_{celltype}_{chrom}.bam"
    output:
        "results/genotyped/{sample}_{celltype}_{chrom}.single_cell_genotype.tsv"
    params:
        temp_dir = temp("results/temp/genotype_{sample}_{celltype}_{chrom}"),
        ref = lambda wildcards: f"{config['paths']['ref_dir']}/{wildcards.chrom}.fa.gz",
        meta = config["inputs"]["cell_type_annotations"],
        script = os.path.join(config["paths"]["scomatic_dir"], "scripts/SingleCellGenotype/SingleCellGenotype.py")
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 600,
        partition = "spot-compute"
    conda:
        "../envs/scomatic.yaml"
    log:
        "logs/genotype_cells/{sample}_{celltype}_{chrom}.log"
    shell:
        """
        python {params.script} \
            --bam {input.bam} \
            --infile {input.vcf} \
            --nprocs {threads} \
            --meta {params.meta} \
            --outfile {output} \
            --tmp_dir {params.temp_dir} \
            --ref {params.ref} \
            > {log} 2>&1
        """