# Snakefile for PhenomeMapper Graph-Based Association Pipeline

SAMPLES = ["sample1", "sample2"]  # extend with more sample names as needed

REF_FASTA = "data/reference/genome.fa"
SAMPLE_ASM = expand("data/assemblies/{sample}.fa", sample=SAMPLES)
READS = expand("data/reads/{sample}.fastq.gz", sample=SAMPLES)
VCFS = expand("results/genotyping/{sample}.vcf", sample=SAMPLES)
PHENO = "data/phenotypes/phenotypes.tsv"
COVAR = "data/phenotypes/covariates.tsv"

rule all:
    input:
        "results/annotation/annotated_loci.tsv"

##########################################################################################

rule build_pangenome_graph:
    input:
        ref=REF_FASTA,
        assemblies=SAMPLE_ASM
    output:
        gfa="results/graphs/pangenome.gfa"
    shell:
        """
        minigraph -xggs {input.ref} {input.assemblies} > {output.gfa}
        """

rule gfa_to_xg:
    input:
        gfa="results/graphs/pangenome.gfa"
    output:
        xg="results/graphs/pangenome.xg"
    shell:
        """
        vg convert -x {input.gfa} > {output.xg}
        """

rule map_reads_graph:
    input:
        xg="results/graphs/pangenome.xg",
        reads="data/reads/{sample}.fastq.gz"
    output:
        gam="results/mapping/{sample}.gam"
    shell:
        """
        vg giraffe -x {input.xg} \
            -H <(vg gbwt --xg {input.xg}) \
            -m <(vg minimizer --xg {input.xg}) \
            -d <(vg distance -x {input.xg}) \
            -f {input.reads} > {output.gam}
        """

rule pack_gam:
    input:
        xg="results/graphs/pangenome.xg",
        gam="results/mapping/{sample}.gam"
    output:
        pack="results/genotyping/{sample}.pack"
    shell:
        """
        vg pack -x {input.xg} -g {input.gam} -o {output.pack}
        """

rule call_variants:
    input:
        xg="results/graphs/pangenome.xg",
        pack="results/genotyping/{sample}.pack"
    output:
        vcf="results/genotyping/{sample}.vcf"
    shell:
        """
        vg call {input.xg} -k {input.pack} > {output.vcf}
        """

rule merge_vcfs:
    input:
        vcfs=VCFS
    output:
        merged_vcf="results/variants/merged.vcf.gz"
    shell:
        """
        bcftools merge {input.vcfs} -Oz -o {output.merged_vcf}
        tabix -p vcf {output.merged_vcf}
        """

rule vcf_to_matrix:
    input:
        merged_vcf="results/variants/merged.vcf.gz"
    output:
        "results/matrix/genotype_matrix.tsv"
    shell:
        """
        bcftools query -f '[%SAMPLE\\t%GT\\n]' {input.merged_vcf} > {output}
        """

rule check_phenotypes:
    input:
        genotypes="results/matrix/genotype_matrix.tsv",
        phenotypes=PHENO,
        covariates=COVAR
    output:
        ready="results/phenotypes/inputs_ready.txt"
    shell:
        """
        python scripts/check_inputs.py --geno {input.genotypes} --pheno {input.phenotypes} --covar {input.covariates} --okflag {output.ready}
        """

rule run_association:
    input:
        genotypes="results/matrix/genotype_matrix.tsv",
        phenotypes=PHENO,
        covariates=COVAR,
        ready="results/phenotypes/inputs_ready.txt"
    output:
        results="results/association/assoc_results.tsv"
    shell:
        """
        Rscript scripts/run_gwas.R --geno {input.genotypes} --pheno {input.phenotypes} --covar {input.covariates} --out {output.results}
        """

rule annotate_loci:
    input:
        results="results/association/assoc_results.tsv"
    output:
        "results/annotation/annotated_loci.tsv"
    shell:
        """
        python scripts/annotate_results.py --assoc {input.results} --gene-anno data/annotation/genes.bed --out {output}
        """

##########################################################################################

# Required Scripts:
# - scripts/check_inputs.py: ensures matrix and phenotype files match and reports problems.
# - scripts/run_gwas.R: runs GWAS/association tests for loci vs. phenotypes (can use plink2, regenie, or your choice).
# - scripts/annotate_results.py: annotates association result loci with genomic regions or gene overlaps from BED/GFF.

# Required Data:
# - data/reference/genome.fa: reference fasta
# - data/assemblies/*.fa: sample genome assemblies
# - data/reads/*.fastq.gz: short or long reads for mapping
# - data/phenotypes/phenotypes.tsv and covariates.tsv
# - data/annotation/genes.bed: gene model annotation


