## Phenome-Mapper

                        To find association between multiple graph variation loci and (sub)phenotypes
  Detailed step-by-step workflow to find associations between multiple graph variation loci (such as SV/VNTR) and (sub)phenotypes:

<div align="center">
<img width="600" alt="Phenome-Mapper workflow" src="https://github.com/user-attachments/assets/4cb6269f-efbc-4a5c-b744-ec03993f7223" />
</div>

1)**Build a sequence graph:** Use pangenome graph construction tools (e.g., minigraph, vg, APAV toolkit) to integrate all samples, capturing all structural variations and VNTRs.
**Step 1: Build a Sequence Graph**
    Goal: Integrate all known and sample-specific genome variations into a comprehensive graph.
    Tools: minigraph, vg, pggb, or similar.
    Input: Reference genome (FASTA), sample assemblies/reads, population VCFs.
    Output: Sequence graph (usually GFA/VG format) capturing SVs, VNTRs, and allelic paths.
    Example: asuming you have created your work space, installed minigraph, spades for assemblies and dowonloaded reference fasta files.
```
    bash
        spades.py -1 sample_R1.fastq -2 sample_R2.fastq -o spades_output
            Output: Assembled contigs in spades_output/contigs.fasta
        quast.py spades_output/contigs.fasta -o quast_report/
            # Optionally add: -r reference.fa for reference-based metricsquast flye_assembly/assembly.fasta -o quast_report
            # use your fasta contigs as input to minigraph
        minigraph -x asm -t8 -l -o graph.gfa reference.fa spades_output/contigs.fasta
            # can add more assemblies fasta files minigraph -x asm -t8 -l -o graph.gfa reference.fa sample1.fa sample2.fa ...
            # Or: vg construct -r ref.fa -v cohort.vcf.gz > graph.vg
```

Index the Sequence Graph for Mapping
minigraph outputs a GFA that may require conversion for some mapping/analysis pipelines.
If using vg, create XG indexes for rapid mapping and variant calling.

```
    bash
            # Convert GFA to vg format (if needed)
        vg convert -g graph.gfa > graph.vg
            # Index for Giraffe (vg)
        vg index -x graph.xg -g graph.gg graph.vg


```

2)**Align reads or assemblies:** Map sequence data back to the constructed graph using graph-aware aligners (GraphAligner, vg giraffe, Paragraph).
```
   bash
        vg giraffe -x graph.xg -g graph.gg -m graph.min -d graph.dist \
        -f sample_R1.fastq -f sample_R2.fastq > sample.gam
```
   Call Variants and Genotypes at Graph Loci
   Generate genotype calls or repeat counts for every sample at each locus in the graph.
   Tools:(vg pack and vg call)
   Example:
 ```
    bash
        vg pack -x graph.xg -g sample.gam -o sample.pack
        vg call graph.xg -k sample.pack > sample.vcf
```

3)**Create sample-by-locus matrix:** For each sample, list allelic state or genotype at each graph variation locus.
*Build Sample-by-Locus Genotype Matrix*
  Convert your VCF(s) for all samples into a matrix or table, with rows as samples and columns as SV/VNTR loci (with genotype state, repeat lengths, etc.).
    Tools:(bcftools, vcftools, pandas, or R.)
    Prepare/Harmonize Phenotype and Covariate Data
    Ensure your phenotype table matches sample IDs.
    Include relevant covariates (age, sex, batch effects, ancestry PCs).
  
    ```
      bash
        bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' input.vcf.gz > genotype_matrix.tsv  #single sample/multisample)
     ```
    
    #Explanation:
                  %CHROM, %POS, %ID: Outputs chromosome, position, and variant ID as locus identifiers.
                  [\t%GT]: For each sample, outputs the genotype (0/0, 0/1, etc.), separated by tabs.
                  Each line is a locus (SV/VNTR/SNP), and columns after the first three are per sample.
                  The output is a tab-delimited file ideal for pandas, R, or spreadsheet analysis.

    To include SV or VNTR-specific information (e.g., repeat lengths, SV type):
    
    ```
      bash
        bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVLEN\t%INFO/REPTYPE[\t%GT]\n' input.vcf.gz > genotype_matrix_with_info.tsv
    ```    

4)**Annotate associated loci:**For loci showing significant association, examine overlap with genes, regulatory regions, and known disease markers


**Future direction **
Perform Association Analysis
Test for statistical association (linear regression/ linear mixed model) between each graph locus (SV/VNTR genotype) and phenotypes/subphenotypes.
Tools:
plink2, REGENIE, SAIGE, rvtests, R, or Python analysis packages.

