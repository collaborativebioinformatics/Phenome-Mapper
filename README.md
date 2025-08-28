# Phenome-Mapper

                        To find association between multiple graph variation loci and (sub)phenotypes
Detailed step-by-step workflow to find associations between multiple graph variation loci (such as SV/VNTR) and (sub)phenotypes. 
1)Build a sequence graph: Use pangenome graph construction tools (e.g., minigraph, vg, APAV toolkit) to integrate all samples, capturing all structural variations and VNTRs.

2)Align reads or assemblies: Map sequence data back to the constructed graph using graph-aware aligners (GraphAligner, vg giraffe, Paragraph).
3)Create sample-by-locus matrix: For each sample, list allelic state or genotype at each graph variation locus.
4)Annotate associated loci: For loci showing significant association, examine overlap with genes, regulatory regions, and known disease markers
    <img width="1880" height="1063" alt="image" src="https://github.com/user-attachments/assets/4cb6269f-efbc-4a5c-b744-ec03993f7223" />
