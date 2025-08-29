#!/usr/bin/env Rscript

#This is an R script for very basic GWAS/association. It uses linear regression as an example: substitute in REGENIE/SAIGE/PLINK2/GWAS as needed.

library(optparse)
option_list = list(
  make_option(c("--geno"), type="character", help="Genotype matrix TSV"),
  make_option(c("--pheno"), type="character", help="Phenotype file TSV"),
  make_option(c("--covar"), type="character", help="Covariate file TSV"),
  make_option(c("--out"), type="character", help="Output associations TSV")
)
opt = parse_args(OptionParser(option_list=option_list))

genos = read.table(opt$geno, header=TRUE, sep="\t", check.names=FALSE)
pheno = read.table(opt$pheno, header=TRUE, sep="\t")
covar = read.table(opt$covar, header=TRUE, sep="\t")

samples = Reduce(intersect, list(colnames(genos)[-1], as.character(pheno[,1]), as.character(covar[,1])))

assoc_results = data.frame(Locus=genos$Locus)
for (i in 2:ncol(genos)) {
  gtype = genos[,i]
  df = data.frame(Geno=gtype, Pheno=pheno[match(colnames(genos)[i], pheno[,1]),2], Covar=covar[match(colnames(genos)[i], covar[,1]),-1])
  fit = lm(Pheno ~ Geno + ., data=df)
  p = summary(fit)$coefficients["Geno", "Pr(>|t|)"]
  assoc_results[i-1, "pvalue"] = p
}
write.table(assoc_results, file=opt$out, sep="\t", row.names=FALSE, quote=FALSE)
