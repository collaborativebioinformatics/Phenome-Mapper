#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Check input file sample matching")
parser.add_argument("--geno", required=True, help="Genotype matrix TSV")
parser.add_argument("--pheno", required=True, help="Phenotype table TSV")
parser.add_argument("--covar", required=True, help="Covariate table TSV")
parser.add_argument("--okflag", required=True, help="OK flag output file")
args = parser.parse_args()

'''
This script checks whether the sample IDs are consistent across genotype, phenotype, and covariate files. It creates a flag file if the check passes.
'''


# Read sample IDs
geno_samples = pd.read_csv(args.geno, sep="\t", nrows=0).columns.tolist()[1:]
pheno = pd.read_csv(args.pheno, sep="\t")
covar = pd.read_csv(args.covar, sep="\t")

pheno_samples = pheno.iloc[:, 0].astype(str).tolist()
covar_samples = covar.iloc[:, 0].astype(str).tolist()

if set(geno_samples) == set(pheno_samples) == set(covar_samples):
    with open(args.okflag, "w") as f:
        f.write("OK\n")
    sys.exit(0)
else:
    sys.stderr.write("Sample IDs do not match across files!\n")
    sys.exit(1)

