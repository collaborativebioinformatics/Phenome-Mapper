#!/usr/bin/env python3

'''
This script annotates loci with overlapping genes from a BED file.
'''



import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Annotate association results with gene overlaps")
parser.add_argument("--assoc", required=True, help="Association results TSV")
parser.add_argument("--gene-anno", required=True, help="Gene annotation BED")
parser.add_argument("--out", required=True, help="Annotated output TSV")
args = parser.parse_args()

assoc = pd.read_csv(args.assoc, sep="\t")
genes = pd.read_csv(args.gene_anno, sep="\t", header=None, names=["chrom", "start", "end", "gene"])

def annotate_locus(row):
    # Example locus: "chr1:12345"
    chrom, pos = row["Locus"].split(":")
    pos = int(pos)
    overlapping = genes[(genes["chrom"] == chrom) & (genes["start"] <= pos) & (genes["end"] >= pos)]
    return ";".join(overlapping["gene"].tolist()) if not overlapping.empty else "NA"

assoc["Genes"] = assoc.apply(annotate_locus, axis=1)
assoc.to_csv(args.out, sep="\t", index=False)
