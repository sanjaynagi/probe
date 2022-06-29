#!/usr/bin/env python
# coding: utf-8

"""
Karyotype Anopheles gambiae s.l
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import probetools as probe
import numpy as np
import pandas as pd
import allel


cloud = snakemake.params['cloud']
dataset = snakemake.params['dataset']
ag3_sample_sets = snakemake.params['ag3_sample_sets']

# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3(pre=True)
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)

    genotypePath = []
    haplotypePath = []
    positionsPath = []
    siteFilterPath = []
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")
    
    genotypePath = snakemake.input['genotypes']
    positionsPath = snakemake.input['positions']
    siteFilterPath = snakemake.input['siteFilters']


invDict = {}

for inversion in probe.inversionDict.keys():
    chrom = probe.inversionDict[inversion][0]

    if cloud:
        snps = allel.GenotypeDaskArray(ag3.snp_genotypes(region=chrom, sample_sets=ag3_sample_sets))
        pos = allel.SortedIndex(ag3.snp_sites(region=chrom, field='POS'))
    else:
        snps, pos = probe.loadZarrArrays(genotypePath, positionsPath, siteFilterPath=siteFilterPath, haplotypes=False, cloud=cloud, contig=contig)

    callset = {'geno':snps, 'pos':pos, 'chrom':chrom}
    print(f"--- Running CompKaryo {inversion}--- ")
    av_gts, total_sites, num_0, num_1, num_2 = probe.compkaryo(callset, inversion)
    total_sites = total_sites.compute()
    invDict[inversion] = pd.DataFrame({'partner_sample_id': metadata['sample_id'], 
                           'inversion':inversion, 
                           'mean_genotype': av_gts, 
                           'total_snp_tags':total_sites})
    
karyotypes = pd.concat(invDict)
karyotypes.to_csv(f"results/{dataset}_karyotypes.tsv", sep="\t")