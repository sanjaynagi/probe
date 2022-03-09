#!/usr/bin/env python
# coding: utf-8

"""

"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import loadZarrArrays, log
from pathlib import Path
import numpy as np
import pandas as pd
import allel
import dask.array as da
import seaborn as sns
import matplotlib.pyplot as plt


cloud = snakemake.params['cloud']
ag3_sample_sets = snakemake.params['ag3_sample_sets']
contig = snakemake.wildcards['contig']
genotypePath = snakemake.params['genotypes'] 
positionsPath = snakemake.params['positions']


# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3()
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")


snps = {}
pos = {}

# Load Arrays
snps[contig], pos[contig] = loadZarrArrays(genotypePath=genotypePath, 
                                            positionsPath=positionsPath,
                                            siteFilterPath=None,
                                            cloud=cloud)
ac = snps[contig].count_alleles()
seg = ac.is_segregating()
snps[contig] = snps[contig].compress(seg, axis=0).compute(numworkers=12)
pos[contig] = pos[contig][seg]

### Load doubletons
dblton = pd.read_csv(snakemake.input['f2variantPairs'], sep="\t")
dblton = dblton.query("contig == @contig")

with open(f"results/f2HapLengths_{contig}.tsv", mode='w') as file:
    file.close()

for idx, row in dblton.iterrows():

    # take one doubleton at a time, get contig, pos, individual1 and 2
    contig, dbltonpos, idx1, idx2 = row[['contig', 'pos', 'idx1', 'idx2']]
    ## subset to individuals, we need whole chrom atm

    # subset genotypes to these individuals
    geno = snps[contig].take([idx1,idx2], axis=1)
    # get boolean of dblton idx
    dblton_idx = np.where(pos[contig] == dbltonpos)[0]
    # subset genotypes to this position, because we need somethign to start the while loop?
    gn1 = geno[dblton_idx, 0]
    gn2 = geno[dblton_idx, 1]

    # we will add to and subtract from the dblton index 
    upperBreakpoint = dblton_idx.copy()
    lowerBreakpoint = dblton_idx.copy()

    # Scan right along genome, as long as two inds are not both homozygous but different
    while not (gn1.is_hom() & gn2.is_hom() & (gn1 != gn2)).all():

        upperBreakpoint += 1
        if upperBreakpoint == pos[contig].shape[0]: # limit the upper breakpoint at end of the contig
            break

        gn1 = geno[upperBreakpoint, 0]
        gn2 = geno[upperBreakpoint, 1]

    # Reset position for next loop
    gn1 = geno[dblton_idx, 0]
    gn2 = geno[dblton_idx, 1]

    # Scan left along genome
    while not (gn1.is_hom() & gn2.is_hom() & (gn1 != gn2)).all():

        lowerBreakpoint -= 1
        if lowerBreakpoint == 0: # limit lower breakpoint at zero, start of contig
            break

        gn1 = geno[lowerBreakpoint, 0]
        gn2 = geno[lowerBreakpoint, 1]

    start = pos[contig][lowerBreakpoint]
    end = pos[contig][upperBreakpoint]
    size = end-start

    with open(f"results/f2HapLengths_{contig}.tsv", mode='a') as file:
        file.write(f"{contig}\t{str(dbltonpos)}\t{str(size[0])}\t{str(start[0])}\t{str(end[0])}\n")
