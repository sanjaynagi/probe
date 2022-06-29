#!/usr/bin/env python
# coding: utf-8

"""

"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import probetools as probe
from pathlib import Path
import numpy as np
import pandas as pd
import allel
import dask.array as da
import seaborn as sns
import matplotlib.pyplot as plt
from numba import njit

cloud = snakemake.params['cloud']
ag3_sample_sets = snakemake.params['ag3_sample_sets']
contig = snakemake.wildcards['contig']
genotypePath = snakemake.params['genotypes'] 
positionsPath = snakemake.params['positions']


# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3(pre=True)
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")


@njit()
def scanRight(geno1, geno2, upperBreakpoint, max): 
    # Reset position for next loop
    gn1 = geno1[upperBreakpoint]
    gn2 = geno2[upperBreakpoint]

    # Scan right along genome, as long as two inds are not both homozygous but different
    while not (gn1[0] == gn1[1]) & (gn2[0] == gn2[1]) & ((gn1 != gn2).all()) & (-1 not in gn1 and -1 not in gn2):

        upperBreakpoint += 1
        if upperBreakpoint == max-1: # limit the upper breakpoint at end of the contig
            return(upperBreakpoint)

        gn1 = geno1[upperBreakpoint]
        gn2 = geno2[upperBreakpoint]
        
    return(upperBreakpoint)

@njit()
def scanLeft(geno1, geno2, lowerBreakpoint): 

    # subset genotypes to this position, because we need somethign to start the while loop?
    gn1 = geno1[lowerBreakpoint]
    gn2 = geno2[lowerBreakpoint]

    # Scan left along genome
    while not (gn1[0] == gn1[1]) & (gn2[0] == gn2[1]) & ((gn1 != gn2).all()) & (-1 not in gn1 and -1 not in gn2):

        lowerBreakpoint -= 1
        if lowerBreakpoint == 0: # limit lower breakpoint at zero, start of contig
            return(lowerBreakpoint)

        gn1 = geno1[lowerBreakpoint]
        gn2 = geno2[lowerBreakpoint]
        
    return(lowerBreakpoint)


@njit()
def f2scans(dblton_arr, snps, pos):

    starts = []# np.empty((len(dblton_arr)),dtype='uint8')
    ends = []# np.empty((len(dblton_arr)),dtype='uint8')
    dbltonpos = []

    for idx in range(0, len(dblton_arr)):

        geno1 = snps[:, dblton_arr[idx][0]]
        geno2 = snps[:, dblton_arr[idx][1]]
        # get boolean of dblton idx
        #dblton_idx = bisect.bisect_left(pos, dblton_arr[idx][2])
        dblton_idx = np.searchsorted(pos, dblton_arr[idx][2])

        upperBreakpoint = scanRight(geno1, geno2, dblton_idx, len(pos))
        # Scan left along genome
        lowerBreakpoint = scanLeft(geno1, geno2, dblton_idx)

        starts.append(pos[lowerBreakpoint])
        ends.append(pos[upperBreakpoint])
        dbltonpos.append(dblton_arr[idx][2])
    
    return(np.array(starts), np.array(ends), np.array(dbltonpos))



########## main #################


snps = {}
pos = {}

# Load Arrays
snps, pos = probe.loadZarrArrays(genotypePath=genotypePath, 
                                            positionsPath=positionsPath,
                                            siteFilterPath=None,
                                            cloud=cloud)
ac = snps.count_alleles()
seg = ac.is_segregating()
snps = snps.compress(seg, axis=0).compute(numworkers=12)
pos = pos[seg]

### Load doubletons
dblton = pd.read_csv(snakemake.input['f2variantPairs'], sep="\t")
dblton = dblton.query("contig == @contig")
dblton_arr = dblton[['idx1', 'idx2', 'pos']].to_numpy()

# extract np array 
snps = snps.values
# Run F2 hap length scans
starts, ends, dbltonpos = f2scans(dblton_arr, snps, pos)


