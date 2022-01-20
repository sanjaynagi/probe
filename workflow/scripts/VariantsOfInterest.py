#!/usr/bin/env python
# coding: utf-8

"""
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import loadZarrArrays, getCohorts, plotRectangular, log
import numpy as np
import pandas as pd
import allel
import dask.array as da
import seaborn as sns


#Selection Scans # 
contigs = ['2L', '2R', '3R', '3L', 'X']
genotypePath = snakemake.params['genotypePath']
positionsPath = snakemake.params['positionPath']


## Read VOI data
vois = pd.read_csv(snakemake.input['variants'], sep="\t")

## separate chrom and pos data and sort 
vois['chrom'] = vois['Location'].str.split(":").str.get(0)
vois['pos'] = vois['Location'].str.split(":").str.get(1).str.split("-").str.get(0)
vois = vois.sort_values(['chrom', 'pos'])


# Read metadata 
metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")

# Load cohorts
cohorts = getCohorts(metadata=metadata, 
                    columns=snakemake.params.columns, 
                    minPopSize=snakemake.params.minPopSize)



snps = {}
pos = {}
for chrom in contigs:
    # Load Arrays
    snps[chrom], pos[chrom] = loadZarrArrays(genotypePath=genotypePath.format(chrom = chrom), 
                                             positionsPath=positionsPath.format(chrom = chrom),
                                             siteFilterPath=None, 
                                             haplotypes=False)
    

Dict = {}
allCohorts = {}

for idx, cohort in cohorts.iterrows():
    
    for i, row in vois.iterrows():
        name = row['Name']
        chrom = row['chrom']
        voiPos = int(row['pos'])
        longName = chrom + ":"+ str(voiPos) + "  " + row['Gene'] + " | " + row['Name']

        VariantsOfInterest = pd.DataFrame([{'chrom':chrom, 'pos':voiPos, 'variant':name, 'name':longName}])

        bool_ = pos[chrom][:] == voiPos
        
        geno = snps[chrom].compress(bool_, axis=0).take(cohort['indices'], axis=1)
        ac = geno.count_alleles().compute()
        # if there are no ALTs lets add the zeros for the ALTs otherwise only REF count returned 
        aclength = ac.shape[1]
        acneeded = 4-aclength
        ac = np.append(ac, np.repeat(0, acneeded))
        #get frequency and round
        freqs = pd.DataFrame([ac/ac.sum().round(2)])
        df2 = freqs.apply(pd.Series).rename(columns={0:'REF', 1:'ALT1', 2:'ALT2', 3:'ALT3'})
        VariantsOfInterest[cohort['cohortText']] = df2.drop(columns=['REF']).sum(axis=1).round(2)
        
        Dict[name] = VariantsOfInterest
    
    allCohorts[idx] = pd.concat(Dict)

# Concatenated the table and write table to TSV
VariantsOfInterest = pd.concat(allCohorts, axis=1).T.drop_duplicates().T.droplevel(level=0, axis=1)
VariantsOfInterest.to_csv("results/variantsOfInterest/VOI.frequencies.tsv", sep="\t")

#Drop unnecessary columns for plotting as heatmap
VariantsOfInterest = VariantsOfInterest.drop(columns=['chrom', 'pos', 'variant']).set_index('name').astype("float64").round(2)
plotRectangular(VariantsOfInterest, path="results/variantsOfInterest/VOI.heatmap.png", figsize=[14,14], xlab='cohort')