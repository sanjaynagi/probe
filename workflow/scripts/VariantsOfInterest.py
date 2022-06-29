#!/usr/bin/env python
# coding: utf-8

"""
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import probetools as probe
import numpy as np
import pandas as pd
import allel
import dask.array as da
import seaborn as sns


#Selection Scans #
cloud = snakemake.params['cloud']
ag3_sample_sets = snakemake.params['ag3_sample_sets']
contigs = ['2L', '2R', '3R', '3L', 'X']
genotypePath = snakemake.params['genotypePath'] if not cloud else "placeholder_{contig}"
positionsPath = snakemake.params['positionPath'] if not cloud else "placeholder2_{contig}"
dataset = snakemake.params['dataset']

## Read VOI data
vois = pd.read_csv(snakemake.input['variants'], sep="\t")

## separate chrom and pos data and sort 
vois['contig'] = vois['Location'].str.split(":").str.get(0)
vois['pos'] = vois['Location'].str.split(":").str.get(1).str.split("-").str.get(0)
vois = vois.sort_values(['contig', 'pos'])


# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3(pre=True)
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")

# Load cohorts
cohorts = probe.getCohorts(metadata=metadata, 
                    columns=snakemake.params.columns, 
                    minPopSize=snakemake.params.minPopSize)



snps = {}
pos = {}
for contig in contigs:

    probe.log(f"Loading arrays for {contig}")
    # Load Arrays
    snps[contig], pos[contig] = probe.loadZarrArrays(genotypePath=genotypePath.format(contig = contig), 
                                             positionsPath=positionsPath.format(contig = contig),
                                             siteFilterPath=None, 
                                             cloud=cloud,
                                             contig=contig,
                                             haplotypes=False)
    

Dict = {}
allCohorts = {}

for idx, cohort in cohorts.iterrows():
    
    for i, row in vois.iterrows():
        name = row['Name']
        contig = row['contig']
        voiPos = int(row['pos'])
        longName = contig + ":"+ str(voiPos) + "  " + row['Gene'] + " | " + row['Name']

        VariantsOfInterest = pd.DataFrame([{'contig':contig, 'pos':voiPos, 'variant':name, 'name':longName}])

        bool_ = pos[contig][:] == voiPos
        
        geno = snps[contig].compress(bool_, axis=0).take(cohort['indices'], axis=1)
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
VariantsOfInterest.to_csv(f"results/variantsOfInterest/VOI.{dataset}.frequencies.tsv", sep="\t")

#Drop unnecessary columns for plotting as heatmap
VariantsOfInterest = VariantsOfInterest.drop(columns=['contig', 'pos', 'variant']).set_index('name').astype("float64").round(2)
probe.plotRectangular(VariantsOfInterest, path=f"results/variantsOfInterest/VOI.{dataset}.heatmap.png", figsize=[14,14], xlab='cohort')
