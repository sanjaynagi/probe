#!/usr/bin/env python
# coding: utf-8

"""

"""

from curses import meta
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
genotypePath = snakemake.params['genotypes'] 
positionsPath = snakemake.params['positions']
siteFiltersPath = snakemake.params['sitefilters']
contigs = snakemake.params['']


# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3(pre=True)
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")


def is_variant01(gn, allele):
    return((gn == allele).any())

def checkSampleID(x, metadata=metadata):
    name = metadata.loc[x,'partner_sample_id']
    return(name)


snps = {}
pos = {}

for contig in contigs:
    # Load Arrays
    snps[contig], pos[contig] = loadZarrArrays(genotypePath=genotypePath, 
                                             positionsPath=positionsPath,
                                             siteFilterPath=siteFiltersPath)
                                        

    ac = snps[contig].count_alleles()
    seg = ac.is_segregating()
    snps[contig] = snps[contig].compress(seg, axis=0)
    pos[contig] = pos[contig][seg]

pos_dbltons = {}
inds_dbltons = {}

for contig in contigs:

    ac = snps[contig].count_alleles()
    doubletons_bool = ac.is_doubleton(allele=1).compute() # Need to do for each allele
    geno = snps[contig].compress(doubletons_bool, axis=0)
    log("Recorded dblton loc")
    pos_dbltons[contig] = pos[contig][doubletons_bool]
    
    n_doubletons = doubletons_bool.sum()
    log(f"There are {n_doubletons.shape} on {contig}")
    log("locating dblton sharers")

    # get 1 hets and 1 homs to each ind of the genotype data
    res = geno.is_het(1).compute()
    res2 = geno.is_hom(1).compute()
    dbhets = res.sum()
    dbhoms = res2.sum()
    totdbinds = dbhets + (2*dbhoms)
    assert (n_doubletons*2) == totdbinds, "Unequal individual samples v n_dbltons!!!"
    
    res = np.logical_or(res, res2)
    # and use np where to get indices. Pandas apply is fast.
    pairs = pd.DataFrame(res).apply(np.where, axis=1).apply(np.asarray).apply(lambda x: x.flatten())
    hom_filter = pairs.apply(len) == 2
    pos_dbltons[contig] = pos_dbltons[contig][hom_filter]
    pairs = pairs[hom_filter]
    
    #make 1d array into two column pd df
    log("organising arrays")
    idxs = pd.DataFrame(np.vstack(pairs.values))
    dblton = pd.DataFrame(np.vstack(pairs.values), columns=['partner_sample_id','partner_sample_id2'])
    # shouldnt be any but remove hom/homs
    dblton = dblton.query("partner_sample_id != partner_sample_id2").reset_index(drop=True)
    dblton = dblton.applymap(checkSampleID) #store WA-XXXX ID
    inds_dbltons[contig] = pd.concat([idxs, dblton], axis=1)
    inds_dbltons[contig]['pos'] = pos_dbltons[contig]

dblton = pd.concat(inds_dbltons).reset_index().drop(columns=['level_1']).rename(columns={'level_0':'contig', 0:'idx1', 1:'idx2'})

dblton = dblton.merge(metadata[['partner_sample_id', 'latitude', 'longitude']])
dblton = dblton.merge(metadata.rename(columns={'partner_sample_id':'partner_sample_id2', 
                                      'latitude': 'latitude2', 
                                      'longitude':  'longitude2'})[['partner_sample_id2', 'latitude2', 'longitude2']])
dblton = dblton.sort_values(by=['contig', 'pos']).reset_index(drop=True)

dblton.to_csv("results/f2variantPairs.tsv", sep="\t", index=None)