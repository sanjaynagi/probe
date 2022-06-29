#!/usr/bin/env python
# coding: utf-8

"""
pca 
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import probetools as probe
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import dask.array as da
import seaborn as sns
import matplotlib.pyplot as plt


# Garuds Selection Scans # 
cloud = snakemake.params['cloud']
ag3_sample_sets = snakemake.params['ag3_sample_sets']
pcaColumn = snakemake.params['cohortColumn']
contig = snakemake.wildcards['contig']
dataset = snakemake.params['dataset']
genotypePath = snakemake.input['genotypes'] if not cloud else []
positionsPath = snakemake.input['positions'] if not cloud else []
siteFilterPath = snakemake.input['siteFilters'] if not cloud else []

results_dir = snakemake.params['data']

# Load metadata 
if cloud:
    import malariagen_data
    ag3 = malariagen_data.Ag3(pre=True)
    metadata = ag3.sample_metadata(sample_sets=ag3_sample_sets)
else:
    metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")

# Load Arrays
snps, pos = probe.loadZarrArrays(genotypePath, positionsPath, siteFilterPath=siteFilterPath, cloud=cloud, haplotypes=False, contig=contig)

# Determine cohorts
cohorts = probe.getCohorts(metadata, columns=pcaColumn)


# choose colours for species
species_palette = px.colors.qualitative.Plotly
species_color_map = {
    'gambiae': species_palette[0],
    'coluzzii': species_palette[1],
    'arabiensis': species_palette[2],
    'intermediate_gambiae_coluzzii': species_palette[3],
    'intermediate_arabiensis_gambiae': species_palette[4],
}


# Run PCA on whole dataset together
data, evr = probe.run_pca(contig=contig, gt=snps, pos=pos, df_samples=metadata,
    sample_sets=dataset, results_dir=results_dir
)
evr = evr.astype("float").round(4) # round decimals for variance explained % 

probe.plot_coords(data, evr, title=f" PCA | {dataset} | {contig}", filename=f"results/PCA/{dataset}.{contig}.html")

fig = plt.figure(figsize=(10, 10))
fig = sns.scatterplot('PC1','PC2', data=data, hue="species_gambiae_coluzzii")
fig.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title(f"PCA | {dataset} | {contig}", fontsize=14)
plt.xlabel(f"PC1 ({evr[0]*100} % variance explained)", fontdict={"fontsize":14})
plt.ylabel(f"PC2 ({evr[1]*100} % variance explained)", fontdict={"fontsize":14})
plt.savefig(f"results/PCA/{dataset}.{contig}.png")



# Loop through each cohort, manipulate genotype arrays and calculate chosen Garuds Statistic
for idx, cohort in cohorts.iterrows():

    # filter to correct loc, year, species individuals
    gt_cohort = snps.take(cohort['indices'], axis=1)
    meta = metadata.take(cohort['indices'])
    
    data, evr = probe.run_pca(contig=contig, gt=gt_cohort, pos=pos, df_samples=meta,
        sample_sets=cohort['cohortNoSpaceText'], results_dir=results_dir
    )
    evr = evr.astype("float").round(4)

    probe.plot_coords(data, evr, title=f" PCA | {cohort['cohortText']} | {contig}", filename=f"results/PCA/{cohort['cohortNoSpaceText']}.{contig}.html")

    fig = plt.figure(figsize=(10, 10))
    fig = sns.scatterplot('PC1','PC2', data=data, hue='location')
    fig.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(f"PCA | {cohort['cohortText']} | {contig}", fontsize=14)
    plt.xlabel(f"PC1 ({evr[0]*100} % variance explained)", fontdict={"fontsize":14})
    plt.ylabel(f"PC2 ({evr[1]*100} % variance explained)", fontdict={"fontsize":14})
    plt.savefig(f"results/PCA/{cohort['cohortNoSpaceText']}.{contig}.png")

