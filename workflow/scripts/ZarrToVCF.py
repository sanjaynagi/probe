#!/usr/bin/env python
# coding: utf-8

"""
Zarr to VCF. Needs providing with REF and ALT and allpositions files.
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import loadZarrArrays, getCohorts, log
from pathlib import Path
import numpy as np
import zarr
import pandas as pd
import allel
import dask.array as da


# Garuds Selection Scans # 
chrom = snakemake.wildcards['chrom']
dataset = snakemake.params['dataset']
genotypePath = snakemake.input['genotypes'] if stat in ['G12', 'G123'] else []
positionsPath = snakemake.input['positions']
siteFilterPath = snakemake.input['siteFilters']
#refPath =
#altPath = 
#allPosPath = 

results_dir = snakemake.params['data']

# Read metadata 
metadata = pd.read_csv(snakemake.params['metadata'], sep=",")
metadata['location'] = metadata['location'].str.split(".").str.get(0)

def write_vcf_header(vcf_file, chrom):
    """
    Writes a VCF header.
    """
    
    print('##fileformat=VCFv4.1', file=vcf_file)
    # write today's date
    today = date.today().strftime('%Y%m%d')
    print('##fileDate=%s' % today, file=vcf_file)
    # write source
    print('##source=scikit-allel-%s + ZarrToVCF.py' % allel.__version__, file=vcf_file)
    #write refs and contigs 
    print('##reference=resources/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa', file=vcf_file)
    print('##contig=<ID=2R,length=61545105>', file=vcf_file) if chrom == '2R' else None
    print('##contig=<ID=3R,length=53200684>', file=vcf_file) if chrom == '3R' else None 
    print('##contig=<ID=2L,length=49364325>', file=vcf_file) if chrom == '2L' else None
    print('##contig=<ID=3L,length=41963435>', file=vcf_file) if chrom == '3L' else None
    print('##contig=<ID=X,length=24393108>', file=vcf_file) if chrom == 'X' else None
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=vcf_file)


def ZarrToPandasToVCF(vcf_file, genotypePath, positionsPath, siteFilterPath, chrom, nchunks=50, snpfilter = "segregating"):
    
    """
    Converts genotype and POS arrays to vcf, using pd dataframes in chunks. 
    Segregating sites only. Needs REF and ALT arrays.
    """
    
    #if file exists ignore and skip
    myfile = Path(f"{vcf_file}.gz")
    if myfile.is_file():
        print(f"File {vcf_file}.gz Exists...")
        return
    
    log(f"Loading array for {chrom}...")

    geno, pos = loadZarrArrays(genotypePath, positionsPath, siteFilterPath=siteFilterPath, haplotypes=False)
    allpos = allel.SortedIndex(zarr.open_array(f"resources/snp_genotypes/all/sites/{chrom}/variants/POS/")[:])
    
    ref_alt_filter = allpos.locate_intersection(pos)[0]
    
    refs = zarr.open_array(f"../resources/snp_genotypes/all/sites/{chrom}/variants/REF/")[:][ref_alt_filter]
    alts = zarr.open_array(f"../resources/snp_genotypes/all/sites/{chrom}/variants/ALT/")[:][ref_alt_filter]
    
    if snpfilter == "segregating":
        log("Find segregating sites...")
        flt = geno.count_alleles().is_segregating()
        geno = geno.compress(flt, axis=0)
        positions = pos[flt]
        refs = refs[flt].astype(str)
        alts = [a +"," + b + "," + c for a,b,c in alts[flt].astype(str)]
    elif snpfilter == 'biallelic01':
        log("Finding biallelic01 sites...")
        flt = geno.count_alleles().is_biallelic_01()
        geno = geno.compress(flt, axis=0)
        positions = pos[flt]
        refs = refs[flt].astype(str)
        alts = [a for a,b,c in alts[flt].astype(str)]
    else:
        assert np.isin(snpfilter, ['segregating', "biallelic01"]).any(), "incorrect snpfilter value"
  
    log("calculating chunks sizes...")
    chunks = np.round(np.arange(0, geno.shape[0], geno.shape[0]/nchunks)).astype(int).tolist()
    chunks.append(geno.shape[0])

    for idx, chunk in enumerate(chunks[:-1]):

        gn = geno[chunks[idx]:chunks[idx+1]].compute()
        pos = positions[chunks[idx]:chunks[idx+1]]
        ref = refs[chunks[idx]:chunks[idx+1]]
        alt = alts[chunks[idx]:chunks[idx+1]]
        
        # Contruct SNP info DF
        vcf_df = pd.DataFrame({'#CHROM': chrom,
                 'POS': pos,
                 'ID': '.',
                 'REF': ref,
                 'ALT': alt,
                 'QUAL': '.',
                 'FILTER': '.',
                 'INFO':'.',
                'FORMAT': 'GT'})

        log(f"Pandas SNP info DataFrame constructed...{idx}")

        # Geno to VCF
        vcf = pd.DataFrame(gn.to_gt().astype(str), columns=metadata['partner_sample_id'])
        log("Concatenating info and genotype dataframes...")
        vcf = pd.concat([vcf_df, vcf], axis=1)

        log(f"Pandas Genotype data constructed...{idx}")

        if (idx==0) is True:
            with open(f"{vcf_file}", 'w') as vcfheader:
                    write_vcf_header(vcfheader, chrom)

        log("Writing to .vcf")

        vcf.to_csv(vcf_file, 
                   sep="\t", 
                   index=False,
                   mode='a',
                  header=(idx==0), 
                  line_terminator="\n")



### MAIN ####

chroms = ['2L', '2R', '3L', '3R', 'X']

ZarrToPandasToVCF(f"../resources/vcfs/ag3_gaardian_{chrom}.multiallelic.vcf", genotypePath, positionsPath, siteFilterPath, chrom, snpfilter="segregating")


ZarrToPandasToVCF(f"../resources/vcfs/ag3_gaardian_{chrom}.biallelic.vcf", genotypePath, positionsPath, siteFilterPath, chrom, snpfilter="biallelic01")