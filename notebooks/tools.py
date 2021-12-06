import pandas as pd
import allel
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import json
import seaborn as sns
import zarr
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import os
import bisect
import hashlib
# quieten dask warnings about large chunks
import plotly
import plotly.express as px
import bokeh.plotting
import bokeh.models
import bokeh.layouts
import bokeh.io
import bokeh.palettes
import itertools
import geopy.distance

inversionDict = {"2La" : ("2L", "karyotype/2La_targets.txt"),
                 "2Rj" : ("2R", "karyotype/2Rj_targets.txt"),
                 "2Rb" : ("2R", "karyotype/2Rb_targets.txt"),
                 "2Rc_col" : ("2R", "karyotype/2Rc_col_targets.txt"),
                 "2Rc_gam" : ("2R", "karyotype/2Rc_gam_targets.txt"),
                 "2Rd" : ("2R", "karyotype/2Rd_targets.txt"),
                 "2Ru" : ("2R", "karyotype/2Ru_targets.txt")}


def import_inversion(inversion):

    '''Load the tag SNPs for the appropriate inversion.'''

    path = inversionDict[inversion][1]
    targets = pd.read_csv(path, header=None).rename(columns={0:'pos'})    
    return targets

def calculate_genotype_at_concordant_sites(callset, targets):

    '''Calculate the average number of alternate alleles for each specimen at
    each tag SNP.'''

    bool_ =  np.isin(callset['pos'], targets['pos'])
    genos = callset["geno"].compress(bool_, axis=0)
    alt_count = genos.to_n_alt()
    is_called = genos.is_called()

    av_gts = np.mean(np.ma.MaskedArray(alt_count,
                                       mask=~is_called), axis=0).data
    match_dict = {0: None, 1: None, 2: None}

    for value in [0, 1, 2]:
        n_matches = np.sum(np.ma.MaskedArray(alt_count,
                                             mask=~is_called) == value,
                            axis=0).data
        match_dict[value] = n_matches

    total_sites = np.sum(is_called, axis=0)

    return av_gts, total_sites, match_dict[0], match_dict[1], match_dict[2]
    

def compkaryo(callset, inversion):

    '''Extract tag SNPs and desired specimens and calculate the average
    number of alternate alleles.'''

    targets = import_inversion(inversion)

    av_gts, total_sites, num_0, num_1, num_2 = calculate_genotype_at_concordant_sites(callset, targets)

    return av_gts, total_sites, num_0, num_1, num_2
   

def get_colour_dict(populations, palette="Set1"):
    """
    This function creates a colour palette for a provided list, returning a dict
    """

    cmap = plt.get_cmap(palette, len(np.unique(populations)))    # PiYG
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    pop_colours = {A: B for A, B in zip(np.unique(populations), colors)}
    return(pop_colours)


def legend_without_duplicate_labels(ax, **kwargs):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), **kwargs)

def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population, samples, pop_colours):
        sns.despine(ax=ax, offset=5)
        x = coords[:, pc1]
        y = coords[:, pc2]
        for pop in sample_population:
            treatment = samples[samples['partner_sample_id'] == pop]['location2'].values[0]
            flt = (sample_population == pop)
            ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[treatment], 
                    label=treatment, markersize=14, mec='k', mew=.5)

        texts = [plt.text(x[i], y[i], pop, fontsize='small', ha='center', va='center') for i, pop in enumerate(sample_population)]
        adjust_text(texts)

        ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
        ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def fig_pca(coords, model, title, path, samples, pop_colours,sample_population=None):
        if sample_population is None:
            sample_population = samples['location2'].values
        # plot coords for PCs 1 vs 2, 3 vs 4
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 2, 1)
        plot_pca_coords(coords, model, 0, 1, ax, sample_population, samples, pop_colours)
        ax = fig.add_subplot(1, 2, 2)
        plot_pca_coords(coords, model, 2, 3, ax, sample_population, samples, pop_colours)
        legend_without_duplicate_labels(ax)
        fig.suptitle(title, y=1.02)
        
        fig.savefig(path, bbox_inches='tight', dpi=300)

def pca(geno, chrom, ploidy, dataset, populations, samples, pop_colours, prune=True, scaler=None):
    if prune is True:
        if ploidy > 1:
            geno = geno.to_n_alt()
        geno = ld_prune(geno, size=500, step=200,threshold=0.2)
    else:
        if ploidy > 1:
            geno = geno.to_n_alt()
        
    coords1, model1 = allel.pca(geno, n_components=10, scaler=scaler)

    fig_pca(coords1, model1, f"PCA {chrom} {dataset}", f"gaardian/PCA/PCA-{chrom}-{dataset}", samples, pop_colours, sample_population=populations)



def ld_prune(gn, size, step, threshold=.1, n_iter=1):
    """
    Performs LD pruning, originally from Alistair Miles' blog. 
    """
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn





def load_arrays_and_metadata(local_path, _set, chrom, sitefilter='gamb_colu'):
        
    """
    This function reads genotype arrays and applys the site filter depending on 
    whether the samples are arabiensis (arab) or gambiae/coluzzii (gamb_colu)
    """
   
    snps = zarr.open_array(f"{local_path}/snp_genotypes/all/{_set}/{chrom}/calldata/GT/", mode = 'r')
    filters = zarr.open(f"{local_path}/site_filters/dt_20200416/{sitefilter}/{chrom}/variants/filter_pass", mode="r")
    positions = zarr.open_array(f"{local_path}/snp_genotypes/all/sites/{chrom}/variants/POS/", mode='r')
    positions = positions[:][filters[:]]    
    snps = allel.GenotypeDaskArray(snps).compress(filters, axis=0)
    
    metadata = pd.read_csv(f"{local_path}/metadata/samples.meta.csv")
    species = pd.read_csv(f"{local_path}/metadata/species_calls.csv")[['sample_id', 'species_gambiae_coluzzii']]
    metadata = pd.merge(metadata, species)
    
    return(snps, positions, metadata)




#####  PCA function from alistair ####

def hash_params(*args, **kwargs):
    """Helper function to hash analysis parameters."""
    o = {
        'args': args,
        'kwargs': kwargs
    }
    s = json.dumps(o, sort_keys=True).encode()
    h = hashlib.md5(s).hexdigest()
    return h


def run_pca(contig,
            gt,
            pos,
            df_samples,
    region_start=None, 
    region_stop=None,
    sample_sets="v3.2",
    sample_query=None,
    #site_mask="gamb_colu_arab",
    site_mask="gamb_colu",
    min_minor_ac=3,
    max_an_missing=0,
    n_snps=100000,
    snp_offset=0,
    n_components=10,
    results_dir=""):
    """Main function to run a PCA.
    
    Parameters
    ----------
    contig : str
        Chromosome arm, e.g., '3L'.
    region_start : int, optional
        Start position of contig region to use.
    region_stop : int, optional
        Stop position of contig region to use.
    sample_sets : str or list of str, optional
        Sample sets to analyse.
    sample_query : str, optional
        A pandas query string to select specific samples.
    site_mask : {'gamb_colu_arab', 'gamb_colu', 'arab'}
        Which site mask to apply.
    min_minor_ac : int
        Minimum minor allele count.
    max_an_missing : int
        Maximum number of missing allele calls.
    n_snps : int
        Approximate number of SNPs to use.
    snp_offset : int
        Offset when thinning SNPs.
    n_components : int
        Number of PCA components to retain.

    Returns
    -------
    data : pandas DataFrame
        Data frame with one row per sample, including columns "PC1", "PC2", etc.
    evr : numpy array
        Explained variance ratio per principal component.
    
    """
    
    # construct a key to save the results under
    results_key = hash_params(
        contig=contig,
        region_start=region_start, 
        region_stop=region_stop,
        sample_sets=sample_sets,
        sample_query=sample_query,
        site_mask=site_mask,
        min_minor_ac=min_minor_ac,
        max_an_missing=max_an_missing,
        n_snps=n_snps,
        snp_offset=snp_offset,
        n_components=n_components
    )

    # define paths for results files
    data_path = f'{results_dir}/{results_key}-data.csv'
    evr_path = f'{results_dir}/{results_key}-evr.npy'

    try:
        # try to load previously generated results
        data = pd.read_csv(data_path)
        evr = np.load(evr_path)
        return data, evr
    except FileNotFoundError:
        # no previous results available, need to run analysis
        print(f'running analysis: {results_key}')
    
    print('setting up inputs')

    if region_start or region_stop:
        # locate region within contig
        loc_region = slice(
            bisect.bisect_left(pos, region_start) if region_start else None,
            bisect.bisect_right(pos, region_stop) if region_stop else None,
        )
        gt = gt[loc_region]
    
    if sample_query:
        # locate selected samples
        loc_samples = df_samples.eval(sample_query).values
        df_samples = df_samples.loc[loc_samples, :]
        gt = da.compress(loc_samples, gt, axis=1)
        
    print('locating segregating sites within desired frequency range')

    # perform allele count
    ac = allel.GenotypeDaskArray(gt).count_alleles(max_allele=3).compute()
    
    # calculate some convenience variables
    n_chroms = gt.shape[1] * 2
    an_called = ac.sum(axis=1)
    an_missing = n_chroms - an_called
    min_ref_ac = min_minor_ac
    max_ref_ac = n_chroms - min_minor_ac

    # here we choose biallelic sites involving the reference allele
    loc_seg = np.nonzero(ac.is_biallelic() & 
                         (ac[:, 0] >= min_ref_ac) & 
                         (ac[:, 0] <= max_ref_ac) & 
                         (an_missing <= max_an_missing))[0]
    
    print('preparing PCA input data')

    # thin SNPs to approximately the desired number
    snp_step = loc_seg.shape[0] // n_snps
    loc_seg_ds = loc_seg[snp_offset::snp_step]

    # subset genotypes to selected sites
    gt_seg = da.take(gt, loc_seg_ds, axis=0)
    
    # convert to genotype alt counts
    gn_seg = allel.GenotypeDaskArray(gt_seg).to_n_alt().compute()
    
    # remove any edge-case variants where all genotypes are identical
    loc_var = np.any(gn_seg != gn_seg[:, 0, np.newaxis], axis=1)
    gn_var = np.compress(loc_var, gn_seg, axis=0)

    print('running PCA')

    # run the PCA
    coords, model = allel.pca(gn_var, n_components=n_components)
    
    # add PCs to dataframe
    data = df_samples.copy()
    for i in range(n_components):
        data[f'PC{i+1}'] = coords[:, i]
    
    # save results
    evr = model.explained_variance_ratio_
    data.to_csv(data_path, index=False)
    np.save(evr_path, evr)
    print(f'saved results: {results_key}')
    
    return data, evr
    

def plot_variance(evr, **kwargs):
    """Plot a bar chart showing variance explained by each principal
    component."""
    
    # prepare variables
    y = evr * 100
    x = [str(i+1) for i in range(len(y))]
    
    # setup plotting options
    plot_kwargs = dict(
        labels={
            'x': 'Principal component',
            'y': 'Explained variance (%)',
        },
        width=600,
        height=400
    )
    # apply any user overrides
    plot_kwargs.update(kwargs)

    # make a bar plot
    fig = px.bar(x=x, y=y, **plot_kwargs)
    fig.show()
    

def jitter(a, f):
    r = a.max() - a.min()
    return a + f * np.random.uniform(-r, r, a.shape)


def plot_coords(
    data,
    x='PC1',
    y='PC2',
    filename='foo',
    jitter_frac=0.02,
    random_seed=42,
    **kwargs,
    ):

    # setup data
    data = data.copy()
    
    # apply jitter if desired - helps spread out points when tightly clustered
    if jitter_frac:
        np.random.seed(random_seed)
        data[x] = jitter(data[x], jitter_frac)
        data[y] = jitter(data[y], jitter_frac)
            
    # convenience variables
    data['country_location'] = data['country'] + ' - ' + data['location']
    data['size'] = 1  # hack to allow us to control marker size
    
    # setup plotting options
    plot_kwargs = dict(
        width=700,
        height=500,
        hover_name='sample_id',
        hover_data=[
            'partner_sample_id',
            #'species', 
            'country', 
            'location',
            # 'insecticide',
            # 'phenotype',
            # 'exposure_time',
            # 'concentration',
            'year', 
        ],
        size='size',
        size_max=8,
        opacity=0.9,
        render_mode='svg',
    )
    # apply any user overrides
    plot_kwargs.update(kwargs)

    # 2D scatter plot
    fig = px.scatter(data, x=x, y=y, **plot_kwargs)
    #fig.write_image(f"{filename}.png", format=".png")
    fig.show()

