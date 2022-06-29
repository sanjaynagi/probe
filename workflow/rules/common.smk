
#from tools import get_colour_dict
#from tools import getCohorts
#import matplotlib.pyplot as plt
#import matplotlib
#import numpy as np

rule set_kernel:
    input:
        srcdir('env/pythonGenomics.yml')
    output:
        touch(f"resources/.kernel.set")
    conda: 'env/pythonGenomics.yml'
    shell: 
        """
        python -m ipykernel install --user --name probe
        """

def getCohorts(metadata, columns=['species_gambiae_coluzzii', 'location'], comparatorColumn=None, minPopSize=15):
    
    # subset metadata dataFrame and find combinations of variables with more than minPopSize individuals
    cohorts = metadata[columns]
    cohorts = cohorts.groupby(columns).size().reset_index().rename(columns={0:'size'})
    cohorts = cohorts[cohorts['size'] > minPopSize][columns]
    
    
    if comparatorColumn != None:
        cols = [i for i in columns if i != comparatorColumn]
    else:
        cols = columns

    idxs = []
    for _, row in cohorts.iterrows():   
        # create the pandas metadata query for each cohort
        mycohortQuery = " & ".join([col + " == " + "'" + row.astype(str)[col] + "'" for col in cohorts.columns])
        # get indices of individuals in each cohort
        idxs.append(metadata.query(mycohortQuery).index.tolist())
    
    cohorts['indices'] = idxs
    cohorts['cohortText'] = cohorts[cols].agg(' | '.join, axis=1)
    cohorts['cohortNoSpaceText'] = cohorts['cohortText'].str.replace("|", ".", regex=False).str.replace(" ", "",regex=False)
    #colours = get_colour_dict(cohorts['species_gambiae_coluzzii'], palette="Set1")
    #cohorts['colour'] = cohorts['species_gambiae_coluzzii'].map(colours)
    if comparatorColumn != None: 
        cols = cols + ['cohortText', 'cohortNoSpaceText']
        cohorts = cohorts.pivot(index=cols, columns=comparatorColumn)
        return(cohorts.reset_index())

    return(cohorts.reset_index(drop=True))


def getZarrArray(type_="Genotype", all_contigs=False, cloud=False):

    if cloud is False:
        if config['Zarr']['activate'] == True:
            Array = config['Zarr'][type_]
            if all_contigs == True:
                Array = Array.replace("{contig}", "{{contig}}")
        elif config['Zarr']['activate'] == False:
            if type_ == "Genotype":
                Array = "resources/Zarr/{dataset}/{contig}/calldata/GT" 
            if type_ == 'Haplotype':
                Array == "resources/Zarr/{dataset}/{contig}/calldata/GT" 
            elif type_ == "Positions":
                Array = "resources/Zarr/{dataset}/{contig}/variants/POS"
            elif type_ == 'SiteFilters' and siteFilters is not None:
                Array = "resources/Zarr/{dataset}/{contig}/variants/siteFilter"

            if all_contigs == True:
                Array = Array.replace("{contig}", "{{contig}}")
        else:
            Array = []
    else:
        Array = [] # if using VObs cloud just return an empty list

    return(Array)


def getVCFs(gz=True, allelism = 'biallelic', bothAllelisms=False, allcontigs=False, allcontigsseparately=False):

    if allcontigs == False:
        if config['VCF']['activate'] == True:
            genotypes = config['VCF'][allelism]
        elif gz == True:
            genotypes = expand("resources/vcfs/{dataset}_{{contig}}.{{allelism}}.vcf.gz", dataset=config['dataset']) if bothAllelisms == True else expand("resources/vcfs/{dataset}_{{contig}}.{allelism}.vcf.gz", dataset=config['dataset'], allelism=allelism)
        elif gz == False:
            genotypes = expand("resources/vcfs/{dataset}_{{contig}}.{{allelism}}.vcf", dataset=config['dataset']) if bothAllelisms == True else expand("resources/vcfs/{dataset}_{{contig}}.{allelism}.vcf", dataset=config['dataset'], allelism=allelism)
    elif allcontigs == True:
        if config['VCF']['activate'] == True:
            genotypes = expand("resources/vcfs/wholegenome/{dataset}.{{allelism}}.vcf.gz", dataset=config['dataset']) if bothAllelisms == True else expand("resources/vcfs/wholegenome/{dataset}.{allelism}.vcf.gz", dataset=config['dataset'], allelism=allelism)
        elif gz == True:
            genotypes = expand("resources/vcfs/wholegenome/{dataset}.{{allelism}}.vcf.gz", dataset=config['dataset']) if bothAllelisms == True else expand("resources/vcfs/wholegenome/{dataset}.{allelism}.vcf.gz", dataset=config['dataset'], allelism=allelism)
        elif gz == False:
            genotypes = expand("resources/vcfs/wholegenome/{dataset}.{{allelism}}.vcf", dataset=config['dataset']) if bothAllelisms == True else expand("resources/vcfs/wholegenome/{dataset}.{allelism}.vcf", dataset=config['dataset'], allelism=allelism)
    
    if allcontigsseparately:
        genotypes = expand(genotypes, contig=contigs)

    return(genotypes)



def getSelectedOutputs(wildcards):

    """
    Function that returns a list of the desired outputs for the rule all, depending on the config.yaml
    configuration file.
    """

    selected_input = []
   
   # selected_input.extend(expand("resources/vcfs/{dataset}_{contig}.{allelism}.vcf.gz", dataset=config['dataset'], contig=config['contigs'], allelism=['biallelic']))

    if config['PopulationStructure']['Relatedness']['activate']:
        selected_input.extend(
            expand(
                "results/relatedness/ngsRelate.{dataset}",
            dataset=config['dataset'],
            )
        )

    if config['PopulationStructure']['PCA']['activate']:
        selected_input.extend(
            expand(
                [
                "results/PCA/{cohort}.{contig}.html",
                ],
                contig=config['contigs'], 
                cohort=PCAcohorts['cohortNoSpaceText'],
            )
        )

    if config['Selection']['VariantsOfInterest']['activate']:
        selected_input.extend(
            expand(
                [
                    "results/variantsOfInterest/VOI.{dataset}.frequencies.tsv",
                     "results/variantsOfInterest/VOI.{dataset}.heatmap.png"
                ],    
                dataset=config['dataset'],
            )
        )


    for stat in ['H12', 'G12','G123']:
        if config['Selection'][stat]['activate']:
            selected_input.extend(
                expand(
                [
                    "results/selection/{stat}/{stat}_{cohort}.{contig}.png"
                ],    
                    contig=config['contigs'], 
                    cohort=cohorts['cohortNoSpaceText'],
                    stat=stat
                )
            )

    if config['Selection']['PBS']['activate']:
        selected_input.extend(
            expand(
            [
                "results/selection/PBS/PBS_{cohort}.{contig}.png"
            ],    
                contig=config['contigs'], 
                cohort=PBScohorts['cohortNoSpaceText'],
                stat=stat
            )
        )


    if config['PopulationStructure']['f2_variants']['activate']:
        selected_input.extend(
            expand(
            [
                "results/f2HapLengths_{contig}.tsv",
            ],    
                contig=config['contigs'], 
            )
        )

    if config['karyotype']['activate']:
        selected_input.extend(
            expand(
            [
                "results/{dataset}_karyotypes.tsv",
            ],    
                dataset=config['dataset'], 
            )
        )

    return(selected_input)



def singleTrue(iterable):
    iterator = iter(iterable)
#    # consume from "i" until first true or it's exhausted
    has_true = any(iterator) 
#    # carry on consuming until another true value / exhausted
    has_another_true = any(iterator) 
#    # True if exactly one true found
    return has_true and not has_another_true