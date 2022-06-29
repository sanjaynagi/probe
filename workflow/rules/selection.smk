
## Optionally specifify PBS
## G123/

rule VariantsOfInterest:
    """
    This rule reports and plots allele frequencies of Variants of Interest specified in VariantsOfInterest.tsv.
    """
    input:
        #genotypes = getZarrArray(type_="Genotypes", all_contigs=True),
        #positions = getZarrArray(type_='Positions', all_contigs=True),
        genotypes = expand(config['Zarr']['Genotypes'], contig = contigs) if not cloud else [],
        positions = expand(config['Zarr']['Positions'], contig = contigs) if not cloud else [],
        variants = config['Selection']['VariantsOfInterest']['path']
    output:
        "results/variantsOfInterest/VOI.{dataset}.heatmap.png",
        "results/variantsOfInterest/VOI.{dataset}.frequencies.tsv"
    log:
        "logs/variantsOfInterest_{dataset}.log"
    priority: 50
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        genotypePath = lambda wildcards: config['Zarr']['Genotypes'],
        positionPath = lambda wildcards: config['Zarr']['Positions'],
        dataset= config['dataset'],
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        minPopSize = config['Selection']['VariantsOfInterest']['minPopSize'],
        cloud = cloud,
        ag3_sample_sets = ag3_sample_sets
    script:
        "../scripts/VariantsOfInterest.py"


rule G12:
    """
    This rule performs G12 selection scans on each specified population
    """
    input:
        nb = "workflow/notebooks/GarudsStatistics.ipynb",
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        siteFilters = getZarrArray(type_ = "SiteFilters", cloud=cloud)
    output:
        nb = "results/notebooks/GarudsStatistics_G12_{contig}.ipynb",
        html = "results/notebooks/GarudsStatistics_G12_{contig}.html",
        plot = expand("results/selection/G12/G12_{cohort}.{{contig}}.png", cohort=cohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/G12/G12_{cohort}.{{contig}}.tsv", cohort=cohorts['cohortNoSpaceText'])
    log:
        "logs/selection/G12.{contig}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        cloud = cloud,
        ag3_sample_sets = ag3_sample_sets,
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        GarudsStat = "G12",
        windowSize = config['Selection']['G12']['windowSize'],
        windowStep = config['Selection']['G12']['windowStep'],
        cutHeight = config['Selection']['G12']['cutHeight'],
        minPopSize = 15,
    shell:
        """
        papermill {input.nb} {output.nb} -k probe -p cloud {params.cloud} -p ag3_sample_sets {params.ag3_sample_sets} \
        -p contig {wildcards.contig} -p stat G12 -p windowSize {params.windowSize} -p windowStep {params.windowStep} \ 
        -p cutHeight {params.CutHeight} -p metaColumns {params.columns} -p minPopSize {params.minPopSize}

        python -m nbconvert {output.nb} --to html --stdout --no-input \
             --ExecutePreprocessor.kernel_name=probe > {output.html}
        """

rule G123:
    """
    This rule performs G123 selection scans on each specified population
    """
    input:
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        siteFilters = getZarrArray(type_ = "SiteFilters", cloud=cloud)
    output:
        plot = expand("results/selection/G123/G123_{cohort}.{{contig}}.png", cohort=cohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/G123/G123_{cohort}.{{contig}}.tsv", cohort=cohorts['cohortNoSpaceText'])
    log:
        "logs/selection/G123.{contig}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        cloud=cloud, 
        ag3_sample_sets = ag3_sample_sets,
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        GarudsStat = "G123",
        windowSize = config['Selection']['G123']['windowSize'],
        windowStep = config['Selection']['G123']['windowStep'],
        cutHeight = config['Selection']['G123']['cutHeight'],
        minPopSize = 15
    script:
        "../scripts/GarudsStatistics.py"


rule H12:
    """
    This rule performs H12 selection scans on each specified population
    """
    input:
        nb = "workflow/notebooks/GarudsStatistics.ipynb",
        haplotypes = getZarrArray(type_="Haplotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        siteFilters = getZarrArray(type_ = "SiteFilters", cloud=cloud)
    output:
        nb = "results/notebooks/GarudsStatistics_G12_{contig}.ipynb",
        html = "results/notebooks/GarudsStatistics_G12_{contig}.html",
        plot = expand("results/selection/H12/H12_{cohort}.{{contig}}.png", cohort=cohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/H12/H12_{cohort}.{{contig}}.tsv", cohort=cohorts['cohortNoSpaceText'])
    log:
        "logs/selection/H12.{contig}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        cloud = cloud, 
        ag3_sample_sets = ag3_sample_sets,
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        GarudsStat = 'H12',
        windowSize = config['Selection']['H12']['windowSize'],
        windowStep = config['Selection']['H12']['windowStep'],
        minPopSize = 15,
    shell:
        """
        papermill {input.nb} {output.nb} -k probe -p cloud {params.cloud} -p ag3_sample_sets {params.ag3_sample_sets} \
        -p contig {wildcards.contig} -p stat G12 -p windowSize {params.windowSize} -p windowStep {params.windowStep} \ 
        -p cutHeight {params.CutHeight} -p metaColumns {params.columns} -p minPopSize {params.minPopSize}

        python -m nbconvert {output.nb} --to html --stdout --no-input \
             --ExecutePreprocessor.kernel_name=probe > {output.html}
        """


rule PopulationBranchStatistic:
    """
    This rule performs PBS selection scans on each specified population
    """
    input:
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
        siteFilters = getZarrArray(type_ = "SiteFilters", cloud=cloud),
        outgroupPath = "/home/sanj/ag1000g/data/phase3/snp_genotypes/all/AG1000G-ML-B/{contig}/calldata/GT/",
        outgroupMetaPath = "/home/sanj/ag1000g/data/phase3/metadata/general/AG1000G-ML-B/samples.meta.csv"
    output:
        plot = expand("results/selection/PBS/PBS_{cohort}.{{contig}}.png", cohort=PBScohorts['cohortNoSpaceText']),
        tsv = expand("results/selection/PBS/PBS_{cohort}.{{contig}}.tsv", cohort=PBScohorts['cohortNoSpaceText'])
    log:
        "logs/selection/PBS.{contig}.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        cloud = cloud, 
        ag3_sample_sets=ag3_sample_sets,
        metadata = config['metadata'],
        columns = config['metadataCohortColumns'],
        comparatorColumn = config['Selection']['PBS']['metadataComparatorColumn'],
        windowSize = config['Selection']['PBS']['windowSize'],
        windowStep = config['Selection']['PBS']['windowStep'],
        minPopSize = 15
    script:
        "../scripts/PopulationBranchStatistic.py"

