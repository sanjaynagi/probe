
rule karyotype:
    """
    Perform karyotype analysis
    """
    input:
        genotypes = getZarrArray(type_="Genotypes", cloud=cloud),
        positions = getZarrArray(type_='Positions', cloud=cloud),
    output:
        "results/{dataset}_karyotypes.tsv",
    log:
        log = "logs/{dataset}_karyotypes.log"
    conda:
        "../envs/pythonGenomics.yaml"
    params:
        metadata = config['metadata'],
        dataset = config['dataset'],
        cloud = cloud,
        ag3_sample_sets = ag3_sample_sets
    script:
        "../scripts/Karyotype.py"