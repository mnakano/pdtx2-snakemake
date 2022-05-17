from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]
preprocessed=config["preprocessed"]
annotation=config["annotation"]

rule get_pset:
    input:
        S3.remote(prefix + "annotation/drugs_with_ids.csv"),
        S3.remote(prefix + "preprocessed_data/cell.csv"),
        S3.remote(prefix + "preprocessed_data/cell_annotation_all.csv"),
        S3.remote(prefix + "preprocessed_data/final_eset.Rda"),
        S3.remote(prefix + "preprocessed_data/raw_drug.csv"),
        S3.remote(prefix + "drug_sensitivity/info.rds"),
        S3.remote(prefix + "drug_sensitivity/dose_viability.rds"),
        S3.remote(prefix + "drug_sensitivity/raw.rds"),
        S3.remote(prefix + "drug_sensitivity/profiles.rds")
    output:
        S3.remote(prefix + filename)
    shell:
        """
        Rscript scripts/getPDTX.R \
        {prefix}preprocessed_data/ \
        {prefix}drug_sensitivity/ \
        {prefix}annotation/ \
        {prefix} \
        {filename}
        """

rule normalize_and_compute_sens:
    input:
        S3.remote(prefix + "published_data/RawDataDrugsSingleAgents.txt"),
        S3.remote(prefix + "published_data/DrugResponsesAUCSamples.txt")
    output:
        S3.remote(prefix + "drug_sensitivity/info.rds"),
        S3.remote(prefix + "drug_sensitivity/dose_viability.rds"),
        S3.remote(prefix + "drug_sensitivity/raw.rds"),
        S3.remote(prefix + "drug_sensitivity/profiles.rds")
    shell:
        """
        Rscript scripts/sensitivity.R \
        {prefix}published_data/ \
        {prefix}drug_sensitivity/
        """

rule get_published_data:
    output:
        S3.remote(prefix + "published_data/RawDataDrugsSingleAgents.txt"),
        S3.remote(prefix + "published_data/DrugResponsesAUCSamples.txt")
    shell:
        """
        wget -O {prefix}published_data/16 https://figshare.com/ndownloader/articles/2069274/versions/16
        unzip -o {prefix}published_data/16 -d {prefix}published_data
        find {prefix}published_data ! -name 'RawDataDrugsSingleAgents.txt' ! -name 'DrugResponsesAUCSamples.txt' -type f -delete
        """

rule get_preprocessed_data:
    output:
        S3.remote(prefix + "preprocessed_data/cell.csv"),
        S3.remote(prefix + "preprocessed_data/cell_annotation_all.csv"),
        S3.remote(prefix + "preprocessed_data/final_eset.Rda"),
        S3.remote(prefix + "preprocessed_data/raw_drug.csv")
    shell:
        """
        wget -O {prefix}preprocessed_data/cell.csv {preprocessed}cell.csv
        wget -O {prefix}preprocessed_data/cell_annotation_all.csv {preprocessed}cell_annotation_all.csv
        wget -O {prefix}preprocessed_data/final_eset.Rda {preprocessed}final_eset.Rda
        wget -O {prefix}preprocessed_data/raw_drug.csv {preprocessed}raw_drug.csv
        """

rule get_annotation:
    output:
        S3.remote(prefix + "annotation/drugs_with_ids.csv")
    shell:
        """
        wget -O {prefix}annotation/drugs_with_ids.csv {annotation}drugs_with_ids.csv
        """
