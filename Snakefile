from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]

rule get_MultiAssayExp:
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv"),
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    output:
        S3.remote(prefix + filename)
    resources:
        mem_mb=3000
    shell:
        """
        Rscript -e \
        '
        load(paste0("{prefix}", "annotation/Gencode.v40.annotation.RData"))
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Riaz", input_dir = paste0("{prefix}", "processed"), expr_with_counts_isoforms=TRUE), 
            "{prefix}{filename}"
        );
        '
        """

rule format_snv:
    input:
        S3.remote(prefix + "download/SNV.txt.gz"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    output:
        S3.remote(prefix + "processed/SNV.csv")
    resources:
        mem_mb=3000
    shell:
        """
        Rscript scripts/Format_SNV.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_expr:
    input:
        S3.remote(prefix + "download/expr_list.rds")
    output:
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv")
    resources:
        mem_mb=3000
    shell:
        """
        Rscript scripts/Format_EXPR.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_clin:
    input:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    output:
        S3.remote(prefix + "processed/CLIN.csv")
    resources:
        mem_mb=1000
    shell:
        """
        Rscript scripts/Format_CLIN.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation
        """

rule format_cased_sequenced:
    input:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/expr_list.rds")
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=1000
    shell:
        """
        Rscript scripts/Format_cased_sequenced.R \
        {prefix}download \
        {prefix}processed
        """

rule format_downloaded_data:
    input:
        S3.remote(prefix + "download/1-s2.0-S0092867417311224-mmc2.xlsx"),
        S3.remote(prefix + "download/1-s2.0-S0092867417311224-mmc3.xlsx"),
        S3.remote(prefix + "download/rnaseq_samples.tsv"),
        S3.remote(prefix + "download/Riaz_kallisto1.zip"),
        S3.remote(prefix + "download/Riaz_kallisto2.zip"),
        S3.remote(prefix + "download/Riaz_kallisto3.zip"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    output:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/expr_list.rds"),
        S3.remote(prefix + "download/SNV.txt.gz")
    shell:
        """
        Rscript scripts/format_downloaded_data.R \
        {prefix}download \
        {prefix}annotation     
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    shell:
        """
        wget -O {prefix}annotation/Gencode.v40.annotation.RData https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_drug.csv -O {prefix}annotation/curation_drug.csv
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_tissue.csv -O {prefix}annotation/curation_tissue.csv 
        """

rule download_data:
    output:
        S3.remote(prefix + "download/1-s2.0-S0092867417311224-mmc2.xlsx"),
        S3.remote(prefix + "download/1-s2.0-S0092867417311224-mmc3.xlsx"),
        S3.remote(prefix + "download/rnaseq_samples.tsv"),
        S3.remote(prefix + "download/Riaz_kallisto1.zip"),
        S3.remote(prefix + "download/Riaz_kallisto2.zip"),
        S3.remote(prefix + "download/Riaz_kallisto3.zip")
    resources:
        mem_mb=1000
    shell:
        """
        wget -O {prefix}download/Riaz_kallisto1.zip https://zenodo.org/record/6968453/files/Riaz_kallisto1.zip?download=1
        wget -O {prefix}download/Riaz_kallisto2.zip https://zenodo.org/record/6968453/files/Riaz_kallisto2.zip?download=1
        wget -O {prefix}download/Riaz_kallisto3.zip https://zenodo.org/record/6968453/files/Riaz_kallisto3.zip?download=1
        wget -O {prefix}download/1-s2.0-S0092867417311224-mmc2.xlsx https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311224-mmc2.xlsx
        wget -O {prefix}download/1-s2.0-S0092867417311224-mmc3.xlsx https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311224-mmc3.xlsx
        wget -O {prefix}download/rnaseq_samples.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA356761&result=read_run&fields=run_accession,sample_title&format=tsv&download=true&limit=0"
        """ 
