library(data.table)
library(readxl) 
library(stringr)
library(tximport)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
annot_dir <- args[2]

# CLIN.txt
clin <- read_excel(file.path(work_dir, '1-s2.0-S0092867417311224-mmc2.xlsx'), sheet='Table S2')
colnames(clin) <- str_replace_all(clin[2, ], '\\W', '.')
clin <- clin[-c(1:2), ]
clin[clin == 'NA'] <- NA
numcols <- c('Time.to.Death...weeks.', 'Mutation.Load', 'Neo.antigen.Load', 'Neo.peptide.Load', 'Cytolytic.Score')
clin[, numcols] <- sapply(clin[, numcols], as.numeric)
clin$Dead.Alive...Dead...True. <- sapply(clin$Dead.Alive...Dead...True., as.logical)
write.table(clin, file.path(work_dir, 'CLIN.txt'), col.names=TRUE, sep='\t')

# expr_list.rds
source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')
load(file.path(annot_dir, "Gencode.v40.annotation.RData"))

dir.create(file.path(work_dir, 'rnaseq'))
zipfiles <- c('Riaz_kallisto1.zip', 'Riaz_kallisto2.zip', 'Riaz_kallisto3.zip')
for(zipfile in zipfiles){
  unzip(file.path(work_dir, zipfile), exdir=file.path(work_dir, 'rnaseq'))
}
unlink(file.path(work_dir, 'rnaseq', '__MACOSX'), recursive = TRUE)

process_kallisto_output(work_dir, tx2gene)

rnaseq_samples <- read.table(file.path(work_dir, 'rnaseq_samples.tsv'), header=TRUE, sep='\t')
rnaseq_samples <- rnaseq_samples[str_detect(rnaseq_samples$sample_title, 'Pre_'), ]
rnaseq_samples$sample_title <- str_replace(rnaseq_samples$sample_title, '_Pre.+', '')
rnaseq_samples <- rnaseq_samples[!rnaseq_samples$sample_title %in% c("Pt37","Pt3","Pt5","Pt78","Pt92"), ]
expr_list <- readRDS(file.path(work_dir, 'expr_list.rds'))

for(assay_name in names(expr_list)){
  expr <- data.frame(expr_list[[assay_name]])
  expr <- expr[, rnaseq_samples$run_accession]
  colnames(expr) <- unlist(lapply(colnames(expr), function(col){
    return(
      rnaseq_samples$sample_title[rnaseq_samples$run_accession == col]
    )
  }))
  expr_list[[assay_name]] <- expr
}

saveRDS(expr_list, file.path(work_dir, 'expr_list.rds'))

# SNV.txt.gz
snv <- read_excel(file.path(work_dir, '1-s2.0-S0092867417311224-mmc3.xlsx'), sheet='Table S3')
colnames(snv) <- snv[3, ]
snv <- snv[-c(1:3), ]
numcols <- c('Start', 'End', 'Tcov', 'Tac', 'Taf')
snv[, numcols] <- sapply(snv[, numcols], as.numeric)
gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
write.table( snv , file=gz , quote=FALSE , sep=";" , col.names=TRUE, row.names = FALSE )
close(gz)

