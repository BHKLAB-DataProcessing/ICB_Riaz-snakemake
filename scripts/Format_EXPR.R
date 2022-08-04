library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr_list <- readRDS(file.path(input_dir, 'expr_list.rds'))
for(assay_name in names(expr_list)){
  write.table( 
    expr_list[[assay_name]], 
    file= file.path(output_dir, paste0('EXPR_', str_replace(assay_name, 'expr_', ''), '.csv')), 
    quote=FALSE, 
    sep=";", 
    col.names=TRUE, 
    row.names=TRUE 
  )
}

