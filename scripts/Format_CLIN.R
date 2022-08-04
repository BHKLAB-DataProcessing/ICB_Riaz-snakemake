library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")

clin_original = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
selected_cols <- c( "Patient","Response","Dead.Alive...Dead...True.","Time.to.Death...weeks.","Subtype","M.Stage" , "Cohort" )
clin = cbind( clin_original[ , selected_cols ] , "Melanoma", "PD-1/PD-L1", NA, NA , NA , NA, NA , NA, NA , NA, NA )
colnames(clin) = c( "patient" , "recist" , "os" , "t.os"  ,"histo" , "stage" , "Cohort" ,"primary" , "drug_type" , "pfs" , "t.pfs" , "sex", "age" , "dna" , "rna" , "response.other.info" , "response" )

clin$stage = ifelse( str_detect(clin$stage, 'M0'), "III" , 
				ifelse( str_detect(clin$stage, 'M1'), "IV" , NA )) 

clin$recist[ clin$recist %in% "NE" ] = NA 

clin$response = Get_Response( data=clin )

clin$os = ifelse(clin$os %in% "TRUE" , 1 , 0)
clin$t.os = clin$t.os / 4

clin$drug_type = ifelse( clin$Cohort %in% "NIV3-PROG" , "Combo" , "PD-1/PD-L1" )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
clin$rna[ clin$patient %in% case[ case$expr %in% 1 , ]$patient ] = "tpm"
clin$dna[ clin$patient %in% case[ case$snv %in% 1 , ]$patient ] = "wes"

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

clin <- format_clin_data(clin_original, 'Patient', selected_cols, clin)

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

