library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

rna = as.matrix( fread( file.path(input_dir, "EXPR_count.txt.gz") , stringsAsFactors=FALSE  , sep="\t" ) )
rna = rna[ , grep( "_Pre" , colnames(rna) ) & !colnames(rna) %in% c("Pt37_Pre","Pt3_Pre","Pt5_Pre","Pt78_Pre","Pt92_Pre") ]
rna = sort( unique( sapply( colnames(rna) , function(x){ unlist( strsplit( x , "_" , fixed=TRUE ) )[1] } ) ) )[-1]

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
rownames(clin) = clin$Patient
patient = sort( unique( clin$Patient ) )

case = as.data.frame( cbind( patient , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) ) )
colnames(case) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = patient

case$snv = as.numeric( as.character( case$snv ) )
case$cna = as.numeric( as.character( case$cna ) )
case$expr = as.numeric( as.character( case$expr ) )

for( i in 1:nrow(case)){
	if( !is.na( clin[ rownames(case)[i] , ]$Mutation.Load ) ){
		case$snv[i] = 1
		case$cna[i] = 1
	}
	if( rownames(case)[i] %in% rna ){
		case$expr[i] = 1
	}
}

write.table( case , file=file.path(output_dir, "cased_sequenced.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
