library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

snv = as.data.frame( fread( file.path(input_dir, "SNV.txt.gz") , stringsAsFactors=FALSE , sep="\t" ))

data = cbind( snv[ , c("Start" , "Patient" , "Hugo Symbol", "Variant Classification"  ) ] ,
				sapply( snv[ , "Chromosome" ] , function(x){ paste( "chr" , x , sep="" ) } ) ,
				sapply( snv[ , "HGVS_c" ] , function(x){ z = unlist( strsplit( x , ">" , fixed=TRUE ) )[1] ; substr( z , nchar(z) ,  nchar(z) ) } ) ,
				sapply( snv[ , "HGVS_c" ] , function(x){ unlist( strsplit( x , ">" , fixed=TRUE ) )[2] } )
			)
colnames(data) = c( "Pos" , "Sample" , "Gene" , "Effect" , "Chr" , "Ref" , "Alt"  )

data$Ref = ifelse( data$Ref %in% "-" , "" , data$Ref )
data$Alt = ifelse( data$Alt %in% "-" , "" , data$Alt )

data = cbind( data ,
				apply( data[ , c( "Ref", "Alt" ) ] , 1 , function(x){ ifelse( nchar(x[1]) != nchar(x[2]) , "INDEL", "SNV") } )
			)

colnames(data) = c( "Pos" , "Sample" , "Gene" , "Effect" , "Chr" , "Ref" , "Alt" , "MutType"  )


case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
data = data[ data$Sample %in% case[ case$snv %in% 1 , ]$patient , c( "Sample" , "Gene" , "Chr" , "Pos" , "Ref" , "Alt" , "Effect" , "MutType" ) ]


data$Effect = ifelse( data$Effect %in% c( "missense_variant","missense_variant&splice_region_variant" ) , "Missense_Mutation" , data$Effect )
data$Effect = ifelse( data$Effect %in% c( "splice_acceptor_variant&intron_variant" , "splice_acceptor_variant&splice_donor_variant&intron_variant" , "splice_donor_variant&intron_variant" , "start_lost&splice_region_variant" ,"stop_gained&splice_region_variant" ) , "Splice_Site" , data$Effect )
data$Effect = ifelse( data$Effect %in% "stop_gained" , "Nonsense_Mutation" , data$Effect )
data$Effect = ifelse( data$Effect %in% "stop_lost" , "Stop_Codon_Del" , data$Effect )
data$Effect = ifelse( data$Effect %in% "start_lost" , "Start_Codon_Ins" , data$Effect )                

write.table( data , file=file.path(output_dir, "SNV.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
