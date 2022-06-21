library(data.table)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

expr = as.matrix( fread( file.path(input_dir, "EXPR_count.txt.gz") , stringsAsFactors=FALSE  , sep="\t" ) )
rownames(expr) = expr[,"HUGO"]
expr = expr[ !is.na(rownames(expr)) , !colnames(expr) %in% "HUGO" ] 

rid = rownames(expr)
cid = colnames(expr)
expr = apply(apply(expr,2,as.character),2,as.numeric)
colnames(expr) = cid
rownames(expr) = rid


##################################
## Get Gene Length
genes <- rownames(expr)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol",values=genes, mart=human)
size=gene_coords$end_position - gene_coords$start_position
names(size) = gene_coords[,"hgnc_symbol"]

##################################
##Remove Duplicated Genes
t_uniq <- size[ !( names(size) %in% names(size)[ duplicated(names(size)) ]) ]
t_dup <- size[ ( names(size) %in% names(size)[ duplicated(names(size)) ]) ]

t_dup <- t_dup[order(names(t_dup))]
id <- unique(names(t_dup))

t.dup.rm <- NULL
for(j in 1:length(id)){
	tmp <- t_dup[which(names(t_dup)%in%id[j])]
	tmp = mean(tmp, na.rm=T)
	t.dup.rm <- c(t.dup.rm,tmp)			
}
size <- c(t_uniq,t.dup.rm)
names(size) <- c(names(t_uniq),id)

expr = expr[names(size),]

##################################
## Compute TPM data

GetTPM <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

TPM = log2( GetTPM(expr,size) + 1 )

TPM = TPM[ , grep( "_Pre" , colnames(TPM) ) ]
TPM = TPM[ , !colnames(TPM) %in% c("Pt37_Pre","Pt3_Pre","Pt5_Pre","Pt78_Pre","Pt92_Pre") ]

colnames(TPM) = sapply( colnames( TPM ) , function(x){ unlist( strsplit( x , "_" , fixed=TRUE ) )[1] } )

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
TPM = TPM[ ,  colnames( TPM ) %in% case[ case$expr %in% 1 , ]$patient ]

write.table( TPM , file=file.path(output_dir, "EXPR.csv"), quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )
