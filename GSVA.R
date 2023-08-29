library(GSVA)
library(data.table)
library(org.Hs.eg.db)
library(dplyr)
library(optparse)

##### Input #####
option_list <- list(
  make_option(c("-g", "--geneset"), type="character",
              help="Path to geneset file"),
  make_option(c("-e", "--expression"), type="character",
              help="Path to expression file"),
  make_option(c("-o", "--output"), type="character",
              help="Path to output file (.tsv) (GSVA results)")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)

args <- parse_args(parser)

if(length(args)==0){
  stop("All argument must be supplied ex) -g geneset.txt -e expression.tsv -o output.tsv",call.=FALSE)
}

geneset_pth <- args$geneset 
exp_pth <- args$expression 
output_pth <- args$output 

##### Load geneset #####
geneset <- fread(geneset_pth) # Read geneset file
geneset_sub <- subset(geneset, tissue == "colon") # Subset tissue of interest
# Generate gmt file
gmt <- list()
for (i in unique(geneset_sub$TF)){
  gmt[[i]] <- subset(geneset_sub, TF == i)$target
}


##### RNA-seq #####
exp <- fread(exp_pth) # Read expression file
# Convert ENSEMBL to SYMBOL (since geneset file is annnotated with SYMBOL)
ens <- exp$V1
ens_1st <- sapply(strsplit(ens, '\\.'), '[', 1)
exp$V1 <- ens_1st
annots_ <- AnnotationDbi::select(org.Hs.eg.db, keys = ens_1st, columns = "SYMBOL", keytype = "ENSEMBL")
annots <- annots_[-which(is.na(annots_$SYMBOL)),]
annots <- annots[!duplicated(annots$ENSEMBL),]
annots <- annots[!duplicated(annots$SYMBOL),]
exp_symbol <- data.frame(merge(exp, annots, by.x = "V1", by.y = "ENSEMBL"))
exp_symbol <- data.frame(exp_symbol[,-c(1, ncol(exp_symbol))], row.names = exp_symbol$SYMBOL)


##### GSVA #####
# Perform GSVA analysis using log-transformed expresson table and gmt
GSVA <- gsva(as.matrix(log2(exp_symbol+1)), gmt, min.sz = 2, method = "gsva")
GSVA <- data.frame(GSVA)

fwrite(GSVA, output_pth, sep = '\t', row.names = T, quote = F)