# PROCESSING RMSK FILE FOR scTE-----------------------------

# libs
library(data.table)
library(biomaRt)
library(fuzzyjoin)

# File download 
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
destination <- "C:/Users/qaedi/Documents/rmsk.txt.gz"

# downloading file 
download.file(url, destination) 

# reading  in R 
rmsk_file <- fread("~/rmsk.txt.gz")

# Keeping intrested columns
rmsk <- rmsk_file[, c(6:8, 11:13, 10, 17)]

# Keeping just TE elements
te <- rmsk[rmsk$V12 %in% c("LINE", "SINE", "DNA", "LTR", "Retroposon"), ]

# Keeping only intergenic elements

# retriving gene coordinates from biomart

mart <- useMart(biomart="ensembl", 
                dataset="hsapiens_gene_ensembl")

# Get attributes and filters for genes
attributes <- c("chromosome_name", "start_position", "end_position", "ensembl_gene_id")
filters <- chromosome_name %in% c(1:22, X, Y)

# Retrieve gene information
genes <- getBM(attributes = attributes,
               values = "ensembl_gene_id",
               mart = mart)

# filter genes
genes <- genes[genes$chromosome_name %in% c(1:22, "X", "Y"),]
colnames(genes) <- NULL
# add chr to the chromosome column
genes[,1] <- paste0("chr",genes[,1])

# set names
names(te)[1:3] <- c("chromosome", "start", "end")
names(genes)[1:3] <- c("chromosome", "start", "end")

# perform intersection : RETAIN INTERGENIC TEs

int_te <- genome_anti_join(te, genes, by = c("chromosome", "start", "end"))

#
fwrite(int_te, "~/hg38_intergenic_TEs.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

