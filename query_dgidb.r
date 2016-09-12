##########################################################
# Query DGIDB API for gene drug interactions and parse results.
# Script is based on rDGIdb R/Bioconductor package
# API: http://dgidb.genome.wustl.edu/api
# 
# Usage: Rscript Query_DGIdb_Script.R path/to/input path/to/output min-db-support
#
# Input: - Input file name containing a list of genes
#        - Ouput file name
#
# Output: Three files are written to disk:
#         - Condensed gene interaction table with scores in brackets
#         - Complete gene interaction table with individual DBs and interaction types listed
#         - list of genes and their categories, where categories for each gene are comma-separated
#
# Version: 1.0
#
##########################################################
# CHANGE HISTORY:
#
# 2016-08-08 -- version 0.1
# Fixed bug while parsing file names. File ending is not removed anymore
# Handling of empty gene list (without stop call)
#
# 
#
# 2016-08-01 -- version 0.0
# Initial version
##########################################################


#######################################
### Prepare session, load package
#######################################
library(rDGIdb)

### Check input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript Query_DGIDB.R path/to/input path/to/output #DBs")
} else {
    srcFile <- args[1]
    outFile <- args[2]
    MIN_DB_SUPPORT <- as.integer(args[3])
}

### Read input file
if (!file.exists(srcFile))
{
    stop("Input file does not exist!")
} else {
    if (file.info(srcFile)$size == 0) {
		# write empty files, necessary for snakemake pipeline
		interactionData_file=paste(outFile, 'CompleteTable', 'txt', sep='.')
		categories_file=paste(outFile, 'GeneCategories', 'txt', sep='.')
		interactionData_cols = c('Gene','Drug','Score','Type')
		write(interactionData_cols, file=interactionData_file, ncolumns=length(interactionData_cols), sep='\t')

		categories_cols = c("Gene","GeneName","Category")
		write(categories_cols, file=categories_file, ncolumns=length(categories_cols), sep='\t')

		output_cols = c('Gene', 'Drugs')
		write(output_cols, file=outFile, ncolumns=length(output_cols), sep='\t')
        quit(save = "default", status = 0)
    } else {
        input <- read.table(srcFile, sep='\t', header = FALSE, stringsAsFactors = FALSE)
    }
}

### Map to unique official gene symbols (limma package)
if (ncol(input) != 1) stop("Wrong input format, single column gene list expected!")
genes <- unique(input$V1)
#library(limma)
#genes <- alias2Symbol(genes, species = "Hs", expand.symbols = FALSE)

#######################################
### Query DGIdb using rDGIdb package
#######################################

geneCategories <- c("clinically actionable", "druggable genome", "tumor suppressor")

result <- queryDGIdb(genes, geneCategories = geneCategories)

emptyQuery = FALSE
if(length(result@matchedTerms) == 0) emptyQuery = TRUE


if(emptyQuery == TRUE)
{	
	# write empty files, necessary for snakemake pipeline
	interactionData_file=paste(outFile, 'CompleteTable', 'txt', sep='.')
	categories_file=paste(outFile, 'GeneCategories', 'txt', sep='.')
	interactionData_cols = c('Gene','Drug','Score','Type')
    write(interactionData_cols, file=interactionData_file, ncolumns=length(interactionData_cols), sep='\t')

    categories_cols = c("Gene","GeneName","Category")
    write(categories_cols, file=categories_file, ncolumns=length(categories_cols), sep='\t')

    output_cols = c('Gene', 'Drugs')
    write(output_cols, file=outFile, ncolumns=length(output_cols), sep='\t')
	
	quit(save = "default", status = 0)
}

#######################################
### Prepare and write results
######################################

### Complete interaction table
interactionData <- resultSummary(result)
interactionData$Type <- apply(interactionData, 1, function(x, details) {
    idx <- which(x['Gene'] == details$Gene & x['Drug'] == details$Drug)
    paste(sort(unique(details$InteractionType[idx])), collapse = ",")
}, result@detailedResults)
interactionData <- interactionData[interactionData$Score >= MIN_DB_SUPPORT,]

write.table(x = interactionData,
            file = paste(outFile, 'CompleteTable', 'txt', sep='.'),
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

### Gene categories file
categories <- byGene(result)
categories <- categories[,c('Gene','GeneName','DruggableGeneCategories')]
# Remove genes that dropped out from MIN_DB_SUPPORT test
categories <- categories[categories$Gene %in% interactionData$Gene,]
write.table(x = categories,
            file = paste(outFile, 'GeneCategories', 'txt', sep='.'),
            sep="\t", row.names = FALSE, col.names = c("Gene","GeneName","Category"), quote = FALSE)


### Write summary interaction table
collapseInteractionTable <- function(gene, interactionData) {
    geneDrug <- c(gene, paste(interactionData$Drug[interactionData$Gene == gene], collapse = ","))
}
interactionData$Drug <- paste(interactionData$Drug, " ", "(", interactionData$Score, ")", sep = "")
geneDrugTable <- t(sapply(unique(interactionData$Gene), collapseInteractionTable, interactionData))
colnames(geneDrugTable) <- c("Gene", "Drug")

write.table(x = geneDrugTable, file = outFile, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
