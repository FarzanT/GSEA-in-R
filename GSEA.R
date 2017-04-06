# ======== Gene set enrichment analysis (GSEA) with Bioconductor packages ======

# "Bioconductor offers a host of GSEA packages. Which ones are good? Compare,
# contrast, and for the three or more interesting ones provide sample scripts
# that work with our subnetworks."

# First, install Bioconductor, its dependencies and recommended packages from
# the Bioconductor website
source("https://bioconductor.org/biocLite.R")
biocLite()
# Select all (a) when prompted to update currently installed recommended
# packages (and update previous ones)

# Note that biocLite() can be given package names as string arguments to
# download specific packages, instead of default of 'all'
# For a list of currently availabe packages go to:
# https://bioconductor.org/packages/release/BiocViews.html#___Software
# For a list of 'Bioconductor workflows' go to:
# https://bioconductor.org/help/workflows/

# Gene Set Enrichment Analysis uses two different types of data: a gene
# expression dataset and a list of gene sets.
# For the list of gene sets, we use MsigDB:

# ==== Download MsigDB =========================================================
# The complete curated gene sets can be downloaded from
# http://software.broadinstitute.org/gsea/msigdb/
# Note that registration is required for download

# If the 'all gene sets with gene symbols' file (tab delimited) is available in
# the current working directory, each package can load it with their own
# functions. We can also have a look at the contents by reading the file line
# by line: .gmt is a tab separated file format where each line contains a gene
# set.
MsigDB <- readLines(con = "inst/extdata/msigdb.v5.2.symbols.gmt.txt")
MsigDB <- strsplit(MsigDB, split = "\t")
MsigDB[1]
# The First line here is the gene set ID, the second line is the URL to a
# description of the gene set on MsigDB's website.
# The rest of the lines are the genes that are in the set, usually up to 256.
# ==============================================================================

# "The NCBI Gene Expression Omnibus (GEO) is a public repository of microarray
# data. Given the rich and varied nature of this resource, it is only natural
# to want to apply BioConductor tools to these data. GEOquery is the bridge
# between GEO and BioConductor."
biocLite("GEOquery")
library(GEOquery)
?getGEOfile
# Let's download gene expression data for "MicroRNA-135b overexpression effect
# on prostate cancer cell line: time course"
GDS <- getGEOfile("GDS6100", AnnotGPL = TRUE)

# Parse the time course data
# Note: The temp path may vary on differernt machines
GDS <- getGEO(filename = GDS, AnnotGPL = TRUE)

getGSEDataTables(GSE = GSE)
# We can convert a GDS class data structure to BioConductor data structure
# using GDS2eSet() (if the your file of interest has such format)
?GDS2eSet
eSet <- GDS2eSet(GDS = GDS)

# Here we use another BioConductor Package, BioBase (which should be installed
# when biocLite() was called)
library(Biobase)

# Get information about the parsed GSE (Platform should usually be indicated in
# annotation, not true for all sets)
experimentData(eSet)

# We can see that the platform for this assay is GPL10558 (?). Consequently,
# the platform specific IDs are used in the expression data:
head(featureNames(eSet))
# Which need to be converted to gene names to be useful in GSEA.
# A search for GPL10558 reveals that the platform name is
# HumanHT-12 v4 Expression BeadChip

# We can use the biomaRt package to get the gene symbols associated with these
# IDs:
biocLite("biomaRt")
library(biomaRt)
?getBM

# Get the Ensembl IDs for homo sapiens (or your organism)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://www.ebi.ac.uk")

# We want to match HGNC symbols with HumanHT-12 v4 Expression BeadChip, a quick
# search using listAttributes reveals:

listAttributes(ensembl, what = "name")
# Can grep for Illumina
grep(listAttributes(ensembl, what = "name"), pattern = "ill", value = TRUE)
# The attribute's name is illumina_humanht_12_v4

# Now for filters:
listFilters(ensembl, what = "name")
grep(listFilters(ensembl, what = "name"), pattern = "hgnc", value = TRUE)
# "hgnc_symbol" is the filter we are looking for

# Get all the genes available in the MsigDB
allGenes <- unlist(unique(MsigDB))
listEnsembl()
# Create mapping
geneMaps <- getBM(attributes = c("hgnc_symbol",
                                 "illumina_humanht_12_v4"),
                  filters = "hgnc_symbol", values = allGenes, mart = ensembl)
head(geneMaps)


ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")



# Delete genes that could not be mapped
geneMaps <- geneMaps[!geneMaps[, 2] == "", ]

# Get expression data
eSetxp <- exprs(eSet)


# Clear genes that are NA
# eSetxp <- eSetxp[!is.na(eSetxp[,1]), ]





# # For most methods, a functionaly annotated database is needed, and Bioconductor
# # offers genome wide annotation for Human, primarily based on mapping using
# # Entrez Gene identifiers (66 MB):
# biocLite("org.Hs.eg.db")
# library(org.Hs.eg.db)
# library(Annotationhub)
# annotationd
# ?org.Hs.eg.db
# columns(org.Hs.eg.db)

# There are a varierty of GSEA related packages offered on Bioconductor (74, as
# of April 2017), but only a few are related to analysis of subnetworks found
# by the rete workflow.
myGenes <-
    c(
        "ADAM9",
        "ITGAV",
        "ITGA6",
        "ITGA3",
        "ITGB5",
        "LIMK1",
        "FGFR2",
        "DLST",
        "UMPS",
        "PAK4",
        "GATAD2A",
        "C2orf65",
        "DOK1",
        "DQX1",
        "LOXL3",
        "SEMA4F"
    )


# ==== clusterProfiler (top 5%) ====
biocLite("clusterProfiler")
library(clusterProfiler)
?clusterProfiler

# This package offers multiple GSEA methods as well as DAVID-style enrichment
# analysis.
# Methods include:

# ==== gseGO ====
# Gene Set Enrichment Analysis of Gene Ontology
# The gene list give to the function below is assumed to be ranked.
clusterProfiler::gseGO(geneList = myGenes, OrgDb = "HGNC", keytype = "HGNC")

# ==== gseKEGG ====

# ==== gseMKEGGG ====

# =========================== fgsea (top 5%) ===================================
BiocInstaller::biocLite("fgsea")
library(fgsea)
?fgsea
# Load example pathways and ranks
data("examplePathways")
data("exampleRanks")

# Read MsigDB file as GMT pathways
MsigDBpathways <- gmtPathways(gmt.file = "inst/extdata/msigdb.v5.2.symbols.gmt.txt")
# Must pass a named numeric vector with gene-level statistics
myStats <- runif(length(myGenes), 50, 100)
names(myStats) <- myGenes
# Set everything else to a value a lot less than the hot genes

# Calculate GSEA stat:
# Usually the statistic is based on the observations from microarray data, but
# one can use 'heat scores' from rete's graph as such statistics
calcGseaStat(myStats, selectedStats = c(1,3,5))

# Calculate gene set enrichment in MsigDB curated gene sets
fgseaRes <- fgsea(MsigDBpathways, exampleRanks, nperm = 10000, maxSize = 500)


