# ======== Gene set enrichment analysis (GSEA) with Bioconductor packages ======

# "Bioconductor offers a host of GSEA packages. Which ones are good? Compare,
# contrast, and for the three or more interesting ones provide sample scripts
# that work with our subnetworks."



# First, install Bioconductor, its dependencies and recommended packages from
# the Bioconductor website
source("https://bioconductor.org/biocLite.R")
biocLite()
# Select all (a) when prompted to update currently installed recommended
# packages

# Note that biocLite() can be given package names as string arguments to
# download specific packages, instead of default of 'all'
# For a list of currently availabe packages go to:
# https://bioconductor.org/packages/release/BiocViews.html#___Software
# For a list of 'Bioconductor workflows' go to:
# https://bioconductor.org/help/workflows/

x <- system.file("extdata","fastMapUniProt.rds", package = "rete", mustWork = TRUE)
load(file = x)

# Before continuing any further it is important to distinguish the following
# two approaches:
#
# "DAVID" determines overlaps between user-supplied gene lists and the curated
# databases, looking for overlaps that are bigger than that expected by random
# chance. You can improve the accuracy of the algorithm by providing a
# background file that contains all genes that were considered/detected in the
# experiment.

# "GSEA" is a tool that uses every datapoint in its statistical algorithm. In the
# "classical" method genes are ranked by from most up-regulated to most
# down-regulated. The rank metric itself varies, but two valid methods are to
# use signed p-value, or lower 90% confidence interval of the fold change. The
# basis of the test is to assess whether members of a gene set appear enriched
# at one end of the profile. To test the enrichment, GSEA performs permutations
# of the profile, calculating the enrichment of the gene set a thousand or more
# times to estimate p-values empirically.

# For most methods, a functionaly annotated database is needed, and Bioconductor
# offers genome wide annotation for Human, primarily based on mapping using
# Entrez Gene identifiers (66 MB):
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
library(Annotationhub)
annotationd
?org.Hs.eg.db
columns(org.Hs.eg.db)
data("inst/extdata/fastMapUniProt.rds")
data("inst/extdata/fastMapENSP.rds")
# There are a varierty of GSEA related packages offered on Bioconductor (74, as
# of April 2017), but only a few are related to analysis of subnetworks found
# by the rete workflow.

# Let's use the miniDatasets for validation to test the packages below.
# Load devPPI
miniPPI <- read.delim("inst/extdata/devPPI.txt")
# Read the interactor uniProtID and refseq IDs
interactorA <- gsub(pattern = ".*:(.*)", replacement = "\\1", x = miniPPI[,1])
interactorB <- gsub(pattern = ".*:(.*)", replacement = "\\1", x = miniPPI[,2])
allInteractors <- unique(c(interactorA, interactorB))

# Get HGNC gene symbols for the uniProtIDs (refseq unavailable)
# load("inst/extdata/fastMapUniProt.rds")
# fastMap(ID = allInteractors, , )

myGenes <- c(
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
