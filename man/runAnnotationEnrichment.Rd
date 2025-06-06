% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFLOMICS-Methods_06_annot_enrichment.R
\name{runAnnotationEnrichment}
\alias{runAnnotationEnrichment}
\alias{runAnnotationEnrichment,RflomicsSE-method}
\alias{runAnnotationEnrichment,RflomicsMAE-method}
\alias{plotClusterProfiler}
\alias{plotClusterProfiler,RflomicsSE-method}
\alias{plotEnrichComp}
\alias{plotEnrichComp,RflomicsSE-method}
\alias{getEnrichRes}
\alias{getEnrichRes,RflomicsSE-method}
\alias{getEnrichRes,RflomicsMAE-method}
\alias{sumORA}
\alias{sumORA,RflomicsSE-method}
\alias{getEnrichSettings}
\alias{getEnrichSettings,RflomicsSE-method}
\alias{getAnnotAnalysesSummary}
\alias{getAnnotAnalysesSummary,RflomicsMAE-method}
\title{Run Gene Enrichment Analysis and process results}
\usage{
\S4method{runAnnotationEnrichment}{RflomicsSE}(
  object,
  featureList = NULL,
  from = "DiffExp",
  universe = NULL,
  database = "custom",
  domain = "no-domain",
  annotation = NULL,
  OrgDb = NULL,
  organism = NULL,
  keyType = NULL,
  pvalueCutoff = 0.05,
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  ...
)

\S4method{runAnnotationEnrichment}{RflomicsMAE}(
  object,
  SE.name,
  featureList = NULL,
  from = "DiffExp",
  universe = NULL,
  database = "custom",
  domain = "no-domain",
  annotation = NULL,
  OrgDb = NULL,
  organism = NULL,
  keyType = NULL,
  pvalueCutoff = 0.05,
  qvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  ...
)

\S4method{plotClusterProfiler}{RflomicsSE}(
  object,
  featureListName = NULL,
  database = NULL,
  domain = "no-domain",
  plotType = "dotplot",
  showCategory = 15,
  searchExpr = "",
  nodeLabel = "all",
  p.adj.cutoff = NULL,
  ...
)

\S4method{plotEnrichComp}{RflomicsSE}(
  object,
  from = "DiffExp",
  database = NULL,
  domain = "no-domain",
  matrixType = "FC",
  nClust = NULL,
  ...
)

\S4method{getEnrichRes}{RflomicsSE}(
  object,
  featureListName = NULL,
  from = "DiffExp",
  database = "GO",
  domain = NULL
)

\S4method{getEnrichRes}{RflomicsMAE}(
  object,
  experiment,
  featureListName = NULL,
  from = "DiffExp",
  database = "GO",
  domain = NULL
)

\S4method{sumORA}{RflomicsSE}(object, from = "DiffExp", database = NULL, featureListName = NULL)

\S4method{getEnrichSettings}{RflomicsSE}(object, from = "DiffExp", database = "GO")

\S4method{getAnnotAnalysesSummary}{RflomicsMAE}(object, from = "DiffExp", matrixType = "presence", ...)
}
\arguments{
\item{object}{An object of class \link{RflomicsSE} or
class \link{RflomicsMAE-class}}

\item{featureList}{name of contrasts (tags or names) from which to extract DE
genes if from is DiffExpAnal.}

\item{from}{indicates if the enrichment results are taken from differential
analysis results (DiffExp) or from the co-expression analysis
results (CoExp)}

\item{universe}{description}

\item{database}{which database (GO, KEGG, custom...)}

\item{domain}{the subonology or subdomain for the database (eg CC, MF or
BP for GO.)}

\item{annotation}{for custom annotation, a data frame of the annotation.
The data frame must contains at least two columns: gene and term, with the
omics name and the associated term id respectively. A column name
can be added with the full name of the term (if  term is not the full name
already). The column domain can be used to indicate either different databases
(grouped analyses of kegg and go for example) or different domains for
a single database (CC, MF and BP) for GO.}

\item{OrgDb}{OrgDb (with enrichGO)}

\item{organism}{supported organism listed in
'https://www.genome.jp/kegg/catalog/org_list.html' (with enrichKEGG)}

\item{keyType}{keytype of input gene with enrichGO
(one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot' with enrichKEGG)}

\item{pvalueCutoff}{adjusted pvalue cutoff on enrichment tests to report}

\item{qvalueCutoff}{qvalue cutoff on enrichment tests to report as significant.
Tests must pass
i) pvalueCutoff on unadjusted pvalues,
ii) pvalueCutoff on adjusted pvalues and
iii) qvalueCutoff on qvalues to be reported.}

\item{minGSSize}{minimal size of genes annotated by Ontology term for testing.}

\item{maxGSSize}{maximal size of genes annotated for testing}

\item{...}{more arguments for ComplexHeatmap::Heatmap.}

\item{SE.name}{SE.name the name of the dataset if the input object
is a \link{RflomicsMAE-class}}

\item{featureListName}{the contrastName or clusterName to retrieve
the results from. If NULL, all results are returned.}

\item{plotType}{type of plot. Define the function used inside.
One of dotplot, heatplot or cnetplot.}

\item{showCategory}{max number of terms to show.}

\item{searchExpr}{expression to search in the showCategory terms.}

\item{nodeLabel}{same as in enrichplot::cnetplot function, defines}

\item{p.adj.cutoff}{pvalueCutoff to define the enrichment threshold.}

\item{matrixType}{Heatmap matrix to plot, one of GeneRatio, p.adjust
or presence.}

\item{nClust}{number of separate cluster to plot on the heatmap, based o
n the clustering.}

\item{experiment}{if the object is a RflomicsMAE, then experiment is the
name of the RflomicsSE to look for.}
}
\value{
A RflomicsMAE or a RflomicsSE, depending on the class of object
parameter. The enrichment results are added to the metadata slot, either
in DiffExpEnrichAnal or CoExpEnrichAnal.

a list of tables or a table

a list with all settings
}
\description{
This function performs over representation analysis (ORA) using
clusterprofiler functions. It can be used with custom annotation file
(via enricher), GO (enrichGO) or KEGG (enrichKEGG) annotations.
}
\section{Accessors}{

A set of getters and setters generic functions to access and
modify objects of the slot metadata of a \link{RflomicsMAE-class} object or
a \link{RflomicsMAE-class} object.


\itemize{
   \item getEnrichRes: get a particular enrichment result.
   return enrichment results given in the form of lists of clusterprofiler
   results.
}


\itemize{
   \item sumORA: Get summary tables from ORA analyses -
   once an enrichment has been conducted.
}


\itemize{
   \item getEnrichSettings:
   get the settings of an enrichment analysis.
}


\itemize{
   \item getAnnotAnalysesSummary:
   return A list of heatmaps, one for each ontology/domain.
}
}

\section{Plots}{

A collection of functions for plotting results from omics analysis steps.


\itemize{
   \item plotClusterProfiler: Plot a dotplot, a cnetplot or an heatplot, using enrichplot
package. It is a wrapper method destined for the RflomicsSE class..
}


\itemize{
   \item plotEnrichComp: plot an heatmap of all the enriched term found for a given
database and a given source (differential analysis or coexpression clusters).
Allow for the comparison of several enrichment results.
}
}

\examples{
# load ecoseed data
library(RFLOMICS)
data(ecoseed.mae)

factorInfo <- data.frame(
  "factorName"   = c("Repeat", "temperature", "imbibition"),
  "factorType"   = c("batch", "Bio", "Bio")
)

# create rflomicsMAE object with ecoseed data
MAE <- createRflomicsMAE(
  projectName = "Tests",
  omicsData   = ecoseed.mae,
  omicsTypes  = c("RNAseq","proteomics","metabolomics"),
  factorInfo  = factorInfo)
#names(MAE) <- c("RNAtest", "protetest", "metatest")

# Set the statistical model and contrasts to test
formulae <- generateModelFormulae(MAE)
MAE <- setModelFormula(MAE, formulae[[1]])  

# Get the contrasts List and choose the first 3 contrasts of type averaged
contrastList <- generateExpressionContrast(MAE, "averaged")

MAE <- setSelectedContrasts(MAE, contrastList = contrastList[c(1, 2, 3),])

# Run the data preprocessing and perform the differential analysis 
MAE <- runDataProcessing(MAE, SE.name = "protetest",  
                         transformMethod = "log2",
                         normMethod = "median")

MAE <- runDiffAnalysis(MAE, 
                       SE.name = "protetest", 
                       method = "limmalmFit")

# Run GO annotation (enrichGO)
# Not run: need org.At.tair.db package
# MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest",
#                                list_args = list(OrgDb = "org.At.tair.db",
#                                                 keyType = "TAIR",
#                                                 pvalueCutoff = 0.05),
#                                from = "DiffExp", database = "GO",
#                                domain = "CC")

# Run KEGG annotation (enrichKEGG)
# need internet connection
# MAE <- runAnnotationEnrichment(MAE, SE.name = "protetest",
#                                list_args = list(organism = "ath",
#                                                 keyType = "kegg",
#                                                 pvalueCutoff = 0.05),
#                                from = "DiffExp", database = "KEGG")


# Search for the pvalue cutoff:
# sumORA(MAE[["protetest"]], from = "DiffExp", database = "KEGG")

# need internet connection
# pathview package
# plotKEGG(MAE[["protetest"]], pathway_id = "ath00710", species = "ath",
#          contrastName = "cluster.4", from = "Coexp")

# From differential analysis proteins lists:
# plotClusterProfiler(MAE[["protetest"]],
#                     contrastName = "(temperatureElevated - temperatureMedium) in mean",
#                     database = "KEGG", from = "DiffExp",
#                     plotType = "heatplot", p.adj.cutoff = 0.05,
#                     domain = "no-domain")
# 
# plotEnrichComp(MAE[["protetest"]], from = "DiffExp",
#                database = "KEGG", matrixType = "FC")

# Get all results from KEGG on differential expression lists:
# results <- getEnrichRes(MAE[["protetest"]],
#                         from = "diffexp", database = "KEGG")

# Search for the pvalue cutoff:
# usedPvalue <- 
#     getEnrichPvalue(MAE[["protetest"]], from = "diffexp", database = "KEGG")
# settings <- 
#     getEnrichSettings(MAE[["protetest"]], from = "diffexp", database = "KEGG")

}
