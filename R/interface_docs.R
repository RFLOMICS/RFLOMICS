
# ---- 03 data exploratory ----
data_explor_docs <- list(
  p(
    "Very important and non trivial steps, crucial for single-omics analysis and 
    for the integration of different omics data. Default settings have been expertized."
  ),
  h4(tags$span("Explore:", style = "color:orange;font-weight:bold")),
  p(
    "Identify the noise and the source of technical and biological variability 
    thanks to the", tags$b("PCA"), "."
  ),
  h4(tags$span("Filter:", style = "color:orange;font-weight:bold")),
  p(
    "Remove outliers samples and low expressed entities which introduce noise 
    in the data."
  ),
  h4(
    tags$span("Transform and normalize:",style = "color:orange;font-weight:bold")
  ),
  p(
    "Transform data to linearize and make it more Gaussian-like and normalize 
    to identify and correct technical biases and make the data comparable across 
    samples. Depending on the omics type, the pre-processing steps will be different: "
  ),
  h5(tags$span("Transcriptomics (RNAseq counts data):", style = "color:blue")),
  p(
    "", tags$b("Details on filtering step:"), ""),
  p(
    "By default, non expressed/non detected genes are removed"),
  p(
    "By default, low expressed genes are removed according to their CPM. 
    By default, genes with a cpm >= 1 in at least min(NbReplicates) samples are 
    kept. The cpm threshold can be changed."
  ),
  p(
    "NB: you can choose the other strategy, which is to remove genes according 
    to a cpm >= 1 in at least NbConditions samples. The cpm threshold can be changed."
  ),
  p(
    "", tags$b("Details on normalization:"), ""),
  p(
    "TMM (from edgeR::calcnormfactors) is the default method.
        It was found to be the best method (Dillies et al., 2013) for counts
        RNA-seq data."
  ),
  h5(
    tags$span("Proteomics and metabolomics data", style = "color:blue")
  ),
  p("", tags$b("Details on transformation:"), ""),
  p(
    "Log2 is the default method for proteomics and metabolomics data 
    transformation (Efstathiou et al, 2017). A small quantity (10^-10) is added 
    to the data before tranformation."
  ),
  p("", tags$b("Normalization:"), ""),
  p(
    "Median is the the default method for proteomics and metabolomics data 
    normalization. All samples will have the same median."
  ))


# ---- 04 diff analysis ----
diff_analysis_docs <- list(
p("Differential expression analysis is performed for each contrast. 
          There are just two options to set: the adjusted-pvalue cut-off and the |logFC| cut-off.
          The results will appear in blocks with the contrast's name and statistics (one per contrast), 
          each block offering a tab panel with several outputs:"),
p("- The graph of Pvalue's distribution: Distribution of pvalue's which has to be check to validate results. The most desirable shape is a pick of p-values at 0 following by a uniform distribution. "),
p("- The MA plot which gives the logFC across the mean of the expression/abundance"),
p("- The Volcano plot: implemented in the EnhancedVolcano R-package (Blighe K, Rana S, Lewis M (2022). EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and labeling R package version 1.12.0.))"),
p("- A dataframe with the statistical results of the differential analysis per DE entities"),
p("- A Heatmap plot: implemented in the ComplexHeatmap package (Gu Z (2022). Complex Heatmap Visualization. iMeta. doi:10.1002/imt2.43.)"),
p("- A PCA on DE entities"),
p("- Boxplot of DE : boxplot showing the expression/abundance profile of a selected DE entity across experimental factors")
)


# ---- 05 coExp analysis ----
