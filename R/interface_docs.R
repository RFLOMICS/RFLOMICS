
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
