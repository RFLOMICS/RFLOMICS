### ============================================================================
### [instructions]
### ----------------------------------------------------------------------------

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
coexp_analysis_docs <- list(
  
  div(
    p("Analyses in this module are conducted using the 
                      coseq R-package. If you have more questions or
                      interest in this package, 
              please check the associated paper or the online vignette
              at https://bioconductor.org/packages/release/bioc/vignettes/coseq/inst/doc/coseq.html"),
    h4(tags$span("Parameters set up:", style = "color:orange")),
    p("You have first to choose between the ",tags$b("union"),
      " or ",tags$b("intersection")," of your contrasts lists 
                according to your biological question."),
    p("All the default parameters have been expertised according 
                to each omics."),
    p("It is then recommanded to do a ",tags$b("first run"),
      " with a large number of K with few iterations.
            If there is a K (Kbest different from Kmin and Kmax) for which 
            the ICL is minimum (check the first graph obtained),
            then a second run has to be done with a larger experiment: 
            the window of K can be centered around the Kbest and the number
            of technical replicates has to be increased to at least 20
            iterations. For this larger experiment, it is recommanded to send
              analysis to remote ressources (see cluster option)"),
    h4(tags$span("Successful analysis:", style = "color:orange")),
    p("For a given K, the result will be considered as successful
              when at least the half of the iteration
              have run. Co-expression analysis will be considered as 
              successful if there is at least a result for more than the
              half of K. In case of unsuccessful results, a detailed table 
              of errors will appear."),
    # h4(tags$span("Cluster option", style = "color:orange")),
    # p("If you have a cluster account, you can configure a remote 
    #   access to it
    #   (", a("see config_file", href="install_clustermq.txt"),")",
    #   "and speed up results obtention."),
  )
)

# ---- 06 annot analysis ----
annot_analysis_docs <- list(
div(
  p(
    "Analyses in this module are conducted using the clusterprofiler
            R-package.
              If you have more questions or interest in this package,
              please check the associated paper or the online vignette at
              https://yulab-smu.top/biomedical-knowledge-mining-book/index.html.
            "
  ),
  p(""),
  h4(tags$span("Parameters set up:", style = "color:orange")),
  p(
    "Choose the lists of omics features you want to run the
            enrichment for.
              Default option selects all the available lists
              (contrasts or co-expression clusters)."
  ),
  p(
    "Then choose the ontology you want to refer to for the analysis.
              Multiple choices are not allowed.
            If you select custom, you'll have to enter an annotation file
            with at least two columns :
              the names of the entity (same as the rownames of your dataset)
            and an id for an ontology term (eg: GO:0030198).
              It can also contains a column for a more explicit name
            for the term (eg: extracellular matrix organization)
              and a column for specifying the domain (eg: MF, BP or CC).
              An enrichment analysis will be ran on each specified domain."
  ),
  p(
    "You will have to specify the names of the columns after
              validating the annotation file."
  ),
  p(
    "If you choose GO, the three GO:BP, GO:MF and GO:CC will
              be analyzed.
              You can chose to only analyze one of them by selecting
              the wanted ontology domain).
              It requires to indicate an R-package for the database,
              in the form of org.*db.
              (eg: org.At.tair.db for Arabidopsis thaliana),
              and to specify which type of identifier is used in the data
              (eg: TAIR for arabidopsis)."
  ),
  p(
    "For KEGG analysis, only four identifiers are possible.
              Please check if the rownames correspond to one of them
              (eg: TAIR is also kegg identifiers)."
  ),
  p(
    "KEGG analysis uses an API to search for pathways online,
              it requires to have access to an internet connection."
  ),
  p(
    "Set the adjusted pvalue threshold.
              Only results below this threshold will be displayed."
  ),
  
)
)