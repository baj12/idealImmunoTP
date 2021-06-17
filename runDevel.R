# library(ideal)
# BiocManager::install("apeglm")
# install.packages("ashr")

# library(ashr)
library(airway)
library(apeglm)
library(DESeq2) 
library(SummarizedExperiment) 
library(GenomicRanges) 
library(IRanges)
library(S4Vectors) 
library(ggplot2) 
library(heatmaply) 
library(plotly)
library(pcaExplorer) 
library(IHW) 
library(gplots) 
library(UpSetR) 
library(goseq) 
library(stringr) 
library(dplyr)
library(limma) 
library(GOstats) 
library(GO.db) 
library(AnnotationDbi) 
library(shiny)
library(shinydashboard) 
library(shinyBS) 
library(DT) 
library(rentrez) 
library(rintrojs) 
library(ggrepel) 
library(knitr)
library(rmarkdown) 
library(shinyAce) 
library(BiocParallel) 
library(grDevices) 
library(base64enc)
library(methods)
library(testthat) 
library(BiocStyle) 
library(airway) 
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(DEFormats) 
library(edgeR)

source("R/helpers.R")

# data(airway)
# dds_airway <- DESeqDataSet(airway,design= ~ cell + dex)
# dds_airway <- DESeq(dds_airway)
# res_airway <- results(dds_airway, contrast = c("dex","trt","untrt"),alpha = 0.05)
# library(org.Hs.eg.db)
# genenames_airway <- mapIds(org.Hs.eg.db,keys = rownames(dds_airway),column = "SYMBOL",keytype="ENSEMBL")
# annotation_airway <- data.frame(gene_id = rownames(dds_airway),
#                                 gene_name = genenames_airway,
#                                 row.names = rownames(dds_airway),
#                                 stringsAsFactors = FALSE)
# save(file = "testData.RData", list = c("genenames_airway", "annotation_airway", "dds_airway", "res_airway"))

load(file = "testData.RData")
ideal_ui <- NULL
ideal_server <- NULL
library(RColorBrewer)
library(topGO)
source("R/helpers.R")
source("R/table-enhancers.R")
source("R/plot_ma.R")
source("R/plot_volcano.R")
source("R/ggplotCounts.R")
# source("R/")  
source("R/res2tbl.R")
source("R/plotCoefficients.R")
source("R/genesignatures.R")
source("R/goseqTable.R")
source("R/ideal.R", local = T)
source("R/iSEE_plug.R")
dds_obj = dds_airway
res_obj = res_airway
annotation_obj = annotation_airway
dds_obj = NULL
res_obj = NULL
annotation_obj = NULL
countmatrix = NULL
  expdesign = NULL
  gene_signatures = NULL

  # app <- shinyApp(ui = ideal_ui, server = ideal_server)
  # runApp(app)
# ideal(dds_obj = dds_airway, annotation_obj = annotation_airway, res_obj = res_airway)
ideal()

# cp = load(file = "test.Rdata")
# dds_obj = ideal_env$ideal_values_20210526_152858$dds_obj 
# res_obj = ideal_env$ideal_values_20210526_152858$res_obj
# annotation_obj = ideal_env$ideal_values_20210526_152858$annotation_obj
# countmatrix = ideal_env$ideal_values_20210526_152858$countmatrix
# expdesign = ideal_env$ideal_values_20210526_152858$expdesign
# gene_signatures = ideal_env$ideal_values_20210526_152858$gene_signatures
# 
# 
# options(browser = "/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome")
# 
# source("R/ideal.R", local = T);
# ideal(dds_obj = dds_obj, 
#       res_obj = res_obj, 
#       annotation_obj = annotation_obj,
#       countmatrix = countmatrix,
#       expdesign = expdesign, 
#       gene_signatures = gene_signatures
#     )
# 
