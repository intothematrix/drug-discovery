load.libraries<-function(){
  #if (!require("dplyr")) install.packages("dplyr")
  #if (!require("ggplot2")) install.packages("ggplot2")
  #if (!require("viridis")) install.packages("viridis")
  #if (!require("tibble")) install.packages("tibble")
  #if (!require("gridExtra")) install.packages("gridExtra")
  #if (!require("stringr")) install.packages("stringr")
  #if (!require("depmap")) install.packages("depmap")
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #  +   install.packages("BiocManager")
  #BiocManager::install("ExperimentHub")
  library("ggridges")
  library("dplyr")
  library("ggplot2")
  library("viridis")
  library("tibble")
  library("gridExtra")
  library("stringr")
  #library("depmap")
  library("ExperimentHub")
  library("tidyr")
  library("lattice")
  library("forcats")
  library("hrbrthemes")
  library("cluster")
  library("fpc")
  library("GGally")
  library("data.table")
}

load.data<-function(){
  #load depmap
  eh <- ExperimentHub()
  query(eh, "depmap")
  
  # get rnai and copynumber
  rnai <- eh[["EH3080"]]
  crispr <- eh[["EH3081"]]
  copyNumber <- eh[["EH3082"]]
  TPM <- eh[["EH3084"]]
  drug_sensitivity <- eh[["EH3087"]]
  mutationCalls <- eh[["EH3085"]]
  metadata <- eh[["EH3086"]]
  
  #get a list of genes from rnai
  rnai %>% dplyr::select(gene_name) -> gene.list
  gene.list<-unique(gene.list)
  rnai %>% dplyr::select(cell_line) -> cl.list
  cl.list<-unique(cl.list)
  gl.list<-gene.list$gene_name
  cl.list<-cl.list$cell_line
  rm(gene.list)
  
  #get rnai matrix
  rnai %>% dplyr::select(cell_line, gene_name, dependency) %>%
    dplyr::filter(cell_line %in% cl.list) %>%
    dplyr::filter(gene_name %in% gl.list) -> df.rnai
  rm(rnai)
  
  #get copynumber matrix
  copyNumber %>% dplyr::select(gene_name,cell_line,log_copy_number) %>%
    dplyr::filter(cell_line %in% cl.list) %>%
    dplyr::filter(gene_name %in% gl.list) -> df.cn
  rm(copyNumber)
  
  #join tables
  dplyr::inner_join(df.rnai,df.cn,by=c("gene_name"="gene_name","cell_line" = "cell_line")) %>%
    dplyr::filter(!is.na(dependency)) %>%
    dplyr::filter(!is.na(log_copy_number))->df
  rm(df.rnai)
  rm(df.cn)
  return(df)
}

build.features<-function(){
  #basic features: min, max, median, range, mean, median-mean
  #corr features all cell_lines: pearson
    #rnai vs copyNumber
    #rnai vs TPM
    #rnai vs crispr
    #rnai vs drug_sensitivity
    #copyNumber vs TPM
    #copyNumber vs crispr
    #copyNumber vs drug_sensitivity
    #TPM vs crispr
    #TPM vs drug_sensitivity
    #crispr vs drug_sensitivity
  #corr features all primary disease using median: pearson
    #rnai vs copyNumber
    #rnai vs TPM
    #rnai vs crispr
    #rnai vs drug_sensitivity
    #copyNumber vs TPM
    #copyNumber vs crispr
    #copyNumber vs drug_sensitivity
    #TPM vs crispr
    #TPM vs drug_sensitivity
    #crispr vs drug_sensitivity
  #corr features all primary disease using mean: pearson
    #rnai vs copyNumber
    #rnai vs TPM
    #rnai vs crispr
    #rnai vs drug_sensitivity
    #copyNumber vs TPM
    #copyNumber vs crispr
    #copyNumber vs drug_sensitivity
    #TPM vs crispr
    #TPM vs drug_sensitivity
    #crispr vs drug_sensitivity
  #corr features all primary disease using max: pearson
    #rnai vs copyNumber
    #rnai vs TPM
    #rnai vs crispr
    #rnai vs drug_sensitivity
    #copyNumber vs TPM
    #copyNumber vs crispr
    #copyNumber vs drug_sensitivity
    #TPM vs crispr
    #TPM vs drug_sensitivity
    #crispr vs drug_sensitivity
  #ranking data: separate H M L using 1/3 and countM, countL, countH
}
