minimal_test <- function() {
    library(citrus)
    setwd("C:/Users/fgherardini/temp/standalone_citrus/data")
    tab <- read.table("Patient20_diseased_unstim.fcs.clustered.txt", header = T, 
        sep = "\t", stringsAsFactors = F, check.names = F)
    metadata.tab <- read.table("samples_tab.txt", header = T, sep = "\t", stringsAsFactors = F, 
        check.names = F)
    features.names <- c("FunctionalMarker1", "FunctionalMarker2", "popsize")
    
    metadata.tab <- metadata.tab[metadata.tab$condition == "unstim", ]
    kk <- grappolo::get_cluster_features(tab, metadata.tab, features.names, "condition", 
        "sample", out.format = "table")
    
    
    ff <- convert_to_citrus_featureset(kk, "sample")
    sam.model <- run_citrus_analysis(ff, metadata.tab$label, "./", "sam")
    glmnet.model <- run_citrus_analysis(ff, metadata.tab$label, "./", "glmnet")
    
    message("#########################")
    
    pam.model <- run_citrus_analysis(ff, metadata.tab$label, "./", "pamr")
    
}

minimal_test2 <- function() {
    library(citrus)
    setwd("C:/Users/fgherardini/temp/standalone_citrus/science")
    tab <- read.table("all.pooled.clustered.txt", header = T, 
                      sep = "\t", stringsAsFactors = F, check.names = F)
    metadata.tab <- read.table("metadata_tab_citrus.txt", header = T, sep = "\t", stringsAsFactors = F, 
                               check.names = F)
    features.names <- c("popsize")
    
    metadata.tab <- metadata.tab[metadata.tab$tissue %in% c("blood", "bone_marrow"), ]
    kk <- grappolo::get_cluster_features(tab, metadata.tab, features.names, predictors = "condition", 
                                         endpoint.grouping = "file", out.format = "table")
    
    
    ff <- convert_to_citrus_featureset(kk, "file")
    sam.model <- run_citrus_analysis(ff, metadata.tab$tissue, "./", "sam")
    glmnet.model <- run_citrus_analysis(ff, metadata.tab$tissue, "./", "glmnet")
    
    message("#########################")
    
    pam.model <- run_citrus_analysis(ff, metadata.tab$tissue, "./", "pamr")
    
}
