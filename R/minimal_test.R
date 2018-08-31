minimal_test <- function() {
    library(citrus)
    setwd("C:/Users/fgherardini/temp/standalone_citrus/data")
    tab <- read.table("Patient20_diseased_unstim.fcs.clustered.txt", header = T, sep = "\t", stringsAsFactors = F, 
        check.names = F)
    metadata.tab <- read.table("samples_tab.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)
    features.names <- c("FunctionalMarker1", "FunctionalMarker2", "popsize")
    
    metadata.tab <- metadata.tab[metadata.tab$condition == "unstim", ]
    kk <- grappolo::get_cluster_features(tab, metadata.tab, features.names, "condition", "sample", out.format = "table")
    
    
    ff <- convert_to_citrus_featureset(kk, "sample")
    run_citrus_analysis(ff, metadata.tab$label, "./", "sam")
    
    message("#########################")
    
    run_citrus_analysis(ff, metadata.tab$label, "./", "pamr")
    
}
