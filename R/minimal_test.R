minimal_test <- function() {
    tab <- read.table("Patient20_diseased_unstim.fcs.clustered.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)
    metadata.tab <- read.table("samples_tab.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)
    features.names <- c("FunctionalMarker1", "FunctionalMarker2", "popsize")

    metadata.tab <- metadata.tab[metadata.tab$condition == "unstim",]
    kk <- scaffold:::calculate_cluster_features(tab, metadata.tab, features.names, "condition", "sample")
    ff <- scaffold:::convert_to_citrus_featureset(kk, "sample")
    scaffold:::run_citrus_analysis(ff, metadata.tab$label, "./", "sam")
    
    message("#########################")
    
    scaffold:::run_citrus_analysis(ff, metadata.tab$label, "./", "pamr")

}