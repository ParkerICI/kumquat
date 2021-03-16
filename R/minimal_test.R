minimal_test3 <- function() {
    library(kumquat)
    setwd("C:/Users/pfgherardini/temp/standalone_citrus")
    
    tab <- read.table("out_pooled_downsample_40000/group1.clustered.txt", header = T,
                      sep = "\t", stringsAsFactors = F, check.names = F)
    
    metadata.tab <- readRDS("metadata_with_filename.rds")
    metadata.tab$sample.id <- gsub(".clustered.txt", "", metadata.tab$filename)
    metadata.tab <- metadata.tab[!metadata.tab$timepoint.id %in% c("BL", "EOT"), ]
    metadata.tab <- metadata.tab[metadata.tab$timepoint.id %in% c("C1D1", "C1D5"), ]
    
    
    features.names <- c("popsize")
    
    features <- grappolo::melt_cluster_results(tab, features.names)
    
    #features.normalized <- multistep_normalize(features, list(timepoint.id = "C1D1"), subject.var = "subject.id")
    

    
    features.citrus <- get_cluster_features(
      features, 
      metadata.tab = metadata.tab,
      sample.col = "sample.id",
      predictors = "timepoint.id",
      endpoint.grouping = "subject.id"
    )
    
    
    endpoint <- metadata.tab[!duplicated(metadata.tab$subject.id), ]
    row.names(endpoint) <- endpoint$subject.id
    
    pam.model <- run_analysis(features.citrus, endpoint$psa.response, "./", "pamr")
    sam.model <- run_analysis(features.citrus, endpoint$psa.response, "./", "pamr")
    glmnet.model <- run_analysis(features.citrus, endpoint$psa.response, "./", "glmnet")
    
}




minimal_test <- function() {
    library(kumquat)
    setwd("C:/Users/fgherardini/temp/standalone_citrus/data")
    tab <- read.table("Patient20_diseased_unstim.fcs.clustered.txt", header = T, 
        sep = "\t", stringsAsFactors = F, check.names = F)
    metadata.tab <- read.table("samples_tab.txt", header = T, sep = "\t", stringsAsFactors = F, 
        check.names = F)
    features.names <- c("FunctionalMarker1", "FunctionalMarker2", "popsize")
    
    metadata.tab <- metadata.tab[metadata.tab$condition == "unstim", ]
    kk <- grappolo::get_cluster_features(tab, metadata.tab, features.names, "condition", 
        "sample", out.format = "table")
    
    
    ff <- convert_to_citrus_featureset(kk)
    sam.model <- run_citrus_analysis(ff, metadata.tab$label, "./", "sam")
    glmnet.model <- run_citrus_analysis(ff, metadata.tab$label, "./", "glmnet")
    
    message("#########################")
    
    pam.model <- run_citrus_analysis(ff, metadata.tab$label, "./", "pamr")
    
}

minimal_test2 <- function() {
    library(kumquat)
    setwd("C:/Users/fgherardini/temp/standalone_citrus/science")
    tab <- read.table("all.pooled.clustered.txt", header = T, 
                      sep = "\t", stringsAsFactors = F, check.names = F)
    metadata.tab <- read.table("metadata_tab_citrus.txt", header = T, sep = "\t", stringsAsFactors = F, 
                               check.names = F)
    features.names <- c("popsize")
    
    metadata.tab <- metadata.tab[metadata.tab$tissue %in% c("blood", "bone_marrow"), ]
    kk <- grappolo::get_cluster_features(tab, metadata.tab, features.names, predictors = "condition", 
                                         endpoint.grouping = "file", out.format = "table")
    
    clusters.data <- readRDS("all.pooled.clustered.all_events.rds")
    
    
    ff <- convert_to_citrus_featureset(kk)
    #sam.model <- get_model(ff, metadata.tab$tissue, "sam")
    #pam.model <- get_model(ff, metadata.tab$tissue, "pamr")
    
    sam.model <- run_analysis(ff, metadata.tab$tissue, "./", "sam")
    pam.model <- run_analysis(ff, metadata.tab$tissue, "./", "pamr")
    plot_stratifying_features(pam.model, "pamr_results", by.cluster = TRUE, all.features = TRUE)
    
    
    ######################################
    
    
    glmnet.model <- run_citrus_analysis(ff, metadata.tab$tissue, "./", "glmnet")
    
    message("#########################")
    
    pam.model <- run_citrus_analysis(ff, metadata.tab$tissue, "./", "pamr")
    
    tab$pamr.cv1se <- 0
    tab$pamr.cv1se[pam.model$differentialFeatures$cv.1se$clusters] <- 1
    
    
    tab$pamr.constrained <- 0
    tab$pamr.constrained[pam.model$differentialFeatures$cv.fdr.constrained$clusters] <- 1
    
    
    sam.res <- get_signifcant_features_matrix(sam.model$finalModel$model)
    scores <- select_score_by_cluster(sam.res,  "Score(d)", max)
    tab[scores$cluster, "sam.qvalue"] <- scores$"q-value"
    tab[scores$cluster, "sam.score"] <- scores$"Score(d)"
    
    
    write.table(tab, "all.pooled.clustered.txt", row.names = F, col.names = T, sep = "\t", quote = F)
    
    #scaffold analysis 
    
    input.files <- list.files(pattern = "*.clustered.txt$", full.names = TRUE)
    
    col.names <- c("CD45.2", "Ly6G", "IgD", "CD11c", "F480", "CD3", "NKp46", "CD23", "CD34", "CD115", 
                   "CD19", "120g8", "CD8", "Ly6C", "CD4", "CD11b", "CD27", "CD16_32", "SiglecF", "Foxp3", "B220", 
                   "CD5", "FceR1a", "TCRgd", "CCR7", "Sca1", "CD49b", "cKit", "CD150", "CD25", "TCRb", "CD43", "CD64",
                   "CD138", "CD103", "IgM", "CD44", "MHCII")
    

    landmarks.data <- vite::load_landmarks_from_dir("landmarks/", asinh.cofactor = 5, transform.data = T)
   
    vite::run_scaffold_analysis(input.files, ref.file = input.files[[1]], landmarks.data, col.names)
    
    
    
    
    
}
