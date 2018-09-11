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
    
    
    ff <- convert_to_citrus_featureset(kk)
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
    
    clusters.data <- readRDS("all.pooled.clustered.all_events.rds")
    
    
    ff <- convert_to_citrus_featureset(kk)
    sam.model <- get_model(ff, metadata.tab$tissue, "./", "sam")

    
    
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
