#' @rdname citrus.buildEndpointModel
#' @name citrus.buildEndpointModel
#' 
citrus.buildModel.classification <- function(features, labels, type, regularizationThresholds, 
    ...) {
    
    addtlArgs <- list(...)
    alpha <- 1
    if ("alpha" %in% names(addtlArgs)) {
        alpha <- addtlArgs[["alpha"]]
    }
    standardize <- T
    if ("standardize" %in% names(addtlArgs)) {
        standardize <- addtlArgs[["standardize"]]
    }
    
    if (all(c("thisFoldIndex", "finalModelIndex") %in% names(addtlArgs))) {
        if ((type == "sam") && (addtlArgs[["thisFoldIndex"]] != addtlArgs[["finalModelIndex"]])) {
            return(NULL)
        }
    }
    
    if (type == "pamr") {
        pamrData <- list(x = t(features), y = labels)
        model <- pamr::pamr.train(data = pamrData, threshold = regularizationThresholds, 
            remove.zeros = F)
    } else if (type == "glmnet") {
        # NOTE THAT THIS IS BINOMIAL EXPLICITLY. DOES MULTINOMIAL WORK THE SAME, IF ONLY
        # 2 CLASSES PROVIDED?
        family <- "binomial"
        if (length(unique(labels)) > 2) {
            family <- "multinomial"
        }
        model <- glmnet::glmnet(x = features, y = labels, family = family, lambda = regularizationThresholds, 
            alpha = alpha, standardize = standardize)
    } else if (type == "sam") {
        if ("family" %in% names(addtlArgs)){
            family <- addtlArgs[["family"]] 
        } else if (!("family" %in% addtlArgs)){
            family <- "Two class unpaired"
            if (length(unique(labels)) > 2) {
                family <- "Multiclass"
            }
        }
        if ("testStatistic" %in% names(addtlArgs)){ 
            testStatistic <- addtlArgs[["testStatistic"]]
        } else {
            testStatistic <- "wilcoxon"
        }
        print(paste("The test statisitc is:", testStatistic))
        noVarianceFeatures <- apply(features, 2, var) == 0
        model <- samr::SAM(x = t(features[, !noVarianceFeatures]), y = labels, resp.type = family, 
            genenames = colnames(features[, !noVarianceFeatures]), nperms = 10000, testStatistic = testStatistic)
    } else {
        stop(paste("Type:", type, "not yet implemented"))
    }
    return(model)
}

#' @rdname citrus.thresholdCVs
#' @name citrus.thresholdCVs
#' 
citrus.thresholdCVs.quick.classification <- function(modelType, features, labels, 
    regularizationThresholds, nCVFolds = 10, ...) {
    
    errorRates <- list()
    errorRates$threshold <- regularizationThresholds
    
    if (modelType == "pamr") {
        pamrData <- list(x = t(features), y = labels)
        pamrModel <- pamr::pamr.train(data = pamrData, threshold = regularizationThresholds, 
            remove.zeros = F)
        pamrCVModel <- pamr::pamr.cv(fit = pamrModel, data = pamrData, nfold = nCVFolds)
        errorRates$cvm <- pamrCVModel$error
        
        cvmSD <- as.vector(apply(sapply(pamrCVModel$folds, function(foldIndices, 
            y, yhat) {
            apply(apply(yhat[foldIndices, ], 2, "==", y[foldIndices]), 2, sum)/length(foldIndices)
        }, y = labels, yhat = pamrCVModel$yhat), 1, sd))
        
        errorRates$cvsd <- cvmSD/sqrt(length(pamrCVModel$folds))
        errorRates$fdr <- pamr.fdr.new(pamrModel, data = pamrData, nperms = 1000)$results[, 
            "Median FDR"]
    } else if (modelType == "glmnet") {
        glmnetFamily <- "binomial"
        if (length(unique(labels)) > 2) {
            glmnetFamily <- "multinomial"
        }
        addtlArgs <- list(...)
        alpha <- 1
        if ("alpha" %in% names(addtlArgs)) {
            alpha <- addtlArgs[["alpha"]]
        }
        standardize <- T
        if ("standardize" %in% names(addtlArgs)) {
            standardize <- addtlArgs[["standardize"]]
        }
        glmnetModel <- glmnet::cv.glmnet(x = features, y = labels, family = glmnetFamily, 
            lambda = regularizationThresholds, type.measure = "class", alpha = alpha, 
            standardize = standardize)
        errorRates$cvm <- glmnetModel$cvm
        errorRates$cvsd <- glmnetModel$cvsd
    } else if (modelType == "sam") {
        warning("No thresholds for SAM. This is Normal.")
        return(NA)
    } else {
        stop(paste("CV for Model type", modelType, "not implemented"))
    }
    
    return(data.frame(errorRates))
}


foldPredict.classification <- function(index, models, features) {
    citrus.predict.classification(models[[index]], features[[index]])
}

foldScore.classification <- function(index, folds, predictions, labels) {
    return(predictions[[index]] != labels[folds[[index]]])
}

#' @rdname citrus.predict
#' @name citrus.predict
#' 
citrus.predict.classification <- function(citrus.endpointModel, newFeatures) {
    if (citrus.endpointModel$type == "glmnet") {
        predictions <- predict(citrus.endpointModel$model, newx = newFeatures, type = "class")
    } else if (citrus.endpointModel$type == "pamr") {
        predictions <- pamr::pamr.predictmany(fit = citrus.endpointModel$model, newx = t(newFeatures))$predclass
    } else {
        stop(paste("don't know how to predict for class", citrus.endpointModel$type))
    }
    rownames(predictions) <- rownames(newFeatures)
    return(predictions)
}

#' @rdname citrus.generateRegularizationThresholds
#' @name citrus.generateRegularizationThresholds
#' 
citrus.generateRegularizationThresholds.classification <- function(features, labels, 
    modelType, n = 100, ...) {
    addtlArgs <- list(...)
    alpha <- 1
    if ("alpha" %in% names(addtlArgs)) {
        alpha <- addtlArgs[["alpha"]]
    }
    standardize <- T
    if ("standardize" %in% names(addtlArgs)) {
        standardize <- addtlArgs[["standardize"]]
    }
    
    if (modelType == "pamr") {
        return(rev(pamr::pamr.train(data = list(x = t(features), y = labels), n.threshold = n)$threshold))
    }
    
    if (modelType == "glmnet") {
        if (length(unique(labels)) == 2) {
            glmfamily <- "binomial"
        } else {
            glmfamily <- "multinomial"
        }
        return(glmnet::glmnet(x = features, y = labels, family = glmfamily, alpha = alpha, 
            nlambda = n, standardize = standardize)$lambda)
    }
    
}


.calculateTypeFDRRate <- function(foldModels, foldFeatures, labels, modelType) {
    if (modelType == "pamr") {
        # calculate FDR Rates for each individual fold model
        foldFDRRates <- mcmapply(FUN = function(foldModel, foldFeatures) {
            pamr.fdr.new(foldModel$model, data = list(x = t(foldFeatures), y = foldModel$model$y), 
                nperms = 1000)$results[, "Median FDR"]
        }, foldModel = foldModels, foldFeatures = foldFeatures)
        
        # Average FDR Rates across all folds for each regularization threshold
        averageFDRRate <- apply(foldFDRRates, 1, mean)
        
        return(averageFDRRate)
    } else {
        return(NULL)
    }
}
