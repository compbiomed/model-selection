#########################
# Author: Yue (Jason) Zhao
# Github: https://github.com/jasonzhao0307




############### Load packages ##############

require(caret)
require(mice)
require(mlbench)
require(DESeq2)
require(ggplot2)
require(plotROC)
require(caretEnsemble)
require(pROC)
library(caTools)
library(gmodels)
require(glmnet)
require(ROCR)
require(e1071)

############### functions define ##############

# get threshold of prediction probability for calculating best sensitivit + specificity
Get_Classification_Threshold <- function(predict, response) {
    perf <- ROCR::performance(ROCR::prediction(predict, response), "sens", "spec")
    df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@y.values[[1]], spec = perf@x.values[[1]])
    df[which.max(df$sens + df$spec), "cut"]
}



### Sensitivity Maximize threshold at 0.9
Get_Classification_Threshold_Sen <- function(predict, response, sen_min = 0.9, spe_min = 0.5) {
    perf <- ROCR::performance(ROCR::prediction(predict, response), "sens", "spec")
    df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@y.values[[1]], spec = perf@x.values[[1]])
    df <- df[-1,]
    df_tmp <- df[which(df$sens >= sen_min),]
    df_tmp <- df_tmp[which(df_tmp$spec >= spe_min),]
    if (nrow(df_tmp) == 0){
        print("Can't meet both sensitivity and specificity target! Will use max(sen+spe) instead")
        cut_off_report <- df[which.max(df$sens + df$spec), "cut"]
    }else{
        cut_off_report <- df_tmp[which.max(df_tmp$sens + df_tmp$spec), "cut"]
    }
    cut_off_report

}



# input a vector of predicted probability and true label, output sensitivity and specificity at best cut-off.
Get_Sensitivity_And_Specificity <- function(predict, response, provide.cutoff = FALSE,
                                            cutoff.value = 0.5,
                                            prevalence = 0.1){
  if (provide.cutoff == FALSE){
    prob.cutoff <- Get_Classification_Threshold(predict, response)
  } else{
    prob.cutoff <- cutoff.value
  }

  pred <- ifelse(predict >= prob.cutoff, 1, 0)
  # Must specify which class is positive !!!
  # also, in this case, we need to use character.
  # provide prevalence!!
  cm = confusionMatrix(as.factor(pred), as.factor(response), positive = "1", prevalence = prevalence)
  sensi <- cm$byClass[1]
  speci <- cm$byClass[2]
  ppv <- cm$byClass[3]
  npv <- cm$byClass[4]

  return(list(sens=sensi,spec=speci,ppv=ppv,npv=npv))
}


# input a vector of predicted probability and true label, output sensitivity and specificity at best cut-off.
Get_Sensitivity_And_Specificity_Sen <- function(predict, response, provide.cutoff = FALSE,
                                            cutoff.value = 0.5,
                                            prevalence = 0.1,
                                            sen_min = sen_min,
                                            spe_min = spe_min){
  if (provide.cutoff == FALSE){
    prob.cutoff <- Get_Classification_Threshold_Sen(predict, response, sen_min, spe_min)
  } else{
    prob.cutoff <- cutoff.value
  }

  pred <- as.character(ifelse(predict >= prob.cutoff, 1, 0))
  # Must specify which class is positive !!!
  # also, in this case, we need to use character.
  # provide prevalence!!
  cm = confusionMatrix(as.factor(pred), as.factor(response), positive = "1", prevalence = prevalence)
  sensi <- cm$byClass[1]
  speci <- cm$byClass[2]
  ppv <- cm$byClass[3]
  npv <- cm$byClass[4]

  return(list(sens=sensi,spec=speci,ppv=ppv,npv=npv))
}

MaxFilter <- function(df, max.value = 10){
  df.filtered <- df[which(apply(df,1,max) >= max.value),]
  return(df.filtered)
}


deseq2_norm_rle <- function(dataFrame){
# RLE normalization: relative log expression
scaling.dataFrame <- estimateSizeFactorsForMatrix(dataFrame)
dataFrame.scaled <- dataFrame
for(i in 1:ncol(dataFrame)){
  dataFrame.scaled[,i] <- dataFrame[,i]/scaling.dataFrame[i]
}

return(dataFrame.scaled)
}



### centering and scaling data
normData <- function(training, testing){
  preProcValues <- preProcess(training, method = c("center", "scale"))
  trainTransformed <- predict(preProcValues, training)
  testTransformed <- predict(preProcValues, testing)
  return(list(trainTransformed = trainTransformed,
              testTransformed = testTransformed))
}



# prepare the dataset: user needs to provide which class in y is positive
prepareData <- function(df.training, df.testing, y.train, y.test, positive){
  if (length(unique(y.train)) == 2  & sum(unique(y.train) %in% c("P","N")) == 2){
    training <- cbind(df.training, y.train)
    testing <- cbind(df.testing, y.test)
  }else{
    y.train.tmp <- as.character(y.train)
    y.test.tmp <- as.character(y.test)
    y.train <- rep("N", length(y.train))
    y.test <- rep("N", length(y.test))
    y.train[y.train.tmp == as.character(positive)] <- "P"
    y.test[y.test.tmp == as.character(positive)] <- "P"
    training <- cbind(df.training, y.train)
    testing <- cbind(df.testing, y.test)
  }

  #training$y.train <- as.factor(training$y.train)
  #testing$y.test <- as.factor(testing$y.test)
  training$y.train <- as.character(training$y.train)
  testing$y.test <- as.character(testing$y.test)
  colnames(training)[ncol(training)] <- "Class"
  colnames(testing)[ncol(testing)] <- "Class"

  return(list(training = training, testing = testing))
}


### Final model selection function
modelSelection <- function(df.training,
                           model.names,
                           num.cv.fold = 10,
                           num.cv.repeat = 10,
                           num.param.tune = 12,
                           seed.use = 100,
                           show.some.models = FALSE){
# const define
model.names.example <- c("lda", "naive_bayes", "gbm", "glmnet", "ranger", "svmLinear",
                 "svmRadial", "xgbLinear", "xgbTree")

if (show.some.models == TRUE){
  print(paste0("Here are some models you could try: ", paste(model.names.example,
                                                             collapse = ",")))
}
# Model training general parameters setting
fitControl <- trainControl(method = "repeatedcv",
                           number = num.cv.fold,
                           repeats = num.cv.repeat,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           savePredictions = TRUE,
                           ## parallel in windows:
                           allowParallel = TRUE,
                           ## Evaluate performance using
                           ## the following function
                           summaryFunction = twoClassSummary)

# train multiple models
list.model.fit <- list()
for (i in 1:length(model.names)){
  print(paste0("Now training using: ", model.names[i]))
 list.model.fit[[i]] <- caretModelTraining(df.training,
                                           model.names[i],
                                           fitControl,
                                           num.param.tune,
                                           seed.use)
}
names(list.model.fit) <- model.names

# add ensemble model



tmp.list <- list.model.fit
class(tmp.list) <- "caretList"

### correlation between models
roc.cor <- modelCor(resamples(tmp.list))

set.seed(seed.use)
greedy_ensemble <- caretEnsemble(
  tmp.list,
  metric="ROC",
  trControl=fitControl)


list.model.fit[[(length(list.model.fit) + 1)]] <- greedy_ensemble$ens_model
names(list.model.fit)[length(list.model.fit)] <- "Ensemble"



# compare models using their best params
resamps <- resamples(list.model.fit)
#print(summary(resamps))
model.summary <- (summary(resamps))$statistics
model.summary <- lapply(model.summary, function(x) {x[,3:4]})
df.model.summary <- as.data.frame(model.summary[[1]])
for (i in 2:length(model.summary)){
  df.model.summary <- cbind(df.model.summary, model.summary[[i]])
}
colnames(df.model.summary) <- c("AUC-median", "AUC-mean",
                                "Sens-median", "Sens-mean",
                                "Spec-median", "Spec-mean")


# Visualize the model comparison in AUC, sens and spec
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
plot.handle <- bwplot(resamps, layout = c(3, 1))




print("Finished.")

return(list(list.model.fit = list.model.fit,
            df.model.summary = df.model.summary,
            comparison.plot = plot.handle,
            model.ensemble = greedy_ensemble,
            roc.cor = roc.cor))
}



### Final model selection function, maximizing sensitivity
modelSelectionSen <- function(df.training,
                           model.names,
                           num.cv.fold = 10,
                           num.cv.repeat = 10,
                           num.param.tune = 12,
                           seed.use = 100,
                           show.some.models = FALSE){
# const define
model.names.example <- c("lda", "naive_bayes", "gbm", "glmnet", "ranger", "svmLinear",
                 "svmRadial", "xgbLinear", "xgbTree")

if (show.some.models == TRUE){
  print(paste0("Here are some models you could try: ", paste(model.names.example,
                                                             collapse = ",")))
}
# Model training general parameters setting
fitControl <- trainControl(method = "repeatedcv",
                           number = num.cv.fold,
                           repeats = num.cv.repeat,
                           ## Estimate class probabilities
                           classProbs = TRUE,
                           savePredictions = TRUE,
                           ## parallel in windows:
                           allowParallel = TRUE,
                           ## Evaluate performance using
                           ## the following function
                           summaryFunction = twoClassSummary)

# train multiple models
list.model.fit <- list()
for (i in 1:length(model.names)){
  print(paste0("Now training using: ", model.names[i]))
 list.model.fit[[i]] <- caretModelTrainingSen(df.training,
                                           model.names[i],
                                           fitControl,
                                           num.param.tune,
                                           seed.use)
}
names(list.model.fit) <- model.names

# add ensemble model



tmp.list <- list.model.fit
class(tmp.list) <- "caretList"

### correlation between models
roc.cor <- modelCor(resamples(tmp.list))

set.seed(seed.use)
greedy_ensemble <- caretEnsemble(
  tmp.list,
  metric="Sens",
  trControl=fitControl)


list.model.fit[[(length(list.model.fit) + 1)]] <- greedy_ensemble$ens_model
names(list.model.fit)[length(list.model.fit)] <- "Ensemble"



# compare models using their best params
resamps <- resamples(list.model.fit)
#print(summary(resamps))
model.summary <- (summary(resamps))$statistics
model.summary <- lapply(model.summary, function(x) {x[,3:4]})
df.model.summary <- as.data.frame(model.summary[[1]])
for (i in 2:length(model.summary)){
  df.model.summary <- cbind(df.model.summary, model.summary[[i]])
}
colnames(df.model.summary) <- c("AUC-median", "AUC-mean",
                                "Sens-median", "Sens-mean",
                                "Spec-median", "Spec-mean")


# Visualize the model comparison in AUC, sens and spec
theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
plot.handle <- bwplot(resamps, layout = c(3, 1))




print("Finished.")

return(list(list.model.fit = list.model.fit,
            df.model.summary = df.model.summary,
            comparison.plot = plot.handle,
            model.ensemble = greedy_ensemble,
            roc.cor = roc.cor))
}



### Individual model training
# train classification model, with the training matrix being:
# dim = m*n, m = #sample, n = (#features + 1)
# The last column is the label factor with colname: "Class"
caretModelTraining <- function(df.training, model.name, fitControl, num.param.tune, seed.use){
  set.seed(seed.use)
  if (model.name %in% c("glm", "glmnet")){
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "ROC")
  } else{
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "ROC",
                 verbose = FALSE)
  }

  return(modelFit)
}




### Individual model training
# train classification model, with the training matrix being:
# dim = m*n, m = #sample, n = (#features + 1)
# The last column is the label factor with colname: "Class"
caretModelTrainingSen <- function(df.training, model.name, fitControl, num.param.tune, seed.use){
  set.seed(seed.use)
  if (model.name %in% c("glm", "glmnet")){
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "Sens",
                 maximize = TRUE)
  } else{
      modelFit <- train(Class ~ ., data = df.training,
                 method = model.name,
                 trControl = fitControl,
                 tuneLength = num.param.tune,
                 metric = "Sens",
                 maximize = TRUE,
                 verbose = FALSE)
  }

  return(modelFit)
}



### predict and evaluate functions
predictAndEvaluation <- function(model.best,
                                 test.data,
                                 prevalence,
                                 is.ensemble = FALSE){
    # create a numerical class var
    class.num <- as.character(test.data$Class)
    class.num[class.num == "P"] <- 1
    class.num[class.num == "N"] <- 0
    class.num <- as.numeric(class.num)

    ## predict
    test.prediction <- predict(model.best, newdata = test.data)
    #print(test.prediction)
    test.prediction.prob <- predict(model.best, newdata = test.data, type = "prob")
    #return(test.prediction.prob)
    ### ROC and AUC
    if (is.ensemble == FALSE){
        g <- ggplot(test.prediction.prob, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))

        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob[,2]
    } else{
        prob.df <- data.frame(P = test.prediction.prob, N = 1-test.prediction.prob)
        g <- ggplot(prob.df, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob
        if (auc.value < 0.5){
          prob.df <- data.frame(P = 1-test.prediction.prob, N = test.prediction.prob)
          g <- ggplot(prob.df, aes(m=P, d=class.num)) +
          geom_roc(n.cuts=0) +
          coord_equal() +
          style_roc()
          auc.value <- round((calc_auc(g))$AUC, 4)
          #### get probabilities for POSITIVE
          model.prob = 1-test.prediction.prob
        }
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))


    }

    #return(list(predict = model.prob, response = class.num))
    ## get sens and spec, ppv and npv (require prevalence!)
    test.prediction.evaluation <- Get_Sensitivity_And_Specificity(predict = model.prob,
                                    response = class.num,
                                    prevalence = prevalence)

    ## combine all metrics
    test.prediction.evaluation$auc <- auc.value
    evaluation.names <- names(test.prediction.evaluation)
    test.prediction.evaluation <- sapply(test.prediction.evaluation, function(x) x)
    names(test.prediction.evaluation) <- evaluation.names

    ## return the results
    return(list(positive.prob = model.prob,
                test.prediction.evaluation = test.prediction.evaluation,
                roc = g))

}



### predict and evaluate functions
predictAndEvaluationSen <- function(model.best,
                                 test.data,
                                 prevalence,
                                 is.ensemble = FALSE,
                                 sen_min = 0.9,
                                 spe_min = 0.2){
    # create a numerical class var
    class.num <- as.character(test.data$Class)
    class.num[class.num == "P"] <- 1
    class.num[class.num == "N"] <- 0
    class.num <- as.numeric(class.num)

    ## predict
    test.prediction <- predict(model.best, newdata = test.data)
    #print(test.prediction)
    test.prediction.prob <- predict(model.best, newdata = test.data, type = "prob")
    #return(test.prediction.prob)
    ### ROC and AUC
    if (is.ensemble == FALSE){
        g <- ggplot(test.prediction.prob, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))

        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob[,2]
    } else{
        prob.df <- data.frame(P = test.prediction.prob, N = 1-test.prediction.prob)
        g <- ggplot(prob.df, aes(m=P, d=class.num)) +
        geom_roc(n.cuts=0) +
        coord_equal() +
        style_roc()
        auc.value <- round((calc_auc(g))$AUC, 4)
        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob
        if (auc.value < 0.5){
          prob.df <- data.frame(P = 1-test.prediction.prob, N = test.prediction.prob)
          g <- ggplot(prob.df, aes(m=P, d=class.num)) +
          geom_roc(n.cuts=0) +
          coord_equal() +
          style_roc()
          auc.value <- round((calc_auc(g))$AUC, 4)
          #### get probabilities for POSITIVE
          model.prob = 1-test.prediction.prob
        }
        g <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", auc.value))


    }

    #return(list(predict = model.prob, response = class.num))
    ## get sens and spec, ppv and npv (require prevalence!)
    test.prediction.evaluation <- Get_Sensitivity_And_Specificity_Sen(predict = model.prob,
                                    response = class.num,
                                    prevalence = prevalence,
                                    sen_min = sen_min,
                                    spe_min = spe_min)

    ## combine all metrics
    test.prediction.evaluation$auc <- auc.value
    evaluation.names <- names(test.prediction.evaluation)
    test.prediction.evaluation <- sapply(test.prediction.evaluation, function(x) x)
    names(test.prediction.evaluation) <- evaluation.names

    ## return the results
    return(list(positive.prob = model.prob,
                test.prediction.evaluation = test.prediction.evaluation,
                roc = g))

}


predictAndEvaluationSenBootstrap <- function(model_use,
                                             data_test,
                                             sen_min,
                                             spe_min,
                                             seed_use = 1,
                                             n_boot = 50,
                                             prevalence = 0.1
                                             ){
set.seed(seed_use)
sen.vec <- c()
spe.vec <- c()
ppv.vec <- c()
npv.vec <- c()
auc.vec <- c()

for (i in 1:n_boot){
    index <- sample(1:nrow(data_test),size = nrow(data_test),replace = TRUE)
    print(paste0("Round: ", i))
    result_tmp <- predictAndEvaluationSen(model.best = model_use,
                                    test.data = data_test[index,],
                                    prevalence = prevalence,
                                    is.ensemble = FALSE,sen_min = sen_min,spe_min = spe_min)

    sen.vec <- c(sen.vec, as.numeric(result_tmp$test.prediction.evaluation[1]))
    spe.vec <- c(spe.vec, as.numeric(result_tmp$test.prediction.evaluation[2]))
    ppv.vec <- c(ppv.vec, as.numeric(result_tmp$test.prediction.evaluation[3]))
    npv.vec <- c(npv.vec, as.numeric(result_tmp$test.prediction.evaluation[4]))
    auc.vec <- c(auc.vec, as.numeric(result_tmp$test.prediction.evaluation[5]))

}

print("Finished!")

test_set_prediction_evaluation <- NULL
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(sen.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(spe.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(ppv.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(npv.vec),4)))
test_set_prediction_evaluation <- rbind(test_set_prediction_evaluation, suppressWarnings(round(ci(auc.vec),4)))

rownames(test_set_prediction_evaluation) <- c("Sensitivity",
                                              "Specificity",
                                              "PPV",
                                              "NPV",
                                              "AUC")


return(test_set_prediction_evaluation)

}

### predict and evaluate functions
predictAndReturnProb <- function(model.best,
                                 test.data,
                                 prevalence,
                                 is.ensemble = FALSE){
    # create a numerical class var
    class.num <- as.character(test.data$Class)
    class.num[class.num == "P"] <- 1
    class.num[class.num == "N"] <- 0
    class.num <- as.numeric(class.num)

    ## predict
    test.prediction <- predict(model.best, newdata = test.data)
    #print(test.prediction)
    test.prediction.prob <- predict(model.best, newdata = test.data, type = "prob")
    #return(test.prediction.prob)
    ### ROC and AUC
    if (is.ensemble == FALSE){
        #### get probabilities for POSITIVE
        model.prob = test.prediction.prob[,2]
    } else{
        model.prob = test.prediction.prob
    }

    ## return the results
    return(list(positive.prob = model.prob))
}








###### Biomarker section #####

# https://github.com/jasonzhao0307/FeatureSelection
require(FeatureSelection)


# get signature by using Lasso
getSignatureFromMultipleGlmnet <- function(dataFrame,
                                           targetVec,
                                           nfolds = 10,
                                           logisticRegression = FALSE,
                                           nRun=100,
                                           alpha = 1){
x <- t(as.matrix(dataFrame))
featureDict <- list()
featureNum <- c()
weights.vec <- rep(0,nrow(dataFrame))
for (i in seq(1,nRun,1)) {
  if (logisticRegression == FALSE){
    fit2 <- cv.glmnet(x, targetVec, alpha = alpha, nfolds = nfolds)
  } else{
    targetVec <- as.factor(targetVec)
    fit2 <- cv.glmnet(x, targetVec, alpha = alpha, family = "binomial", type.measure = "auc", nfolds = nfolds)
  }
  weights.mat <- as.matrix(coef(fit2, s = "lambda.min"))
  weights.vec <- (weights.vec + weights.mat[-1,1])

  tmp_vec <- as.vector((coef(fit2, s="lambda.min") != 0))
  featureFromGlmnet <- rownames(dataFrame)[tmp_vec]
  featureNum <- c(featureNum, length(featureFromGlmnet))
  for (k in seq(1,length(featureFromGlmnet),1)){
    gene <- featureFromGlmnet[k]
    if (gene %in% names(featureDict)){
      featureDict[[gene]] <- featureDict[[gene]] + 1
    }
    else{
      if (is.na(gene) == FALSE){
        featureDict[[gene]] <- 1
      }
    }
  }
}
#print(featureDict)
featureSelectionComplete <- names(featureDict)
numFloor <- floor(mean(featureNum))
featureDictInverse <- list()
for (i in seq(1,length(featureDict),1)){
  numTmp <- featureDict[[i]]
  #print(numTmp)
  numTmpChr <- as.character(numTmp)
  if (numTmp %in% names(featureDictInverse)){
    featureDictInverse[[numTmpChr]] <- c(featureDictInverse[[numTmpChr]], names(featureDict)[i])
  }
  else {
    featureDictInverse[[numTmpChr]] <- c(names(featureDict)[i])
  }
}
numIndex <- sort(as.numeric(names(featureDictInverse)), decreasing = TRUE)
featureSelectionFloor <- c()
for (i in seq(1,length(numIndex),1)){
  numTmp <- numIndex[i]
  numTmpChr <- as.character(numTmp)
  featureSelectionFloor <- c(featureSelectionFloor, featureDictInverse[[numTmpChr]])
  if (length(featureSelectionFloor) > numFloor) {
    break
  }
}
return.list <- list()
return.list[["feature"]] <- featureSelectionFloor
return.list[["counts"]] <- featureDict
return.list[["counts.inverse"]] <- featureDictInverse
return.list[["weights"]] <- weights.vec
# Here we get the features:
return(return.list)
}



# Feature selection by applying Feature_Selection_Wrapper to get a smaller set.
# Also, run multiple times with different seeds.
Feature_Selection_Wrapper_Multiple_Seeds_One_Round <- function(df.train, label.num.vec, seed.vec = 1:100){
  signature.list <- list()
  for(i in 1:length(seed.vec)){
    print(paste0("Start Round ", i))
    fs.round.1 <- Feature_Selection_Wrapper(df.train, label.num.vec, seed = seed.vec[i])
    signature.list[[i]] <- fs.round.1
  }
  return(signature.list)
}



# use rf, lasso and xgboost for fs
Feature_Selection_Wrapper <- function(df.train, label.num.vec, seed = 10){
  set.seed(seed)
  X_train <- as.data.frame(t(as.matrix(df.train)))
  y_train <- label.num.vec
  params_glmnet = list(alpha = 1, family = 'binomial', nfolds = 5, parallel = TRUE)

  params_xgboost = list( params = list("objective" = "binary:logistic", "bst:eta" = 0.001,
                                     "subsample" = 0.75, "max_depth" = 5, "colsample_bytree" = 0.75,
                                     "nthread" = 6),nrounds = 1000, print.every.n = 250, maximize = FALSE)
  params_ranger = list(dependent.variable.name = 'y', probability = FALSE, num.trees = 1000, verbose = TRUE,
                     mtry = 5, min.node.size = 10, num.threads = 6, classification = TRUE,
                     importance = 'permutation')
  params_features = list(keep_number_feat = NULL, union = TRUE)
  feat = wrapper_feat_select(X = X_train, y = y_train, params_glmnet = params_glmnet,
                           params_xgboost = params_xgboost, params_ranger = params_ranger,
                           xgb_sort = 'Gain', CV_folds = 5, stratified_regr = FALSE,
                           scale_coefs_glmnet = FALSE, cores_glmnet = 5,
                           params_features = params_features, verbose = TRUE)
  fs.union <- feat$union_feat
  fs.union$feature <- as.character(fs.union$feature)
  dfList.loocv <- list()
  loocv.importance.vec <- seq(0.5,0.8,0.05)
  for (i in 1:length(loocv.importance.vec)){
    fs.union.frequent.tmp <- fs.union$feature[fs.union$Frequency >= 2 &
                                                fs.union$importance >= loocv.importance.vec[i]]
    # if fs.union.frequent.tmp is less than two genes:
    if (length(fs.union.frequent.tmp) <= 2){
      fs.union.frequent.tmp <- fs.union$feature[fs.union$Frequency >= 2]
    }


    dfList.loocv[[i]] <- df.train[which(rownames(df.train) %in% fs.union.frequent.tmp), ]
  }
  loocv.auc <- LOOAUC_simple_multiple_noplot(dfList.loocv, label.num.vec)
  #plot(loocv.auc, x = loocv.importance.vec, main = "Feature Importance Cut-off LOOCV", xlab = "Ensemble model importance", ylab = "AUC")
  importance.threshold <- max(loocv.importance.vec[which(loocv.auc == max(loocv.auc))])

  fs.union.frequent <- fs.union$feature[fs.union$Frequency >= 2 & fs.union$importance >= importance.threshold]
  return(fs.union.frequent)
}







# Bootstrap LOOCV with logistic regression
Bootstrap_LOOCV_LR_AUC <- function(df, target.vec, nboot){
  output.auc.vec <- c()
	output.other.df <- NULL
  for (i in 1:nboot){
    index.boot <- sample(1:ncol(df), ncol(df), replace = T)
    df.tmp <- df[,index.boot]
    loo.output.list <- LOOAUC_simple_multiple_noplot_one_df(df.tmp, target.vec[index.boot])
    output.auc.vec[i] <- loo.output.list[[1]]
		output.other.df <- rbind(output.other.df, loo.output.list[[2]])
  }

	output.list <- list()
	output.list[[1]] <- output.auc.vec
	output.list[[2]] <- as.data.frame(output.other.df)
	names(output.list) <- c("auc", "other")
	return(output.list)
}




# Use logistic regression and Bootstrap LOOCV for signature evaluation.
SignatureQuantitative <- function(df.input,
                                  target.vec.num,
                                  signature.list,
                                  signature.name.vec,
                                  num.boot = 100){

  df.list <- list()
  for (i in 1:length(signature.list)){
    df.list[[i]] <- df.input[signature.list[[i]],]
  }

  auc.result <- list()
  auc.result.ci <- list()
  sensitivity.ci <- list()
  specificity.ci <- list()
  for(i in 1:length(df.list)){
    boot.output.list <- Bootstrap_LOOCV_LR_AUC(df.list[[i]],
                                               target.vec.num,
                                               nboot = num.boot)
    #auc
    auc.result[[i]] <- boot.output.list[[1]]
    auc.result.ci[[i]] <- ci(auc.result[[i]])
    names(auc.result)[i] <- signature.name.vec[i]
    names(auc.result.ci)[i] <- signature.name.vec[i]
    # sensitivity
    sensitivity.ci[[i]] <- ci(boot.output.list[[2]]$Sensitivity)
    names(sensitivity.ci)[i] <- signature.name.vec[i]
    # specificity
    specificity.ci[[i]] <- ci(boot.output.list[[2]]$Specificity)
    names(specificity.ci)[i] <- signature.name.vec[i]
  }

pdf("boxplot.pdf", width = 13, height = 8)
boxplot(auc.result, main = "Quantitative evaluation of signatures: AUC")
dev.off()

# output df instead of list
df.auc.ci <- data.frame(matrix(unlist(auc.result.ci),
                                       nrow=length(auc.result.ci), byrow=T))
colnames(df.auc.ci) <- names(auc.result.ci[[1]])
rownames(df.auc.ci) <- signature.name.vec

df.sensitivity.ci <- data.frame(matrix(unlist(sensitivity.ci),
                                       nrow=length(sensitivity.ci), byrow=T))
colnames(df.sensitivity.ci) <- names(sensitivity.ci[[1]])
rownames(df.sensitivity.ci) <- signature.name.vec

df.specificity.ci <- data.frame(matrix(unlist(specificity.ci),
                                       nrow=length(specificity.ci), byrow=T))
colnames(df.specificity.ci) <- names(specificity.ci[[1]])
rownames(df.specificity.ci) <- signature.name.vec

return(list(df.auc.ci = df.auc.ci,
            df.sensitivity.ci = df.sensitivity.ci,
            df.specificity.ci = df.specificity.ci))
}



# Get Average_Top_Features from a list
Get_Average_Top_Features <- function(signature.list){
      # mean signature length
    num.single.best.feature <- ceiling(mean(sapply(signature.list, length)))

    # get the count for each gene
    signature.count <- sort(table(unlist(signature.list)), decreasing = TRUE)

    # get the top N genes.
    single.best.feature <- names(signature.count)[1:num.single.best.feature]
    return(single.best.feature)

}




# filter
rowsumFilter <- function(df.input, row.sum.min = 20){
  df.input <- df.input[rowSums(df.input) > row.sum.min,]
  return(df.input)
}



# normalization
deseq2_norm_rle <- function(dataFrame){
# RLE normalization: relative log expression
scaling.dataFrame <- estimateSizeFactorsForMatrix(dataFrame)
dataFrame.scaled <- dataFrame
for(i in 1:ncol(dataFrame)){
  dataFrame.scaled[,i] <- dataFrame[,i]/scaling.dataFrame[i]
}
return(dataFrame.scaled)
}



# get signature by using Lasso
getSignatureFromMultipleGlmnet <- function(dataFrame,
                                           targetVec,
                                           nfolds = 10,
                                           logisticRegression = FALSE,
                                           nRun=100,
                                           alpha = 1){
x <- t(as.matrix(dataFrame))
featureDict <- list()
featureNum <- c()
weights.vec <- rep(0,nrow(dataFrame))
for (i in seq(1,nRun,1)) {
  if (logisticRegression == FALSE){
    fit2 <- cv.glmnet(x, targetVec, alpha = alpha, nfolds = nfolds)
  } else{
    targetVec <- as.factor(targetVec)
    fit2 <- cv.glmnet(x, targetVec, alpha = alpha, family = "binomial", type.measure = "auc", nfolds = nfolds)
  }
  weights.mat <- as.matrix(coef(fit2, s = "lambda.min"))
  weights.vec <- (weights.vec + weights.mat[-1,1])

  tmp_vec <- as.vector((coef(fit2, s="lambda.min") != 0))
  featureFromGlmnet <- rownames(dataFrame)[tmp_vec]
  featureNum <- c(featureNum, length(featureFromGlmnet))
  for (k in seq(1,length(featureFromGlmnet),1)){
    gene <- featureFromGlmnet[k]
    if (gene %in% names(featureDict)){
      featureDict[[gene]] <- featureDict[[gene]] + 1
    }
    else{
      if (is.na(gene) == FALSE){
        featureDict[[gene]] <- 1
      }
    }
  }
}
#print(featureDict)
featureSelectionComplete <- names(featureDict)
numFloor <- floor(mean(featureNum))
featureDictInverse <- list()
for (i in seq(1,length(featureDict),1)){
  numTmp <- featureDict[[i]]
  #print(numTmp)
  numTmpChr <- as.character(numTmp)
  if (numTmp %in% names(featureDictInverse)){
    featureDictInverse[[numTmpChr]] <- c(featureDictInverse[[numTmpChr]], names(featureDict)[i])
  }
  else {
    featureDictInverse[[numTmpChr]] <- c(names(featureDict)[i])
  }
}
numIndex <- sort(as.numeric(names(featureDictInverse)), decreasing = TRUE)
featureSelectionFloor <- c()
for (i in seq(1,length(numIndex),1)){
  numTmp <- numIndex[i]
  numTmpChr <- as.character(numTmp)
  featureSelectionFloor <- c(featureSelectionFloor, featureDictInverse[[numTmpChr]])
  if (length(featureSelectionFloor) > numFloor) {
    break
  }
}
return.list <- list()
return.list[["feature"]] <- featureSelectionFloor
return.list[["counts"]] <- featureDict
return.list[["counts.inverse"]] <- featureDictInverse
return.list[["weights"]] <- weights.vec
# Here we get the features:
return(return.list)
}






# LOOCV with logistic regression
LOOAUC_simple_multiple_noplot_one_df <- function(df, targetVec){
auc.vec <- c()
	nSample <- ncol(df)
	testPredictionClassVec <- c()
  testPredictionProbVec <- c()
	for (j in 1:nSample){
		train = t(as.matrix(df[,-j]))
		test = t(as.matrix(df[,j]))
  	 	 fit <- glmnet(train, targetVec[-j], family = "binomial")
			 testPredictionClassVec[j] <- predict(fit,type="class", newx = test, s = 0)
       testPredictionProbVec[j] <- predict(fit,type="response", newx = test, s = 0)
	}
	loo.pred = prediction(testPredictionProbVec, targetVec)
	loo.perf = performance(loo.pred,"tpr","fpr")
	auc <- performance(loo.pred,"auc")
	auc <- unlist(slot(auc, "y.values"))
	aucRound <- round(auc,3)
	auc.vec <- c(auc.vec, aucRound)
	# for other metric
	testPredictionClassVec <- as.numeric(testPredictionClassVec)
	cm = confusionMatrix(as.factor(testPredictionClassVec), as.factor(targetVec))
 	output.list <- list()
	output.list[[1]] <- auc.vec
	output.list[[2]] <- cm$byClass
  output.list[[3]] <- testPredictionProbVec
	names(output.list) <- c("auc", "other", "prob")
  return(output.list)
}



# Bootstrap LOOCV with logistic regression
Bootstrap_LOOCV_LR_AUC <- function(df, target.vec, nboot){
  output.auc.vec <- c()
	output.other.df <- NULL
  for (i in 1:nboot){
    index.boot <- sample(1:ncol(df), ncol(df), replace = T)
    df.tmp <- df[,index.boot]
    loo.output.list <- LOOAUC_simple_multiple_noplot_one_df(df.tmp, target.vec[index.boot])
    output.auc.vec[i] <- loo.output.list[[1]]
		output.other.df <- rbind(output.other.df, loo.output.list[[2]])
  }

	output.list <- list()
	output.list[[1]] <- output.auc.vec
	output.list[[2]] <- as.data.frame(output.other.df)
	names(output.list) <- c("auc", "other")
	return(output.list)
}




# Use logistic regression and Bootstrap LOOCV for signature evaluation.
SignatureQuantitative <- function(df.input,
                                  target.vec.num,
                                  signature.list,
                                  signature.name.vec,
                                  num.boot = 100){

  df.list <- list()
  for (i in 1:length(signature.list)){
    df.list[[i]] <- df.input[signature.list[[i]],]
  }

  auc.result <- list()
  auc.result.ci <- list()
  sensitivity.ci <- list()
  specificity.ci <- list()
  for(i in 1:length(df.list)){
    boot.output.list <- Bootstrap_LOOCV_LR_AUC(df.list[[i]],
                                               target.vec.num,
                                               nboot = num.boot)
    #auc
    auc.result[[i]] <- boot.output.list[[1]]
    auc.result.ci[[i]] <- ci(auc.result[[i]])
    names(auc.result)[i] <- signature.name.vec[i]
    names(auc.result.ci)[i] <- signature.name.vec[i]
    # sensitivity
    sensitivity.ci[[i]] <- ci(boot.output.list[[2]]$Sensitivity)
    names(sensitivity.ci)[i] <- signature.name.vec[i]
    # specificity
    specificity.ci[[i]] <- ci(boot.output.list[[2]]$Specificity)
    names(specificity.ci)[i] <- signature.name.vec[i]
  }

pdf("boxplot.pdf", width = 13, height = 8)
boxplot(auc.result, main = "Quantitative evaluation of signatures: AUC")
dev.off()

# output df instead of list
df.auc.ci <- data.frame(matrix(unlist(auc.result.ci),
                                       nrow=length(auc.result.ci), byrow=T))
colnames(df.auc.ci) <- names(auc.result.ci[[1]])
rownames(df.auc.ci) <- signature.name.vec

df.sensitivity.ci <- data.frame(matrix(unlist(sensitivity.ci),
                                       nrow=length(sensitivity.ci), byrow=T))
colnames(df.sensitivity.ci) <- names(sensitivity.ci[[1]])
rownames(df.sensitivity.ci) <- signature.name.vec

df.specificity.ci <- data.frame(matrix(unlist(specificity.ci),
                                       nrow=length(specificity.ci), byrow=T))
colnames(df.specificity.ci) <- names(specificity.ci[[1]])
rownames(df.specificity.ci) <- signature.name.vec

return(list(df.auc.ci = df.auc.ci,
            df.sensitivity.ci = df.sensitivity.ci,
            df.specificity.ci = df.specificity.ci))
}





# Fisher Exact Test for signature enrichment in DE gene list
FisherTestGeneListEnrichment <- function(signature.vec, DE.vec, all.vec){
  len.signature.DE <- length(intersect(signature.vec, DE.vec))
  len.signature.notDE <- length(signature.vec) - len.signature.DE
  len.notSignature.DE <- length(DE.vec) - len.signature.DE
  len.notSignature.notDE <- length(all.vec) - len.signature.DE - len.notSignature.DE - len.signature.notDE
  fisher.mat <- matrix(c(len.signature.DE, len.signature.notDE, len.notSignature.DE, len.notSignature.notDE),
                nrow = 2,
                dimnames = list(DE = c("DE", "Not DE"), SIGNATURE = c("Signature", "Not signature")))
  tmp <- fisher.test(fisher.mat, alternative = "greater")
  print(tmp)
  return(tmp)
}





# LOOCV with logistic regression
LOOAUC_simple_multiple_noplot_one_df <- function(df, targetVec){
auc.vec <- c()
	nSample <- ncol(df)
	testPredictionClassVec <- c()
  testPredictionProbVec <- c()
	for (j in 1:nSample){
		train = t(as.matrix(df[,-j]))
		test = t(as.matrix(df[,j]))
  	 	 fit <- glmnet(train, targetVec[-j], family = "binomial")
			 testPredictionClassVec[j] <- predict(fit,type="class", newx = test, s = 0)
       testPredictionProbVec[j] <- predict(fit,type="response", newx = test, s = 0)
	}
	loo.pred = prediction(testPredictionProbVec, targetVec)
	loo.perf = performance(loo.pred,"tpr","fpr")
	auc <- performance(loo.pred,"auc")
	auc <- unlist(slot(auc, "y.values"))
	aucRound <- round(auc,3)
	auc.vec <- c(auc.vec, aucRound)
	# for other metric
	testPredictionClassVec <- as.numeric(testPredictionClassVec)
	cm = confusionMatrix(as.factor(testPredictionClassVec), as.factor(targetVec))
 	output.list <- list()
	output.list[[1]] <- auc.vec
	output.list[[2]] <- cm$byClass
  output.list[[3]] <- testPredictionProbVec
	names(output.list) <- c("auc", "other", "prob")
  return(output.list)
}



LOOAUC_simple_multiple_noplot <- function(dfList, targetVec){
auc.vec <- c()
for (i in 1:length(dfList)){
	df <- dfList[[i]]
	nSample <- ncol(df)
	vecProbTmp <- c()
	for (j in 1:nSample){
		train = t(as.matrix(df[,-j]))
		test = t(as.matrix(df[,j]))
  	 	 fit <- glmnet(train, targetVec[-j], family = "binomial")
  	 	 testProb <- predict(fit,type="response", newx = test, s = 0)
  	 	 vecProbTmp <- c(vecProbTmp, testProb)
	}
	loo.pred = prediction(vecProbTmp, targetVec)
	loo.perf = performance(loo.pred,"tpr","fpr")
	auc <- performance(loo.pred,"auc")
	auc <- unlist(slot(auc, "y.values"))
	aucRound <- round(auc,3)
	auc.vec <- c(auc.vec, aucRound)

}
  return(auc.vec)
}
