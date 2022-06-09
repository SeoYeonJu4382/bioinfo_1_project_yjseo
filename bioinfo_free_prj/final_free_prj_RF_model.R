# 02.build LIN28A binding motif predictor 
#library 불러오기
library(readr)
library(dplyr)
library(caret)
library(randomForest)
library(ComplexHeatmap)
library(ROCR)
library("Biostrings")


# 함수
# AUC & ROC 그리는 함수
plot.roc.curve <- function(model, data, actual, title.text){
  predict_prob <- predict(model, newdata = data,type = 'prob')
  predictions <- prediction(predict_prob[,2], actual)
  
  perf <- performance(predictions, "tpr", "fpr")
  plot(perf,col="black",lty=1, lwd=2,
       main=title.text, cex.main=0.6, cex.lab=0.8,xaxs="i", yaxs="i")
  abline(0,1, col="red")
  auc <- performance(predictions,"auc")
  auc <- unlist(slot(auc, "y.values"))
  #auc <- round(auc,2)
  legend(0.4,0.4,legend=c(paste0("AUC: ",auc)),cex=0.6,bty = "n",box.col = "white")
}

# 1. data load
data = readDNAStringSet("/home/yjseo/bioinfo_1_project_yjseo/final_free_prj_data/chr1.fa_len_norm.out") #41357
seq_name = names(data)
sequence = paste(data)
df <- data.frame(seq_name, sequence)

pre_data = readDNAStringSet("/home/yjseo/bioinfo_1_project_yjseo/final_free_prj_data/chr1.fa.out") #41357


null_data = read.table("/home/yjseo/bioinfo_1_project_yjseo/final_free_prj_data/null3_seq.txt")
null_df  <- data.frame(seq_name, null_data$V1)


# 2. state 만들기 (motif 있으면 = 1, 없으면 = 0) 
# 2-1) motif : AAGNNG, AAGNGN,NUGUGN
df$state <- ifelse((c(grepl("AAG..G", df$sequence))==TRUE) |(c(grepl("AAG.G.", df$sequence))==TRUE) |(c(grepl(".UGUG.", df$sequence))==TRUE)  ,1, 0)
null_df$state <- ifelse((c(grepl("AAG..G", null_df$sequence))==TRUE) |(c(grepl("AAG.G.", null_df$sequence))==TRUE) |(c(grepl(".UGUG.", null_df$sequence))==TRUE)  ,1, 0)
df$state <- as.factor(df$state)
null_df$state <- as.factor(null_data$state)
summary(df$state) 
summary(null_df$state) 

# 2-2) motif : AAGNNG,
df$state <- ifelse((c(grepl("AAG..G", df$sequence))==TRUE)  ,1, 0)
null_df$state <- ifelse((c(grepl("AAG..G", null_df$null_data.V1))==TRUE) ,1, 0)
df$state <- as.factor(df$state)
null_df$state <- as.factor(null_df$state)
summary(df$state) # 0: 21376 1: 19981
summary(null_df$state)  # 0: 28543 1: 12814 

# index 부여하고 train data와 test data 나누기(ulms 고려해서)
set.seed(2022)
train.idx <- createDataPartition(df$state, p = 0.8)[[1]]
data.train <- df[train.idx,] 
data.test <- df[-train.idx,]


summary(data.train$state ) # 0:30180  1: 2907
summary(data.test$state ) #0: 7544  1: 726 
df$state <- as.factor(df$state)


write.csv(data.train, '/home/yjseo/bioinfo_1_project_yjseo/final_free_prj_data/train_data4.csv', row.names=FALSE)
write.csv(data.test, '/home/yjseo/bioinfo_1_project_yjseo/final_free_prj_data/test_data4.csv', row.names=FALSE)
write.csv(null_df, '/home/yjseo/bioinfo_1_project_yjseo/final_free_prj_data/null_data4.csv', row.names=FALSE)



#데이터 준비
data.train_extract <- data.train %>% dplyr::select(-seq_name)
data.test_extract <- data.test %>% dplyr::select(-seq_name) 
null_data_extract <- null_df %>% dplyr::select(-seq_name)



# randomforest model 학습
#data.train_extract <- data.train_extract[sample(NROW(data.train_extract), NROW(data.train_extract)),]

#################################################################################################################
set.seed(7777)
rf8.fit = randomForest(state ~ .
                       , data=data.train_extract,ntree = 500, importance = T )

#################################################################################################################
rf1.fit = randomForest(state ~ .
                      , data=data.train_extract[1:600000,],  ntree = 500, norm.votes=FALSE, do.trace=10, importance = T,proximity = FALSE)

rf2.fit = randomForest(state ~ .
                      , data=data.train_extract[600001:1200001,],  ntree = 500, norm.votes=FALSE, do.trace=10
                      , importance = T,proximity = FALSE,na.action=na.exclude)

rf3.fit = randomForest(state ~ .
                       , data=data.train_extract[1200002:1800002,],  ntree = 500, norm.votes=FALSE, do.trace=10
                       , importance = T,proximity = FALSE,na.action=na.exclude)

rf4.fit = randomForest(state ~ .
                       , data=data.train_extract[1800003:2070687,],  ntree = 500, norm.votes=FALSE, do.trace=10
                       , importance = T,proximity = FALSE,na.action=na.exclude)


my_combine <- function (...) 
{
  pad0 <- function(x, len) c(x, rep(0, len - length(x)))
  padm0 <- function(x, len) rbind(x, matrix(0, nrow = len - 
                                              nrow(x), ncol = ncol(x)))
  rflist <- list(...)
  areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
  if (any(!areForest)) 
    stop("Argument must be a list of randomForest objects")
  rf <- rflist[[1]]
  classRF <- rf$type == "classification"
  trees <- sapply(rflist, function(x) x$ntree)
  ntree <- sum(trees)
  rf$ntree <- ntree
  nforest <- length(rflist)
  haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
  vlist <- lapply(rflist, function(x) rownames(importance(x)))
  numvars <- sapply(vlist, length)
  if (!all(numvars[1] == numvars[-1])) 
    stop("Unequal number of predictor variables in the randomForest objects.")
  for (i in seq_along(vlist)) {
    if (!all(vlist[[i]] == vlist[[1]])) 
      stop("Predictor variables are different in the randomForest objects.")
  }
  haveForest <- sapply(rflist, function(x) !is.null(x$forest))
  if (all(haveForest)) {
    nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
    rf$forest$nrnodes <- nrnodes
    rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
    rf$forest$nodestatus <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$nodestatus, nrnodes)))
    rf$forest$bestvar <- do.call("cbind", lapply(rflist, 
                                                 function(x) padm0(x$forest$bestvar, nrnodes)))
    rf$forest$xbestsplit <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$xbestsplit, nrnodes)))
    rf$forest$nodepred <- do.call("cbind", lapply(rflist, 
                                                  function(x) padm0(x$forest$nodepred, nrnodes)))
    tree.dim <- dim(rf$forest$treemap)
    if (classRF) {
      rf$forest$treemap <- array(unlist(lapply(rflist, 
                                               function(x) apply(x$forest$treemap, 2:3, pad0, 
                                                                 nrnodes))), c(nrnodes, 2, ntree))
    }
    else {
      rf$forest$leftDaughter <- do.call("cbind", lapply(rflist, 
                                                        function(x) padm0(x$forest$leftDaughter, nrnodes)))
      rf$forest$rightDaughter <- do.call("cbind", lapply(rflist, 
                                                         function(x) padm0(x$forest$rightDaughter, nrnodes)))
    }
    rf$forest$ntree <- ntree
    if (classRF) 
      rf$forest$cutoff <- rflist[[1]]$forest$cutoff
  }
  else {
    rf$forest <- NULL
  }
  #
  #Tons of stuff removed here...
  #
  if (classRF) {
    rf$confusion <- NULL
    rf$err.rate <- NULL
    if (haveTest) {
      rf$test$confusion <- NULL
      rf$err.rate <- NULL
    }
  }
  else {
    rf$mse <- rf$rsq <- NULL
    if (haveTest) 
      rf$test$mse <- rf$test$rsq <- NULL
  }
  rf
}
rf.all <- my_combine(rf1.fit,rf2.fit,rf3.fit,rf4.fit)



# 결과 확인 1
# train data 결과
confusionMatrix(
  predict(rf.all, newdata = data.train_extract),
  data.train_extract$state
)


plot.roc.curve( model = rf.all, 
                data = data.train_extract,
                actual = data.train_extract$state, 
                title.text = "Train data ROC Curve")

# test data 결과
confusionMatrix(
  predict(rf.all, newdata = data.test_extract ),
  data.test_extract$state
)


plot.roc.curve( model = rf.all, 
                data = data.test_extract,
                actual = data.test_extract$state, 
                title.text = "Test data ROC Curve")





# null_data 결과
confusionMatrix(
  predict(rf.all, newdata = null_data_extract ),
  null_data_extract$state
)


plot.roc.curve( model = rf.all, 
                data = null_data_extract,
                actual = null_data_extract$state, 
                title.text = "Null data ROC Curve")


# 결과 확인 2
# train data 결과
confusionMatrix(
  predict(rf8.fit, newdata = data.train_extract),
  data.train_extract$state
)


plot.roc.curve( model = rf8.fit, 
                data = data.train_extract,
                actual = data.train_extract$state, 
                title.text = "Train data ROC Curve")

# test data 결과
confusionMatrix(
  predict(rf8.fit, newdata = data.test_extract ),
  data.test_extract$state
)


plot.roc.curve( model = rf8.fit, 
                data = data.test_extract,
                actual = data.test_extract$state, 
                title.text = "Test data ROC Curve")



# null_data 결과
confusionMatrix(
  predict(rf8.fit, newdata = null_data_extract ),
  null_data_extract$state
)


plot.roc.curve( model = rf8.fit, 
                data = null_data_extract,
                actual = null_data_extract$state, 
                title.text = "Null data ROC Curve")
