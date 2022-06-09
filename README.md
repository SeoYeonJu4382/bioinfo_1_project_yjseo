## bioinfo_1_project_yjseo
* * *
# **Final free mission project**

**< Introduction>** 
* **1. Thema** :
Random forest를 이용하여 LIN28A binding motif predictor 만들기
* **2. Tools** :
Running from Mac and using python and R 
* **3. Object** :

  -  CLIP-seq로 LIN28A-RNA complex를 구현 함으로써 LIN28A의 binding motif에 대한 특징을 파악할 수 있다.
  -  LIN28A-bound sequence들로 기계학습으로 예측 모델을 만들고 평가해 봄으로써 key sequence들을 파악할 수 있다.
  -  LIN28A binding motif predictor를 Random forest 모델을 적용하여 생성하고 학습시킬 수 있다.
  -  Random forest 모델링을 진행해 봄으로써 모델에 대한 이해와 방법론에 대한 원리에 대한 이해를 넓힐 수 있다.
 

**< Methods>** 

* **1.LIN28A-bound sequence** 

  ![스크린샷 2022-05-27 오후 3 05 38](https://user-images.githubusercontent.com/64352388/170640313-999877d4-cfd0-4344-8e84-71c889907be4.png) 
  
     - (1) LIN28A-bound sequence 얻기 
     - (2) Input 값 데이터로 정제 후 추출





* **2.Model** : Random forest 

    **<Random forest이란?>**
    * 다수의 의사결정나무모델에 의한 예측을 종합하는 앙상블 방법
    * 일반적으로 하나의 의사결정나무모델 보다 높은 예측 정확성을 보여줌
    * 관측치 수에 비해 변수의 수가 많은 고차원 데이터에서 중요 변수 선택 기법으로 널리 활용됨
     
     
       <img width="319" alt="스크린샷 2022-05-27 오후 2 42 23" src="https://user-images.githubusercontent.com/64352388/170637482-42acd62e-efbe-4d93-957d-51b34164f82e.png">

    **<핵심 아이디어>**
    * Diversity , Random 확보
    * 의사결정나무모델 구축 시 변수 무작위로 선택 
    
    
    
    
    
* **3.Code**
* - python input 데이터 전처리 만들기 코드 : **Final_Free_Proj_2022_24764.ipynb** 업로드 함 (2022.06.03 updated) 
* - R code 
      
      *코드 업로드(수정중)*
     
      
      #2022/05/20 yjseo
      # 01.make input data : LIN28A binding motif sequence 
       -   - **Final_Free_Proj_2022_24764.ipynb** 업로드 함
     - <img width="936" alt="스크린샷 2022-06-03 오후 10 35 38" src="https://user-images.githubusercontent.com/64352388/171864863-cb7fe18c-65b0-4895-8e53-8919697b59c3.png">
     - <img width="898" alt="스크린샷 2022-06-03 오후 10 35 56" src="https://user-images.githubusercontent.com/64352388/171864951-3a3ed185-c919-4631-a96a-8560bc35cdfd.png">


      
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



**< Results>** 

       (1) 잠재적인 LIN28A 결합 부위 주변의 패턴을 분석 -> 자주 돌연변이되는 G 앞에는 A 또는 U에 대한 선호도가 높은 2개의 염기가 오고 그 뒤에 G 또는 A를 선호하는 3개의 염기가 옴을 확인.
     
       (2) LIN28A와 상호작용 하는 AAGGAG, AAGAG 및 UGUG 요소를 포함하는 key RNA sequence를 확인.
     
       (3) 2)단계에서 찾은 key RNA sequence (AAGGAG, AAGAG 및 UGUG)를 target sequence로 잡고 Random forest로 제작한 예측 모델을 제시.
     
       (4) Acurracy, specificity, sensitivity ,Blanced Acurracy 로 모델의 성능을 평가한 결과 제시. (Confusion matrix도 추가로 제시)



**< Discussion>** 
     
     - 모델의 성능이 Acurracy기준 0.9, Blanced Acurracy기준 0.8 이하일 시 해당 모델의 성능을 높이기 위한 추가 방법론에 대한 제시 및 논의
     - 해당 LIN28A binding motif predictor를 기반으로 추후 유사한 연구의 예측모델 적용 방안에 대한 제시 및 논의 
      


* * *
student ID : 2022-27464

Date: 2022.05.20 ~ 2022.06.10

Reference: 
      
   - Cho, J., Chang, H., Kwon, S. C., Kim, B., Kim, Y., Choe, J., ... & Kim, V. N. (2012). LIN28A is a suppressor of ER-associated translation in embryonic stem cells. Cell, 151(4), 765-777.
   -  Random forest 모델 , https://www.youtube.com/watch?v=lIT5-piVtRw&list=PLpIPLT0Pf7IoTxTCi2MEQ94MZnHaxrP0j&index=19
