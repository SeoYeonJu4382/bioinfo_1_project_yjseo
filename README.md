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
      library 불러오기
      library(readr)
      library(dplyr)
      library(caret)
      library(randomForest)
      library(ComplexHeatmap)
      library(ROCR)
      
      
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


      data <- read.csv(*생성 예정*)
      data$state <- as.factor(data$state)

      data.train <- data %>% filter(*훈련 데이터셋 index지정해서 추출 할 예정*) %>% select(-index)
      data.test <- data %>% filter(*테스트 데이터셋 index지정해서 추출 할 예정*) %>% select(-index) 

      
      #데이터 준비
      data.train_extract <- data.train %>% select(-sample)
      data.test_extract <- data.test %>% select(-sample) 

      # randomforest model 학습
      data.train_extract <- data.train_extract[sample(NROW(data.train_extract), NROW(data.train_extract)),]
      rf.fit = randomForest(state ~ .
                      , data=data.train_extract, mtry = floor(sqrt(500)), ntree = 500, importance = T)

      # 결과 확인
      # train data 결과
      confusionMatrix(
        predict(rf.fit, newdata = data.train_extract),
        data.train_extract$state
      )

      
      plot.roc.curve( model = rf.fit, 
                data = data.train_extract,
                actual = data.train_extract$state, 
                title.text = "Train data ROC Curve")

      # test data 결과
      confusionMatrix(
        predict(rf.fit, newdata = data.test_extract ),
        data.test_extract$state
      )

      plot.roc.curve( model = rf.fit, 
                data = data.test_extract,
                actual = data.test_extract$state, 
                title.text = "Test data ROC Curve")



**< Results>** 

       (1) 잠재적인 LIN28A 결합 부위 주변의 패턴을 분석 -> 자주 돌연변이되는 G 앞에는 A 또는 U에 대한 선호도가 높은 2개의 염기가 오고 그 뒤에 G 또는 A를 선호하는 3개의 염기가 옴을 확인.
     
       (2) LIN28A와 상호작용 하는 AAGGAG, AAGAG 및 UGUG 요소를 포함하는 key RNA sequence를 확인.
     
       (3) 2)단계에서 찾은 key RNA sequence (AAGGAG, AAGAG 및 UGUG)를 target sequence로 잡고 Random forest로 제작한 예측 모델을 제시.
     
       (4) Acurracy, specificity, sensitivity ,Blanced Acurracy 로 모델의 성능을 평가한 결과 제시. (Confusion matrix도 추가로 제시)



**< Discussion>** 
     
     - 모델의 성능이 Acurracy기준 0.9, Blanced Acurracy기준 0.8 이하일 시 해당 모델의 성능을 높이기 위한 추가 방법론에 대한 제시 및 논의
     - 해당 LIN28A binding motif predictor를 기반으로 추후 유사한 연구의 예측모델 적용 방안에 대한 제시 및 논의 
       *(프로젝트 결과에 따라 변동 존재)*
   


* * *
student ID : 2022-27464

Date: 2022.05.20 ~ 2022.06.10

Reference: 
      
   - Cho, J., Chang, H., Kwon, S. C., Kim, B., Kim, Y., Choe, J., ... & Kim, V. N. (2012). LIN28A is a suppressor of ER-associated translation in embryonic stem cells. Cell, 151(4), 765-777.
   -  Random forest 모델 , https://www.youtube.com/watch?v=lIT5-piVtRw&list=PLpIPLT0Pf7IoTxTCi2MEQ94MZnHaxrP0j&index=19
