---
title: 'Projet: TCGA HNSC LUSC ML - Partie #4'
author: "Benoit Fiset"
date: "27/04/2021"
output:
  html_document: 
    keep_md: yes
---

<style type="text/css">
/*https://stackoverflow.com/questions/38367392/override-rmarkdown-theme-in-order-to-change-html-page-width/38373846*/
body .main-container {
  max-width: 1100px !important;
  width: 1100px !important;
}
body {
  margin: auto;
  max-width: 1100px !important;
}
</style>



# Classification par apprentissage automatique du type de cancer à partir de données d'expression génétique de séquençage RNA-Seq

***

### But de la section : Entrainement des modèles de machine learning. 
### Pour le projet les algorithmes de prediction **SVM** et **randomForest** seront utilisés.

### Sur serveur HPC:
```
library(knitr)      # kable
library(doParallel)
library(caret)      # caret (Classification And REgression Training) - creating predictive models
library(randomForest)
library(kernlab)    # SVM
```

```r
# Important pour tjrs avoir les memes resultats

set.seed(1234) 
```
***
# Notes:
```
J'ai fait le choix du package R caret (Classification And REgression Training) pour les 
modèles de prédictions pour ce projet. Il en existe d’autres comme e1071 (interface R à libsvm)
et kernlab (kvsm) pour des SVM. 

L’avantage de caret et qu’il fait déjà interface « under the hood » avec les package kernlab, e1071,
randomForest et 236 autre outils de régression et classification et le tout configurable dépendant
du modèle et méthode choisi. C'est vraiment un "couteau suisse" pour créer des modèles de prédictions.
```
```
Une des forces du package caret et de pouvoir directement influencer l’entrainement du modèle
avec des paramètres d‘optimisation dès le début.  Ceci veut dire que qu’au lieu de débuter 
avec les paramètres de base pour le premier entrainement et plus tard essayer d’optimiser les
performances, nous pouvons déjà dès le début dire au model d’essayer de trouver le meilleur
modèle de prédiction en utilisant des variances des métriques d’influences.

La fonction de caret qui permet cela est trainControl() qui génère les paramètres qui 
permettent de contrôler comment les modèles optimisés sont créés.  Des options possibles 
dont des « bootstrap », « K-fold cross-validation », « Leave One Out cross-validation , 
« Repeated K-fold cross-validation », etc..
```

***

#### Initialisation des « worker services » 
```
Ceci permet l’entrainement des modèles pour bénéficier du parallélisme des serveurs haute performance.
L’utilisation de 32 CPU aide énormément à réduire le temps nécessaire pour entrainer les modèles. 

Exemple, avec 32 CPU, l’entrainement de Random Forest prends presque 2 heures. Lorsque essayé avec ordinateur 
personnel, en utilisant l’équivalent de 5 CPU, l’entrainement n’était pas fini après 18 heures de calculs !!
```


```r
set.seed(1234)
cl <- makeCluster(32, type='PSOCK', outfile="OutCaret.txt")
registerDoParallel(cl)
```

```
Dans le projet j’ai choisi de faire l’optimisation avec « CV » qui est « K-fold cross-validation » 
avec un nombre de « 10 » qui control le nombre de « folds ».
```
#### 

```r
# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10, allowParallel = TRUE, verboseIter = TRUE)
metric <- "Accuracy"
```

## Support Vector Machines (SVM) avec "linear kernel"


```r
#Le "seed" du nombre aléatoire garantit que les résultats sont directement comparables.
set.seed(1234)

fitsvmLinear <- caret::train(Type~., data=trainingDataset.df , method='svmLinear', metric=metric, trControl=control)
```

***

## Random Forest (RF)

```r
set.seed(1234)
fitrf <- caret::train(Type~., data=trainingDataset.df , method="rf", metric=metric, trControl=control, verbose = TRUE)
```

## Bonus: Comparaison des performances de 6 algorithmes de ML avec les différentes normalisations des données

```
J'ai fait des test de comparaison avec ces 6 models 

```


```r
control <- trainControl(method="cv", number=10, allowParallel = TRUE, verboseIter = TRUE)
metric <- "Accuracy"

# a) linear algorithms
set.seed(1234)
fitlda     <- caret::train(Type~., data=trainingDataset.df , method="lda", metric=metric, trControl=control)
fitlda_Log <- caret::train(Type~., data=trainingDataset.df_Log , method="lda", metric=metric, trControl=control)
fitlda_VST <- caret::train(Type~., data=trainingDataset.df_VST , method="lda", metric=metric, trControl=control)

# b) nonlinear algorithms
# CART
set.seed(1234)
fitcart     <- caret::train(Type~., data=trainingDataset.df , method="rpart", metric=metric, trControl=control)
fitcart_Log <- caret::train(Type~., data=trainingDataset.df_Log , method="rpart", metric=metric, trControl=control)
fitcart_VST <- caret::train(Type~., data=trainingDataset.df_VST , method="rpart", metric=metric, trControl=control)

# kNN
set.seed(1234)
fitknn     <- caret::train(Type~., data=trainingDataset.df , method="knn", metric=metric, trControl=control)
fitknn_Log <- caret::train(Type~., data=trainingDataset.df_Log , method="knn", metric=metric, trControl=control)
fitknn_VST <- caret::train(Type~., data=trainingDataset.df_VST , method="knn", metric=metric, trControl=control)

# c) advanced algorithms
# SVM
set.seed(1234)
fitsvmRadial     <- caret::train(Type~., data=trainingDataset.df , method="svmRadial", metric=metric, trControl=control)
fitsvmRadial_Log <- caret::train(Type~., data=trainingDataset.df_Log , method="svmRadial", metric=metric, trControl=control)
fitsvmRadial_VST <- caret::train(Type~., data=trainingDataset.df_VST , method="svmRadial", metric=metric, trControl=control)

# SVM
set.seed(1234)
fitsvmLinear     <- caret::train(Type~., data=trainingDataset.df , method='svmLinear', metric=metric, trControl=control)
fitsvmLinear_Log <- caret::train(Type~., data=trainingDataset.df_Log , method='svmLinear', metric=metric, trControl=control)
fitsvmLinear_VST <- caret::train(Type~., data=trainingDataset.df_VST , method='svmLinear', metric=metric, trControl=control)

# Random Forest
set.seed(1234)
fitrf      <- caret::train(Type~., data=trainingDataset.df , method="rf", metric=metric, trControl=control, verbose = TRUE)
fitrf_Log  <- caret::train(Type~., data=trainingDataset.df_Log , method="rf", metric=metric, trControl=control, verbose = TRUE)
fitrf_VST  <- caret::train(Type~., data=trainingDataset.df_VST , method="rf", metric=metric, trControl=control, verbose = TRUE)

# The End
```


```r
results <- resamples(list(lda=fitlda,lda_Log=fitlda_Log,lda_VST=fitlda_VST, cart=fitcart, cart_Log=fitcart_Log,
          cart_VST=fitcart_VST, knn=fitknn,knn_Log=fitknn_Log,knn_VST=fitknn_VST,svmRadial=fitsvmRadial,
          svmRadial_Log=fitsvmRadial_Log,svmRadial_VST=fitsvmRadial_VST, svmLinear=fitsvmLinear,
          svmLinear_Log=fitsvmLinear_Log,svmLinear_VST=fitsvmLinear_VST, rf=fitrf,rf_Log=fitrf_Log, 
          rf_VST=fitrf_VST))

summary(results)
```
```
Call:
summary.resamples(object = results)

Models: lda, lda_Log, lda_VST, cart, cart_Log, cart_VST, knn, knn_Log, knn_VST, svmRadial, svmRadial_Log, svmRadial_VST, svmLinear, svmLinear_Log, svmLinear_VST, rf, rf_Log, rf_VST 
Number of resamples: 10 

Accuracy 
                   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
lda           0.9125000 0.9412037 0.9685918 0.9563264 0.9750000 0.9875000    0
lda_Log       0.9625000 0.9780854 0.9875000 0.9837650 0.9876157 1.0000000    0
lda_VST       0.9625000 0.9750000 0.9875000 0.9824992 0.9876157 1.0000000    0
cart          0.7375000 0.8255401 0.8323302 0.8352422 0.8495253 0.9125000    0
cart_Log      0.7375000 0.8255401 0.8323302 0.8352422 0.8495253 0.9125000    0
cart_VST      0.8227848 0.8406250 0.8750000 0.8638063 0.8843750 0.8888889    0
knn           0.8000000 0.8485759 0.8687500 0.8700107 0.8885417 0.9506173    0
knn_Log       0.9000000 0.9127701 0.9245253 0.9262940 0.9343750 0.9753086    0
knn_VST       0.9000000 0.9250000 0.9316358 0.9338102 0.9498418 0.9629630    0
svmRadial     0.9125000 0.9369066 0.9437500 0.9474826 0.9597222 0.9876543    0
svmRadial_Log 0.9500000 0.9530063 0.9625000 0.9662488 0.9752315 0.9876543    0
svmRadial_VST 0.9500000 0.9621440 0.9689043 0.9712488 0.9875000 0.9876543    0
svmLinear     0.9500000 0.9531250 0.9625000 0.9662805 0.9752315 0.9876543    0
svmLinear_Log 0.9625000 0.9750000 0.9875000 0.9824992 0.9876157 1.0000000    0
svmLinear_VST 0.9620253 0.9781250 0.9875000 0.9824834 0.9876157 1.0000000    0
rf            0.9000000 0.9250000 0.9317130 0.9363727 0.9500000 0.9746835    0
rf_Log        0.9000000 0.9250000 0.9317130 0.9363260 0.9500000 0.9753086    0
rf_VST        0.9125000 0.9187500 0.9378858 0.9400760 0.9623813 0.9629630    0

Kappa 
                   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
lda           0.8245614 0.8825499 0.9371938 0.9126176 0.9500312 0.9749687    0
lda_Log       0.9250000 0.9561322 0.9749844 0.9675078 0.9751890 1.0000000    0
lda_VST       0.9250000 0.9499687 0.9749844 0.9649750 0.9751890 1.0000000    0
cart          0.4716981 0.6501503 0.6652746 0.6695214 0.6977937 0.8243413    0
cart_Log      0.4716981 0.6501503 0.6652746 0.6695214 0.6977937 0.8243413    0
cart_VST      0.6429955 0.6804040 0.7492163 0.7265439 0.7679974 0.7772686    0
knn           0.5977373 0.6962417 0.7367309 0.7393602 0.7768661 0.9014599    0
knn_Log       0.7998749 0.8251914 0.8485060 0.8522640 0.8684205 0.9504587    0
knn_VST       0.7998749 0.8497182 0.8630833 0.8673280 0.8993139 0.9258920    0
svmRadial     0.8243413 0.8731380 0.8872650 0.8946446 0.9193721 0.9752521    0
svmRadial_Log 0.8998121 0.9058833 0.9248591 0.9323767 0.9504564 0.9752521    0
svmRadial_VST 0.8998121 0.9239935 0.9377730 0.9424018 0.9749687 0.9752521    0
svmLinear     0.8998121 0.9061796 0.9250000 0.9324843 0.9503886 0.9752973    0
svmLinear_Log 0.9250000 0.9499217 0.9749844 0.9649687 0.9751890 1.0000000    0
svmLinear_VST 0.9238677 0.9562187 0.9749844 0.9649323 0.9751890 1.0000000    0
rf            0.7991212 0.8495769 0.8631675 0.8724097 0.8999061 0.9491961    0
rf_Log        0.7991212 0.8495298 0.8631675 0.8723097 0.8999061 0.9505495    0
rf_VST        0.8241206 0.8372574 0.8757434 0.8798831 0.9245392 0.9257562    0
```


```r
results <- resamples(list(lda=fitlda, cart=fitcart, knn=fitknn, svmRadial=fitsvmRadial, svmLinear=fitsvmLinear, 
            rf=fitrf, lda_Log=fitlda_Log, cart_Log=fitcart_Log, knn_Log=fitknn_Log, svmRadial_Log=fitsvmRadial_Log,
            svmLinear_Log=fitsvmLinear_Log, rf_Log=fitrf_Log, lda_VST=fitlda_VST, cart_VST=fitcart_VST, 
            knn_VST=fitknn_VST, svmRadial_VST=fitsvmRadial_VST, svmLinear_VST=fitsvmLinear_VST, rf_VST=fitrf_VST))
 
summary(results)
```

```
Call:
summary.resamples(object = results)

Models: lda, cart, knn, svmRadial, svmLinear, rf, lda_Log, cart_Log, knn_Log, svmRadial_Log, svmLinear_Log, rf_Log, lda_VST, cart_VST, knn_VST, svmRadial_VST, svmLinear_VST, rf_VST 
Number of resamples: 10 

Accuracy 
                   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
lda           0.9125000 0.9412037 0.9685918 0.9563264 0.9750000 0.9875000    0
cart          0.7375000 0.8255401 0.8323302 0.8352422 0.8495253 0.9125000    0
knn           0.8000000 0.8485759 0.8687500 0.8700107 0.8885417 0.9506173    0
svmRadial     0.9125000 0.9369066 0.9437500 0.9474826 0.9597222 0.9876543    0
svmLinear     0.9500000 0.9531250 0.9625000 0.9662805 0.9752315 0.9876543    0
rf            0.9000000 0.9250000 0.9317130 0.9363727 0.9500000 0.9746835    0
lda_Log       0.9625000 0.9780854 0.9875000 0.9837650 0.9876157 1.0000000    0
cart_Log      0.7375000 0.8255401 0.8323302 0.8352422 0.8495253 0.9125000    0
knn_Log       0.9000000 0.9127701 0.9245253 0.9262940 0.9343750 0.9753086    0
svmRadial_Log 0.9500000 0.9530063 0.9625000 0.9662488 0.9752315 0.9876543    0
svmLinear_Log 0.9625000 0.9750000 0.9875000 0.9824992 0.9876157 1.0000000    0
rf_Log        0.9000000 0.9250000 0.9317130 0.9363260 0.9500000 0.9753086    0
lda_VST       0.9625000 0.9750000 0.9875000 0.9824992 0.9876157 1.0000000    0
cart_VST      0.8227848 0.8406250 0.8750000 0.8638063 0.8843750 0.8888889    0
knn_VST       0.9000000 0.9250000 0.9316358 0.9338102 0.9498418 0.9629630    0
svmRadial_VST 0.9500000 0.9621440 0.9689043 0.9712488 0.9875000 0.9876543    0
svmLinear_VST 0.9620253 0.9781250 0.9875000 0.9824834 0.9876157 1.0000000    0
rf_VST        0.9125000 0.9187500 0.9378858 0.9400760 0.9623813 0.9629630    0

Kappa 
                   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
lda           0.8245614 0.8825499 0.9371938 0.9126176 0.9500312 0.9749687    0
cart          0.4716981 0.6501503 0.6652746 0.6695214 0.6977937 0.8243413    0
knn           0.5977373 0.6962417 0.7367309 0.7393602 0.7768661 0.9014599    0
svmRadial     0.8243413 0.8731380 0.8872650 0.8946446 0.9193721 0.9752521    0
svmLinear     0.8998121 0.9061796 0.9250000 0.9324843 0.9503886 0.9752973    0
rf            0.7991212 0.8495769 0.8631675 0.8724097 0.8999061 0.9491961    0
lda_Log       0.9250000 0.9561322 0.9749844 0.9675078 0.9751890 1.0000000    0
cart_Log      0.4716981 0.6501503 0.6652746 0.6695214 0.6977937 0.8243413    0
knn_Log       0.7998749 0.8251914 0.8485060 0.8522640 0.8684205 0.9504587    0
svmRadial_Log 0.8998121 0.9058833 0.9248591 0.9323767 0.9504564 0.9752521    0
svmLinear_Log 0.9250000 0.9499217 0.9749844 0.9649687 0.9751890 1.0000000    0
rf_Log        0.7991212 0.8495298 0.8631675 0.8723097 0.8999061 0.9505495    0
lda_VST       0.9250000 0.9499687 0.9749844 0.9649750 0.9751890 1.0000000    0
cart_VST      0.6429955 0.6804040 0.7492163 0.7265439 0.7679974 0.7772686    0
knn_VST       0.7998749 0.8497182 0.8630833 0.8673280 0.8993139 0.9258920    0
svmRadial_VST 0.8998121 0.9239935 0.9377730 0.9424018 0.9749687 0.9752521    0
svmLinear_VST 0.9238677 0.9562187 0.9749844 0.9649323 0.9751890 1.0000000    0
rf_VST        0.8241206 0.8372574 0.8757434 0.8798831 0.9245392 0.9257562    0
```
