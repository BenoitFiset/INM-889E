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
