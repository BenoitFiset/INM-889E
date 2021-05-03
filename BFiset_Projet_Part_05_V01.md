---
title: 'Projet: TCGA HNSC LUSC ML - Partie #5'
author: "Benoit Fiset"
date: "30/04/2021"
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

## Sur l'ordinateur personnel:

### But de la section: Prédictions et résultats avec les modèles entrainé fitsvmLinear_VST et fitrf_VST avec les jeux de données Test et Normale.


```
testDataset.df_VST <-readRDS("testDataset_df_VST_Normal.rds")
fitsvmLinear_VST <- readRDS("fitsvmLinear.rds")
fitrf_VST  <- readRDS("fitrf.rds")
```

#### Bref retour sur le nombre d’échantillons dans chaque jeu de donnés

![](figures/Final_Count_DataSets_01.png){width=70%}
![](figures/Final_Count_DataSets_Normal_01.png){width=70%}

#### Spécifiquement sur le jeu de Tests


```r
print(dim(testDataset.df_VST))
kable(testDataset.df_VST [43:53,1:6],"rst")
```
```
print(dim(testDataset.df_VST))
[1]   201 11062

kable(testDataset.df_VST [43:53,1:6],"rst")


============================  ====  ==================  ==================  ==================  ==================  ==================
\                             Type  ENSG00000000003.13  ENSG00000000419.11  ENSG00000000457.12  ENSG00000000460.15  ENSG00000000938.11
============================  ====  ==================  ==================  ==================  ==================  ==================
HNSC.CV.6959.01A.11R.1915.07  HNSC             9.80056           11.043547            8.705865            9.189783            9.459919
HNSC.CV.7415.01A.11R.2081.07  HNSC            10.59721           10.644880            8.525398            8.133543            8.707573
HNSC.CV.7407.01A.11R.2081.07  HNSC            10.41610           10.160302           10.085017           10.142132            9.643313
HNSC.UF.A7J9.01A.12R.A34R.07  HNSC            10.78716           10.951942            8.982523            9.340772            9.099182
HNSC.CN.A63W.01A.11R.A30B.07  HNSC            10.90019           11.404775            9.198699            8.631291            7.636378
HNSC.BA.5557.01A.01R.1514.07  HNSC            10.67046           10.896745            8.960470            8.240411           12.008264
HNSC.BA.6869.01A.11R.1873.07  HNSC            11.38521           10.874937            8.914856            8.525548            7.969842
HNSC.UF.A7JH.01A.21R.A34R.07  HNSC            11.16086           11.083761            9.714146            9.564517            9.122768
HNSC.MT.A67A.01A.11R.A30B.07  HNSC            10.55457           11.689214            9.237489            8.203214            8.721378
HNSC.D6.A6EM.01A.21R.A31N.07  HNSC            10.43486           11.705890            8.683907            9.806529            8.970645
HNSC.CN.5360.01A.01R.1436.07  HNSC            10.13117            9.655543            9.931549            8.932524           10.058532
============================  ====  ==================  ==================  ==================  ==================  ==================
```

### C’est maintenant le temps de faire des prédictions avec le modèle fitsvmLinear_VST et le jeu de Tests pour voir comment le modèle est précis avec des données nouvelles.

Pour pas que le modèle « triche » dans ses prédictions, nous devons enlever la colonne Type dand le jeu Tests qui permet l’identification de la classe de prédiction du modèle. Serait trop facile !!!

```r
dim(testDataset.df_VST) 
dim(testDataset.df_VST[,-1]) # Remove the Type colunm from TestDataset
kable(testDataset.df_VST[,-1] [43:53,1:5],"rst")
```
```
> dim(testDataset.df_VST) 
[1]   201 11062

> dim(testDataset.df_VST[,-1]) # Remove the Type colunm from TestDataset
[1]   201 11061

> kable(testDataset.df_VST[,-1] [43:53,1:5],"rst")

============================  ==================  ==================  ==================  ==================  ==================
\                             ENSG00000000003.13  ENSG00000000419.11  ENSG00000000457.12  ENSG00000000460.15  ENSG00000000938.11
============================  ==================  ==================  ==================  ==================  ==================
HNSC.CV.6959.01A.11R.1915.07             9.80056           11.043547            8.705865            9.189783            9.459919
HNSC.CV.7415.01A.11R.2081.07            10.59721           10.644880            8.525398            8.133543            8.707573
HNSC.CV.7407.01A.11R.2081.07            10.41610           10.160302           10.085017           10.142132            9.643313
HNSC.UF.A7J9.01A.12R.A34R.07            10.78716           10.951942            8.982523            9.340772            9.099182
HNSC.CN.A63W.01A.11R.A30B.07            10.90019           11.404775            9.198699            8.631291            7.636378
HNSC.BA.5557.01A.01R.1514.07            10.67046           10.896745            8.960470            8.240411           12.008264
HNSC.BA.6869.01A.11R.1873.07            11.38521           10.874937            8.914856            8.525548            7.969842
HNSC.UF.A7JH.01A.21R.A34R.07            11.16086           11.083761            9.714146            9.564517            9.122768
HNSC.MT.A67A.01A.11R.A30B.07            10.55457           11.689214            9.237489            8.203214            8.721378
HNSC.D6.A6EM.01A.21R.A31N.07            10.43486           11.705890            8.683907            9.806529            8.970645
HNSC.CN.5360.01A.01R.1436.07            10.13117            9.655543            9.931549            8.932524           10.058532
============================  ==================  ==================  ==================  ==================  ==================
```

***

#### Maintenant faisons les prédictions tant attendues du modèle fitsvmLinear_VST et le jeu de Tests


```r
predictions <- predict(fitsvmLinear_VST, testDataset.df_VST[,-1])  # Remove the Type colunm from TestDataset
```

### Voyons ce que la matrice de confusion (résultats de prédiction du modèle) donne.


```r
confusionMatrix(predictions, testDataset.df_VST$Type)
```
```
Confusion Matrix and Statistics

          Reference
Prediction HNSC LUSC
      HNSC  108    0
      LUSC    3   90
                                         
               Accuracy : 0.9851         
                 95% CI : (0.957, 0.9969)
    No Information Rate : 0.5522         
    P-Value [Acc > NIR] : <2e-16         
                                         
                  Kappa : 0.9699         
                                         
 Mcnemar's Test P-Value : 0.2482         
                                         
            Sensitivity : 0.9730         
            Specificity : 1.0000         
         Pos Pred Value : 1.0000         
         Neg Pred Value : 0.9677         
             Prevalence : 0.5522         
         Detection Rate : 0.5373         
   Detection Prevalence : 0.5373         
      Balanced Accuracy : 0.9865         
                                         
       'Positive' Class : HNSC
```
Nous pouvons voir de ces résultats que le modèle à un très haut taux de précision «Accuracy» =  98.51%. Il est aussi important de regarder les valeurs de « Sensitivity » = (True positive rate) = 97.30% et « Specificity » = (True negative rate) = 100% et la valeur de Kappa = 96.99%

La valeur de précision «Accuracy »  est le pourcentage que le modèle classifie correctement une valeur de tous les valeurs.

La valeur de Kappa (Cohen’s Kappa) est comme la précision mais ajoute un élément de correction de probabilité de classification dû  « au hasard » par le modèle.  Plus le Kappa est haut plus que la classification est parfaite. 

Ici le Kappa à 96.99% est tres haut... pouvons dire que pas de classifications dû au hasard ici.

Ce qui dit la matrice de confusion :

* 108 « True Positive - TP »  -  0  « False Positive - FP»
* 3   « False Negative - FN » -  90 « True Negative - TN »

Donc en d’autres termes il y a eu 3 échantillons LUSC identifiés faussement comme HNSC (« Sensitivity » = (True positive rate) = 97.30% ) et qu’il y a eu aucun échantillon HNSC identifié comme LUSC (« Specificity » = (True negative rate) = 100%).

Les calculs que font la fonction confusionMatrix() sont le suivants:

* « Accuracy »    = (TP + TN) / (TP + FP + FN + TN) 
* « Sensitivity » = (True positive rate - TPR ): TPR = TP / P = TP / (TP+FN)
* « Specificity » = (True negative rate - TNR) : TN / N = TN / (TN+FP)

***


#### Maintenant faisons les prédictions du modèle fitrf_VST (Random Forest) et le jeu de Tests


```r
predictions <- predict(fitrf_VST, testDataset.df_VST[,-1])  # Remove the Type colunm from TestDataset
```

### Voyons ce que la matrice de confusion (résultats de prédiction du modèle) donne.


```r
confusionMatrix(predictions, testDataset.df_VST$Type)
```
```
Confusion Matrix and Statistics

          Reference
Prediction HNSC LUSC
      HNSC   94    3
      LUSC   17   87
                                          
               Accuracy : 0.9005          
                 95% CI : (0.8505, 0.9382)
    No Information Rate : 0.5522          
    P-Value [Acc > NIR] : < 2e-16         
                                          
                  Kappa : 0.8017          
                                          
 Mcnemar's Test P-Value : 0.00365         
                                          
            Sensitivity : 0.8468          
            Specificity : 0.9667          
         Pos Pred Value : 0.9691          
         Neg Pred Value : 0.8365          
             Prevalence : 0.5522          
         Detection Rate : 0.4677          
   Detection Prevalence : 0.4826          
      Balanced Accuracy : 0.9068          
                                          
       'Positive' Class : HNSC  

```

Nous pouvons voir de ces résultats que le modèle à une précision «Accuracy» = 90.05%. Il est aussi important de regarder les valeurs de « Sensitivity » = (True positive rate) = 84.68% et « Specificity » = (True negative rate) = 96.67% et la valeur de Kappa = 80.17%

Ici le Kappa à 80.17% est assez haut... pas beaucoup de classification dû au hasard.

Ce qui dit la matrice de confusion :

* 94 « True Positive - TP »    -  3  « False Positive - FP»
* 17   « False Negative - FN » -  87 « True Negative - TN »

Donc en d’autres termes il y a eu 17 échantillons LUSC identifiés faussement comme HNSC (« Sensitivity » = (True positive rate) = 84.68% ) et qu’il y a eu 3 échantillon HNSC identifié comme LUSC (« Specificity » = (True negative rate) = 96.67%).

Nous pouvons voir clairement que le model SVM fitsvmLinear_VST est plus précis avec sa classification que Random Forest fitrf_VST.


```r
#Prediect with Normal

fitsvmLinear <- readRDS("fitsvmLinear.rds")
normalDataset.df_VST  <- readRDS("NormalDataset_df_VST_Normal.rds")

dim(normalDataset.df )
kable(normalDataset.df  [43:53,1:6],"rst")
normalDataset.df  <- normalDataset.df [,-1] # Remove the Type colunm from TestDataset
dim(normalDataset.df )
kable(normalDataset.df  [43:53,1:6],"rst")

TrainCols <- colnames(fitsvmLinear_VST$trainingData)[-1] # Colnames of model and remove first value "Outcome"
NormalCols <- colnames(normalDataset.df) # Colnames of dataset

NormalNotInTrain <- subset(TrainCols, !(TrainCols %in% NormalCols))  
#TrainNotInCols2Normal <- subset(NormalCols, !(NormalCols %in% TrainCols))

normalDataset.df[NormalNotInTrain] <- 0 # Create Missing Columns and put value to 0

predictions <- predict(fitsvmLinear_VST, newdata=normalDataset.df)
normalDataset.df <-readRDS("NormalDataset_df_VST_Normal.rds") # reload to get Type Colunm back
confusionMatrix(predictions, normalDataset.df$Type)

predictions <- predict(fitrf_VST, newdata=normalDataset.df)
normalDataset.df <-readRDS("NormalDataset_df_VST_Normal.rds") # reload to get Type Colunm back
confusionMatrix(predictions, normalDataset.df$Type)
```
