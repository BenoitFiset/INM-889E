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

## Maintenant faisons les prédictions tant attendues du modèle fitsvmLinear_VST et le jeu de Tests


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

## Maintenant faisons les prédictions du modèle fitrf_VST (Random Forest) et le jeu de Tests


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

### Nous pouvons voir clairement que le model SVM fitsvmLinear_VST est plus précis avec sa classification que Random Forest fitrf_VST.

***

## Maintenant faisons les prédictions des modèles fitsvmLinear_VST et fitrf_VST (Random Forest) et le jeu Normale

Il y aura un mimimum de sorties... allons se concentre sur les résultats


```r
# normalDataset.df_VST  <- readRDS("NormalDataset_df_VST_Normal.rds")

dim(normalDataset.df_VST )
kable(normalDataset.df_VST  [43:53,1:6],"rst")
```
```
 dim(normalDataset.df_VST )
[1]    93 17427

kable(normalDataset.df_VST  [43:53,1:6],"rst")

============================  ====  ==================  ==================  ==================  ==================  ==================
\                             Type  ENSG00000000003.13  ENSG00000000419.11  ENSG00000000457.12  ENSG00000000460.15  ENSG00000000938.11
============================  ====  ==================  ==================  ==================  ==================  ==================
HNSC.CV.7423.11A.01R.2081.07  HNSC           10.363988           10.254031            7.277398            6.154062            6.730621
HNSC.CV.7416.11A.01R.2081.07  HNSC           10.074296           10.618110            6.754814            6.101676            6.063347
LUSC.51.4079.11A.01R.1758.07  LUSC            8.603351            8.154606            7.070855            4.398018           10.302034
LUSC.56.7582.11A.01R.2045.07  LUSC            9.170808            8.478738            7.512434            5.312511           10.315651
LUSC.85.7710.11A.01R.2125.07  LUSC            8.108887            8.971991            8.025179            5.596142           11.611160
LUSC.56.8083.11A.01R.2247.07  LUSC            8.392597            8.378330            7.883732            4.595641           11.461918
LUSC.43.6647.11A.01R.1820.07  LUSC            8.059541            8.271135            7.189719            4.535879           12.269652
LUSC.56.7579.11A.01R.2045.07  LUSC            9.507510            9.134634            7.860343            4.894959           11.730293
LUSC.43.7658.11A.01R.2125.07  LUSC            9.835925            8.850212            8.112907            4.789651           11.254981
LUSC.77.7138.11A.01R.2045.07  LUSC            8.318903            8.777017            7.797314            4.556016           11.362884
LUSC.56.8309.11A.01R.2296.07  LUSC           10.813934            8.682587            7.417639            4.962103           12.689825
============================  ====  ==================  ==================  ==================  ==================  ==================
```

Il y a une différence du nombre de gènes post filtrages pour le jeu de données Normale. Il y a au finale 17427 gènes restant vs 11062 pour les jeux « Training » et « Test ». 

Les gènes qui sont inclus dans les modèles entrainés fitsvmLinear_VST et fitrf_VST **DOIVENT** se retrouver dans le jeu de données Normale sinon il y aura erreur de la fonction predict() qui dira pas possible de prédire car il manque tel ou tel gène….

Cela dit, il faut faire une étape de validation pour etre certain que tous les gènes des modèles entrainés fitsvmLinear_VST et fitrf_VST se retrouvent dans les jeu de données Normale. 


```r
TrainCols <- colnames(fitsvmLinear_VST$trainingData)[-1] # Colnames of model and remove first value "Outcome"
NormalCols <- colnames(normalDataset.df_VST) # Colnames of dataset
```

Pouvons voir la taille différente des jeux de données .
```
print(length(TrainCols))
[1] 11057

print(length(NormalCols))
[1] 17427
```

C’est ici que nous validons la concordance des gènes entre le modèle et le jeu de données Normale. S’il y en a une, ajout de colonne du/des gènes manquants en mettant un compte de 0 pour ce/ces gènes. 

```r
NormalNotInTrain <- subset(TrainCols, !(TrainCols %in% NormalCols))  # Check missing Genes between Train and Normal
normalDataset.df_VST[NormalNotInTrain] <- 0 # Create Missing Columns and put value to 0

print(NormalNotInTrain) # See what is missing...
```
Sommes chanceux, il ne manquait qu'un gène
```
[1] "ENSG00000135094.9"
```

```r
print(normalDataset.df_VST[,"ENSG00000135094.9"])
```
La colonne est bien ajoutée avec des valeurs de 0
```
 [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
[77] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

## Maintenant faisons les prédictions du modèle fitrf_VST (Random Forest) et le jeu Normale


```r
predictions <- predict(fitsvmLinear_VST, newdata=normalDataset.df_VST[,-1]) # Remove the Type colunm from normalDataset
confusionMatrix(predictions, normalDataset.df_VST$Type)
```
```
Confusion Matrix and Statistics

          Reference
Prediction HNSC LUSC
      HNSC   44    0
      LUSC    0   49
                                     
               Accuracy : 1          
                 95% CI : (0.9611, 1)
    No Information Rate : 0.5269     
    P-Value [Acc > NIR] : < 2.2e-16  
                                     
                  Kappa : 1          
                                     
 Mcnemar's Test P-Value : NA         
                                     
            Sensitivity : 1.0000     
            Specificity : 1.0000     
         Pos Pred Value : 1.0000     
         Neg Pred Value : 1.0000     
             Prevalence : 0.4731     
         Detection Rate : 0.4731     
   Detection Prevalence : 0.4731     
      Balanced Accuracy : 1.0000     
                                     
       'Positive' Class : HNSC  
```

Nous pouvons voir de ces résultats que le modèle à une précision «Accuracy» = 100%. Le reste des valeurs «Sensitivity», «Specificity» et Kappa sont aussi 100%

Ce qui dit la matrice de confusion : pas grand-chose sauf que le modèle a prédit tous les échantillons correctement.

***

## Maintenant faisons les prédictions du modèle fitrf_VST (Random Forest) et le jeu Normale

```r
predictions <- predict(fitrf_VST, newdata=normalDataset.df_VST[,-1]) # Remove the Type colunm from normalDataset
confusionMatrix(predictions, normalDataset.df_VST$Type)
```
```
Confusion Matrix and Statistics

          Reference
Prediction HNSC LUSC
      HNSC   42   43
      LUSC    2    6
                                          
               Accuracy : 0.5161          
                 95% CI : (0.4101, 0.6211)
    No Information Rate : 0.5269          
    P-Value [Acc > NIR] : 0.6228          
                                          
                  Kappa : 0.0735          
                                          
 Mcnemar's Test P-Value : 2.479e-09       
                                          
            Sensitivity : 0.9545          
            Specificity : 0.1224          
         Pos Pred Value : 0.4941          
         Neg Pred Value : 0.7500          
             Prevalence : 0.4731          
         Detection Rate : 0.4516          
   Detection Prevalence : 0.9140          
      Balanced Accuracy : 0.5385          
                                          
       'Positive' Class : HNSC    
```

Nous pouvons voir de ces résultats, avec un précision «Accuracy» = 51.61%,  que le modèle à de la difficulté à prédire avec ce jeu de données surtout la classe LUSC. 

Ceci se reflète avec la valeur de « Sensitivity » = (True positive rate) = 95.45 % et « Specificity » = (True negative rate) = 12.24%. La valeur de Kappa = 7.35%

Ici le Kappa à 7.35% est très très bas... nous serions mieux de faire des prédictions à la main au hasard !!!

Ce qui dit la matrice de confusion :

* 42 « True Positive   - TP » -  43 « False Positive - FP»
* 2  « False Negative  - FN » -  6  « True Negative - TN »

Donc en d’autres termes il y a eu:

* 42 échantillons HNSC bien identifié comme HNSC
* 2  échantillons HNSC identifiés faussement comme LUSC (« Sensitivity » = (True positive rate) = 95.45%
* 43 échantillons LUSC identifiés faussement comme HNSC (« Specificity » = (True negative rate) = 12.24%)
* 6  échantillons LUSC bien identifié comme LUSC


