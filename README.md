# INM-889E
Repo pour le cours INM 889E - Projet de ML

***

# Classification par apprentissage automatique du type de cancer à partir de données d'expression génétique de séquençage RNA-Seq

* Vu la taille des fichiers et le temps de calcul nécessaire une partie du traitement des fichiers ont été effectués sur des serveurs HPC de Calcul Québec (Béluga et Graham):
  + Le téléchargement des fichiers incluant le regroupement et décompression
  + Le regroupement des fichiers individuels en une matrice unifié de comptes - Tumeurs
  + Le regroupement des fichiers individuels en une matrice unifié de comptes - Normale
  + Entrainement des modèles (Algo. SVM et Random Forest)
  + “Tuning" des modèles
  + Prediction des modèles 

* Sur le l’ordinateur local:
  + Le filtrage des genes qui ont comportent plus de 95% de comptes de valeur 0
  + Le filtrage des genes qui ont moins de 1% de comptes
  + La creation de le colonne qui sera la “Classe” du projet. Se nomme “Type”
  + Le découpage du jeu de données “Training”  80% et “Test” 20%. 
  + Normalisation des données “Training” avec VST de DESeq2
  + Normalisation des données “Test” avec VST de DESeq2
  + Normalisation des données “Normale” avec VST de DESeq2
  + Filtrer les genes qui ont un indice de correlation de plus de 98% (Tumeur et Normale)

***

```
Le projet à été fait en R avec les libraires suivantes:
    library(caret)        # caret (Classification And REgression Training) - creating predictive models
    library(DESeq2)       # variance stabilizing transformation (VST) Normalisation  
    library(doParallel)
    library(FactoMineR)   # PCA
    library(kernlab)      # SVM
    library(knitr)        # kable
    library(randomForest)
    library(reshape2)     # melt
    library(stringr)      # str_sub and str_replace
    library(tidyverse)    # purrr, reader, tibble, ...
```
***

# 5 Sections pour ce projet
  + Section 1 - [Sélection et Téléchargement des fichiers TCGA](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_01_V01.md)
  + Section 2 - [Pré-Traitement des fichiers (Regroupement des fichiers comptes)](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_02_V01.md)
  + Section 3 - [Pré-Traitement des fichiers (Filtrage, Découpe, Normalisation, Corrélation) ](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md)
    + Lien Rapide: [Filtrage des données](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#filtrage-des-g%C3%A8nes-qui-ont-plus-de-95-de-0-comme-compte)
    + Lien Rapide: [Colonne d'idenfication du Type de l’échantillon](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#ajouter-une-colonne-pour-identifier-le-type-de-l%C3%A9chantillon)
    + Lien Rapide: [Découpage Training / Test](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#d%C3%A9coupage-des-donn%C3%A9es-en-80-training-et-20-test)
    + Lien Rapide: [Normalisation VST](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#normalisation-avec-vst-de-deseq2)
    + Lien Rapide: [Filtrage des ind. corrélation](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#filtrage-des-%C3%A9chantillons-qui-ont-un-coefficient-de-corr%C3%A9lation-de-plus-de-98)
    + Lien Rapide: [Compte final des jeux de donnés](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#compte-final)
    + Lien Rapide: [Filtrage et normalisation VST du jeu de test normale](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#filtrage-et-normalisation-vst-du-jeu-de-test-normale)
  + Section 4 - [Entrainement](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_04_V01.md)
    + Lien Rapide: [Comparaison des performances de 6 algorithmes](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_04_V01.md#comparaison-des-performances-de-6-algorithmes-de-ml-avec-des-diff%C3%A9rentes-normalisations-des-donn%C3%A9es)
    + Lien Rapide: [BoxPlot des resultats de test de performance et choix](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_04_V01.md#un-beau-boxplot-plus-gentil-pour-loeil-que-des-colonnes-de-chiffres-pour-les-resultats)
    + Lien Rapide: [Test pour tenter d'améliorer la précision](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_04_V01.md#essayons-dam%C3%A9liorer-la-pr%C3%A9cision-du-mod%C3%A8le-fitrf_vst-pour-avoir-une-pr%C3%A9cision-de-plus-que-94-avec-les-param%C3%A8tres-traincontrol-de-caret)
    + Lien Rapide: [Resultats et conclusion de ces tests d'amélioration de la précision](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_04_V01.md#agr%C3%A9gation-des-r%C3%A9sultats-pour-comparaison-des-2-tests-dessai-daugmentation-de-pr%C3%A9cision)
  + Section 5 - [Prédictions et Résultats](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_05_V01.md)
    + Lien Rapide: [Prédictions avec le modèle SVM fitsvmLinear_VST et le jeu de Tests](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_05_V01.md#maintenant-faisons-les-pr%C3%A9dictions-tant-attendues-du-mod%C3%A8le-svm-fitsvmlinear_vst-et-le-jeu-de-tests)
    + Lien Rapide: [Prédictions du modèle fitrf_VST (Random Forest) et le jeu de Tests](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_05_V01.md#maintenant-faisons-les-pr%C3%A9dictions-du-mod%C3%A8le-fitrf_vst-random-forest-et-le-jeu-de-tests)
    + Lien Rapide: [Prédictions avec le modèle SVM fitsvmLinear_VST et le jeu de donnés Normal](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_05_V01.md#maintenant-faisons-les-pr%C3%A9dictions-du-mod%C3%A8le-svm-fitsvmlinear_vst-et-le-jeu-de-donn%C3%A9es-normal)
    + Lien Rapide: [Prédictions du modèle fitrf_VST (Random Forest) et le jeu de donnés Normal](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_05_V01.md#maintenant-faisons-les-pr%C3%A9dictions-du-mod%C3%A8le-fitrf_vst-random-forest-et-le-jeu-de-donn%C3%A9es-normal)

### Note:  Le jeu de données d’échantillons Normaux sera utilisé pour tester les prédictions des modèles. Ici le modèle fera des prédictions sur le type de cancer HNSC vs LUSC et non “Tumeur” vs “Normale”. Le jeu de données Normales comporte le type de cancer HNSC et LUSC dans son metadata.**
