# INM-889E
Repo pour le cours INM 889E - Projet de ML

***

# Classification par apprentissage automatique du type de cancer à partir de données d'expression génétique de séquençage RNA-Seq

**Note:  Le jeu de donnés d’échantillons Normale sera utilisé pour tester les prédictions des modèle. Ici le modèle fera des prédictions sur le type de cancer HNSC vs LUSC et non “Tumeur” vs “Normale”. Le jeu de donnés Normale comporte le type de cancer HNSC et LUSC dans son metadata.**

* Vu la taille des fichier et le temps de calcul nécessaire une partie du traitement des fichiers ont été effectués sur des serveurs HPC de Calcul Québec (Béluga):
  + Le téléchargement des fichiers incluant le regroupement et décompression
  + Le regroupement des fichiers individuels en une matrice unifié de comptes - Tumeurs
  + Le regroupement des fichiers individuels en une matrice unifié de comptes - Normale
  + Entrainement des modèles
  + “Tuning" des modèles
  + Prediction des modèles 

* Sur le l’ordinateur local:
  + Le filtrage des genes qui ont comportent plus de 95% de comptes de valeur 0
  + Le filtrage des genes qui ont moins de 1% de comptes
  + La creation de le colonne qui sera la “Classe” du projet. Se nomme “Type”
  + Le découpage du jeu de donnes “Training”  80% et “Test” 20%. 
  + Normalisation des données “Training” avec VST de DESeq2
  + Normalisation des données “Test” avec VST de DESeq2
  + Normalisation des données “Normale” avec VST de DESeq2
  + Filtrer les genes qui ont un indice de correlation de plus de 98% (Tumeur et Normale)

***

# 6 Sections pour ce projet
  + Section 1 - [Sélection et Téléchargement des fichiers TCGA](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_01_V01.md)
  + Section 2 - [Pré-Traitement des fichiers (Regroupement des fichiers comptes)](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_02_V01.md)
  + Section 3 - [Pré-Traitement des fichiers (Filtrage, Découpe, Normalisation, Corrélation) ](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md)
    + Lien Rapide: [Filtrage des données](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#filtrage-des-g%C3%A8nes-qui-on-plus-de-95-de-0-comme-compte)
    + Lien Rapide: [Colonne d'idenfication du Type de l’échantillon](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#ajouter-dune-colonne-pour-identifier-le-type-de-l%C3%A9chantillon)
    + Lien Rapide: [Découpage Training / Test](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#decoupage-de-des-donnees-en-80-training-et-20-test)
    + Lien Rapide: [Normalisation VST](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#normalisation-avec-vst-de-deseq2)
    + Lien Rapide: [Filtrage des ind. corrélation](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#filtrage-des-%C3%A9chantillons-qui-ont-un-coefficient-de-corr%C3%A9lation-de-plus-de-98)
    + Lien Rapide: [Compte final des jeu de donnés](https://github.com/BenoitFiset/INM-889E/blob/main/BFiset_Projet_Part_03_V01.md#voici-le-compte-final-des-%C3%A9chantillons-qui-serviront-pour-lentrainement-et-les-tests-des-mod%C3%A8les-du-projet-post-tous-les-filtrages)
    + Bonus - [Autres Normalistions: Aucune et ln()]
  + Section 4 - [Entrainement]
    + Bonus - [Comparaison des performances de 6 algorithmes de ML avec les différentes normalisations des données]
  + Section 5 - [Prédictions]
  + Section 6 - [Résultats]

***

# Autre Section
