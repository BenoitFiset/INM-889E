---
title: 'Projet: TCGA HNSC LUSC ML - Partie #2'
author: "Benoit Fiset"
date: "27/04/2021"
output:
  html_document: 
    keep_md: yes
---

<!-- if you want to keep only markdown and get GitHub-flavored markdown. In that case, make your YAML look like this: output: github_document -->




<center>
# Classification par apprentissage automatique du type de cancer à partir de données d'expression génétique de séquençage RNA-Seq
</center>

**Note:  Le jeu de donnés d’échantillons Normale sera utilisé pour tester les prédictions des modèle. Ici le modèle fera des prédictions sur le type de cancer HNSC vs LUSC et non “Tumeur” vs “Normale”. Le jeu de donnés Normale comporte le type de cancer HNSC et LUSC dans son metadata.**

* Vu la taille des fichier et le temps de calcul nécessaire une partie du traitement des fichiers ont été effectués sur des serveurs HPC de Calcul Québec (Béluga):
  + Le téléchargement des fichiers incluant le regroupement et décompression
  + Le regroupement des fichiers individuels en une matrice unifié de comptes - Tumeurs
  + Le regroupement des fichiers individuels en une matrice unifié de comptes - Normale
  + Entrainement des modèles
  + “Tuning: des modèles
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
## Sur le serveur HPC:


```r
set.seed(1234) # Important pour tjrs avoir les memes resultats
```

## Le regroupement des fichiers individuels en une matrice unifié de comptes - Tumeurs er Normale


#### Ouverture du ficher "Sample Sheet" qui contient les métas données des fichiers télécharges:

```r
TCGA_SampleDict <- read.table("gdc_sample_sheet.2021-04-26.tsv", sep ="\t", header=TRUE, stringsAsFactors = FALSE )
```

#### Nous pouvons voir les compte des échantillions par type:

```r
table(TCGA_SampleDict$Sample.Type)
```

```
## 
##          Metastatic       Primary Tumor Solid Tissue Normal 
##                   2                1002                  93
```

#### Création d’une nouvelle colonne « BarCode » qui contiendra la concaténation du type de cancer et le nom de l’échantillon.

```r
# Create new column "BarCode" Replace TCGA in Sample.ID with TCGA Project name in new column
TCGA_SampleDict$BarCode <- str_replace(TCGA_SampleDict$Sample.ID,"TCGA",str_sub(TCGA_SampleDict$Project.ID,-4))
# make name Column name data.frame friendly
TCGA_SampleDict$BarCode <- str_replace_all(TCGA_SampleDict$BarCode,"-",".")
```

#### Filtrer les fichier par type d’échantillon:

```r
# Extract the data subsets by type.
normalSubset <-subset(TCGA_SampleDict[,c("File.Name","BarCode")],TCGA_SampleDict$Sample.Type == "Solid Tissue Normal")
tumorSubset <-subset(TCGA_SampleDict[,c("File.Name","BarCode")],TCGA_SampleDict$Sample.Type == "Primary Tumor")
allDataset  <- TCGA_SampleDict[,c("File.Name","BarCode")]
```

#### Exemple de contenu suite au filtre:

```r
kable(tumorSubset[1:5,])
```
```
File.Name	                                            BarCode
d1489e83-a47f-4d2e-b330-c0698b1ccc27.htseq.counts.gz	HNSC.CV.A45V.01A
f38a09f0-50e3-45c9-92a2-18e7639bfd6f.htseq.counts.gz	HNSC.H7.7774.01A
352b68d8-a9a2-4210-9f28-071b777654c0.htseq.counts.gz	HNSC.CQ.A4CG.01A
537d05c7-a7d0-48ab-aae0-90fde055a277.htseq.counts.gz	HNSC.CV.5970.01A
3bdbc229-cbc8-401b-836c-90e76ff3866b.htseq.counts.gz	LUSC.94.7943.01A
```


***
#### Sauvegarde des fichiers d’échantillons par type:


```r
write.table(normalSubset,sep=",","TCGA_Barcode_Dict_Normal.tsv", row.names = FALSE, quote = FALSE)
write.table(tumorSubset,sep=",","TCGA_Barcode_Dict_Tumor.tsv", row.names = FALSE, quote = FALSE) 
write.table(allDataset,sep=",","TCGA_Barcode_Dict_All2Lung.tsv", row.names = FALSE, quote = FALSE) 
```

#### Lecture du ficher Tumeur pour la suite des opérations de regroupement des fichiers individuels Tumeur en une matrice unifié de comptes Tumeur

```r
TCGA_Barcodes <- read_csv("TCGA_Barcode_Dict_Tumor.tsv",col_names = c("Filename","BarCode"),skip =1 )
```

#### Lecture et faire un liste de tous les fichiers possiblement à traiter:

```r
locationOfFiles <- "./UnzippedFiles"

f_files<- list.files(paste(locationOfFiles,sep = "/"), pattern = "*.tsv", full.names = T)
```

```r
head(f_files)
```
```
[1] "./UnzippedFiles/00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.tsv"
[2] "./UnzippedFiles/0056ef0b-ecb9-443e-a694-1c0797a21e5c.htseq.counts.tsv"
[3] "./UnzippedFiles/0063c2da-4f47-4974-afea-03b4cdb21e33.htseq.counts.tsv"
[4] "./UnzippedFiles/007d074e-3b7c-484b-b757-1a60e1172506.htseq.counts.tsv"
[5] "./UnzippedFiles/00bed789-1573-4f1a-bc93-4ace31a32a6a.htseq.counts.tsv"
[6] "./UnzippedFiles/00f32c8e-2a91-485b-99c9-6e2e22addec1.htseq.counts.tsv"
```

#### Les fichier de comptes génétique sont en format .tsv et le nom du fichier dans le « Sample Sheet » est en format .gz donc pour pouvoir retrouver le ficher dans le dictionnaire de référence « File.Name et BarCode » doit créer un référence des fichiers lu sur le disque et avoir la même référence .tsv - .gz

```r
originalFilename <- data.frame(FilenameOnDisk=f_files,FileName=str_replace(str_replace(f_files,locationOfFiles,""),".tsv",".gz"))
```

```
kable(originalFilename[1:5,1:2],"rst")

FilenameOnDisk                                                          FileName                                              
----------------------------------------------------------------------  ------------------------------------------------------
./UnzippedFiles/00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.tsv   /00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.gz 
./UnzippedFiles/0056ef0b-ecb9-443e-a694-1c0797a21e5c.htseq.counts.tsv   /0056ef0b-ecb9-443e-a694-1c0797a21e5c.htseq.counts.gz 
./UnzippedFiles/0063c2da-4f47-4974-afea-03b4cdb21e33.htseq.counts.tsv   /0063c2da-4f47-4974-afea-03b4cdb21e33.htseq.counts.gz 
./UnzippedFiles/007d074e-3b7c-484b-b757-1a60e1172506.htseq.counts.tsv   /007d074e-3b7c-484b-b757-1a60e1172506.htseq.counts.gz 
./UnzippedFiles/00bed789-1573-4f1a-bc93-4ace31a32a6a.htseq.counts.tsv   /00bed789-1573-4f1a-bc93-4ace31a32a6a.htseq.counts.gz 
```


<!-- # ```{r pressure, echo=FALSE} -->
<!-- # plot(pressure) -->
<!-- # ``` -->
<!-- Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot. -->
