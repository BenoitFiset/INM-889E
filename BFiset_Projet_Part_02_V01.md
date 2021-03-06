---
title: 'Projet: TCGA HNSC LUSC ML - Partie #2'
author: "Benoit Fiset"
date: "27/04/2021"
---

# Classification par apprentissage automatique du type de cancer à partir de données d'expression génétique de séquençage RNA-Seq

***

## Sur le serveur HPC:

### But de la section: Regroupement des fichiers individuels en une matrice unifié de comptes - Tumeurs et Normaux

***
```
library(knitr)      #kable
library(stringr)    #str_sub and str_replace
library(tidyverse)  #purrr, reader, tibble, ...
```

```r
set.seed(1234) # Important pour toujours avoir les memes resultats
```

#### Ouverture du ficher "Sample Sheet" qui contient les métas données des fichiers télécharges:

```r
TCGA_SampleDict <- read.table("gdc_sample_sheet.2021-04-26.tsv", sep ="\t", header=TRUE, stringsAsFactors = FALSE )

head(TCGA_SampleDict)
```

![](figures/TCGA_SampleDict_01.png)

#### Nous pouvons voir les comptes des échantillions par type:

```r
print(table(TCGA_SampleDict$Sample.Type))
```
```
         Metastatic       Primary Tumor Solid Tissue Normal 
                  2                1002                  93 
```

#### Création d’une nouvelle colonne « BarCode » qui contiendra la concaténation du type de cancer et le nom de l’échantillon.

```r
# Create new column "BarCode" Replace TCGA in Sample.ID with TCGA Project name in new column
TCGA_SampleDict$BarCode <- str_replace(TCGA_SampleDict$Sample.ID,"TCGA",str_sub(TCGA_SampleDict$Project.ID,-4))
# make name Column name data.frame friendly
TCGA_SampleDict$BarCode <- str_replace_all(TCGA_SampleDict$BarCode,"-",".")
```

#### Filtrer les fichiers par type d’échantillon. Ne garde pas les 2 "Metastatic"

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
```
> TCGA_Barcodes
# A tibble: 1,002 x 2
   Filename                                             BarCode         
   <chr>                                                <chr>           
 1 d1489e83-a47f-4d2e-b330-c0698b1ccc27.htseq.counts.gz HNSC.CV.A45V.01A
 2 f38a09f0-50e3-45c9-92a2-18e7639bfd6f.htseq.counts.gz HNSC.H7.7774.01A
 3 352b68d8-a9a2-4210-9f28-071b777654c0.htseq.counts.gz HNSC.CQ.A4CG.01A
 4 537d05c7-a7d0-48ab-aae0-90fde055a277.htseq.counts.gz HNSC.CV.5970.01A
 5 3bdbc229-cbc8-401b-836c-90e76ff3866b.htseq.counts.gz LUSC.94.7943.01A
 6 db2ba932-0f99-4e43-ab09-355bed4651e7.htseq.counts.gz LUSC.68.8251.01A
 7 211f4551-8720-4f0d-9c59-88108fcc037f.htseq.counts.gz HNSC.HD.8635.01A
 8 1aebe8a0-7fcb-407d-aa21-55e386f7297e.htseq.counts.gz HNSC.CV.7250.01A
 9 728d636b-1d0c-4835-9a06-9d84c6ae21a8.htseq.counts.gz HNSC.CR.6487.01A
10 84742003-9470-44dd-8fbf-c1a73468b8d8.htseq.counts.gz HNSC.CN.6018.01A
# … with 992 more rows
```

#### Lecture et faire une liste de tous les fichiers possiblement à traiter:

```r
locationOfFiles <- "./UnzippedFiles"

f_files<- list.files(locationOfFiles, pattern = "*.tsv", full.names = T)
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

#### Les fichiers de comptes génétique sont en format .tsv et le nom du fichier dans le « Sample Sheet » est en format .gz donc pour pouvoir retrouver le ficher dans le dictionnaire de référence « File.Name et BarCode » doit créer une référence des fichiers lus sur le disque et avoir la même référence .tsv - .gz

```r
originalFilename <- data.frame(FilenameOnDisk=f_files,FileName=str_replace(str_replace(f_files,paste(locationOfFiles,"/",sep=""),""),".tsv",".gz"))
```

```
kable(originalFilename[1:5,1:2],"simple")

FilenameOnDisk                                                          FileName                                              
----------------------------------------------------------------------  ------------------------------------------------------
./UnzippedFiles/00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.tsv   00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.gz 
./UnzippedFiles/0056ef0b-ecb9-443e-a694-1c0797a21e5c.htseq.counts.tsv   0056ef0b-ecb9-443e-a694-1c0797a21e5c.htseq.counts.gz 
./UnzippedFiles/0063c2da-4f47-4974-afea-03b4cdb21e33.htseq.counts.tsv   0063c2da-4f47-4974-afea-03b4cdb21e33.htseq.counts.gz 
./UnzippedFiles/007d074e-3b7c-484b-b757-1a60e1172506.htseq.counts.tsv   007d074e-3b7c-484b-b757-1a60e1172506.htseq.counts.gz 
./UnzippedFiles/00bed789-1573-4f1a-bc93-4ace31a32a6a.htseq.counts.tsv   00bed789-1573-4f1a-bc93-4ace31a32a6a.htseq.counts.gz 
```


```r
# Filter the files that only the files that are to be used are selected. 
# Meaning use only files in Barcodes VS all files on disk.
fileListFiltered <- subset(originalFilename, (originalFilename$FileName %in% TCGA_Barcodes$Filename))

# Factors may have empty levels after subsetting; unused levels are not automatically removed so ...
fileListFiltered  <- droplevels(fileListFiltered )

# New filtered file list to be used as input to read .tsv files.
fileListInput <- as.vector(fileListFiltered[["FilenameOnDisk"]])
```

#### Le filtre a enlevé 95 fichiers (Normale et Métastases) … donc Il y 1002 fichiers de Tumeurs pour travailler avec les modèles de ML de ce projet.
```
str(originalFilename)
'data.frame':	1097 obs. of  2 variables:
 $ FilenameOnDisk: Factor w/ 1097 levels "./UnzippedFiles/00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.tsv",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ FileName      : Factor w/ 1097 levels "00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.gz",..: 1 2 3 4 5 6 7 8 9 10 ...

str(fileListFiltered)
'data.frame':	1002 obs. of  2 variables:
 $ FilenameOnDisk: Factor w/ 1002 levels "./UnzippedFiles/00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.tsv",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ FileName      : Factor w/ 1002 levels "00068002-f4f0-4610-bfe6-67169c760d21.htseq.counts.gz",..: 1 2 3 4 5 6 7 8 9 10 ...
```

#### Cette fonction fait la magie d’extraire le « BarCode » qui sera le nom de la colonne pour identifier l’échantillon, lire les fichiers .tsv de comptes de gènes et enlever les lignes non utilisées. Cette fonction retourne un data.frame avec 2 colonnes : Les noms des gènes dans une colonne et BarCode dans l’autre qui contient le compte (nombre) de gènes. 


```r
locationOfFiles <- "./UnzippedFiles"

read_in_feature_counts<- function(file){
  cleanFilename <- str_replace(str_replace(file,paste(locationOfFiles,"/",sep=""),""),".tsv",".gz")
  Barcode_Name <- toString(TCGA_Barcodes[TCGA_Barcodes$Filename==cleanFilename,2])
  # comment parameter = remove last 5 lines of HTSeqCount which stat with --
  cnt<- read_tsv(file, col_names = c("Genes",Barcode_Name), comment = c("__"))  
  return(cnt)
}
```

#### Exécution de la fonction read_in_feature_counts pour tous les fichiers Tumeurs à lire qui seront dans une structure readr. Ensuite purrr::reduce est utilisé pour « joindre » tous les comptes dans un tibble qui aura qu’une colonne pour les noms des gènes et le restant des colonnes, identifiées par le « BarCode », pour TOUS les comptes des gènes. 

```r
# read the Count Files and create a readr "dataframe list" with the read_in_feature_counts function
rawCounts<- map(fileListInput, read_in_feature_counts)

 # Joint all the count file in a single tibble (dataframe) merged on the Gene column
rawCounts.df<- purrr::reduce(rawCounts, left_join, by = "Genes")

head(str(rawCounts))
```

```
List of 1002
 $ : spec_tbl_df[,2] [60,483 × 2] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
  ..$ Genes           : chr [1:60483] "ENSG00000000003.13" "ENSG00000000005.5" "ENSG00000000419.11" "ENSG00000000457.12" ...
  ..$ LUSC.56.7823.01B: num [1:60483] 1988 1 1194 449 748 ...
  ..- attr(*, "spec")=
  .. .. cols(
  .. ..   Genes = col_character(),
  .. ..   LUSC.56.7823.01B = col_double()
  .. .. )
 $ : spec_tbl_df[,2] [60,483 × 2] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
  ..$ Genes           : chr [1:60483] "ENSG00000000003.13" "ENSG00000000005.5" "ENSG00000000419.11" "ENSG00000000457.12" ...
  ..$ LUSC.33.4538.01A: num [1:60483] 4640 2 3645 1189 1348 ...
  ..- attr(*, "spec")=
  .. .. cols(
  .. ..   Genes = col_character(),
  .. ..   LUSC.33.4538.01A = col_double()
  .. .. )
  ```

#### Ici nous pouvons voir que le tibble (data.fram/matrice) de Genes et tous les comptes, par échantillon est créé. La dimension du tibble est de 60483 lignes (Genes) par 1003 colonnes (Mais 1002 Échantillions):

```r
head(rawCounts.df)
```
```
# A tibble: 60,483 x 1,003
   Genes     LUSC.56.7823.01B LUSC.33.4538.01A LUSC.34.2608.01A LUSC.XC.AA0X.01A
   <chr>                <dbl>            <dbl>            <dbl>            <dbl>
 1 ENSG0000…             1988             4640             4641             2193
 2 ENSG0000…                1                2                0                0
 3 ENSG0000…             1194             3645             2478             1702
 4 ENSG0000…              449             1189             1039              565
 5 ENSG0000…              748             1348              480              338
 6 ENSG0000…              194              443             1406             2927
 7 ENSG0000…             1515             1593             6691             2074
 8 ENSG0000…             2312             3488             4718             2663
 9 ENSG0000…            20363            27375            14515            12626
10 ENSG0000…             2346             2236             2028             1337
# … with 60,473 more rows, and 998 more variables: HNSC.BA.A4II.01A <dbl>,
#   HNSC.DQ.5629.01A <dbl>, LUSC.77.8153.01A <dbl>, HNSC.BB.4227.01A <dbl>,
#   HNSC.CN.A63U.01A <dbl>, LUSC.22.4613.01A <dbl>, LUSC.85.8277.01A <dbl>,
#   LUSC.60.2696.01A <dbl>, LUSC.18.3408.01A <dbl>, HNSC.CV.6935.01A <dbl>,
#   HNSC.CN.5366.01A <dbl>, LUSC.58.A46L.01A <dbl>, HNSC.IQ.A6SG.01A <dbl>,
```

#### Il ne reste qu'à sauvegarder ce fichier des comptes de gènes regroupés. Il a une conversion du tibble en data.frame pour sauvgarder en format .tsv

```r
write.table(as.data.frame(rawCounts.df),sep="\t","Tumor_Merged_2Lung_Counts.tsv", row.names = FALSE, quote = FALSE)
```

#### Les mêmes étapes ont été faites pour les fichiers de comptes Normaux (total de 93 échantillions) . Voici le code "directe" sans sortie ni explication.

```r
TCGA_Barcodes <- read_csv("TCGA_Barcode_Dict_Normal.tsv",col_names = c("Filename","BarCode"),skip =1 )
locationOfFiles <- "./UnzippedFiles"
f_files<- list.files(locationOfFiles, pattern = "*.tsv", full.names = T)
originalFilename <- data.frame(FilenameOnDisk=f_files,FileName=str_replace(str_replace(f_files,paste(locationOfFiles,"/",sep=""),""),".tsv",".gz"))
fileListFiltered <- subset(originalFilename, (originalFilename$FileName %in% TCGA_Barcodes$Filename))
fileListFiltered  <- droplevels(fileListFiltered )
fileListInput <- as.vector(fileListFiltered[["FilenameOnDisk"]])
rawCounts<- map(fileListInput, read_in_feature_counts)
rawCounts.df<- purrr::reduce(rawCounts, left_join, by = "Genes")
write.table(as.data.frame(rawCounts.df),sep="\t","Normal_Merged_2Lung_Counts.tsv", row.names = FALSE, quote = FALSE)
```

