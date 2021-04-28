---
title: 'Projet: TCGA HNSC LUSC ML - Partie #3'
author: "Benoit Fiset"
date: "27/04/2021"
output:
  html_document: 
    keep_md: yes
---

# Classification par apprentissage automatique du type de cancer à partir de données d'expression génétique de séquençage RNA-Seq

***

## Sur le ordinateur personnel:

### But de la section: Filtrage, Découpe ("Training / Test sets"), Normalisation, Indice de Corrélation - Tumeurs er Normale


```r
set.seed(1234) # Important pour tjrs avoir les memes resultats
```

#### Ouverture du ficher qui regroupe tous les Tumeurs des comptes de gènes

```r
Filename <- "Tumor_Merged_2Lung_Counts.tsv"

rawCountData <- read.table(Filename, header=TRUE, row.names=1)

head(str(rawCountData))
```

```
'data.frame':	60483 obs. of  1002 variables:
 $ LUSC.56.7823.01B  : int  1988 1 1194 449 748 194 1515 2312 20363 2346 ...
 $ LUSC.33.4538.01A  : int  4640 2 3645 1189 1348 443 1593 3488 27375 2236 ...
 $ LUSC.34.2608.01A  : int  4641 0 2478 1039 480 1406 6691 4718 14515 2028 ...
 $ LUSC.XC.AA0X.01A  : int  2193 0 1702 565 338 2927 2074 2663 12626 1337 ...
 $ HNSC.BA.A4II.01A  : int  1609 1 1474 672 413 350 1180 1965 7509 2171 ...
```


```r
kable(rawCountData[1:5,498:502],"rst")
```

```
==================  ================  ================  ================  ================  ================
\                   HNSC.CV.A6JU.01A  HNSC.CR.5247.01A  HNSC.D6.6827.01A  LUSC.52.7622.01A  LUSC.60.2724.01A
==================  ================  ================  ================  ================  ================
ENSG00000000003.13              2984              1440              5457              4092              5252
ENSG00000000005.5                  2                 0                 0                 0                 1
ENSG00000000419.11              2548              2041              2258              2754              3139
ENSG00000000457.12               228               241              3026              1241               522
ENSG00000000460.15               229               297              3479               820               577
==================  ================  ================  ================  ================  ================
```


```r
percentFilterZero = 95

 # If Gene Count = 0 add 1 them sum all * 100 then / by number of columns
percentZerosCount <- 100*apply(rawCountData == 0, 1, sum) / ncol(rawCountData) 

percentZeros.df <- data.frame(PercentZero=percentZerosCount, TotalZero=apply(rawCountData == 0, 1, sum))

kable(percentZeros.df[order(-percentZeros.df$PercentZero)[2200:2210],],format = "rst")
```

```
=================  ===========  =========
\                  PercentZero  TotalZero
=================  ===========  =========
ENSGR0000277120.3     100.0000       1002
ENSGR0000280767.1     100.0000       1002
ENSGR0000281849.1     100.0000       1002
ENSG00000172288.7      99.9002       1001
ENSG00000183336.7      99.9002       1001
ENSG00000189393.5      99.9002       1001
ENSG00000199231.1      99.9002       1001
ENSG00000199311.1      99.9002       1001
ENSG00000199453.1      99.9002       1001
ENSG00000199782.1      99.9002       1001
ENSG00000199870.1      99.9002       1001
=================  ===========  =========
```


```r
ggplot(data=percentZeros.df,aes(x=percentZerosCount)) +  geom_histogram(binwidth=5,color="darkblue", fill="lightblue") +
  labs(title=paste("Filter",sum(percentZerosCount >= percentFilterZero)," genes with more than",percentFilterZero,"% of zero counts"),
       x="Percent count of samples that are Zero", y = "Number of genes") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = percentFilterZero, color = "red")
```

![](figures/Plot_95_01.png)


```r
# Filter Genes that have less than 1% of counts
geneMinCount = round((0.01 * ncol(rawCountData)))

# Find the Min Gene count of that Gene
minGeneCount <- apply(rawCountData, 1, min) 

# Create a TRUE / FALSE where TRUE is to remove that Gene
toFilterGenes <- (percentZerosCount > percentFilterZero) | (minGeneCount < geneMinCount) 
table(toFilterGenes)
```

```
toFilterGenes
FALSE  TRUE 
11061 49422 
```


```r
# Reverse the TRUE / FALSE to put TRUE to keep the good Genes
filteredGenes <- !toFilterGenes 
table(filteredGenes)
```

```
filteredGenes
FALSE  TRUE 
49422 11061 
```


```r
dim(rawCountData)
```
```
[1] 60483  1002
```


```r
# This Filters the genes taht don't pass the more tnat 95% 0 and have less that 1% of counts.
cleanCountData <- rawCountData[filteredGenes, ]

dim(cleanCountData)
```

```
[1] 11061  1002
```


```r
# Function to return a substring, here the 4 first chars of the string sent as parameter (Ex: HNSC, LUSC)
substrColName = function(x){ substr(x,1,4) }

cleanCountDataMat <- as.matrix(cleanCountData)
kable(cleanCountDataMat[495:505,1:6],"rst")
```
```
==================  ================  ================  ================  ================  ================  ================
\                   LUSC.56.7823.01B  LUSC.33.4538.01A  LUSC.34.2608.01A  LUSC.XC.AA0X.01A  HNSC.BA.A4II.01A  HNSC.DQ.5629.01A
==================  ================  ================  ================  ================  ================  ================
ENSG00000051596.8               1063              1251              1551               810               750               817
ENSG00000051620.9               3886              3901              4355              3162              3095              2947
ENSG00000051825.13              1041              1741               909               664               975               789
ENSG00000052126.13              2672              3568              2285              2296              1529              3881
ENSG00000052723.10               903              1805              2416              1354              3144              1256
ENSG00000052749.12              1395              1850              2159              3150              4810              3735
ENSG00000052795.11               708               870              2507              2325              2472               931
ENSG00000052802.11              2279              4071              6115              3916              7638              3404
ENSG00000052841.13              2298              5045              4052              2555              1744              1955
ENSG00000053254.14              4926              2690              6331              3936              3270              3891
ENSG00000053371.11              1491               613              2348              1485              1817              1516
==================  ================  ================  ================  ================  ================  ================
```


```r
df.data <- data.frame(Type=apply(as.matrix(rownames(t(cleanCountDataMat))),1,substrColName), t(cleanCountDataMat) )
kable(df.data[495:505,1:6],"rst")
```

```
================  ====  ==================  ==================  ==================  ==================  ==================
\                 Type  ENSG00000000003.13  ENSG00000000419.11  ENSG00000000457.12  ENSG00000000460.15  ENSG00000000938.11
================  ====  ==================  ==================  ==================  ==================  ==================
LUSC.22.4604.01A  LUSC                5487                1872                 919                 730                2212
LUSC.22.0940.01A  LUSC                2668                2199                1083                1421                1312
LUSC.58.A46K.01A  LUSC                1436                1974                 495                 448                 207
HNSC.CV.A6JU.01A  HNSC                2984                2548                 228                 229                 313
HNSC.CR.5247.01A  HNSC                1440                2041                 241                 297                  90
HNSC.D6.6827.01A  HNSC                5457                2258                3026                3479                 404
LUSC.52.7622.01A  LUSC                4092                2754                1241                 820                2616
LUSC.60.2724.01A  LUSC                5252                3139                 522                 577                 547
HNSC.CV.A465.01A  HNSC                1295                3827                 372                 299                 129
LUSC.21.1070.01A  LUSC                2974                2882                1097                1123                 846
LUSC.66.2781.01A  LUSC                4255                5514                1106                 697                1273
================  ====  ==================  ==================  ==================  ==================  ==================
```

***

## Decoupage de des donnees en 80% Training et 20% Test


```r
# Random Select percentForTraining for training DataSet
trainingDataSelect <- sample(1:nrow(df.data), size = percentForTraining * nrow(df.data))
```



```r
# Select percentForTraining for training DataSet
trainingDataset.df <- df.data[trainingDataSelect ,] 
dim(trainingDataset.df)
kable(trainingDataset.df[495:505,1:6],"rst")

# Select rest of percentForTraining for test DataSet
testDataset.df <- df.data[-trainingDataSelect ,].
dim(testDataset.df)
kable(testDataset.df[107:117,1:6],"rst")
```


## Normalisation avec VST de DESeq2

```r
condition <- condition <- factor(trainingDataset.df$Type)

# Transpose and Convert Training Data Frame to Matrix without the Type column
tempDDSdftoMat <- as.matrix(t(trainingDataset.df[,-1]))
coldata <- data.frame(row.names=colnames(tempDDSdftoMat), condition)

dds <- DESeqDataSetFromMatrix(countData=tempDDSdftoMat, colData=coldata, design=~condition)

data_VST <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
```

