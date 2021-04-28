---
title: "R Notebook"
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
```

```
str(rawCountData)
'data.frame':	60483 obs. of  1002 variables:
 $ LUSC.56.7823.01B  : int  1988 1 1194 449 748 194 1515 2312 20363 2346 ...
 $ LUSC.33.4538.01A  : int  4640 2 3645 1189 1348 443 1593 3488 27375 2236 ...
 $ LUSC.34.2608.01A  : int  4641 0 2478 1039 480 1406 6691 4718 14515 2028 ...
 $ LUSC.XC.AA0X.01A  : int  2193 0 1702 565 338 2927 2074 2663 12626 1337 ...
 $ HNSC.BA.A4II.01A  : int  1609 1 1474 672 413 350 1180 1965 7509 2171 ...
 
kable(rawCountData[1:5,498:502],"rst")

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
```

```
kable(percentZeros.df[order(-percentZeros.df$PercentZero)[2200:2210],],format = "rst")

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

![](BFiset_Projet_Part_03_V01_files/figure-html/unnamed-chunk-4-1.png)<!-- -->
