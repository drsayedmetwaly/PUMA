<br> <br>

Introduction
============

### What and Why?

**rMetabolomics** is an R Markdown project I made specifically to
perform data processing and statistical analysis of metabolomics data.
Most of the popular statistical software packages perform data
procssing/analysis on a single metabolite at a time with limited
capability to do batch processing/analysis thus they are of limited
benefit when handling metabolomics data. There are some online tools
dedicated for the purpose of metabolomics data analysis but issues with
internet connectivity and service run down during  
maintenance times limit their reliability.

### Requirements

I programmed **rMetabolomics** to be fully automated, the only required
input is the name of the dataset file. Right now only CSV files are
accepted. Dataset must comply with the generally accepted metabolomics
data configuration: metabolites in columns, samples in rows and first 2
columns contain sample ID and class respectively. Currently
rMetabolomics can process 2 classes only (e.g. disease vs control,
survivors vs non- survivors … etc.).

### Structure

The code of **rMetabolomics** is organized in 3 parts: starting, data
processing and data analysis. Each part is composed of sections called
“modules” which carry out specific tasks. Modules are connected and
depend on each other. I added comments to explain what is the task of
the module and how the code works.

### Future plans

I am planning to include analysis of variance (ANOVA), principal
component analysis (PCA), partial least square discriminant analysis
(PLS-DA) and orthogonal partial least square discriminant analysis
(oPLS-DA) are in future versions of rMetabolomics. Importantly, I am
planning to make next versions able to automatically run the optimal
statistical test based on the data entered.

------------------------------------------------------------------------

<br> <br>

Part1: Let’s get started!
=========================

This part does initialization and loads necessary custom functions.

<br>

### 1.1 Initialization module

Here we specify the dataset file. This module does not produce any
output, it just loads the dataset and required packages.

``` r
#Initialization Module: loads dataset and the required packages.

#Clear R environment.
rm(list=ls(all=TRUE))

#load the dataset from a csv file. 
#Note: CSV file layout must comply with generally acceptable metabolomics dataset structure (1st column: samples ID, 2nd column: class/sample type, samples in rows and metabolites in columns) 
rawDataset <-read.csv("C:/Users/metwa/Medicine/Study/Canada/PhD CCM @UCalgary/Courses, Workshops, Seminars/11 R Wizardry/Final Project/example_dataset_oxylipins_raw.csv", header=TRUE, sep=",", check.names=FALSE, stringsAsFactors=TRUE, comment.char="")

#Load required packages
library("ggplot2")
library("reshape2")
library("plyr")
library("scales")
library("gridExtra")
library("grid")
library("gtable")
library("pwr")
library("plot3D")
library("rgl")
library("plot3Drgl")
```

<br>

### 1.2 Graphics module

This module contains 3 custom made graphics functions (matrixPlot,
matrixCompare and drawTable). Functions in this modules are used
throughout the rest of rMetabolomics.

``` r
#Graphics module contains custom  plotting functions used throughout the project.


#-------------------------------------------Function matrixPlot----------------------------------------------- 
# Custom function to do boxploting for an entire matrix
# Inputs: numeric matrix, title, x axis label, y axis label,boxplot orientation
matrixPlot<-function(data_mat , title="",x_lab="", y_lab="" , boxplot_orient="v" ){
 

meltData <- melt(data_mat) #converts data_mat into long format
x<-factor(meltData[,2])
y<-meltData[,3]

p <- ggplot(meltData, aes(x,y, fill=x ))#Initialize the plot


if (boxplot_orient=="h") { #Horizontal orientation of boxplots
  p + geom_boxplot(color="black",position = "dodge" )+ 
  coord_flip()+
  theme(legend.position="none", axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y  = element_text(angle=0, vjust=0.5, size=6,face="bold"))+
  labs(title=title,x=x_lab, y = y_lab) 
  
} else { #Vertical orientation of boxplots
 p + geom_boxplot(color="black",position = "dodge" )+ 
  theme(legend.position="none", axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x  = element_text(angle=90, hjust=0.5, size=6,face="bold"))+
  labs(title=title,x=x_lab, y = y_lab)   
}

}#------------------------------------------------------------------------------------------------------------


#-------------------------------------------Function matrixCompare-------------------------------------------- 
# Custom function to compare an the mean value of 2 matrices
# Inputs: 2 numeric matrices, 2 titles
matrixCompare<-function(mat1 , mat2, mat1_title="",mat2_title="" ){
 
# set graphics window to plot side by side
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))  
par(mfrow = c(1,2), mai=c(1,1,1,0))

#replace all non-numeric data with 0
mat1[which(is.numeric(mat1)==FALSE)] <- 0 
mat2[which(is.numeric(mat2)==FALSE)] <- 0

#Plot the denisty of the mean of each column (metabolite) for the 2 numeric matrices
plot(density(apply(mat1, 2, mean), na.rm = TRUE), main=mat1_title)
plot(density(apply(mat2, 2, mean), na.rm = TRUE), main=mat2_title)
  
}#------------------------------------------------------------------------------------------------------------



#-------------------------------------------Function drawTable------------------------------------------------ 
# Custom function to draw graphic tables
# Inputs: matrix, title, footer, title text size, footer font face
drawTable<-function(d , title_txt="", footer_txt="", title_txt_size=20,footer_fontface="plain" ){

table <- tableGrob(d) #Initialize the table

#Initialize table title and footer
title <- textGrob(title_txt,gp=gpar(fontsize=title_txt_size))
footnote <- textGrob(footer_txt, x=0, hjust=0,
                     gp=gpar( fontface=footer_fontface))

padding <- unit(1,"line")#Sets the distance between elements

#Draw table + title
table <- gtable_add_rows(table, 
                         heights = grobHeight(title) + padding,
                         pos = 0)
#Adds footnote
table <- gtable_add_rows(table, 
                         heights = grobHeight(footnote)+ padding)

#Sets the alignment of title and footer
table <- gtable_add_grob(table, list(title, footnote),
                         t=c(1, nrow(table)), l=c(1,1), 
                         r=ncol(table))

#Plot the complete table
grid.newpage()
grid.draw(table)
}#------------------------------------------------------------------------------------------------------------
```

------------------------------------------------------------------------

<br>

Part2: Data Processing
======================

Data processing consists of plotting, normalization, transformation and
scaling of data.

<br>

### 2.1 Plot Raw Data

After entering the dataset file name rMetabolomics will Load dataset and
plot raw data.

``` r
numericData<-data.matrix(rawDataset[,3:ncol(rawDataset)]) #Excludes the primary ID and sample type columns

matrixPlot(data_mat=numericData , title="Raw Data",x_lab="", y_lab="" , boxplot_orient="v")#Plot the numeric data
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/raw%20Data-1.png)

### 2.2 Normalization

Typically metabolomics samples show high variability in concentration,
and the derived metabolic profiles have a heteroscedastic noise
structure characterized by increasing variance as a function of
increased signal intensity. These sources of experimental and
instrumental noise substantially complicate information recovery when
statistical tools are used. Median fold change normalization is least
compromised by the biologically relevant changes in mixture components
and is thus preferable normalization method.

``` r
#-------------------------------------------Function MDFCN---------------------------------------------------- 
# Custom function to perform median fold change normalisation
# Inputs: numeric matrix
# Priniple: calculates median of the rows --> divide each cell by its row median --> calcuate columns median --> divide each cell by the column median

MDFCN <- function(mat) {

  #Replace any non-numeric values in the dataset with 0
  mat[which(is.numeric(mat)==FALSE)] <- 0 
  
  #Calculates median of the rows of the input matrix and store them in the vector medSam
  medSam <- apply(mat, 1, median)
  medSam[which(medSam==0)] <- 0.0001 #Replace any median 0 with 0.0001 (to avoid divide by zero later on)
  
  #Divide cells by the its respective row median, calculates  columns median and divide each cell by the column median
  mat <- apply(mat, 2, function(mat, medSam){
    medFDiSmpl <- mat/medSam
    vec<-mat/median(medFDiSmpl)
    return(vec)
  }, medSam)

  return (mat)
}#------------------------------------------------------------------------------------------------------------

numericData.normalized<-MDFCN(numericData) #Call MDFCN function to normalize numericData

#Plot the normalized data
matrixPlot(data_mat=numericData.normalized , title="Data After Median Fold Change Normalization (MDFCN)",x_lab="", y_lab="" , boxplot_orient="v")
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/Median%20Fold%20Change%20Normalizatio-1.png)

### 2.3 Data Transformation

Logarithmic transformation of data is usually required given the highly
skewed nature of the metabolomics results. **rMetabolomics** will do
base 10 logarithmic transformation and plot transformed data.

``` r
numericData.log<-log10(numericData.normalized)#perform log transformation to the data

#plot the transformed data
matrixPlot(data_mat=numericData.log , title="Data After Base 10 Log Transformation",x_lab="", y_lab="" , boxplot_orient="v")
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/Log_Transformation-1.png)

### 2.4 Data Scaling

Data scaling aims to centralize the data based on the mean of each
metabolite. It is commonly the final step in data processing before
statistical analysis is performed. **rMetabolomics** will do z scaling
and plot scaled data.

``` r
#----------------------------------------------Function zScale------------------------------------------------
# Calculates z_score, formula: (x-u)/sd
# inputs: a numeric matrix
zScale<-function (num_mat){
  
  #Replace any non-numeric values in the dataset with 0
  num_mat[which(is.numeric(num_mat)==FALSE)] <- 0 
  
  
  #Subtracts the column mean from each cell and divide the result by the column sd
  num_mat <- apply(num_mat, 2, function(num_mat){
  
  #Calculates mean of the columns of the input matrix and store them in the vector meanCol
  meanCol<-mean(num_mat)
  
  #Calculates sd of the columns of the input matrix and store them in the vector sdCol
  sdCol <- sd(num_mat)
  sdCol[which(sdCol==0)] <- 0.0001 #Replace any mean 0 with 0.0001 (to avoid divide by zero later on)
  
    num_mat <- num_mat-meanCol
    num_mat<-num_mat/sdCol
    return(num_mat)
  })

return(num_mat)
}#------------------------------------------------------------------------------------------------------------

numericData.z<-zScale(numericData.log)#Call zScale function to normalize numericData.log

#Plot z-scaled data
matrixPlot(data_mat=numericData.z , title="Data After z Scaling",x_lab="", y_lab="" , boxplot_orient="v")
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/zScalling-1.png)

``` r
processedDataset <- cbind(rawDataset[,1:2],numericData.z)#saves the final processed data: composed of first 2 columns of the raw data (sample ID and Class) and processed numeric matrix
```

<br>

### 2.5 Data Processing Summary

As a final step in data processing, **rMetabolomics** will provide an
overal visual reprsentation of the data before and after processing.

``` r
#Calls matrixCompare function to compare between raw data and processed data
matrixCompare(numericData , numericData.z, mat1_title="Raw Data",mat2_title="Data After Processing" )
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/before_after_plot-1.png)

------------------------------------------------------------------------

<br>

Part3: Statistical Analysis
===========================

In this part **rMetabolomics** will do univariate statistical analysis
of the processed data. Analysis consists of doing correlation matrix,
heatmap, t-test, sample size calculation and **“Power Array”**

<br>

### 3.1 Correlation Matrix

Correlation matrix provides an overview of the correlation between all
metabolites. Spearman correlation is preferentially selected given the
skewed nature of most metabolomics data. I programmed **rMetabolomics**
to automatically do a correlation matrix (Spearman correlation) for all
metabolites and plot the result.

``` r
#----------------------------------------Function corSpearman--------------------------------------------------
#Plots a Spearman Correlation Matrix using using qplot
#Inputs: processed numeric data matrix (normalized, transformed and scaled data)
corSpearman<-function (processed_mat){

#perform spearman correlation and convert the output matrix into the long format
data=melt(cor(processed_mat, method="spearman"))

#Plots the output matrix
qplot(x=data[,1], y=data[,2], data=data, fill=value, geom="tile") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.title.x=element_blank(), axis.title.y=element_blank())+ 
  scale_fill_gradient(low="#000433", high="#62c3f6")+
  guides(fill = guide_legend(title = "Spearman Correlation Coeff." ))
}#-------------------------------------------------------------------------------------------------------------

corSpearman(numericData.z)
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/corMatrix-1.png)

<br>

### 3.2 Heatmap

Heatmap provides an intuitive way to evaluate the dataset. I programmed
**rMetabolomics** to rearrange the samples based on the sample
type/class. Then sample type indicator is plotted in the first column of
the heatmap. As a general rule, the more a metabolite column shows
agreement with the first column (sample type indicator) in the heatmap,
the more it can be used to separate between sample classes.

``` r
#-----------------------------------------------Function ggHeatmap--------------------------------------------
#Plots a heatmap of the dataset
#Inputs: processed numeric data matrix (normalized, transformed and scaled data)
ggHeatmap<-function (processed_dataset){

gg_heatmap_mat<-processed_dataset
gg_heatmap_mat[,2]<-as.numeric(as.factor(gg_heatmap_mat[,2]))#convert sample type into factor then into number

#Arranges sample ID by the class/sample type and plots arranged sample IDs on Y axis 
gg_heatmap_mat[,1]<- reorder(gg_heatmap_mat[,1],gg_heatmap_mat[,2])

gg_heatmap_mat.melt <- melt(gg_heatmap_mat)#converts it into long format

#Rescale all variable (columns) values in the output matrix so that they lie between 0-100
gg_heatmap_mat.melt <- ddply(gg_heatmap_mat.melt, .(variable), transform,rescale = rescale(value, to=c(0,100)))

#Plot the result using geom_tile
(p <- ggplot(gg_heatmap_mat.melt, aes(variable, gg_heatmap_mat.melt[,1])) + 
    geom_tile(aes(fill = rescale),colour = "#000433") + 
    scale_fill_gradient(low = "#000433", high = "#62c3f6"))#sets the gradient colour range
base_size <- 9
p + theme_grey(base_size = base_size) + 
  labs(x = "", y = colnames(gg_heatmap_mat)[1]) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme( axis.text.x = element_text(size = base_size * 0.8, angle = 90, hjust = 1))+
  guides(fill =
           guide_legend(title = "Intensity %",#the title of the legend
             title.theme = element_text(
               angle = 0
             )
           )
  )
  

}#-------------------------------------------------------------------------------------------------------------

ggHeatmap(processedDataset)
```

    ## Using ID as id variables

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/heatmap-1.png)

<br>

### 3.3 Statistical Significance Test (t-test)

Next I programmed **rMetabolomics** to do a series of 2 sample t-tests
comparing the levels of each metabolite in both classes, perform
Bonferroni adjustment of p-values for family wise error rate FWER
correction, plot significant metabolites and prints a table with the
details.

``` r
#Work on the processed dataset (though it converts any non-numeric data into 0)
dat<-processedDataset
dat[which(is.numeric(dat)==FALSE)] <- 0

iniFlag<-FALSE #Flip-flop logical flag to mark the first time a condition occurs ;-))
for (metaCounter in 3:ncol(dat)){ #Here the loop starts from the 3rd column of the data (to avoid working on sample ID or class values)
#A saftey condition to run t-test calculation only if the column does not contain NaN or NA
if (is.nan(sum(dat[,metaCounter]))==FALSE & is.na(sum(dat[,metaCounter]))==FALSE){
t.table<-t.test(dat[,metaCounter] ~ dat[,2], data=dat[,c(2,metaCounter)])#performs t-test on column metaCounter of data (contains the metabolite values) by the group in column 2 (contains sample type/classes)
pval<-t.table$p.value
}else{pval<-1}#If we encounter NaN or NA -> set p-value to 1  
      if (iniFlag==FALSE){iniFlag<-TRUE;  raw_p_value.list<-c(pval)}else{
       raw_p_value.list<-c(raw_p_value.list, pval)#saves the raw p-values for all metabolites
      
    }
    
}

  
#applies Bonferroni family wise error rate FWER correction for multiple comparisons
bonf_p_value.list<-p.adjust(raw_p_value.list, method = "bonferroni", n = length(raw_p_value.list))


iniFlag<-FALSE #Flip-flop logical flag to mark the first time a condition occurs ;-))
for (metaCounter in 1:length(bonf_p_value.list)){

    if (bonf_p_value.list[metaCounter]<0.05){#select adjusted p-value <0.05 (significant p-values)
      if (iniFlag==FALSE){iniFlag<-TRUE;  sig.list<-c(metaCounter)}else{
       sig.list<-c(sig.list, metaCounter)#saves the positions of significant p-values after Bonferroni correction
      
      }
    } 
  
}


dataSubset<-matrix(rep(0,length(sig.list)*3), nrow=length(sig.list),ncol=3)#Initialize matrix
dataSubset[,1]<-colnames(dat)[sig.list]#Retrive names of significant metabolites
dataSubset[,2]<-round(raw_p_value.list[sig.list],7)#Retrive raw p-value of significant metabolites
dataSubset[,3]<-round(bonf_p_value.list[sig.list],7)#Retrive adjusted p-value of significant metabolites
colnames(dataSubset)<-c("Metabolite","Raw p-value","Bonferroni corrected p-value")

#calls darwTable to draw a table of the rearranged the matrix (by the values of the Bonferroni corrected p-values)
drawTable(dataSubset[order(as.numeric(dataSubset[,3])),] , title_txt="Significant Metabolites", footer_txt="Metabolites are arranged by the Bonferroni corrected p-values.", title_txt_size=22,footer_fontface="italic" )
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/t-test-1.png)

``` r
#Retrives the columns of significant metabolites + the sample class
dataSubset<-melt(dat[,c(2,sig.list)])
```

    ## Using Type as id variables

``` r
#Plots boxplots of the significant metabolites comparing the level of metabolite between the classes
ggplot(dataSubset, aes(dataSubset[,1], dataSubset[,3], fill=dataSubset[,2]))+
ggtitle("Significant Metabolities (Bonferroni corrected for FWER multiple comparisons)") + 
geom_boxplot() + 
facet_wrap(~dataSubset[,2], scales="free")+
  labs(x="", y="")+
  guides(fill =
           guide_legend(title = "Metabolite",
             title.theme = element_text(
               face = "bold",
               angle = 0
             )
           )
         )
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/t-test-2.png)

<br>

### 3.4 Sample Size Calculation

Next **rMetabolomics** will Calculate required sample size (n) per group
for the previously identified significant metabolites. I programmed
**rMetabolomics** to calculate sample size using the widely acceptable
standards: significance level (alpha) = 0.05, power = 80%.

``` r
#Define the input data, significance level and power
dat<-processedDataset
alpha_level=0.05
power_level=0.80


dataSubset<-dat[,c(2,sig.list)]#retrives only the significant metabolites

class1<-levels(dataSubset[,1])[1]#Gets the 1st class name (1st level in the sample type/class cloumn)
class2<-levels(dataSubset[,1])[2]#Gets the 2st class name (2nd level in the sample type/class cloumn)

iniFlag<-FALSE #Flip-flop logical flag to mark the first time a condition occurs ;-))
for (metaCounter in 2:ncol(dataSubset)){ #a loop to calculate mean delta, sd and sample size for each metabolite
       means_delta<-mean(dataSubset[dataSubset[,1]==class1, metaCounter])-mean(dataSubset[dataSubset[,1]==class2, metaCounter])
       col_sd<-sd(dataSubset[, metaCounter])
       power_t <- power.t.test(n = NULL, sd = col_sd, delta = means_delta, sig.level=alpha_level, power = power_level, type="two.sample", alternative = "two.sided")
       if (iniFlag==FALSE){
         iniFlag<-TRUE
         n.list<-c(round_any(power_t$n,1, f = ceiling)) #round sample size to the upper value to be statistically more conservative ;-))
         sd_list<-c(col_sd)
         delta_list<-c(means_delta)
       
       }else{
         
         n.list<-c(n.list, round_any(power_t$n,1, f = ceiling))#round sample size to the upper value to be statistically more conservative ;-))
         sd_list<-c(sd_list, col_sd)
         delta_list<-c(delta_list, means_delta)

      }

}


dataSubset<-matrix(rep(0,length(sig.list)*5), nrow=length(sig.list),ncol=5)#initialize dataSubset
dataSubset[,1]<-colnames(dat)[sig.list]#gets the names of significant metabolites
dataSubset[,2]<-n.list #list of sample sizes
dataSubset[,3]<-sig.list #positions of significant metabolites
dataSubset[,4]<-sd_list #list of sd
dataSubset[,5]<-delta_list #list of mean delta
colnames(dataSubset)<-c("Metabolite","Sample size (n per group)","","","")

#rearrange the matrix by the n size (ascending order ie from the smallest sample size to the greatest)
dataSubset<-dataSubset[order(as.numeric(dataSubset[,2])),]

#Calls drawTable to draw a table of the 1st and 2d columns only (metabolites name and n)
drawTable( dataSubset[,1:2], title_txt="Sample size calculation", footer_txt=paste("Significance level (Alpha)=",alpha_level,", Power=",power_level, sep=""), title_txt_size=15,footer_fontface="italic" )
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/powercalc-1.png)

``` r
#rearrange the matrix by the n size - this time in descending order (just for the sake of a better looking 3D plot next ;-))
dataSubset<-dataSubset[order(-as.numeric(dataSubset[,2])),]
sig.list.arranged<-as.numeric(dataSubset[,3])
sd_list.arranged<-as.numeric(dataSubset[,4])
delta_list.arranged<-as.numeric(dataSubset[,5])
```

<br>

### 3.5 Power Array **(Interactive 3D Model)**

Finally I programmed **rMetabolomics** to plot an innovative **“Power
Array”** which is an interactive 3D model. Power Array shows the
statistical power obtained for each significant metabolites at a set of
fixed sample sizes (n of 3,6,10,16,24,40,60,100,150 and 200). Its aim is
to help researchers decide the final sample size and metabolites to be
used in further studies in an easy intuitive way. Power Array 3D model
can be manipulated by the computer mouse.

``` r
testN<-c(3,6,10,16,24,40,60,100,150,200)# sets fixed sample sizes at which we will compare the power of significant metabolites
dataSubset<-matrix(rep(0,length(sig.list.arranged)*10), nrow=length(sig.list.arranged),ncol=10)#initialize matrix
rownames(dataSubset)<-colnames(dat)[sig.list.arranged]# get names of the significant metabolites
colnames(dataSubset)<-testN

#A nested loop to cycle through the significant metabolites and perform power calculation at the preset fixed sample sizes selected above
for (metaCounter in 1:length(sig.list.arranged)){
  for (pwrlvlCounter in 1: length(testN)){
          power_t <- power.t.test(n = testN[pwrlvlCounter], sd = sd_list.arranged[metaCounter], delta = delta_list.arranged[metaCounter], sig.level=alpha_level, power = NULL, type="two.sample", alternative = "two.sided")
       
       dataSubset[metaCounter,pwrlvlCounter ] <-  power_t$power #save the power values to the dataSubset (these are the z values of the 3D bars)
  }
}

#set the initial details
nClmns<-ncol(dataSubset)
nRws<-nrow(dataSubset)
x_label="Metabolomics Technique"
y_label="N patients per group" 
z_label="Power"
grphTitle="Power Array (Interactive 3D Model)"

# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("red",  "orange","yellow" ,"green4") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)


        #Plot graphics
  hist3D (x = 1:nRws, y = 1:nClmns, z = dataSubset,
            bty = "g", phi = 20,  theta = -60,
            xlab = x_label, ylab = y_label, zlab = z_label, main = grphTitle, 
            col=color, border = "black", shade = 0.8,box=FALSE,  axes = FALSE,
            ticktype = "detailed", space = 0.15, d = 2, cex.axis = 1e-9, add=FALSE, plot = FALSE)
    
    # Use text3D to label x axis
  text3D(x = 1:nRws, y = rep(-1, nRws), z = rep(0, nRws),
            labels = rownames(dataSubset),
            add = TRUE, adj = 1, plot = FALSE)
    # Use text3D to label y axis
  text3D(x = rep(0.5, nClmns),   y = 1:nClmns, z = rep(0, nClmns),
            labels  = colnames(dataSubset),
            add = TRUE, adj = 1, plot = FALSE)              

    #Use RGL to create an interactive plot ;-)))
  plotrgl (smooth = TRUE, x = 1:nRws, y = 1:nClmns, z = dataSubset,
            bty = "g", phi = 20,  theta = -60,expand = 0.1,
            xlab = "", ylab = "", zlab = "", main = grphTitle,
            col = color, border = "black", shade = 0.4,box=TRUE, axes = FALSE,
            ticktype = "detailed", space = 0.9, d = 2, cex.axis = 1e-9)

    
    #put Z axiz scale
    axes3d( edges=c("z") )
    
    
    
#That is it, Thank you ;-))))
```

![](Sayed-Metwaly---Final-Project--updated-_files/figure-markdown_github/powergrid%20-1.png)
<script>/*
* Copyright (C) 2009 Apple Inc. All Rights Reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY APPLE INC. ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
* PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL APPLE INC. OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
* OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* Copyright (2016) Duncan Murdoch - fixed CanvasMatrix4.ortho,
* cleaned up.
*/
/*
CanvasMatrix4 class
This class implements a 4x4 matrix. It has functions which
duplicate the functionality of the OpenGL matrix stack and
glut functions.
IDL:
[
Constructor(in CanvasMatrix4 matrix),           // copy passed matrix into new CanvasMatrix4
Constructor(in sequence<float> array)           // create new CanvasMatrix4 with 16 floats (row major)
Constructor()                                   // create new CanvasMatrix4 with identity matrix
]
interface CanvasMatrix4 {
attribute float m11;
attribute float m12;
attribute float m13;
attribute float m14;
attribute float m21;
attribute float m22;
attribute float m23;
attribute float m24;
attribute float m31;
attribute float m32;
attribute float m33;
attribute float m34;
attribute float m41;
attribute float m42;
attribute float m43;
attribute float m44;
void load(in CanvasMatrix4 matrix);                 // copy the values from the passed matrix
void load(in sequence<float> array);                // copy 16 floats into the matrix
sequence<float> getAsArray();                       // return the matrix as an array of 16 floats
WebGLFloatArray getAsCanvasFloatArray();           // return the matrix as a WebGLFloatArray with 16 values
void makeIdentity();                                // replace the matrix with identity
void transpose();                                   // replace the matrix with its transpose
void invert();                                      // replace the matrix with its inverse
void translate(in float x, in float y, in float z); // multiply the matrix by passed translation values on the right
void scale(in float x, in float y, in float z);     // multiply the matrix by passed scale values on the right
void rotate(in float angle,                         // multiply the matrix by passed rotation values on the right
in float x, in float y, in float z);    // (angle is in degrees)
void multRight(in CanvasMatrix matrix);             // multiply the matrix by the passed matrix on the right
void multLeft(in CanvasMatrix matrix);              // multiply the matrix by the passed matrix on the left
void ortho(in float left, in float right,           // multiply the matrix by the passed ortho values on the right
in float bottom, in float top,
in float near, in float far);
void frustum(in float left, in float right,         // multiply the matrix by the passed frustum values on the right
in float bottom, in float top,
in float near, in float far);
void perspective(in float fovy, in float aspect,    // multiply the matrix by the passed perspective values on the right
in float zNear, in float zFar);
void lookat(in float eyex, in float eyey, in float eyez,    // multiply the matrix by the passed lookat
in float ctrx, in float ctry, in float ctrz,    // values on the right
in float upx, in float upy, in float upz);
}
*/
CanvasMatrix4 = function(m)
{
if (typeof m == 'object') {
if ("length" in m && m.length >= 16) {
this.load(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
return;
}
else if (m instanceof CanvasMatrix4) {
this.load(m);
return;
}
}
this.makeIdentity();
};
CanvasMatrix4.prototype.load = function()
{
if (arguments.length == 1 && typeof arguments[0] == 'object') {
var matrix = arguments[0];
if ("length" in matrix && matrix.length == 16) {
this.m11 = matrix[0];
this.m12 = matrix[1];
this.m13 = matrix[2];
this.m14 = matrix[3];
this.m21 = matrix[4];
this.m22 = matrix[5];
this.m23 = matrix[6];
this.m24 = matrix[7];
this.m31 = matrix[8];
this.m32 = matrix[9];
this.m33 = matrix[10];
this.m34 = matrix[11];
this.m41 = matrix[12];
this.m42 = matrix[13];
this.m43 = matrix[14];
this.m44 = matrix[15];
return;
}
if (arguments[0] instanceof CanvasMatrix4) {
this.m11 = matrix.m11;
this.m12 = matrix.m12;
this.m13 = matrix.m13;
this.m14 = matrix.m14;
this.m21 = matrix.m21;
this.m22 = matrix.m22;
this.m23 = matrix.m23;
this.m24 = matrix.m24;
this.m31 = matrix.m31;
this.m32 = matrix.m32;
this.m33 = matrix.m33;
this.m34 = matrix.m34;
this.m41 = matrix.m41;
this.m42 = matrix.m42;
this.m43 = matrix.m43;
this.m44 = matrix.m44;
return;
}
}
this.makeIdentity();
};
CanvasMatrix4.prototype.getAsArray = function()
{
return [
this.m11, this.m12, this.m13, this.m14,
this.m21, this.m22, this.m23, this.m24,
this.m31, this.m32, this.m33, this.m34,
this.m41, this.m42, this.m43, this.m44
];
};
CanvasMatrix4.prototype.getAsWebGLFloatArray = function()
{
return new WebGLFloatArray(this.getAsArray());
};
CanvasMatrix4.prototype.makeIdentity = function()
{
this.m11 = 1;
this.m12 = 0;
this.m13 = 0;
this.m14 = 0;
this.m21 = 0;
this.m22 = 1;
this.m23 = 0;
this.m24 = 0;
this.m31 = 0;
this.m32 = 0;
this.m33 = 1;
this.m34 = 0;
this.m41 = 0;
this.m42 = 0;
this.m43 = 0;
this.m44 = 1;
};
CanvasMatrix4.prototype.transpose = function()
{
var tmp = this.m12;
this.m12 = this.m21;
this.m21 = tmp;
tmp = this.m13;
this.m13 = this.m31;
this.m31 = tmp;
tmp = this.m14;
this.m14 = this.m41;
this.m41 = tmp;
tmp = this.m23;
this.m23 = this.m32;
this.m32 = tmp;
tmp = this.m24;
this.m24 = this.m42;
this.m42 = tmp;
tmp = this.m34;
this.m34 = this.m43;
this.m43 = tmp;
};
CanvasMatrix4.prototype.invert = function()
{
// Calculate the 4x4 determinant
// If the determinant is zero,
// then the inverse matrix is not unique.
var det = this._determinant4x4();
if (Math.abs(det) < 1e-8)
return null;
this._makeAdjoint();
// Scale the adjoint matrix to get the inverse
this.m11 /= det;
this.m12 /= det;
this.m13 /= det;
this.m14 /= det;
this.m21 /= det;
this.m22 /= det;
this.m23 /= det;
this.m24 /= det;
this.m31 /= det;
this.m32 /= det;
this.m33 /= det;
this.m34 /= det;
this.m41 /= det;
this.m42 /= det;
this.m43 /= det;
this.m44 /= det;
};
CanvasMatrix4.prototype.translate = function(x,y,z)
{
if (x === undefined)
x = 0;
if (y === undefined)
y = 0;
if (z === undefined)
z = 0;
var matrix = new CanvasMatrix4();
matrix.m41 = x;
matrix.m42 = y;
matrix.m43 = z;
this.multRight(matrix);
};
CanvasMatrix4.prototype.scale = function(x,y,z)
{
if (x === undefined)
x = 1;
if (z === undefined) {
if (y === undefined) {
y = x;
z = x;
}
else
z = 1;
}
else if (y === undefined)
y = x;
var matrix = new CanvasMatrix4();
matrix.m11 = x;
matrix.m22 = y;
matrix.m33 = z;
this.multRight(matrix);
};
CanvasMatrix4.prototype.rotate = function(angle,x,y,z)
{
// angles are in degrees. Switch to radians
angle = angle / 180 * Math.PI;
angle /= 2;
var sinA = Math.sin(angle);
var cosA = Math.cos(angle);
var sinA2 = sinA * sinA;
// normalize
var length = Math.sqrt(x * x + y * y + z * z);
if (length === 0) {
// bad vector, just use something reasonable
x = 0;
y = 0;
z = 1;
} else if (length != 1) {
x /= length;
y /= length;
z /= length;
}
var mat = new CanvasMatrix4();
// optimize case where axis is along major axis
if (x == 1 && y === 0 && z === 0) {
mat.m11 = 1;
mat.m12 = 0;
mat.m13 = 0;
mat.m21 = 0;
mat.m22 = 1 - 2 * sinA2;
mat.m23 = 2 * sinA * cosA;
mat.m31 = 0;
mat.m32 = -2 * sinA * cosA;
mat.m33 = 1 - 2 * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else if (x === 0 && y == 1 && z === 0) {
mat.m11 = 1 - 2 * sinA2;
mat.m12 = 0;
mat.m13 = -2 * sinA * cosA;
mat.m21 = 0;
mat.m22 = 1;
mat.m23 = 0;
mat.m31 = 2 * sinA * cosA;
mat.m32 = 0;
mat.m33 = 1 - 2 * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else if (x === 0 && y === 0 && z == 1) {
mat.m11 = 1 - 2 * sinA2;
mat.m12 = 2 * sinA * cosA;
mat.m13 = 0;
mat.m21 = -2 * sinA * cosA;
mat.m22 = 1 - 2 * sinA2;
mat.m23 = 0;
mat.m31 = 0;
mat.m32 = 0;
mat.m33 = 1;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
} else {
var x2 = x*x;
var y2 = y*y;
var z2 = z*z;
mat.m11 = 1 - 2 * (y2 + z2) * sinA2;
mat.m12 = 2 * (x * y * sinA2 + z * sinA * cosA);
mat.m13 = 2 * (x * z * sinA2 - y * sinA * cosA);
mat.m21 = 2 * (y * x * sinA2 - z * sinA * cosA);
mat.m22 = 1 - 2 * (z2 + x2) * sinA2;
mat.m23 = 2 * (y * z * sinA2 + x * sinA * cosA);
mat.m31 = 2 * (z * x * sinA2 + y * sinA * cosA);
mat.m32 = 2 * (z * y * sinA2 - x * sinA * cosA);
mat.m33 = 1 - 2 * (x2 + y2) * sinA2;
mat.m14 = mat.m24 = mat.m34 = 0;
mat.m41 = mat.m42 = mat.m43 = 0;
mat.m44 = 1;
}
this.multRight(mat);
};
CanvasMatrix4.prototype.multRight = function(mat)
{
var m11 = (this.m11 * mat.m11 + this.m12 * mat.m21 +
this.m13 * mat.m31 + this.m14 * mat.m41);
var m12 = (this.m11 * mat.m12 + this.m12 * mat.m22 +
this.m13 * mat.m32 + this.m14 * mat.m42);
var m13 = (this.m11 * mat.m13 + this.m12 * mat.m23 +
this.m13 * mat.m33 + this.m14 * mat.m43);
var m14 = (this.m11 * mat.m14 + this.m12 * mat.m24 +
this.m13 * mat.m34 + this.m14 * mat.m44);
var m21 = (this.m21 * mat.m11 + this.m22 * mat.m21 +
this.m23 * mat.m31 + this.m24 * mat.m41);
var m22 = (this.m21 * mat.m12 + this.m22 * mat.m22 +
this.m23 * mat.m32 + this.m24 * mat.m42);
var m23 = (this.m21 * mat.m13 + this.m22 * mat.m23 +
this.m23 * mat.m33 + this.m24 * mat.m43);
var m24 = (this.m21 * mat.m14 + this.m22 * mat.m24 +
this.m23 * mat.m34 + this.m24 * mat.m44);
var m31 = (this.m31 * mat.m11 + this.m32 * mat.m21 +
this.m33 * mat.m31 + this.m34 * mat.m41);
var m32 = (this.m31 * mat.m12 + this.m32 * mat.m22 +
this.m33 * mat.m32 + this.m34 * mat.m42);
var m33 = (this.m31 * mat.m13 + this.m32 * mat.m23 +
this.m33 * mat.m33 + this.m34 * mat.m43);
var m34 = (this.m31 * mat.m14 + this.m32 * mat.m24 +
this.m33 * mat.m34 + this.m34 * mat.m44);
var m41 = (this.m41 * mat.m11 + this.m42 * mat.m21 +
this.m43 * mat.m31 + this.m44 * mat.m41);
var m42 = (this.m41 * mat.m12 + this.m42 * mat.m22 +
this.m43 * mat.m32 + this.m44 * mat.m42);
var m43 = (this.m41 * mat.m13 + this.m42 * mat.m23 +
this.m43 * mat.m33 + this.m44 * mat.m43);
var m44 = (this.m41 * mat.m14 + this.m42 * mat.m24 +
this.m43 * mat.m34 + this.m44 * mat.m44);
this.m11 = m11;
this.m12 = m12;
this.m13 = m13;
this.m14 = m14;
this.m21 = m21;
this.m22 = m22;
this.m23 = m23;
this.m24 = m24;
this.m31 = m31;
this.m32 = m32;
this.m33 = m33;
this.m34 = m34;
this.m41 = m41;
this.m42 = m42;
this.m43 = m43;
this.m44 = m44;
};
CanvasMatrix4.prototype.multLeft = function(mat)
{
var m11 = (mat.m11 * this.m11 + mat.m12 * this.m21 +
mat.m13 * this.m31 + mat.m14 * this.m41);
var m12 = (mat.m11 * this.m12 + mat.m12 * this.m22 +
mat.m13 * this.m32 + mat.m14 * this.m42);
var m13 = (mat.m11 * this.m13 + mat.m12 * this.m23 +
mat.m13 * this.m33 + mat.m14 * this.m43);
var m14 = (mat.m11 * this.m14 + mat.m12 * this.m24 +
mat.m13 * this.m34 + mat.m14 * this.m44);
var m21 = (mat.m21 * this.m11 + mat.m22 * this.m21 +
mat.m23 * this.m31 + mat.m24 * this.m41);
var m22 = (mat.m21 * this.m12 + mat.m22 * this.m22 +
mat.m23 * this.m32 + mat.m24 * this.m42);
var m23 = (mat.m21 * this.m13 + mat.m22 * this.m23 +
mat.m23 * this.m33 + mat.m24 * this.m43);
var m24 = (mat.m21 * this.m14 + mat.m22 * this.m24 +
mat.m23 * this.m34 + mat.m24 * this.m44);
var m31 = (mat.m31 * this.m11 + mat.m32 * this.m21 +
mat.m33 * this.m31 + mat.m34 * this.m41);
var m32 = (mat.m31 * this.m12 + mat.m32 * this.m22 +
mat.m33 * this.m32 + mat.m34 * this.m42);
var m33 = (mat.m31 * this.m13 + mat.m32 * this.m23 +
mat.m33 * this.m33 + mat.m34 * this.m43);
var m34 = (mat.m31 * this.m14 + mat.m32 * this.m24 +
mat.m33 * this.m34 + mat.m34 * this.m44);
var m41 = (mat.m41 * this.m11 + mat.m42 * this.m21 +
mat.m43 * this.m31 + mat.m44 * this.m41);
var m42 = (mat.m41 * this.m12 + mat.m42 * this.m22 +
mat.m43 * this.m32 + mat.m44 * this.m42);
var m43 = (mat.m41 * this.m13 + mat.m42 * this.m23 +
mat.m43 * this.m33 + mat.m44 * this.m43);
var m44 = (mat.m41 * this.m14 + mat.m42 * this.m24 +
mat.m43 * this.m34 + mat.m44 * this.m44);
this.m11 = m11;
this.m12 = m12;
this.m13 = m13;
this.m14 = m14;
this.m21 = m21;
this.m22 = m22;
this.m23 = m23;
this.m24 = m24;
this.m31 = m31;
this.m32 = m32;
this.m33 = m33;
this.m34 = m34;
this.m41 = m41;
this.m42 = m42;
this.m43 = m43;
this.m44 = m44;
};
CanvasMatrix4.prototype.ortho = function(left, right, bottom, top, near, far)
{
var tx = (left + right) / (left - right);
var ty = (top + bottom) / (bottom - top);
var tz = (far + near) / (near - far);
var matrix = new CanvasMatrix4();
matrix.m11 = 2 / (right - left);
matrix.m12 = 0;
matrix.m13 = 0;
matrix.m14 = 0;
matrix.m21 = 0;
matrix.m22 = 2 / (top - bottom);
matrix.m23 = 0;
matrix.m24 = 0;
matrix.m31 = 0;
matrix.m32 = 0;
matrix.m33 = -2 / (far - near);
matrix.m34 = 0;
matrix.m41 = tx;
matrix.m42 = ty;
matrix.m43 = tz;
matrix.m44 = 1;
this.multRight(matrix);
};
CanvasMatrix4.prototype.frustum = function(left, right, bottom, top, near, far)
{
var matrix = new CanvasMatrix4();
var A = (right + left) / (right - left);
var B = (top + bottom) / (top - bottom);
var C = -(far + near) / (far - near);
var D = -(2 * far * near) / (far - near);
matrix.m11 = (2 * near) / (right - left);
matrix.m12 = 0;
matrix.m13 = 0;
matrix.m14 = 0;
matrix.m21 = 0;
matrix.m22 = 2 * near / (top - bottom);
matrix.m23 = 0;
matrix.m24 = 0;
matrix.m31 = A;
matrix.m32 = B;
matrix.m33 = C;
matrix.m34 = -1;
matrix.m41 = 0;
matrix.m42 = 0;
matrix.m43 = D;
matrix.m44 = 0;
this.multRight(matrix);
};
CanvasMatrix4.prototype.perspective = function(fovy, aspect, zNear, zFar)
{
var top = Math.tan(fovy * Math.PI / 360) * zNear;
var bottom = -top;
var left = aspect * bottom;
var right = aspect * top;
this.frustum(left, right, bottom, top, zNear, zFar);
};
CanvasMatrix4.prototype.lookat = function(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz)
{
var matrix = new CanvasMatrix4();
// Make rotation matrix
// Z vector
var zx = eyex - centerx;
var zy = eyey - centery;
var zz = eyez - centerz;
var mag = Math.sqrt(zx * zx + zy * zy + zz * zz);
if (mag) {
zx /= mag;
zy /= mag;
zz /= mag;
}
// Y vector
var yx = upx;
var yy = upy;
var yz = upz;
// X vector = Y cross Z
xx =  yy * zz - yz * zy;
xy = -yx * zz + yz * zx;
xz =  yx * zy - yy * zx;
// Recompute Y = Z cross X
yx = zy * xz - zz * xy;
yy = -zx * xz + zz * xx;
yx = zx * xy - zy * xx;
// cross product gives area of parallelogram, which is < 1.0 for
// non-perpendicular unit-length vectors; so normalize x, y here
mag = Math.sqrt(xx * xx + xy * xy + xz * xz);
if (mag) {
xx /= mag;
xy /= mag;
xz /= mag;
}
mag = Math.sqrt(yx * yx + yy * yy + yz * yz);
if (mag) {
yx /= mag;
yy /= mag;
yz /= mag;
}
matrix.m11 = xx;
matrix.m12 = xy;
matrix.m13 = xz;
matrix.m14 = 0;
matrix.m21 = yx;
matrix.m22 = yy;
matrix.m23 = yz;
matrix.m24 = 0;
matrix.m31 = zx;
matrix.m32 = zy;
matrix.m33 = zz;
matrix.m34 = 0;
matrix.m41 = 0;
matrix.m42 = 0;
matrix.m43 = 0;
matrix.m44 = 1;
matrix.translate(-eyex, -eyey, -eyez);
this.multRight(matrix);
};
// Support functions
CanvasMatrix4.prototype._determinant2x2 = function(a, b, c, d)
{
return a * d - b * c;
};
CanvasMatrix4.prototype._determinant3x3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3)
{
return a1 * this._determinant2x2(b2, b3, c2, c3) -
b1 * this._determinant2x2(a2, a3, c2, c3) +
c1 * this._determinant2x2(a2, a3, b2, b3);
};
CanvasMatrix4.prototype._determinant4x4 = function()
{
var a1 = this.m11;
var b1 = this.m12;
var c1 = this.m13;
var d1 = this.m14;
var a2 = this.m21;
var b2 = this.m22;
var c2 = this.m23;
var d2 = this.m24;
var a3 = this.m31;
var b3 = this.m32;
var c3 = this.m33;
var d3 = this.m34;
var a4 = this.m41;
var b4 = this.m42;
var c4 = this.m43;
var d4 = this.m44;
return a1 * this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4) -
b1 * this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4) +
c1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4) -
d1 * this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
};
CanvasMatrix4.prototype._makeAdjoint = function()
{
var a1 = this.m11;
var b1 = this.m12;
var c1 = this.m13;
var d1 = this.m14;
var a2 = this.m21;
var b2 = this.m22;
var c2 = this.m23;
var d2 = this.m24;
var a3 = this.m31;
var b3 = this.m32;
var c3 = this.m33;
var d3 = this.m34;
var a4 = this.m41;
var b4 = this.m42;
var c4 = this.m43;
var d4 = this.m44;
// Row column labeling reversed since we transpose rows & columns
this.m11  =   this._determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
this.m21  = - this._determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
this.m31  =   this._determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
this.m41  = - this._determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
this.m12  = - this._determinant3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
this.m22  =   this._determinant3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
this.m32  = - this._determinant3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
this.m42  =   this._determinant3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);
this.m13  =   this._determinant3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
this.m23  = - this._determinant3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
this.m33  =   this._determinant3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
this.m43  = - this._determinant3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);
this.m14  = - this._determinant3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
this.m24  =   this._determinant3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
this.m34  = - this._determinant3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
this.m44  =   this._determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
};</script>
<script>// To generate the help pages for this library, use
// jsdoc --destination ../../../doc/rglwidgetClass --template ~/node_modules/jsdoc-baseline rglClass.src.js
// To validate, use
// setwd(".../inst/htmlwidgets/lib/rglClass")
// hints <- js::jshint(readLines("rglClass.src.js"))
// hints[, c("line", "reason")]
/**
* The class of an rgl widget
* @class
*/
rglwidgetClass = function() {
this.canvas = null;
this.userMatrix = new CanvasMatrix4();
this.types = [];
this.prMatrix = new CanvasMatrix4();
this.mvMatrix = new CanvasMatrix4();
this.vp = null;
this.prmvMatrix = null;
this.origs = null;
this.gl = null;
this.scene = null;
this.select = {state: "inactive", subscene: null, region: {p1: {x:0, y:0}, p2: {x:0, y:0}}};
this.drawing = false;
};
/**
* Multiply matrix by vector
* @returns {number[]}
* @param M {number[][]} Left operand
* @param v {number[]} Right operand
*/
rglwidgetClass.prototype.multMV = function(M, v) {
return [ M.m11 * v[0] + M.m12 * v[1] + M.m13 * v[2] + M.m14 * v[3],
M.m21 * v[0] + M.m22 * v[1] + M.m23 * v[2] + M.m24 * v[3],
M.m31 * v[0] + M.m32 * v[1] + M.m33 * v[2] + M.m34 * v[3],
M.m41 * v[0] + M.m42 * v[1] + M.m43 * v[2] + M.m44 * v[3]
];
};
/**
* Multiply row vector by Matrix
* @returns {number[]}
* @param v {number[]} left operand
* @param M {number[][]} right operand
*/
rglwidgetClass.prototype.multVM = function(v, M) {
return [ M.m11 * v[0] + M.m21 * v[1] + M.m31 * v[2] + M.m41 * v[3],
M.m12 * v[0] + M.m22 * v[1] + M.m32 * v[2] + M.m42 * v[3],
M.m13 * v[0] + M.m23 * v[1] + M.m33 * v[2] + M.m43 * v[3],
M.m14 * v[0] + M.m24 * v[1] + M.m34 * v[2] + M.m44 * v[3]
];
};
/**
* Euclidean length of a vector
* @returns {number}
* @param v {number[]}
*/
rglwidgetClass.prototype.vlen = function(v) {
return Math.sqrt(this.dotprod(v, v));
};
/**
* Dot product of two vectors
* @instance rglwidgetClass
* @returns {number}
* @param a {number[]}
* @param b {number[]}
*/
rglwidgetClass.prototype.dotprod = function(a, b) {
return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
};
/**
* Cross product of two vectors
* @returns {number[]}
* @param a {number[]}
* @param b {number[]}
*/
rglwidgetClass.prototype.xprod = function(a, b) {
return [a[1]*b[2] - a[2]*b[1],
a[2]*b[0] - a[0]*b[2],
a[0]*b[1] - a[1]*b[0]];
};
/**
* Bind vectors or matrices by columns
* @returns {number[][]}
* @param a {number[]|number[][]}
* @param b {number[]|number[][]}
*/
rglwidgetClass.prototype.cbind = function(a, b) {
if (b.length < a.length)
b = this.repeatToLen(b, a.length);
else if (a.length < b.length)
a = this.repeatToLen(a, b.length);
return a.map(function(currentValue, index, array) {
return currentValue.concat(b[index]);
});
};
/**
* Swap elements
* @returns {any[]}
* @param a {any[]}
* @param i {number} Element to swap
* @param j {number} Other element to swap
*/
rglwidgetClass.prototype.swap = function(a, i, j) {
var temp = a[i];
a[i] = a[j];
a[j] = temp;
};
/**
* Flatten a matrix into a vector
* @returns {any[]}
* @param a {any[][]}
*/
rglwidgetClass.prototype.flatten = function(arr, result) {
var value;
if (typeof result === "undefined") result = [];
for (var i = 0, length = arr.length; i < length; i++) {
value = arr[i];
if (Array.isArray(value)) {
this.flatten(value, result);
} else {
result.push(value);
}
}
return result;
};
/**
* set element of 1d or 2d array as if it was flattened.
* Column major, zero based!
* @returns {any[]|any[][]}
* @param {any[]|any[][]} a - array
* @param {number} i - element
* @param {any} value
*/
rglwidgetClass.prototype.setElement = function(a, i, value) {
if (Array.isArray(a[0])) {
var dim = a.length,
col = Math.floor(i/dim),
row = i % dim;
a[row][col] = value;
} else {
a[i] = value;
}
};
/**
* Transpose an array
* @returns {any[][]}
* @param {any[][]} a
*/
rglwidgetClass.prototype.transpose = function(a) {
var newArray = [],
n = a.length,
m = a[0].length,
i;
for(i = 0; i < m; i++){
newArray.push([]);
}
for(i = 0; i < n; i++){
for(var j = 0; j < m; j++){
newArray[j].push(a[i][j]);
}
}
return newArray;
};
/**
* Calculate sum of squares of a numeric vector
* @returns {number}
* @param {number[]} x
*/
rglwidgetClass.prototype.sumsq = function(x) {
var result = 0, i;
for (i=0; i < x.length; i++)
result += x[i]*x[i];
return result;
};
/**
* Convert a matrix to a CanvasMatrix4
* @returns {CanvasMatrix4}
* @param {number[][]|number[]} mat
*/
rglwidgetClass.prototype.toCanvasMatrix4 = function(mat) {
if (mat instanceof CanvasMatrix4)
return mat;
var result = new CanvasMatrix4();
mat = this.flatten(this.transpose(mat));
result.load(mat);
return result;
};
/**
* Convert an R-style numeric colour string to an rgb vector
* @returns {number[]}
* @param {string} s
*/
rglwidgetClass.prototype.stringToRgb = function(s) {
s = s.replace("#", "");
var bigint = parseInt(s, 16);
return [((bigint >> 16) & 255)/255,
((bigint >> 8) & 255)/255,
(bigint & 255)/255];
};
/**
* Take a component-by-component product of two 3 vectors
* @returns {number[]}
* @param {number[]} x
* @param {number[]} y
*/
rglwidgetClass.prototype.componentProduct = function(x, y) {
if (typeof y === "undefined") {
this.alertOnce("Bad arg to componentProduct");
}
var result = new Float32Array(3), i;
for (i = 0; i<3; i++)
result[i] = x[i]*y[i];
return result;
};
/**
* Get next higher power of two
* @returns { number }
* @param { number } value - input value
*/
rglwidgetClass.prototype.getPowerOfTwo = function(value) {
var pow = 1;
while(pow<value) {
pow *= 2;
}
return pow;
};
/**
* Unique entries
* @returns { any[] }
* @param { any[] } arr - An array
*/
rglwidgetClass.prototype.unique = function(arr) {
arr = [].concat(arr);
return arr.filter(function(value, index, self) {
return self.indexOf(value) === index;
});
};
/**
* Shallow compare of arrays
* @returns { boolean }
* @param { any[] } a - An array
* @param { any[] } b - Another array
*/
rglwidgetClass.prototype.equalArrays = function(a, b) {
return a === b || (a && b &&
a.length === b.length &&
a.every(function(v, i) {return v === b[i];}));
};
/**
* Repeat an array to a desired length
* @returns {any[]}
* @param {any | any[]} arr The input array
* @param {number} len The desired output length
*/
rglwidgetClass.prototype.repeatToLen = function(arr, len) {
arr = [].concat(arr);
while (arr.length < len/2)
arr = arr.concat(arr);
return arr.concat(arr.slice(0, len - arr.length));
};
/**
* Give a single alert message, not to be repeated.
* @param {string} msg  The message to give.
*/
rglwidgetClass.prototype.alertOnce = function(msg) {
if (typeof this.alerted !== "undefined")
return;
this.alerted = true;
alert(msg);
};
rglwidgetClass.prototype.f_is_lit = 1;
rglwidgetClass.prototype.f_is_smooth = 2;
rglwidgetClass.prototype.f_has_texture = 4;
rglwidgetClass.prototype.f_depth_sort = 8;
rglwidgetClass.prototype.f_fixed_quads = 16;
rglwidgetClass.prototype.f_is_transparent = 32;
rglwidgetClass.prototype.f_is_lines = 64;
rglwidgetClass.prototype.f_sprites_3d = 128;
rglwidgetClass.prototype.f_sprite_3d = 256;
rglwidgetClass.prototype.f_is_subscene = 512;
rglwidgetClass.prototype.f_is_clipplanes = 1024;
rglwidgetClass.prototype.f_fixed_size = 2048;
rglwidgetClass.prototype.f_is_points = 4096;
rglwidgetClass.prototype.f_is_twosided = 8192;
rglwidgetClass.prototype.f_fat_lines = 16384;
rglwidgetClass.prototype.f_is_brush = 32768;
/**
* Which list does a particular id come from?
* @returns { string }
* @param {number} id The id to look up.
*/
rglwidgetClass.prototype.whichList = function(id) {
var obj = this.getObj(id),
flags = obj.flags;
if (obj.type === "light")
return "lights";
if (flags & this.f_is_subscene)
return "subscenes";
if (flags & this.f_is_clipplanes)
return "clipplanes";
if (flags & this.f_is_transparent)
return "transparent";
return "opaque";
};
/**
* Get an object by id number.
* @returns { Object }
* @param {number} id
*/
rglwidgetClass.prototype.getObj = function(id) {
if (typeof id !== "number") {
this.alertOnce("getObj id is "+typeof id);
}
return this.scene.objects[id];
};
/**
* Get ids of a particular type from a subscene or the whole scene
* @returns { number[] }
* @param {string} type What type of object?
* @param {number} subscene  Which subscene?  If not given, find in the whole scene
*/
rglwidgetClass.prototype.getIdsByType = function(type, subscene) {
var
result = [], i, self = this;
if (typeof subscene === "undefined") {
Object.keys(this.scene.objects).forEach(
function(key) {
key = parseInt(key, 10);
if (self.getObj(key).type === type)
result.push(key);
});
} else {
ids = this.getObj(subscene).objects;
for (i=0; i < ids.length; i++) {
if (this.getObj(ids[i]).type === type) {
result.push(ids[i]);
}
}
}
return result;
};
/**
* Get a particular material property for an id
* @returns { any }
* @param {number} id  Which object?
* @param {string} property Which material property?
*/
rglwidgetClass.prototype.getMaterial = function(id, property) {
var obj = this.getObj(id), mat;
if (typeof obj.material === "undefined")
console.error("material undefined");
mat = obj.material[property];
if (typeof mat === "undefined")
mat = this.scene.material[property];
return mat;
};
/**
* Is a particular id in a subscene?
* @returns { boolean }
* @param {number} id Which id?
* @param {number} subscene Which subscene id?
*/
rglwidgetClass.prototype.inSubscene = function(id, subscene) {
return this.getObj(subscene).objects.indexOf(id) > -1;
};
/**
* Add an id to a subscene.
* @param {number} id Which id?
* @param {number} subscene Which subscene id?
*/
rglwidgetClass.prototype.addToSubscene = function(id, subscene) {
var thelist,
thesub = this.getObj(subscene),
ids = [id],
obj = this.getObj(id), i;
if (typeof obj != "undefined" && typeof (obj.newIds) !== "undefined") {
ids = ids.concat(obj.newIds);
}
thesub.objects = [].concat(thesub.objects);
for (i = 0; i < ids.length; i++) {
id = ids[i];
if (thesub.objects.indexOf(id) == -1) {
thelist = this.whichList(id);
thesub.objects.push(id);
thesub[thelist].push(id);
}
}
};
/**
* Delete an id from a subscene
* @param { number } id - the id to add
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.delFromSubscene = function(id, subscene) {
var thelist,
thesub = this.getObj(subscene),
obj = this.getObj(id),
ids = [id], i;
if (typeof obj !== "undefined" && typeof (obj.newIds) !== "undefined")
ids = ids.concat(obj.newIds);
thesub.objects = [].concat(thesub.objects); // It might be a scalar
for (j=0; j<ids.length;j++) {
id = ids[j];
i = thesub.objects.indexOf(id);
if (i > -1) {
thesub.objects.splice(i, 1);
thelist = this.whichList(id);
i = thesub[thelist].indexOf(id);
thesub[thelist].splice(i, 1);
}
}
};
/**
* Set the ids in a subscene
* @param { number[] } ids - the ids to set
* @param { number } subsceneid - the id of the subscene
*/
rglwidgetClass.prototype.setSubsceneEntries = function(ids, subsceneid) {
var sub = this.getObj(subsceneid);
sub.objects = ids;
this.initSubscene(subsceneid);
};
/**
* Get the ids in a subscene
* @returns {number[]}
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.getSubsceneEntries = function(subscene) {
return this.getObj(subscene).objects;
};
/**
* Get the ids of the subscenes within a subscene
* @returns { number[] }
* @param { number } subscene - the id of the subscene
*/
rglwidgetClass.prototype.getChildSubscenes = function(subscene) {
return this.getObj(subscene).subscenes;
};
/**
* Start drawing
* @returns { boolean } Previous state
*/
rglwidgetClass.prototype.startDrawing = function() {
var value = this.drawing;
this.drawing = true;
return value;
};
/**
* Stop drawing and check for context loss
* @param { boolean } saved - Previous state
*/
rglwidgetClass.prototype.stopDrawing = function(saved) {
this.drawing = saved;
if (!saved && this.gl && this.gl.isContextLost())
this.restartCanvas();
};
/**
* Generate the vertex shader for an object
* @returns {string}
* @param { number } id - Id of object
*/
rglwidgetClass.prototype.getVertexShader = function(id) {
var obj = this.getObj(id),
userShader = obj.userVertexShader,
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
nclipplanes = this.countClipplanes(),
fixed_size = flags & this.f_fixed_size,
is_points = flags & this.f_is_points,
is_twosided = flags & this.f_is_twosided,
fat_lines = flags & this.f_fat_lines,
is_brush = flags & this.f_is_brush,
result;
if (type === "clipplanes" || sprites_3d) return;
if (typeof userShader !== "undefined") return userShader;
result = "  /* ****** "+type+" object "+id+" vertex shader ****** */\n"+
"  attribute vec3 aPos;\n"+
"  attribute vec4 aCol;\n"+
" uniform mat4 mvMatrix;\n"+
" uniform mat4 prMatrix;\n"+
" varying vec4 vCol;\n"+
" varying vec4 vPosition;\n";
if ((is_lit && !fixed_quads && !is_brush) || sprite_3d)
result = result + "  attribute vec3 aNorm;\n"+
" uniform mat4 normMatrix;\n"+
" varying vec3 vNormal;\n";
if (has_texture || type === "text")
result = result + " attribute vec2 aTexcoord;\n"+
" varying vec2 vTexcoord;\n";
if (fixed_size)
result = result + "  uniform vec2 textScale;\n";
if (fixed_quads)
result = result + "  attribute vec2 aOfs;\n";
else if (sprite_3d)
result = result + "  uniform vec3 uOrig;\n"+
"  uniform float uSize;\n"+
"  uniform mat4 usermat;\n";
if (is_twosided)
result = result + "  attribute vec3 aPos1;\n"+
"  attribute vec3 aPos2;\n"+
"  varying float normz;\n";
if (fat_lines) {
result = result +   "  attribute vec3 aNext;\n"+
"  attribute vec2 aPoint;\n"+
"  varying vec2 vPoint;\n"+
"  varying float vLength;\n"+
"  uniform float uAspect;\n"+
"  uniform float uLwd;\n";
}
result = result + "  void main(void) {\n";
if ((nclipplanes || (!fixed_quads && !sprite_3d)) && !is_brush)
result = result + "    vPosition = mvMatrix * vec4(aPos, 1.);\n";
if (!fixed_quads && !sprite_3d && !is_brush)
result = result + "    gl_Position = prMatrix * vPosition;\n";
if (is_points) {
var size = this.getMaterial(id, "size");
result = result + "    gl_PointSize = "+size.toFixed(1)+";\n";
}
result = result + "    vCol = aCol;\n";
if (is_lit && !fixed_quads && !sprite_3d && !is_brush)
result = result + "    vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n";
if (has_texture || type == "text")
result = result + "    vTexcoord = aTexcoord;\n";
if (fixed_size)
result = result + "    vec4 pos = prMatrix * mvMatrix * vec4(aPos, 1.);\n"+
"   pos = pos/pos.w;\n"+
"   gl_Position = pos + vec4(aOfs*textScale, 0.,0.);\n";
if (type == "sprites" && !fixed_size)
result = result + "    vec4 pos = mvMatrix * vec4(aPos, 1.);\n"+
"   pos = pos/pos.w + vec4(aOfs, 0., 0.);\n"+
"   gl_Position = prMatrix*pos;\n";
if (sprite_3d)
result = result + "   vNormal = normalize((normMatrix * vec4(aNorm, 1.)).xyz);\n"+
"   vec4 pos = mvMatrix * vec4(uOrig, 1.);\n"+
"   vPosition = pos/pos.w + vec4(uSize*(vec4(aPos, 1.)*usermat).xyz,0.);\n"+
"   gl_Position = prMatrix * vPosition;\n";
if (is_twosided)
result = result + "   vec4 pos1 = prMatrix*(mvMatrix*vec4(aPos1, 1.));\n"+
"   pos1 = pos1/pos1.w - gl_Position/gl_Position.w;\n"+
"   vec4 pos2 = prMatrix*(mvMatrix*vec4(aPos2, 1.));\n"+
"   pos2 = pos2/pos2.w - gl_Position/gl_Position.w;\n"+
"   normz = pos1.x*pos2.y - pos1.y*pos2.x;\n";
if (fat_lines) 
/* This code was inspired by Matt Deslauriers' code in https://mattdesl.svbtle.com/drawing-lines-is-hard */
result = result + "   vec2 aspectVec = vec2(uAspect, 1.0);\n"+
"   mat4 projViewModel = prMatrix * mvMatrix;\n"+
"   vec4 currentProjected = projViewModel * vec4(aPos, 1.0);\n"+
"   currentProjected = currentProjected/currentProjected.w;\n"+
"   vec4 nextProjected = projViewModel * vec4(aNext, 1.0);\n"+
"   vec2 currentScreen = currentProjected.xy * aspectVec;\n"+
"   vec2 nextScreen = (nextProjected.xy / nextProjected.w) * aspectVec;\n"+
"   float len = uLwd;\n"+
"   vec2 dir = vec2(1.0, 0.0);\n"+
"   vPoint = aPoint;\n"+
"   vLength = length(nextScreen - currentScreen)/2.0;\n"+
"   vLength = vLength/(vLength + len);\n"+
"   if (vLength > 0.0) {\n"+
"     dir = normalize(nextScreen - currentScreen);\n"+
"   }\n"+
"   vec2 normal = vec2(-dir.y, dir.x);\n"+
"   dir.x /= uAspect;\n"+
"   normal.x /= uAspect;\n"+
"   vec4 offset = vec4(len*(normal*aPoint.x*aPoint.y - dir), 0.0, 0.0);\n"+
"   gl_Position = currentProjected + offset;\n";
if (is_brush)
result = result + "   gl_Position = vec4(aPos, 1.);\n";
result = result + "  }\n";
// console.log(result);
return result;
};
/**
* Generate the fragment shader for an object
* @returns {string}
* @param { number } id - Id of object
*/
rglwidgetClass.prototype.getFragmentShader = function(id) {
var obj = this.getObj(id),
userShader = obj.userFragmentShader,
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
sprites_3d = flags & this.f_sprites_3d,
is_twosided = (flags & this.f_is_twosided) > 0,
fat_lines = flags & this.f_fat_lines,
is_transparent = flags & this.f_is_transparent,
nclipplanes = this.countClipplanes(), i,
texture_format, nlights,
result;
if (type === "clipplanes" || sprites_3d) return;
if (typeof userShader !== "undefined") return userShader;
if (has_texture)
texture_format = this.getMaterial(id, "textype");
result = "/* ****** "+type+" object "+id+" fragment shader ****** */\n"+
"#ifdef GL_ES\n"+
"#ifdef GL_FRAGMENT_PRECISION_HIGH\n"+
"  precision highp float;\n"+
"#else\n"+
"  precision mediump float;\n"+
"#endif\n"+
"#endif\n"+
"  varying vec4 vCol; // carries alpha\n"+
"  varying vec4 vPosition;\n";
if (has_texture || type === "text")
result = result + "  varying vec2 vTexcoord;\n"+
" uniform sampler2D uSampler;\n";
if (is_lit && !fixed_quads)
result = result + "  varying vec3 vNormal;\n";
for (i = 0; i < nclipplanes; i++)
result = result + "  uniform vec4 vClipplane"+i+";\n";
if (is_lit) {
nlights = this.countLights();
if (nlights)
result = result + "  uniform mat4 mvMatrix;\n";
else
is_lit = false;
}
if (is_lit) {
result = result + "   uniform vec3 emission;\n"+
"   uniform float shininess;\n";
for (i=0; i < nlights; i++) {
result = result + "   uniform vec3 ambient" + i + ";\n"+
"   uniform vec3 specular" + i +"; // light*material\n"+
"   uniform vec3 diffuse" + i + ";\n"+
"   uniform vec3 lightDir" + i + ";\n"+
"   uniform bool viewpoint" + i + ";\n"+
"   uniform bool finite" + i + ";\n";
}
}
if (is_twosided)
result = result + "   uniform bool front;\n"+
"   varying float normz;\n";
if (fat_lines)
result = result + "   varying vec2 vPoint;\n"+
"   varying float vLength;\n";
result = result + "  void main(void) {\n";
if (fat_lines) {
result = result + "    vec2 point = vPoint;\n"+
"    bool neg = point.y < 0.0;\n"+
"    point.y = neg ? "+
"      (point.y + vLength)/(1.0 - vLength) :\n"+
"     -(point.y - vLength)/(1.0 - vLength);\n";
if (is_transparent && type == "linestrip")
result = result+"    if (neg && length(point) <= 1.0) discard;\n";
result = result + "    point.y = min(point.y, 0.0);\n"+
"    if (length(point) > 1.0) discard;\n";
}
for (i=0; i < nclipplanes;i++)
result = result + "    if (dot(vPosition, vClipplane"+i+") < 0.0) discard;\n";
if (fixed_quads) {
result = result +   "    vec3 n = vec3(0., 0., 1.);\n";
} else if (is_lit) {
result = result +   "    vec3 n = normalize(vNormal);\n";
}
if (is_twosided) {
result = result +   "    if ((normz <= 0.) != front) discard;\n";
}
if (is_lit) {
result = result + "    vec3 eye = normalize(-vPosition.xyz);\n"+
"   vec3 lightdir;\n"+
"   vec4 colDiff;\n"+
"   vec3 halfVec;\n"+
"   vec4 lighteffect = vec4(emission, 0.);\n"+
"   vec3 col;\n"+
"   float nDotL;\n";
if (!fixed_quads) {
result = result +   "   n = -faceforward(n, n, eye);\n";
}
for (i=0; i < nlights; i++) {
result = result + "   colDiff = vec4(vCol.rgb * diffuse" + i + ", vCol.a);\n"+
"   lightdir = lightDir" + i + ";\n"+
"   if (!viewpoint" + i +")\n"+
"     lightdir = (mvMatrix * vec4(lightdir, 1.)).xyz;\n"+
"   if (!finite" + i + ") {\n"+
"     halfVec = normalize(lightdir + eye);\n"+
"   } else {\n"+
"     lightdir = normalize(lightdir - vPosition.xyz);\n"+
"     halfVec = normalize(lightdir + eye);\n"+
"   }\n"+
"    col = ambient" + i + ";\n"+
"   nDotL = dot(n, lightdir);\n"+
"   col = col + max(nDotL, 0.) * colDiff.rgb;\n"+
"   col = col + pow(max(dot(halfVec, n), 0.), shininess) * specular" + i + ";\n"+
"   lighteffect = lighteffect + vec4(col, colDiff.a);\n";
}
} else {
result = result +   "   vec4 colDiff = vCol;\n"+
"    vec4 lighteffect = colDiff;\n";
}
if (type === "text")
result = result +   "    vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n";
if (has_texture) {
result = result + {
rgb:            "   vec4 textureColor = lighteffect*vec4(texture2D(uSampler, vTexcoord).rgb, 1.);\n",
rgba:           "   vec4 textureColor = lighteffect*texture2D(uSampler, vTexcoord);\n",
alpha:          "   vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
"   float luminance = dot(vec3(1.,1.,1.), textureColor.rgb)/3.;\n"+
"   textureColor =  vec4(lighteffect.rgb, lighteffect.a*luminance);\n",
luminance:      "   vec4 textureColor = vec4(lighteffect.rgb*dot(texture2D(uSampler, vTexcoord).rgb, vec3(1.,1.,1.))/3., lighteffect.a);\n",
"luminance.alpha":"    vec4 textureColor = texture2D(uSampler, vTexcoord);\n"+
"   float luminance = dot(vec3(1.,1.,1.),textureColor.rgb)/3.;\n"+
"   textureColor = vec4(lighteffect.rgb*luminance, lighteffect.a*textureColor.a);\n"
}[texture_format]+
"   gl_FragColor = textureColor;\n";
} else if (type === "text") {
result = result +   "    if (textureColor.a < 0.1)\n"+
"     discard;\n"+
"   else\n"+
"     gl_FragColor = textureColor;\n";
} else
result = result +   "   gl_FragColor = lighteffect;\n";
//if (fat_lines)
//  result = result +   "   gl_FragColor = vec4(0.0, abs(point.x), abs(point.y), 1.0);"
result = result + "  }\n";
// console.log(result);
return result;
};
/**
* Call gl functions to create and compile shader
* @returns {Object}
* @param { number } shaderType - gl code for shader type
* @param { string } code - code for the shader
*/
rglwidgetClass.prototype.getShader = function(shaderType, code) {
var gl = this.gl, shader;
shader = gl.createShader(shaderType);
gl.shaderSource(shader, code);
gl.compileShader(shader);
if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS) && !gl.isContextLost())
alert(gl.getShaderInfoLog(shader));
return shader;
};
/**
* Handle a texture after its image has been loaded
* @param { Object } texture - the gl texture object
* @param { Object } textureCanvas - the canvas holding the image
*/
rglwidgetClass.prototype.handleLoadedTexture = function(texture, textureCanvas) {
var gl = this.gl || this.initGL();
gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
gl.bindTexture(gl.TEXTURE_2D, texture);
gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, textureCanvas);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
gl.generateMipmap(gl.TEXTURE_2D);
gl.bindTexture(gl.TEXTURE_2D, null);
};
/**
* Get maximum dimension of texture in current browser.
* @returns {number}
*/
rglwidgetClass.prototype.getMaxTexSize = function() {
var gl = this.gl || this.initGL();  
return Math.min(4096, gl.getParameter(gl.MAX_TEXTURE_SIZE));
};
/**
* Load an image to a texture
* @param { string } uri - The image location
* @param { Object } texture - the gl texture object
*/
rglwidgetClass.prototype.loadImageToTexture = function(uri, texture) {
var canvas = this.textureCanvas,
ctx = canvas.getContext("2d"),
image = new Image(),
self = this;
image.onload = function() {
var w = image.width,
h = image.height,
canvasX = self.getPowerOfTwo(w),
canvasY = self.getPowerOfTwo(h),
gl = self.gl || self.initGL(),
maxTexSize = self.getMaxTexSize();
while (canvasX > 1 && canvasY > 1 && (canvasX > maxTexSize || canvasY > maxTexSize)) {
canvasX /= 2;
canvasY /= 2;
}
canvas.width = canvasX;
canvas.height = canvasY;
ctx.imageSmoothingEnabled = true;
ctx.drawImage(image, 0, 0, canvasX, canvasY);
self.handleLoadedTexture(texture, canvas);
self.drawScene();
};
image.src = uri;
};
/**
* Draw text to the texture canvas
* @returns { Object } object with text measurements
* @param { string } text - the text
* @param { number } cex - expansion
* @param { string } family - font family
* @param { number } font - font number
*/
rglwidgetClass.prototype.drawTextToCanvas = function(text, cex, family, font) {
var canvasX, canvasY,
textY,
scaling = 20,
textColour = "white",
backgroundColour = "rgba(0,0,0,0)",
canvas = this.textureCanvas,
ctx = canvas.getContext("2d"),
i, textHeight = 0, textHeights = [], width, widths = [], 
offsetx, offsety = 0, line, lines = [], offsetsx = [],
offsetsy = [], lineoffsetsy = [], fontStrings = [],
maxTexSize = this.getMaxTexSize(),
getFontString = function(i) {
textHeights[i] = scaling*cex[i];
var fontString = textHeights[i] + "px",
family0 = family[i],
font0 = font[i];
if (family0 === "sans")
family0 = "sans-serif";
else if (family0 === "mono")
family0 = "monospace";
fontString = fontString + " " + family0;
if (font0 === 2 || font0 === 4)
fontString = "bold " + fontString;
if (font0 === 3 || font0 === 4)
fontString = "italic " + fontString;
return fontString;
};
cex = this.repeatToLen(cex, text.length);
family = this.repeatToLen(family, text.length);
font = this.repeatToLen(font, text.length);
canvasX = 1;
line = -1;
offsetx = maxTexSize;
for (i = 0; i < text.length; i++)  {
ctx.font = fontStrings[i] = getFontString(i);
width = widths[i] = ctx.measureText(text[i]).width;
if (offsetx + width > maxTexSize) {
line += 1;
offsety = lineoffsetsy[line] = offsety + 2*textHeight;
if (offsety > maxTexSize)
console.error("Too many strings for texture.");
textHeight = 0;
offsetx = 0;
}
textHeight = Math.max(textHeight, textHeights[i]);
offsetsx[i] = offsetx;
offsetx += width;
canvasX = Math.max(canvasX, offsetx);
lines[i] = line;
}
offsety = lineoffsetsy[line] = offsety + 2*textHeight;
for (i = 0; i < text.length; i++) {
offsetsy[i] = lineoffsetsy[lines[i]];
}
canvasX = this.getPowerOfTwo(canvasX);
canvasY = this.getPowerOfTwo(offsety);
canvas.width = canvasX;
canvas.height = canvasY;
ctx.fillStyle = backgroundColour;
ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);
ctx.textBaseline = "alphabetic";
for(i = 0; i < text.length; i++) {
ctx.font = fontStrings[i];
ctx.fillStyle = textColour;
ctx.textAlign = "left";
ctx.fillText(text[i], offsetsx[i],  offsetsy[i]);
}
return {canvasX:canvasX, canvasY:canvasY,
widths:widths, textHeights:textHeights,
offsetsx:offsetsx, offsetsy:offsetsy};
};
/**
* Set the gl viewport and scissor test
* @param { number } id - id of subscene
*/
rglwidgetClass.prototype.setViewport = function(id) {
var gl = this.gl || this.initGL(),
vp = this.getObj(id).par3d.viewport,
x = vp.x*this.canvas.width,
y = vp.y*this.canvas.height,
width = vp.width*this.canvas.width,
height = vp.height*this.canvas.height;
this.vp = {x:x, y:y, width:width, height:height};
gl.viewport(x, y, width, height);
gl.scissor(x, y, width, height);
gl.enable(gl.SCISSOR_TEST);
};
/**
* Set the projection matrix for a subscene
* @param { number } id - id of subscene
*/
rglwidgetClass.prototype.setprMatrix = function(id) {
var subscene = this.getObj(id),
embedding = subscene.embeddings.projection;
if (embedding === "replace")
this.prMatrix.makeIdentity();
else
this.setprMatrix(subscene.parent);
if (embedding === "inherit")
return;
// This is based on the Frustum::enclose code from geom.cpp
var bbox = subscene.par3d.bbox,
scale = subscene.par3d.scale,
ranges = [(bbox[1]-bbox[0])*scale[0]/2,
(bbox[3]-bbox[2])*scale[1]/2,
(bbox[5]-bbox[4])*scale[2]/2],
radius = Math.sqrt(this.sumsq(ranges))*1.1; // A bit bigger to handle labels
if (radius <= 0) radius = 1;
var observer = subscene.par3d.observer,
distance = observer[2],
FOV = subscene.par3d.FOV, ortho = FOV === 0,
t = ortho ? 1 : Math.tan(FOV*Math.PI/360),
near = distance - radius,
far = distance + radius,
hlen,
aspect = this.vp.width/this.vp.height,
z = subscene.par3d.zoom,
userProjection = subscene.par3d.userProjection;
if (far < 0.0)
far = 1.0;
if (near < far/100.0)
near = far/100.0;
hlen = t*near;
if (ortho) {
if (aspect > 1)
this.prMatrix.ortho(-hlen*aspect*z, hlen*aspect*z,
-hlen*z, hlen*z, near, far);
else
this.prMatrix.ortho(-hlen*z, hlen*z,
-hlen*z/aspect, hlen*z/aspect,
near, far);
} else {
if (aspect > 1)
this.prMatrix.frustum(-hlen*aspect*z, hlen*aspect*z,
-hlen*z, hlen*z, near, far);
else
this.prMatrix.frustum(-hlen*z, hlen*z,
-hlen*z/aspect, hlen*z/aspect,
near, far);
}
this.prMatrix.multRight(userProjection);
};
/**
* Set the model-view matrix for a subscene
* @param { number } id - id of the subscene
*/
rglwidgetClass.prototype.setmvMatrix = function(id) {
var observer = this.getObj(id).par3d.observer;
this.mvMatrix.makeIdentity();
this.setmodelMatrix(id);
this.mvMatrix.translate(-observer[0], -observer[1], -observer[2]);
};
/**
* Set the model matrix for a subscene
* @param { number } id - id of the subscene
*/
rglwidgetClass.prototype.setmodelMatrix = function(id) {
var subscene = this.getObj(id),
embedding = subscene.embeddings.model;
if (embedding !== "inherit") {
var scale = subscene.par3d.scale,
bbox = subscene.par3d.bbox,
center = [(bbox[0]+bbox[1])/2,
(bbox[2]+bbox[3])/2,
(bbox[4]+bbox[5])/2];
this.mvMatrix.translate(-center[0], -center[1], -center[2]);
this.mvMatrix.scale(scale[0], scale[1], scale[2]);
this.mvMatrix.multRight( subscene.par3d.userMatrix );
}
if (embedding !== "replace")
this.setmodelMatrix(subscene.parent);
};
/**
* Set the normals matrix for a subscene
* @param { number } subsceneid - id of the subscene
*/
rglwidgetClass.prototype.setnormMatrix = function(subsceneid) {
var self = this,
recurse = function(id) {
var sub = self.getObj(id),
embedding = sub.embeddings.model;
if (embedding !== "inherit") {
var scale = sub.par3d.scale;
self.normMatrix.scale(1/scale[0], 1/scale[1], 1/scale[2]);
self.normMatrix.multRight(sub.par3d.userMatrix);
}
if (embedding !== "replace")
recurse(sub.parent);
};
self.normMatrix.makeIdentity();
recurse(subsceneid);
};
/**
* Set the combined projection-model-view matrix
*/
rglwidgetClass.prototype.setprmvMatrix = function() {
this.prmvMatrix = new CanvasMatrix4( this.mvMatrix );
this.prmvMatrix.multRight( this.prMatrix );
};
/**
* Count clipping planes in a scene
* @returns {number}
*/
rglwidgetClass.prototype.countClipplanes = function() {
return this.countObjs("clipplanes");
};
/**
* Count lights in a scene
* @returns { number }
*/
rglwidgetClass.prototype.countLights = function() {
return this.countObjs("light");
};
/**
* Count objects of specific type in a scene
* @returns { number }
* @param { string } type - Type of object to count
*/
rglwidgetClass.prototype.countObjs = function(type) {
var self = this,
bound = 0;
Object.keys(this.scene.objects).forEach(
function(key) {
if (self.getObj(parseInt(key, 10)).type === type)
bound = bound + 1;
});
return bound;
};
/**
* Initialize a subscene
* @param { number } id - id of subscene.
*/
rglwidgetClass.prototype.initSubscene = function(id) {
var sub = this.getObj(id),
i, obj;
if (sub.type !== "subscene")
return;
sub.par3d.userMatrix = this.toCanvasMatrix4(sub.par3d.userMatrix);
sub.par3d.userProjection = this.toCanvasMatrix4(sub.par3d.userProjection);
sub.par3d.userProjection.transpose();
sub.par3d.listeners = [].concat(sub.par3d.listeners);
sub.backgroundId = undefined;
sub.subscenes = [];
sub.clipplanes = [];
sub.transparent = [];
sub.opaque = [];
sub.lights = [];
for (i=0; i < sub.objects.length; i++) {
obj = this.getObj(sub.objects[i]);
if (typeof obj === "undefined") {
sub.objects.splice(i, 1);
i--;
} else if (obj.type === "background")
sub.backgroundId = obj.id;
else
sub[this.whichList(obj.id)].push(obj.id);
}
};
/**
* Copy object
* @param { number } id - id of object to copy
* @param { string } reuse - Document id of scene to reuse
*/
rglwidgetClass.prototype.copyObj = function(id, reuse) {
var obj = this.getObj(id),
prev = document.getElementById(reuse);
if (prev !== null) {
prev = prev.rglinstance;
var
prevobj = prev.getObj(id),
fields = ["flags", "type",
"colors", "vertices", "centers",
"normals", "offsets",
"texts", "cex", "family", "font", "adj",
"material",
"radii",
"texcoords",
"userMatrix", "ids",
"dim",
"par3d", "userMatrix",
"viewpoint", "finite",
"pos"],
i;
for (i = 0; i < fields.length; i++) {
if (typeof prevobj[fields[i]] !== "undefined")
obj[fields[i]] = prevobj[fields[i]];
}
} else
console.warn("copyObj failed");
};
/**
* Update the triangles used to display a plane
* @param { number } id - id of the plane
* @param { Object } bbox - bounding box in which to display the plane
*/
rglwidgetClass.prototype.planeUpdateTriangles = function(id, bbox) {
var perms = [[0,0,1], [1,2,2], [2,1,0]],
x, xrow, elem, A, d, nhits, i, j, k, u, v, w, intersect, which, v0, v2, vx, reverse,
face1 = [], face2 = [], normals = [],
obj = this.getObj(id),
nPlanes = obj.normals.length;
obj.bbox = bbox;
obj.vertices = [];
obj.initialized = false;
for (elem = 0; elem < nPlanes; elem++) {
//    Vertex Av = normal.getRecycled(elem);
x = [];
A = obj.normals[elem];
d = obj.offsets[elem][0];
nhits = 0;
for (i=0; i<3; i++)
for (j=0; j<2; j++)
for (k=0; k<2; k++) {
u = perms[0][i];
v = perms[1][i];
w = perms[2][i];
if (A[w] !== 0.0) {
intersect = -(d + A[u]*bbox[j+2*u] + A[v]*bbox[k+2*v])/A[w];
if (bbox[2*w] < intersect && intersect < bbox[1+2*w]) {
xrow = [];
xrow[u] = bbox[j+2*u];
xrow[v] = bbox[k+2*v];
xrow[w] = intersect;
x.push(xrow);
face1[nhits] = j + 2*u;
face2[nhits] = k + 2*v;
nhits++;
}
}
}
if (nhits > 3) {
/* Re-order the intersections so the triangles work */
for (i=0; i<nhits-2; i++) {
which = 0; /* initialize to suppress warning */
for (j=i+1; j<nhits; j++) {
if (face1[i] == face1[j] || face1[i] == face2[j] ||
face2[i] == face1[j] || face2[i] == face2[j] ) {
which = j;
break;
}
}
if (which > i+1) {
this.swap(x, i+1, which);
this.swap(face1, i+1, which);
this.swap(face2, i+1, which);
}
}
}
if (nhits >= 3) {
/* Put in order so that the normal points out the FRONT of the faces */
v0 = [x[0][0] - x[1][0] , x[0][1] - x[1][1], x[0][2] - x[1][2]];
v2 = [x[2][0] - x[1][0] , x[2][1] - x[1][1], x[2][2] - x[1][2]];
/* cross-product */
vx = this.xprod(v0, v2);
reverse = this.dotprod(vx, A) > 0;
for (i=0; i<nhits-2; i++) {
obj.vertices.push(x[0]);
normals.push(A);
for (j=1; j<3; j++) {
obj.vertices.push(x[i + (reverse ? 3-j : j)]);
normals.push(A);
}
}
}
}
obj.pnormals = normals;
};
rglwidgetClass.prototype.getAdj = function (pos, offset, text) {
switch(pos) {
case 1: return [0.5, 1 + offset];
case 2: return [1 + offset/text.length, 0.5];
case 3: return [0.5, -offset];
case 4: return [-offset/text.length, 0.5];
}
}
/**
* Initialize object for display
* @param { number } id - id of object to initialize
*/
rglwidgetClass.prototype.initObj = function(id) {
var obj = this.getObj(id),
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
is_lines = flags & this.f_is_lines,
fat_lines = flags & this.f_fat_lines,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
is_transparent = obj.is_transparent,
depth_sort = flags & this.f_depth_sort,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
fixed_size = flags & this.f_fixed_size,
is_twosided = (flags & this.f_is_twosided) > 0,
is_brush = flags & this.f_is_brush,
gl = this.gl || this.initGL(),
polygon_offset,
texinfo, drawtype, nclipplanes, f, nrows, oldrows,
i,j,v,v1,v2, mat, uri, matobj, pass, passes, pmode,
dim, nx, nz, attr;
if (typeof id !== "number") {
this.alertOnce("initObj id is "+typeof id);
}
obj.initialized = true;
if (type === "bboxdeco" || type === "subscene")
return;
if (type === "light") {
obj.ambient = new Float32Array(obj.colors[0].slice(0,3));
obj.diffuse = new Float32Array(obj.colors[1].slice(0,3));
obj.specular = new Float32Array(obj.colors[2].slice(0,3));
obj.lightDir = new Float32Array(obj.vertices[0]);
return;
}
if (type === "clipplanes") {
obj.vClipplane = this.flatten(this.cbind(obj.normals, obj.offsets));
return;
}
if (type === "background" && typeof obj.ids !== "undefined") {
obj.quad = this.flatten([].concat(obj.ids));
return;
}
polygon_offset = this.getMaterial(id, "polygon_offset");
if (polygon_offset[0] != 0 || polygon_offset[1] != 0)
obj.polygon_offset = polygon_offset;
if (is_transparent) {
depth_sort = ["triangles", "quads", "surface",
"spheres", "sprites", "text"].indexOf(type) >= 0;
}
if (is_brush)
this.initSelection(id);
if (typeof obj.vertices === "undefined")
obj.vertices = [];
v = obj.vertices;
obj.vertexCount = v.length;
if (!obj.vertexCount) return;
if (is_twosided) {
if (typeof obj.userAttributes === "undefined")
obj.userAttributes = {};
v1 = Array(v.length);
v2 = Array(v.length);
if (obj.type == "triangles" || obj.type == "quads") {
if (obj.type == "triangles")
nrow = 3;
else
nrow = 4;
for (i=0; i<Math.floor(v.length/nrow); i++)
for (j=0; j<nrow; j++) {
v1[nrow*i + j] = v[nrow*i + ((j+1) % nrow)];
v2[nrow*i + j] = v[nrow*i + ((j+2) % nrow)];
}
} else if (obj.type == "surface") {
dim = obj.dim[0];
nx = dim[0];
nz = dim[1];
for (j=0; j<nx; j++) {
for (i=0; i<nz; i++) {
if (i+1 < nz && j+1 < nx) {
v2[j + nx*i] = v[j + nx*(i+1)];
v1[j + nx*i] = v[j+1 + nx*(i+1)];
} else if (i+1 < nz) {
v2[j + nx*i] = v[j-1 + nx*i];
v1[j + nx*i] = v[j + nx*(i+1)];
} else {
v2[j + nx*i] = v[j + nx*(i-1)];
v1[j + nx*i] = v[j-1 + nx*(i-1)];
}
}
}
}
obj.userAttributes.aPos1 = v1;
obj.userAttributes.aPos2 = v2;
}
if (!sprites_3d) {
if (gl.isContextLost()) return;
obj.prog = gl.createProgram();
gl.attachShader(obj.prog, this.getShader( gl.VERTEX_SHADER,
this.getVertexShader(id) ));
gl.attachShader(obj.prog, this.getShader( gl.FRAGMENT_SHADER,
this.getFragmentShader(id) ));
//  Force aPos to location 0, aCol to location 1
gl.bindAttribLocation(obj.prog, 0, "aPos");
gl.bindAttribLocation(obj.prog, 1, "aCol");
gl.linkProgram(obj.prog);
var linked = gl.getProgramParameter(obj.prog, gl.LINK_STATUS);
if (!linked) {
// An error occurred while linking
var lastError = gl.getProgramInfoLog(obj.prog);
console.warn("Error in program linking:" + lastError);
gl.deleteProgram(obj.prog);
return;
}
}
if (type === "text") {
texinfo = this.drawTextToCanvas(obj.texts,
this.flatten(obj.cex),
this.flatten(obj.family),
this.flatten(obj.family));
}
if (fixed_quads && !sprites_3d) {
obj.ofsLoc = gl.getAttribLocation(obj.prog, "aOfs");
}
if (sprite_3d) {
obj.origLoc = gl.getUniformLocation(obj.prog, "uOrig");
obj.sizeLoc = gl.getUniformLocation(obj.prog, "uSize");
obj.usermatLoc = gl.getUniformLocation(obj.prog, "usermat");
}
if (has_texture || type == "text") {
if (!obj.texture)
obj.texture = gl.createTexture();
obj.texLoc = gl.getAttribLocation(obj.prog, "aTexcoord");
obj.sampler = gl.getUniformLocation(obj.prog, "uSampler");
}
if (has_texture) {
mat = obj.material;
if (typeof mat.uri !== "undefined")
uri = mat.uri;
else if (typeof mat.uriElementId === "undefined") {
matobj = this.getObj(mat.uriId);
if (typeof matobj !== "undefined") {
uri = matobj.material.uri;
} else {
uri = "";
}
} else
uri = document.getElementById(mat.uriElementId).rglinstance.getObj(mat.uriId).material.uri;
this.loadImageToTexture(uri, obj.texture);
}
if (type === "text") {
this.handleLoadedTexture(obj.texture, this.textureCanvas);
}
var stride = 3, nc, cofs, nofs, radofs, oofs, tofs, vnew, fnew,
nextofs = -1, pointofs = -1, alias, colors, key, selection, filter, adj, pos, offset;
obj.alias = undefined;
colors = obj.colors;
j = this.scene.crosstalk.id.indexOf(id);
if (j >= 0) {
key = this.scene.crosstalk.key[j];
options = this.scene.crosstalk.options[j];
colors = colors.slice(0); 
for (i = 0; i < v.length; i++)
colors[i] = obj.colors[i % obj.colors.length].slice(0);
if ( (selection = this.scene.crosstalk.selection) &&
(selection.length || !options.selectedIgnoreNone) )
for (i = 0; i < v.length; i++) {
if (!selection.includes(key[i])) {
if (options.deselectedColor)
colors[i] = options.deselectedColor.slice(0);
colors[i][3] = colors[i][3]*options.deselectedFade;   /* default: mostly transparent if not selected */
} else if (options.selectedColor)
colors[i] = options.selectedColor.slice(0);
}
if ( (filter = this.scene.crosstalk.filter) )
for (i = 0; i < v.length; i++) 
if (!filter.includes(key[i])) {
if (options.filteredColor)
colors[i] = options.filteredColor.slice(0);
colors[i][3] = colors[i][3]*options.filteredFade;   /* default: completely hidden if filtered */
}
}  
nc = obj.colorCount = colors.length;
if (nc > 1) {
cofs = stride;
stride = stride + 4;
v = this.cbind(v, colors);
} else {
cofs = -1;
obj.onecolor = this.flatten(colors);
}
if (typeof obj.normals !== "undefined") {
nofs = stride;
stride = stride + 3;
v = this.cbind(v, typeof obj.pnormals !== "undefined" ? obj.pnormals : obj.normals);
} else
nofs = -1;
if (typeof obj.radii !== "undefined") {
radofs = stride;
stride = stride + 1;
// FIXME:  always concat the radii?
if (obj.radii.length === v.length) {
v = this.cbind(v, obj.radii);
} else if (obj.radii.length === 1) {
v = v.map(function(row, i, arr) { return row.concat(obj.radii[0]);});
}
} else
radofs = -1;
// Add default indices
f = Array(v.length);
for (i = 0; i < v.length; i++)
f[i] = i;
obj.f = [f,f];
if (type == "sprites" && !sprites_3d) {
tofs = stride;
stride += 2;
oofs = stride;
stride += 2;
vnew = new Array(4*v.length);
fnew = new Array(4*v.length);
alias = new Array(v.length);
var rescale = fixed_size ? 72 : 1,
size = obj.radii, s = rescale*size[0]/2;
last = v.length;
f = obj.f[0];
for (i=0; i < v.length; i++) {
if (size.length > 1)
s = rescale*size[i]/2;
vnew[i]  = v[i].concat([0,0,-s,-s]);
fnew[4*i] = f[i];
vnew[last]= v[i].concat([1,0, s,-s]);
fnew[4*i+1] = last++;
vnew[last]= v[i].concat([1,1, s, s]);
fnew[4*i+2] = last++;
vnew[last]= v[i].concat([0,1,-s, s]);
fnew[4*i+3] = last++;
alias[i] = [last-3, last-2, last-1];
}
v = vnew;
obj.vertexCount = v.length;
obj.f = [fnew, fnew];
} else if (type === "text") {
tofs = stride;
stride += 2;
oofs = stride;
stride += 2;
vnew = new Array(4*v.length);
f = obj.f[0];
fnew = new Array(4*f.length);
alias = new Array(v.length);
last = v.length;
adj = this.flatten(obj.adj);
if (typeof obj.pos !== "undefined") {
pos = this.flatten(obj.pos);
offset = adj[0];
}
for (i=0; i < v.length; i++) {
if (typeof pos !== "undefined")
adj = this.getAdj(pos[i % pos.length], offset, obj.texts[i]);
vnew[i]  = v[i].concat([0,-0.5]).concat(adj);
fnew[4*i] = f[i];
vnew[last] = v[i].concat([1,-0.5]).concat(adj);
fnew[4*i+1] = last++;
vnew[last] = v[i].concat([1, 1.5]).concat(adj);
fnew[4*i+2] = last++;
vnew[last] = v[i].concat([0, 1.5]).concat(adj);
fnew[4*i+3] = last++;
alias[i] = [last-3, last-2, last-1];
for (j=0; j < 4; j++) {
v1 = vnew[fnew[4*i+j]];
v1[tofs+2] = 2*(v1[tofs]-v1[tofs+2])*texinfo.widths[i];
v1[tofs+3] = 2*(v1[tofs+1]-v1[tofs+3])*texinfo.textHeights[i];
v1[tofs] = (texinfo.offsetsx[i] + v1[tofs]*texinfo.widths[i])/texinfo.canvasX;
v1[tofs+1] = 1.0-(texinfo.offsetsy[i] -
v1[tofs+1]*texinfo.textHeights[i])/texinfo.canvasY;
vnew[fnew[4*i+j]] = v1;
}
}
v = vnew;
obj.vertexCount = v.length;
obj.f = [fnew, fnew];
} else if (typeof obj.texcoords !== "undefined") {
tofs = stride;
stride += 2;
oofs = -1;
v = this.cbind(v, obj.texcoords);
} else {
tofs = -1;
oofs = -1;
}
obj.alias = alias;
if (typeof obj.userAttributes !== "undefined") {
obj.userAttribOffsets = {};
obj.userAttribLocations = {};
obj.userAttribSizes = {};
for (attr in obj.userAttributes) {
obj.userAttribLocations[attr] = gl.getAttribLocation(obj.prog, attr);
if (obj.userAttribLocations[attr] >= 0) { // Attribute may not have been used
obj.userAttribOffsets[attr] = stride;
v = this.cbind(v, obj.userAttributes[attr]);
stride = v[0].length;
obj.userAttribSizes[attr] = stride - obj.userAttribOffsets[attr];
}
}
}
if (typeof obj.userUniforms !== "undefined") {
obj.userUniformLocations = {};
for (attr in obj.userUniforms)
obj.userUniformLocations[attr] = gl.getUniformLocation(obj.prog, attr);
}
if (sprites_3d) {
obj.userMatrix = new CanvasMatrix4(obj.userMatrix);
obj.objects = this.flatten([].concat(obj.ids));
is_lit = false;
for (i=0; i < obj.objects.length; i++)
this.initObj(obj.objects[i]);
}
if (is_lit && !fixed_quads) {
obj.normLoc = gl.getAttribLocation(obj.prog, "aNorm");
}
nclipplanes = this.countClipplanes();
if (nclipplanes && !sprites_3d) {
obj.clipLoc = [];
for (i=0; i < nclipplanes; i++)
obj.clipLoc[i] = gl.getUniformLocation(obj.prog,"vClipplane" + i);
}
if (is_lit) {
obj.emissionLoc = gl.getUniformLocation(obj.prog, "emission");
obj.emission = new Float32Array(this.stringToRgb(this.getMaterial(id, "emission")));
obj.shininessLoc = gl.getUniformLocation(obj.prog, "shininess");
obj.shininess = this.getMaterial(id, "shininess");
obj.nlights = this.countLights();
obj.ambientLoc = [];
obj.ambient = new Float32Array(this.stringToRgb(this.getMaterial(id, "ambient")));
obj.specularLoc = [];
obj.specular = new Float32Array(this.stringToRgb(this.getMaterial(id, "specular")));
obj.diffuseLoc = [];
obj.lightDirLoc = [];
obj.viewpointLoc = [];
obj.finiteLoc = [];
for (i=0; i < obj.nlights; i++) {
obj.ambientLoc[i] = gl.getUniformLocation(obj.prog, "ambient" + i);
obj.specularLoc[i] = gl.getUniformLocation(obj.prog, "specular" + i);
obj.diffuseLoc[i] = gl.getUniformLocation(obj.prog, "diffuse" + i);
obj.lightDirLoc[i] = gl.getUniformLocation(obj.prog, "lightDir" + i);
obj.viewpointLoc[i] = gl.getUniformLocation(obj.prog, "viewpoint" + i);
obj.finiteLoc[i] = gl.getUniformLocation(obj.prog, "finite" + i);
}
}
obj.passes = is_twosided + 1;
obj.pmode = new Array(obj.passes);
for (pass = 0; pass < obj.passes; pass++) {
if (type === "triangles" || type === "quads" || type === "surface")
pmode = this.getMaterial(id, (pass === 0) ? "front" : "back");
else pmode = "filled";
obj.pmode[pass] = pmode;
}
obj.f.length = obj.passes;
for (pass = 0; pass < obj.passes; pass++) {
f = fnew = obj.f[pass];
pmode = obj.pmode[pass];
if (pmode === "culled")
f = [];
else if (pmode === "points") {
// stay with default
} else if ((type === "quads" || type === "text" ||
type === "sprites") && !sprites_3d) {
nrows = Math.floor(obj.vertexCount/4);
if (pmode === "filled") {
fnew = Array(6*nrows);
for (i=0; i < nrows; i++) {
fnew[6*i] = f[4*i];
fnew[6*i+1] = f[4*i + 1];
fnew[6*i+2] = f[4*i + 2];
fnew[6*i+3] = f[4*i];
fnew[6*i+4] = f[4*i + 2];
fnew[6*i+5] = f[4*i + 3];
}
} else {
fnew = Array(8*nrows);
for (i=0; i < nrows; i++) {
fnew[8*i] = f[4*i];
fnew[8*i+1] = f[4*i + 1];
fnew[8*i+2] = f[4*i + 1];
fnew[8*i+3] = f[4*i + 2];
fnew[8*i+4] = f[4*i + 2];
fnew[8*i+5] = f[4*i + 3];
fnew[8*i+6] = f[4*i + 3];
fnew[8*i+7] = f[4*i];
}
}
} else if (type === "triangles") {
nrows = Math.floor(obj.vertexCount/3);
if (pmode === "filled") {
fnew = Array(3*nrows);
for (i=0; i < fnew.length; i++) {
fnew[i] = f[i];
}
} else if (pmode === "lines") {
fnew = Array(6*nrows);
for (i=0; i < nrows; i++) {
fnew[6*i] = f[3*i];
fnew[6*i + 1] = f[3*i + 1];
fnew[6*i + 2] = f[3*i + 1];
fnew[6*i + 3] = f[3*i + 2];
fnew[6*i + 4] = f[3*i + 2];
fnew[6*i + 5] = f[3*i];
}
}
} else if (type === "spheres") {
// default
} else if (type === "surface") {
dim = obj.dim[0];
nx = dim[0];
nz = dim[1];
if (pmode === "filled") {
fnew = [];
for (j=0; j<nx-1; j++) {
for (i=0; i<nz-1; i++) {
fnew.push(f[j + nx*i],
f[j + nx*(i+1)],
f[j + 1 + nx*(i+1)],
f[j + nx*i],
f[j + 1 + nx*(i+1)],
f[j + 1 + nx*i]);
}
}
} else if (pmode === "lines") {
fnew = [];
for (j=0; j<nx; j++) {
for (i=0; i<nz; i++) {
if (i+1 < nz)
fnew.push(f[j + nx*i],
f[j + nx*(i+1)]);
if (j+1 < nx)
fnew.push(f[j + nx*i],
f[j+1 + nx*i]);
}
}
}
}
obj.f[pass] = fnew;
if (depth_sort) {
drawtype = "DYNAMIC_DRAW";
} else {
drawtype = "STATIC_DRAW";
}
}
if (fat_lines) {
alias = undefined;
obj.nextLoc = gl.getAttribLocation(obj.prog, "aNext");
obj.pointLoc = gl.getAttribLocation(obj.prog, "aPoint");
obj.aspectLoc = gl.getUniformLocation(obj.prog, "uAspect");
obj.lwdLoc = gl.getUniformLocation(obj.prog, "uLwd");
// Expand vertices to turn each segment into a pair of triangles
for (pass = 0; pass < obj.passes; pass++) {
f = obj.f[pass];    
oldrows = f.length;
if (obj.pmode[pass] === "lines") 
break;
}
if (type === "linestrip") 
nrows = 4*(oldrows - 1); 
else
nrows = 2*oldrows;
vnew = new Array(nrows);
fnew = new Array(1.5*nrows);
var fnext = new Array(nrows),
fpt = new Array(nrows), 
pt, start, gap = type === "linestrip" ? 3 : 1;
// We're going to turn each pair of vertices into 4 new ones, with the "next" and "pt" attributes
// added.
// We do this by copying the originals in the first pass, adding the new attributes, then in a 
// second pass add new vertices at the end.
for (i = 0; i < v.length; i++) {
vnew[i] = v[i].concat([0,0,0,0,0]); 
}
nextofs = stride;
pointofs = stride + 3;
stride = stride + 5;
// Now add the extras
last = v.length - 1;
ind = 0;
alias = new Array(f.length);
for (i = 0; i < f.length; i++)
alias[i] = [];
for (i = 0; i < f.length - 1; i++) {
if (type !== "linestrip" && i % 2 == 1)
continue;
k = ++last;
vnew[k] = vnew[f[i]].slice();
for (j=0; j<3; j++)
vnew[k][nextofs + j] = vnew[f[i+1]][j];
vnew[k][pointofs] = -1;
vnew[k][pointofs+1] = -1;
fnew[ind] = k;
last++;
vnew[last] = vnew[k].slice();
vnew[last][pointofs] = 1;
fnew[ind+1] = last;
alias[f[i]].push(last-1, last);
last++;
k = last;
vnew[k] = vnew[f[i+1]].slice();
for (j=0; j<3; j++)
vnew[k][nextofs + j] = vnew[f[i]][j];
vnew[k][pointofs] = -1;
vnew[k][pointofs+1] = 1;
fnew[ind+2] = k;
fnew[ind+3] = fnew[ind+1];
last++;
vnew[last] = vnew[k].slice();
vnew[last][pointofs] = 1;
fnew[ind+4] = last;
fnew[ind+5] = fnew[ind+2];
ind += 6;
alias[f[i+1]].push(last-1, last);
}
vnew.length = last+1;
v = vnew;
obj.vertexCount = v.length;
if (typeof alias !== "undefined" && typeof obj.alias !== "undefined") {  // Already have aliases from previous section?
var oldalias = obj.alias, newalias = Array(obj.alias.length);
for (i = 0; i < newalias.length; i++) {
newalias[i] = oldalias[i].slice();
for (j = 0; j < oldalias[i].length; j++)
Array.prototype.push.apply(newalias[i], alias[oldalias[j]]); // pushes each element 
}
obj.alias = newalias;
} else
obj.alias = alias;
for (pass = 0; pass < obj.passes; pass++)
if (type === "lines" || type === "linestrip" || obj.pmode[pass] == "lines") {
obj.f[pass] = fnew;
}
if (depth_sort) 
drawtype = "DYNAMIC_DRAW";
else
drawtype = "STATIC_DRAW";
}
for (pass = 0; pass < obj.passes; pass++) {
if (obj.vertexCount > 65535) {
if (this.index_uint) {
obj.f[pass] = new Uint32Array(obj.f[pass]);
obj.index_uint = true;
} else
this.alertOnce("Object has "+obj.vertexCount+" vertices, not supported in this browser.");
} else {
obj.f[pass] = new Uint16Array(obj.f[pass]);
obj.index_uint = false;
}
}
if (stride !== v[0].length) {
this.alertOnce("problem in stride calculation");
}
obj.vOffsets = {vofs:0, cofs:cofs, nofs:nofs, radofs:radofs, oofs:oofs, tofs:tofs,
nextofs:nextofs, pointofs:pointofs, stride:stride};
obj.values = new Float32Array(this.flatten(v));
if (type !== "spheres" && !sprites_3d) {
obj.buf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW); //
obj.ibuf = Array(obj.passes);
obj.ibuf[0] = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[0]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[0], gl[drawtype]);
if (is_twosided) {
obj.ibuf[1] = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[1]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, obj.f[1], gl[drawtype]);
}
}
if (!sprites_3d) {
obj.mvMatLoc = gl.getUniformLocation(obj.prog, "mvMatrix");
obj.prMatLoc = gl.getUniformLocation(obj.prog, "prMatrix");
}
if (fixed_size) {
obj.textScaleLoc = gl.getUniformLocation(obj.prog, "textScale");
}
if (is_lit && !sprites_3d) {
obj.normMatLoc = gl.getUniformLocation(obj.prog, "normMatrix");
}
if (is_twosided) {
obj.frontLoc = gl.getUniformLocation(obj.prog, "front");
}
};
/**
* Set gl depth test based on object's material
* @param { number } id - object to use
*/
rglwidgetClass.prototype.setDepthTest = function(id) {
var gl = this.gl || this.initGL(),
tests = {never: gl.NEVER,
less:  gl.LESS,
equal: gl.EQUAL,
lequal:gl.LEQUAL,
greater: gl.GREATER,
notequal: gl.NOTEQUAL,
gequal: gl.GEQUAL,
always: gl.ALWAYS},
test = tests[this.getMaterial(id, "depth_test")];
gl.depthFunc(test);
};
rglwidgetClass.prototype.mode4type = {points : "POINTS",
linestrip : "LINE_STRIP",
abclines : "LINES",
lines : "LINES",
sprites : "TRIANGLES",
planes : "TRIANGLES",
text : "TRIANGLES",
quads : "TRIANGLES",
surface : "TRIANGLES",
triangles : "TRIANGLES"};
/**
* Sort objects from back to front
* @returns { number[] }
* @param { Object } obj - object to sort
*/
rglwidgetClass.prototype.depthSort = function(obj) {
var n = obj.centers.length,
depths = new Float32Array(n),
result = new Array(n),
compare = function(i,j) { return depths[j] - depths[i]; },
z, w;
for(i=0; i<n; i++) {
z = this.prmvMatrix.m13*obj.centers[i][0] +
this.prmvMatrix.m23*obj.centers[i][1] +
this.prmvMatrix.m33*obj.centers[i][2] +
this.prmvMatrix.m43;
w = this.prmvMatrix.m14*obj.centers[i][0] +
this.prmvMatrix.m24*obj.centers[i][1] +
this.prmvMatrix.m34*obj.centers[i][2] +
this.prmvMatrix.m44;
depths[i] = z/w;
result[i] = i;
}
result.sort(compare);
return result;
};
rglwidgetClass.prototype.disableArrays = function(obj, enabled) {
var gl = this.gl || this.initGL(),
objLocs = ["normLoc", "texLoc", "ofsLoc", "pointLoc", "nextLoc"],
thisLocs = ["posLoc", "colLoc"], i, attr;
for (i = 0; i < objLocs.length; i++) 
if (enabled[objLocs[i]]) gl.disableVertexAttribArray(obj[objLocs[i]]);
for (i = 0; i < thisLocs.length; i++)
if (enabled[thisLocs[i]]) gl.disableVertexAttribArray(this[objLocs[i]]);
if (typeof obj.userAttributes !== "undefined") {
for (attr in obj.userAttribSizes) {  // Not all attributes may have been used
gl.disableVertexAttribArray( obj.userAttribLocations[attr] );
}
}
}
/**
* Draw an object in a subscene
* @param { number } id - object to draw
* @param { number } subsceneid - id of subscene
*/
rglwidgetClass.prototype.drawObj = function(id, subsceneid) {
var obj = this.getObj(id),
subscene = this.getObj(subsceneid),
flags = obj.flags,
type = obj.type,
is_lit = flags & this.f_is_lit,
has_texture = flags & this.f_has_texture,
fixed_quads = flags & this.f_fixed_quads,
is_transparent = flags & this.f_is_transparent,
depth_sort = flags & this.f_depth_sort,
sprites_3d = flags & this.f_sprites_3d,
sprite_3d = flags & this.f_sprite_3d,
is_lines = flags & this.f_is_lines,
fat_lines = flags & this.f_fat_lines,
is_points = flags & this.f_is_points,
fixed_size = flags & this.f_fixed_size,
is_twosided = (flags & this.f_is_twosided) > 0,
gl = this.gl || this.initGL(),
mat,
sphereMV, baseofs, ofs, sscale, i, count, light,
pass, mode, pmode, attr,
enabled = {};
if (typeof id !== "number") {
this.alertOnce("drawObj id is "+typeof id);
}
if (type === "planes") {
if (obj.bbox !== subscene.par3d.bbox || !obj.initialized) {
this.planeUpdateTriangles(id, subscene.par3d.bbox);
}
}
if (!obj.initialized)
this.initObj(id);
if (type === "clipplanes") {
count = obj.offsets.length;
var IMVClip = [];
for (i=0; i < count; i++) {
IMVClip[i] = this.multMV(this.invMatrix, obj.vClipplane.slice(4*i, 4*(i+1)));
}
obj.IMVClip = IMVClip;
return;
}
if (type === "light" || type === "bboxdeco" || !obj.vertexCount)
return;
if (!is_transparent &&
obj.someHidden) {
is_transparent = true;
depth_sort = ["triangles", "quads", "surface",
"spheres", "sprites", "text"].indexOf(type) >= 0;
}        
this.setDepthTest(id);
if (sprites_3d) {
var norigs = obj.vertices.length,
savenorm = new CanvasMatrix4(this.normMatrix);
this.origs = obj.vertices;
this.usermat = new Float32Array(obj.userMatrix.getAsArray());
this.radii = obj.radii;
this.normMatrix = subscene.spriteNormmat;
for (this.iOrig=0; this.iOrig < norigs; this.iOrig++) {
for (i=0; i < obj.objects.length; i++) {
this.drawObj(obj.objects[i], subsceneid);
}
}
this.normMatrix = savenorm;
return;
} else {
gl.useProgram(obj.prog);
}
if (typeof obj.polygon_offset !== "undefined") {
gl.polygonOffset(obj.polygon_offset[0],
obj.polygon_offset[1]);
gl.enable(gl.POLYGON_OFFSET_FILL);
}
if (sprite_3d) {
gl.uniform3fv(obj.origLoc, new Float32Array(this.origs[this.iOrig]));
if (this.radii.length > 1) {
gl.uniform1f(obj.sizeLoc, this.radii[this.iOrig][0]);
} else {
gl.uniform1f(obj.sizeLoc, this.radii[0][0]);
}
gl.uniformMatrix4fv(obj.usermatLoc, false, this.usermat);
}
if (type === "spheres") {
gl.bindBuffer(gl.ARRAY_BUFFER, this.sphere.buf);
} else {
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
}
gl.uniformMatrix4fv( obj.prMatLoc, false, new Float32Array(this.prMatrix.getAsArray()) );
gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(this.mvMatrix.getAsArray()) );
var clipcheck = 0,
clipplaneids = subscene.clipplanes,
clip, j;
for (i=0; i < clipplaneids.length; i++) {
clip = this.getObj(clipplaneids[i]);
for (j=0; j < clip.offsets.length; j++) {
gl.uniform4fv(obj.clipLoc[clipcheck + j], clip.IMVClip[j]);
}
clipcheck += clip.offsets.length;
}
if (typeof obj.clipLoc !== "undefined")
for (i=clipcheck; i < obj.clipLoc.length; i++)
gl.uniform4f(obj.clipLoc[i], 0,0,0,0);
if (is_lit) {
gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(this.normMatrix.getAsArray()) );
gl.uniform3fv( obj.emissionLoc, obj.emission);
gl.uniform1f( obj.shininessLoc, obj.shininess);
for (i=0; i < subscene.lights.length; i++) {
light = this.getObj(subscene.lights[i]);
if (!light.initialized) this.initObj(subscene.lights[i]);
gl.uniform3fv( obj.ambientLoc[i], this.componentProduct(light.ambient, obj.ambient));
gl.uniform3fv( obj.specularLoc[i], this.componentProduct(light.specular, obj.specular));
gl.uniform3fv( obj.diffuseLoc[i], light.diffuse);
gl.uniform3fv( obj.lightDirLoc[i], light.lightDir);
gl.uniform1i( obj.viewpointLoc[i], light.viewpoint);
gl.uniform1i( obj.finiteLoc[i], light.finite);
}
for (i=subscene.lights.length; i < obj.nlights; i++) {
gl.uniform3f( obj.ambientLoc[i], 0,0,0);
gl.uniform3f( obj.specularLoc[i], 0,0,0);
gl.uniform3f( obj.diffuseLoc[i], 0,0,0);
}
}
if (fixed_size) {
gl.uniform2f( obj.textScaleLoc, 0.75/this.vp.width, 0.75/this.vp.height);
}
gl.enableVertexAttribArray( this.posLoc );
enabled.posLoc = true;
var nc = obj.colorCount;
count = obj.vertexCount;
if (type === "spheres") {
subscene = this.getObj(subsceneid);
var scale = subscene.par3d.scale,
scount = count, indices;
gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
gl.enableVertexAttribArray(obj.normLoc );
enabled.normLoc = true;
gl.vertexAttribPointer(obj.normLoc,  3, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,  0);
gl.disableVertexAttribArray( this.colLoc );
var sphereNorm = new CanvasMatrix4();
sphereNorm.scale(scale[0], scale[1], scale[2]);
sphereNorm.multRight(this.normMatrix);
gl.uniformMatrix4fv( obj.normMatLoc, false, new Float32Array(sphereNorm.getAsArray()) );
if (nc == 1) {
gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
}
if (has_texture) {
gl.enableVertexAttribArray( obj.texLoc );
enabled.texLoc = true;
gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*this.sphere.vOffsets.stride,
4*this.sphere.vOffsets.tofs);
gl.activeTexture(gl.TEXTURE0);
gl.bindTexture(gl.TEXTURE_2D, obj.texture);
gl.uniform1i( obj.sampler, 0);
}
if (depth_sort)
indices = this.depthSort(obj);
for (i = 0; i < scount; i++) {
sphereMV = new CanvasMatrix4();
if (depth_sort) {
baseofs = indices[i]*obj.vOffsets.stride;
} else {
baseofs = i*obj.vOffsets.stride;
}
ofs = baseofs + obj.vOffsets.radofs;
sscale = obj.values[ofs];
sphereMV.scale(sscale/scale[0], sscale/scale[1], sscale/scale[2]);
sphereMV.translate(obj.values[baseofs],
obj.values[baseofs+1],
obj.values[baseofs+2]);
sphereMV.multRight(this.mvMatrix);
gl.uniformMatrix4fv( obj.mvMatLoc, false, new Float32Array(sphereMV.getAsArray()) );
if (nc > 1) {
ofs = baseofs + obj.vOffsets.cofs;
gl.vertexAttrib4f( this.colLoc, obj.values[ofs],
obj.values[ofs+1],
obj.values[ofs+2],
obj.values[ofs+3] );
}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.sphere.ibuf);
gl.drawElements(gl.TRIANGLES, this.sphere.sphereCount, gl.UNSIGNED_SHORT, 0);
}
this.disableArrays(obj, enabled);
if (typeof obj.polygon_offset !== "undefined") 
gl.disable(gl.POLYGON_OFFSET_FILL);
return;
} else {
if (obj.colorCount === 1) {
gl.disableVertexAttribArray( this.colLoc );
gl.vertexAttrib4fv( this.colLoc, new Float32Array(obj.onecolor));
} else {
gl.enableVertexAttribArray( this.colLoc );
enabled.colLoc = true;
gl.vertexAttribPointer(this.colLoc, 4, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.cofs);
}
}
if (is_lit && obj.vOffsets.nofs > 0) {
gl.enableVertexAttribArray( obj.normLoc );
enabled.normLoc = true;
gl.vertexAttribPointer(obj.normLoc, 3, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.nofs);
}
if (has_texture || type === "text") {
gl.enableVertexAttribArray( obj.texLoc );
enabled.texLoc = true;
gl.vertexAttribPointer(obj.texLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.tofs);
gl.activeTexture(gl.TEXTURE0);
gl.bindTexture(gl.TEXTURE_2D, obj.texture);
gl.uniform1i( obj.sampler, 0);
}
if (fixed_quads) {
gl.enableVertexAttribArray( obj.ofsLoc );
enabled.ofsLoc = true;
gl.vertexAttribPointer(obj.ofsLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.oofs);
}
if (typeof obj.userAttributes !== "undefined") {
for (attr in obj.userAttribSizes) {  // Not all attributes may have been used
gl.enableVertexAttribArray( obj.userAttribLocations[attr] );
gl.vertexAttribPointer( obj.userAttribLocations[attr], obj.userAttribSizes[attr],
gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.userAttribOffsets[attr]);
}
}
if (typeof obj.userUniforms !== "undefined") {
for (attr in obj.userUniformLocations) {
var loc = obj.userUniformLocations[attr];
if (loc !== null) {
var uniform = obj.userUniforms[attr];
if (typeof uniform.length === "undefined")
gl.uniform1f(loc, uniform);
else if (typeof uniform[0].length === "undefined") {
uniform = new Float32Array(uniform);
switch(uniform.length) {
case 2: gl.uniform2fv(loc, uniform); break;
case 3: gl.uniform3fv(loc, uniform); break;
case 4: gl.uniform4fv(loc, uniform); break;
default: console.warn("bad uniform length");
}
} else if (uniform.length == 4 && uniform[0].length == 4)
gl.uniformMatrix4fv(loc, false, new Float32Array(uniform.getAsArray()));
else
console.warn("unsupported uniform matrix");
}
}
}
for (pass = 0; pass < obj.passes; pass++) {
pmode = obj.pmode[pass];
if (pmode === "culled")
continue;
mode = fat_lines && (is_lines || pmode == "lines") ? "TRIANGLES" : this.mode4type[type];
if (depth_sort && pmode == "filled") {// Don't try depthsorting on wireframe or points
var faces = this.depthSort(obj),
nfaces = faces.length,
frowsize = Math.floor(obj.f[pass].length/nfaces);
if (type !== "spheres") {
var f = obj.index_uint ? new Uint32Array(obj.f[pass].length) : new Uint16Array(obj.f[pass].length);
for (i=0; i<nfaces; i++) {
for (j=0; j<frowsize; j++) {
f[frowsize*i + j] = obj.f[pass][frowsize*faces[i] + j];
}
}
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, f, gl.DYNAMIC_DRAW);
}
}
if (is_twosided)
gl.uniform1i(obj.frontLoc, pass !== 0);
if (type !== "spheres") 
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, obj.ibuf[pass]);
if (type === "sprites" || type === "text" || type === "quads") {
count = count * 6/4;
} else if (type === "surface") {
count = obj.f[pass].length;
}
count = obj.f[pass].length;
if (!is_lines && pmode === "lines" && !fat_lines) {
mode = "LINES";
} else if (pmode === "points") {
mode = "POINTS";
}
if ((is_lines || pmode === "lines") && fat_lines) {
gl.enableVertexAttribArray(obj.pointLoc);
enabled.pointLoc = true;
gl.vertexAttribPointer(obj.pointLoc, 2, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.pointofs);
gl.enableVertexAttribArray(obj.nextLoc );
enabled.nextLoc = true;
gl.vertexAttribPointer(obj.nextLoc, 3, gl.FLOAT, false, 4*obj.vOffsets.stride, 4*obj.vOffsets.nextofs);
gl.uniform1f(obj.aspectLoc, this.vp.width/this.vp.height);
gl.uniform1f(obj.lwdLoc, this.getMaterial(id, "lwd")/this.vp.height);
}
gl.vertexAttribPointer(this.posLoc,  3, gl.FLOAT, false, 4*obj.vOffsets.stride,  4*obj.vOffsets.vofs);
gl.drawElements(gl[mode], count, obj.index_uint ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT, 0);
this.disableArrays(obj, enabled);
}
if (typeof obj.polygon_offset !== "undefined") 
gl.disable(gl.POLYGON_OFFSET_FILL);
};
/**
* Draw the background for a subscene
* @param { number } id - id of background object
* @param { number } subsceneid - id of subscene
*/
rglwidgetClass.prototype.drawBackground = function(id, subsceneid) {
var gl = this.gl || this.initGL(),
obj = this.getObj(id),
bg, i;
if (!obj.initialized)
this.initObj(id);
if (obj.colors.length) {
bg = obj.colors[0];
gl.clearColor(bg[0], bg[1], bg[2], bg[3]);
gl.depthMask(true);
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
}
if (typeof obj.quad !== "undefined") {
this.prMatrix.makeIdentity();
this.mvMatrix.makeIdentity();
gl.disable(gl.BLEND);
gl.disable(gl.DEPTH_TEST);
gl.depthMask(false);
for (i=0; i < obj.quad.length; i++)
this.drawObj(obj.quad[i], subsceneid);
}
};
/**
* Draw a subscene
* @param { number } subsceneid - id of subscene
* @param { boolean } opaquePass - is this the opaque drawing pass?
*/
rglwidgetClass.prototype.drawSubscene = function(subsceneid, opaquePass) {
var gl = this.gl || this.initGL(),
sub = this.getObj(subsceneid),
objects = this.scene.objects,
subids = sub.objects,
subscene_has_faces = false,
subscene_needs_sorting = false,
flags, i, obj;
if (sub.par3d.skipRedraw)
return;
for (i=0; i < subids.length; i++) {
obj = objects[subids[i]];
flags = obj.flags;
if (typeof flags !== "undefined") {
subscene_has_faces |= (flags & this.f_is_lit)
& !(flags & this.f_fixed_quads);
obj.is_transparent = (flags & this.f_is_transparent) || obj.someHidden;
subscene_needs_sorting |= (flags & this.f_depth_sort) || obj.is_transparent;
}
}
this.setViewport(subsceneid);
if (typeof sub.backgroundId !== "undefined" && opaquePass)
this.drawBackground(sub.backgroundId, subsceneid);
if (subids.length) {
this.setprMatrix(subsceneid);
this.setmvMatrix(subsceneid);
if (subscene_has_faces) {
this.setnormMatrix(subsceneid);
if ((sub.flags & this.f_sprites_3d) &&
typeof sub.spriteNormmat === "undefined") {
sub.spriteNormmat = new CanvasMatrix4(this.normMatrix);
}
}
if (subscene_needs_sorting)
this.setprmvMatrix();
var clipids = sub.clipplanes;
if (typeof clipids === "undefined") {
console.warn("bad clipids");
}
if (clipids.length > 0) {
this.invMatrix = new CanvasMatrix4(this.mvMatrix);
this.invMatrix.invert();
for (i = 0; i < clipids.length; i++)
this.drawObj(clipids[i], subsceneid);
}
subids = sub.opaque.concat(sub.transparent);
if (opaquePass) {
gl.enable(gl.DEPTH_TEST);
gl.depthMask(true);
gl.disable(gl.BLEND);
for (i = 0; i < subids.length; i++) {
if (!this.getObj(subids[i]).is_transparent) 
this.drawObj(subids[i], subsceneid);
}
} else {
gl.depthMask(false);
gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA,
gl.ONE, gl.ONE);
gl.enable(gl.BLEND);
for (i = 0; i < subids.length; i++) {
if (this.getObj(subids[i]).is_transparent)
this.drawObj(subids[i], subsceneid);
}
}
subids = sub.subscenes;
for (i = 0; i < subids.length; i++) {
this.drawSubscene(subids[i], opaquePass);
}
}
};
/**
* Respond to brush change
*/
rglwidgetClass.prototype.selectionChanged = function() {
var i, j, k, id, subid = this.select.subscene, subscene,
objids, obj,
p1 = this.select.region.p1, p2 = this.select.region.p2,
filter, selection = [], handle, keys, xmin, x, xmax, ymin, y, ymax, z, v,
someHidden;
if (!subid)
return;
subscene = this.getObj(subid);
objids = subscene.objects;
filter = this.scene.crosstalk.filter;
this.setmvMatrix(subid);
this.setprMatrix(subid);
this.setprmvMatrix();
xmin = Math.min(p1.x, p2.x);
xmax = Math.max(p1.x, p2.x);
ymin = Math.min(p1.y, p2.y);
ymax = Math.max(p1.y, p2.y);
for (i = 0; i < objids.length; i++) {
id = objids[i];
j = this.scene.crosstalk.id.indexOf(id);
if (j >= 0) {
keys = this.scene.crosstalk.key[j];
obj = this.getObj(id);
someHidden = false;
for (k = 0; k < keys.length; k++) {
if (filter && filter.indexOf(keys[k]) < 0) {
someHidden = true;
continue;
}
v = [].concat(obj.vertices[k]).concat(1.0);
v = this.multVM(v, this.prmvMatrix);
x = v[0]/v[3];
y = v[1]/v[3];
z = v[2]/v[3];
if (xmin <= x && x <= xmax && ymin <= y && y <= ymax && -1.0 <= z && z <= 1.0) {
selection.push(keys[k]);
} else
someHidden = true;
}
obj.someHidden = someHidden && (filter || selection.length);
obj.initialized = false;
/* Who should we notify?  Only shared data in the current subscene, or everyone? */
if (!this.equalArrays(selection, this.scene.crosstalk.selection)) {
handle = this.scene.crosstalk.sel_handle[j];
handle.set(selection, {rglSubsceneId: this.select.subscene});
}
}
}
};
/**
* Respond to selection or filter change from crosstalk
* @param { Object } event - crosstalk event
* @param { boolean } filter - filter or selection?
*/
rglwidgetClass.prototype.selection = function(event, filter) {
var i, j, ids, obj, keys, crosstalk = this.scene.crosstalk,
selection, someHidden;
// Record the message and find out if this event makes some objects have mixed values:
crosstalk = this.scene.crosstalk;
if (filter) {
filter = crosstalk.filter = event.value;
selection = crosstalk.selection;
} else {  
selection = crosstalk.selection = event.value;
filter = crosstalk.filter;
}
ids = crosstalk.id;
for (i = 0; i < ids.length ; i++) {
obj = this.getObj(ids[i]);
obj.initialized = false;
keys = crosstalk.key[i];
someHidden = false;
for (j = 0; j < keys.length && !someHidden; j++) {
if ((filter && filter.indexOf(keys[j]) < 0) ||
(selection.length && selection.indexOf(keys[j]) < 0))
someHidden = true;
}
obj.someHidden = someHidden;
}
this.drawScene();
};
/**
* Clear the selection brush
* @param { number } except - Subscene that should ignore this request
*/
rglwidgetClass.prototype.clearBrush = function(except) {
if (this.select.subscene != except) {
this.select.state = "inactive";
this.delFromSubscene(this.scene.brushId, this.select.subscene);
}
this.drawScene();
};
/**
* Compute mouse coordinates relative to current canvas
* @returns { Object }
* @param { Object } event - event object from mouse click
*/
rglwidgetClass.prototype.relMouseCoords = function(event) {
var totalOffsetX = 0,
totalOffsetY = 0,
currentElement = this.canvas;
do {
totalOffsetX += currentElement.offsetLeft;
totalOffsetY += currentElement.offsetTop;
currentElement = currentElement.offsetParent;
}
while(currentElement);
var canvasX = event.pageX - totalOffsetX,
canvasY = event.pageY - totalOffsetY;
return {x:canvasX, y:canvasY};
};
/**
* Set mouse handlers for the scene
*/
rglwidgetClass.prototype.setMouseHandlers = function() {
var self = this, activeSubscene, handler,
handlers = {}, drag = 0;
handlers.rotBase = 0;
this.screenToVector = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height,
radius = Math.max(width, height)/2.0,
cx = width/2.0,
cy = height/2.0,
px = (x-cx)/radius,
py = (y-cy)/radius,
plen = Math.sqrt(px*px+py*py);
if (plen > 1.e-6) {
px = px/plen;
py = py/plen;
}
var angle = (Math.SQRT2 - plen)/Math.SQRT2*Math.PI/2,
z = Math.sin(angle),
zlen = Math.sqrt(1.0 - z*z);
px = px * zlen;
py = py * zlen;
return [px, py, z];
};
handlers.trackballdown = function(x,y) {
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
handlers.rotBase = this.screenToVector(x, y);
this.saveMat = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
}
};
handlers.trackballmove = function(x,y) {
var rotCurrent = this.screenToVector(x,y),
rotBase = handlers.rotBase,
dot = rotBase[0]*rotCurrent[0] +
rotBase[1]*rotCurrent[1] +
rotBase[2]*rotCurrent[2],
angle = Math.acos( dot/this.vlen(rotBase)/this.vlen(rotCurrent) )*180.0/Math.PI,
axis = this.xprod(rotBase, rotCurrent),
objects = this.scene.objects,
activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
l = activeModel.par3d.listeners,
i;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.userMatrix.load(objects[l[i]].saveMat);
activeSub.par3d.userMatrix.rotate(angle, axis[0], axis[1], axis[2]);
}
this.drawScene();
};
handlers.trackballend = 0;
this.clamp = function(x, lo, hi) {
return Math.max(lo, Math.min(x, hi));
};
this.screenToPolar = function(x,y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height,
r = Math.min(width, height)/2,
dx = this.clamp(x - width/2, -r, r),
dy = this.clamp(y - height/2, -r, r);
return [Math.asin(dx/r), Math.asin(-dy/r)];
};
handlers.polardown = function(x,y) {
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
handlers.dragBase = this.screenToPolar(x, y);
this.saveMat = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
activeSub.camBase = [-Math.atan2(activeSub.saveMat.m13, activeSub.saveMat.m11),
Math.atan2(activeSub.saveMat.m32, activeSub.saveMat.m22)];
}
};
handlers.polarmove = function(x,y) {
var dragCurrent = this.screenToPolar(x,y),
activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
objects = this.scene.objects,
l = activeModel.par3d.listeners,
i, changepos = [];
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
for (j=0; j<2; j++)
changepos[j] = -(dragCurrent[j] - handlers.dragBase[j]);
activeSub.par3d.userMatrix.makeIdentity();
activeSub.par3d.userMatrix.rotate(changepos[0]*180/Math.PI, 0,-1,0);
activeSub.par3d.userMatrix.multRight(objects[l[i]].saveMat);
activeSub.par3d.userMatrix.rotate(changepos[1]*180/Math.PI, -1,0,0);
}
this.drawScene();
};
handlers.polarend = 0;
handlers.axisdown = function(x,y) {
handlers.rotBase = this.screenToVector(x, this.canvas.height/2);
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.saveMat = new CanvasMatrix4(activeSub.par3d.userMatrix);
}
};
handlers.axismove = function(x,y) {
var rotCurrent = this.screenToVector(x, this.canvas.height/2),
rotBase = handlers.rotBase,
angle = (rotCurrent[0] - rotBase[0])*180/Math.PI,
rotMat = new CanvasMatrix4();
rotMat.rotate(angle, handlers.axis[0], handlers.axis[1], handlers.axis[2]);
var activeSub = this.getObj(activeSubscene),
activeModel = this.getObj(this.useid(activeSub.id, "model")),
i, l = activeModel.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.userMatrix.load(activeSub.saveMat);
activeSub.par3d.userMatrix.multLeft(rotMat);
}
this.drawScene();
};
handlers.axisend = 0;
handlers.y0zoom = 0;
handlers.zoom0 = 0;
handlers.zoomdown = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
handlers.y0zoom = y;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.zoom0 = Math.log(activeSub.par3d.zoom);
}
};
handlers.zoommove = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.zoom = Math.exp(activeSub.zoom0 + (y-handlers.y0zoom)/this.canvas.height);
}
this.drawScene();
};
handlers.zoomend = 0;
handlers.y0fov = 0;
handlers.fovdown = function(x, y) {
handlers.y0fov = y;
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.fov0 = activeSub.par3d.FOV;
}
};
handlers.fovmove = function(x, y) {
var activeSub = this.getObj(activeSubscene),
activeProjection = this.getObj(this.useid(activeSub.id, "projection")),
i, l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = this.getObj(l[i]);
activeSub.par3d.FOV = Math.max(1, Math.min(179, activeSub.fov0 +
180*(y-handlers.y0fov)/this.canvas.height));
}
this.drawScene();
};
handlers.fovend = 0;
handlers.selectingdown = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height, 
p = {x: 2.0*x/width - 1.0, y: 2.0*y/height - 1.0};
this.select.region = {p1: p, p2: p};
if (this.select.subscene && this.select.subscene != activeSubscene)
this.delFromSubscene(this.scene.brushId, this.select.subscene);
this.select.subscene = activeSubscene;
this.addToSubscene(this.scene.brushId, activeSubscene);
this.select.state = "changing";
if (typeof this.scene.brushId !== "undefined")
this.getObj(this.scene.brushId).initialized = false;
this.drawScene();
};
handlers.selectingmove = function(x, y) {
var viewport = this.getObj(activeSubscene).par3d.viewport,
width = viewport.width*this.canvas.width,
height = viewport.height*this.canvas.height;
if (this.select.state === "inactive") 
return;
this.select.region.p2 = {x: 2.0*x/width - 1.0, y: 2.0*y/height - 1.0};
if (typeof this.scene.brushId !== "undefined")
this.getObj(this.scene.brushId).initialized = false;
this.drawScene();
};
handlers.selectingend = 0;
this.canvas.onmousedown = function ( ev ){
if (!ev.which) // Use w3c defns in preference to MS
switch (ev.button) {
case 0: ev.which = 1; break;
case 1:
case 4: ev.which = 2; break;
case 2: ev.which = 3;
}
drag = ["left", "middle", "right"][ev.which-1];
var coords = self.relMouseCoords(ev);
coords.y = self.canvas.height-coords.y;
activeSubscene = self.whichSubscene(coords);
var sub = self.getObj(activeSubscene), f;
handler = sub.par3d.mouseMode[drag];
switch (handler) {
case "xAxis":
handler = "axis";
handlers.axis = [1.0, 0.0, 0.0];
break;
case "yAxis":
handler = "axis";
handlers.axis = [0.0, 1.0, 0.0];
break;
case "zAxis":
handler = "axis";
handlers.axis = [0.0, 0.0, 1.0];
break;
}
f = handlers[handler + "down"];
if (f) {
coords = self.translateCoords(activeSubscene, coords);
f.call(self, coords.x, coords.y);
ev.preventDefault();
} else
console.warn("Mouse handler '" + handler + "' is not implemented.");
};
this.canvas.onmouseup = function ( ev ){
if ( drag === 0 ) return;
var f = handlers[handler + "end"];
if (f) {
f.call(self);
ev.preventDefault();
}
drag = 0;
};
this.canvas.onmouseout = this.canvas.onmouseup;
this.canvas.onmousemove = function ( ev ) {
if ( drag === 0 ) return;
var f = handlers[handler + "move"];
if (f) {
var coords = self.relMouseCoords(ev);
coords.y = self.canvas.height - coords.y;
coords = self.translateCoords(activeSubscene, coords);
f.call(self, coords.x, coords.y);
}
};
handlers.wheelHandler = function(ev) {
var del = 1.02, i;
if (ev.shiftKey) del = 1.002;
var ds = ((ev.detail || ev.wheelDelta) > 0) ? del : (1 / del);
if (typeof activeSubscene === "undefined")
activeSubscene = self.scene.rootSubscene;
var activeSub = self.getObj(activeSubscene),
activeProjection = self.getObj(self.useid(activeSub.id, "projection")),
l = activeProjection.par3d.listeners;
for (i = 0; i < l.length; i++) {
activeSub = self.getObj(l[i]);
activeSub.par3d.zoom *= ds;
}
self.drawScene();
ev.preventDefault();
};
this.canvas.addEventListener("DOMMouseScroll", handlers.wheelHandler, false);
this.canvas.addEventListener("mousewheel", handlers.wheelHandler, false);
};
/**
* Find a particular subscene by inheritance
* @returns { number } id of subscene to use
* @param { number } subsceneid - child subscene
* @param { string } type - type of inheritance:  "projection" or "model"
*/
rglwidgetClass.prototype.useid = function(subsceneid, type) {
var sub = this.getObj(subsceneid);
if (sub.embeddings[type] === "inherit")
return(this.useid(sub.parent, type));
else
return subsceneid;
};
/**
* Check whether point is in viewport of subscene
* @returns {boolean}
* @param { Object } coords - screen coordinates of point
* @param { number } subsceneid - subscene to check
*/
rglwidgetClass.prototype.inViewport = function(coords, subsceneid) {
var viewport = this.getObj(subsceneid).par3d.viewport,
x0 = coords.x - viewport.x*this.canvas.width,
y0 = coords.y - viewport.y*this.canvas.height;
return 0 <= x0 && x0 <= viewport.width*this.canvas.width &&
0 <= y0 && y0 <= viewport.height*this.canvas.height;
};
/**
* Find which subscene contains a point
* @returns { number } subscene id
* @param { Object } coords - coordinates of point
*/
rglwidgetClass.prototype.whichSubscene = function(coords) {
var self = this,
recurse = function(subsceneid) {
var subscenes = self.getChildSubscenes(subsceneid), i, id;
for (i=0; i < subscenes.length; i++) {
id = recurse(subscenes[i]);
if (typeof(id) !== "undefined")
return(id);
}
if (self.inViewport(coords, subsceneid))
return(subsceneid);
else
return undefined;
},
rootid = this.scene.rootSubscene,
result = recurse(rootid);
if (typeof(result) === "undefined")
result = rootid;
return result;
};
/**
* Translate from window coordinates to viewport coordinates
* @returns { Object } translated coordinates
* @param { number } subsceneid - which subscene to use?
* @param { Object } coords - point to translate
*/
rglwidgetClass.prototype.translateCoords = function(subsceneid, coords) {
var viewport = this.getObj(subsceneid).par3d.viewport;
return {x: coords.x - viewport.x*this.canvas.width,
y: coords.y - viewport.y*this.canvas.height};
};
/**
* Initialize the sphere object
*/
rglwidgetClass.prototype.initSphere = function() {
var verts = this.scene.sphereVerts,
reuse = verts.reuse, result;
if (typeof reuse !== "undefined") {
var prev = document.getElementById(reuse).rglinstance.sphere;
result = {values: prev.values, vOffsets: prev.vOffsets, it: prev.it};
} else
result = {values: new Float32Array(this.flatten(this.cbind(this.transpose(verts.vb),
this.transpose(verts.texcoords)))),
it: new Uint16Array(this.flatten(this.transpose(verts.it))),
vOffsets: {vofs:0, cofs:-1, nofs:-1, radofs:-1, oofs:-1,
tofs:3, nextofs:-1, pointofs:-1, stride:5}};
result.sphereCount = result.it.length;
this.sphere = result;
};
/**
* Set the vertices in the selection box object
*/
rglwidgetClass.prototype.initSelection = function(id) {
if (typeof this.select.region === "undefined")
return;
var obj = this.getObj(id),
width = this.canvas.width,
height = this.canvas.height, 
p1 = this.select.region.p1,
p2 = this.select.region.p2;
obj.vertices = [[p1.x, p1.y, 0.0],
[p2.x, p1.y, 0.0],
[p2.x, p2.y, 0.0],
[p1.x, p2.y, 0.0],
[p1.x, p1.y, 0.0]];
};
/**
* Do the gl part of initializing the sphere
*/
rglwidgetClass.prototype.initSphereGL = function() {
var gl = this.gl || this.initGL(), sphere = this.sphere;
if (gl.isContextLost()) return;
sphere.buf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, sphere.buf);
gl.bufferData(gl.ARRAY_BUFFER, sphere.values, gl.STATIC_DRAW);
sphere.ibuf = gl.createBuffer();
gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, sphere.ibuf);
gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, sphere.it, gl.STATIC_DRAW);
return;
};
/**
* Initialize the DOM object
* @param { Object } el - the DOM object
* @param { Object } x - the scene data sent by JSON from R
*/
rglwidgetClass.prototype.initialize = function(el, x) {
this.textureCanvas = document.createElement("canvas");
this.textureCanvas.style.display = "block";
this.scene = x;
this.normMatrix = new CanvasMatrix4();
this.saveMat = {};
this.distance = null;
this.posLoc = 0;
this.colLoc = 1;
if (el) {
el.rglinstance = this;
this.el = el;
this.webGLoptions = el.rglinstance.scene.webGLoptions;
this.initCanvas();
}
if (typeof Shiny !== "undefined") {
var self = this;
Shiny.addCustomMessageHandler("shinyGetPar3d",
function(message) {
var i, param, 
subscene = self.getObj(message.subscene),
parameters = [].concat(message.parameters),
result = {tag: message.tag, subscene: message.subscene};
if (typeof subscene !== "undefined") {
for (i = 0; i < parameters.length; i++) {
param = parameters[i];
result[param] = subscene.par3d[param];
};
} else {
console.log("subscene "+message.subscene+" undefined.")
}
Shiny.setInputValue("par3d:shinyPar3d", result, {priority: "event"});
});
Shiny.addCustomMessageHandler("shinySetPar3d",
function(message) {
var param = message.parameter, 
subscene = self.getObj(message.subscene);
if (typeof subscene !== "undefined") {
subscene.par3d[param] = message.value;
subscene.initialized = false;
self.drawScene();
} else {
console.log("subscene "+message.subscene+" undefined.")
}
})
}
};
/**
* Restart the WebGL canvas
*/
rglwidgetClass.prototype.restartCanvas = function() {
var newcanvas = document.createElement("canvas"),
self = this;
newcanvas.width = this.el.width;
newcanvas.height = this.el.height;
newcanvas.addEventListener("webglcontextrestored",
this.onContextRestored, false);
newcanvas.addEventListener("webglcontextlost",
this.onContextLost, false);
while (this.el.firstChild) {
this.el.removeChild(this.el.firstChild);
}
this.el.appendChild(newcanvas);
this.canvas = newcanvas;
this.setMouseHandlers();
if (this.gl) 
Object.keys(this.scene.objects).forEach(function(key){
self.getObj(parseInt(key, 10)).texture = undefined; 
});
this.gl = null;
};
/**
* Initialize the WebGL canvas
*/
rglwidgetClass.prototype.initCanvas = function() {
this.restartCanvas();
var objs = this.scene.objects,
self = this;
Object.keys(objs).forEach(function(key){
var id = parseInt(key, 10),
obj = self.getObj(id);
if (typeof obj.reuse !== "undefined")
self.copyObj(id, obj.reuse);
});
Object.keys(objs).forEach(function(key){
self.initSubscene(parseInt(key, 10));
});
this.setMouseHandlers();
this.initSphere();
this.onContextRestored = function(event) {
self.initGL();
self.drawScene();
};
this.onContextLost = function(event) {
if (!self.drawing)
this.gl = null;
event.preventDefault();
};
this.initGL0();
this.lazyLoadScene = function() {
if (typeof self.slide === "undefined")
self.slide = self.getSlide();
if (self.isInBrowserViewport()) {
if (!self.gl || self.gl.isContextLost())
self.initGL();
self.drawScene();
}
};
window.addEventListener("DOMContentLoaded", this.lazyLoadScene, false);
window.addEventListener("load", this.lazyLoadScene, false);
window.addEventListener("resize", this.lazyLoadScene, false);
window.addEventListener("scroll", this.lazyLoadScene, false);
this.slide = this.getSlide();
if (this.slide) {
if (typeof this.slide.rgl === "undefined")
this.slide.rgl = [this];
else
this.slide.rgl.push(this);
if (this.scene.context.rmarkdown) 
if (this.scene.context.rmarkdown === "ioslides_presentation") {
this.slide.setAttribute("slideenter", "this.rgl.forEach(function(scene) { scene.lazyLoadScene.call(window);})");
} else if (this.scene.context.rmarkdown === "slidy_presentation") {
// This method would also work in ioslides, but it gets triggered
// something like 5 times per slide for every slide change, so
// you'd need a quicker function than lazyLoadScene.
var MutationObserver = window.MutationObserver || window.WebKitMutationObserver || window.MozMutationObserver,
observer = new MutationObserver(function(mutations) {
mutations.forEach(function(mutation) {
self.slide.rgl.forEach(function(scene) { scene.lazyLoadScene.call(window); });});});
observer.observe(this.slide, { attributes: true, attributeFilter:["class"] });
}
}
};
/**
* Start the writeWebGL scene. This is only used by writeWebGL; rglwidget has
no debug element and does the drawing in rglwidget.js.
*/
rglwidgetClass.prototype.start = function() {
if (typeof this.prefix !== "undefined") {
this.debugelement = document.getElementById(this.prefix + "debug");
this.debug("");
}
this.drag = 0;
this.drawScene();
};
/**
* Display a debug message
* @param { string } msg - The message to display
* @param { Object } [img] - Image to insert before message
*/
rglwidgetClass.prototype.debug = function(msg, img) {
if (typeof this.debugelement !== "undefined" && this.debugelement !== null) {
this.debugelement.innerHTML = msg;
if (typeof img !== "undefined") {
this.debugelement.insertBefore(img, this.debugelement.firstChild);
}
} else if (msg !== "")
alert(msg);
};
/**
* Get the snapshot image of this scene
* @returns { Object } The img DOM element
*/
rglwidgetClass.prototype.getSnapshot = function() {
var img;
if (typeof this.scene.snapshot !== "undefined") {
img = document.createElement("img");
img.src = this.scene.snapshot;
img.alt = "Snapshot";
}
return img;
};
/**
* Initial test for WebGL
*/
rglwidgetClass.prototype.initGL0 = function() {
if (!window.WebGLRenderingContext){
alert("Your browser does not support WebGL. See http://get.webgl.org");
return;
}
};
/**
* If we are in an ioslides or slidy presentation, get the
* DOM element of the current slide
* @returns { Object }
*/
rglwidgetClass.prototype.getSlide = function() {
var result = this.el, done = false;
while (result && !done && this.scene.context.rmarkdown) {
switch(this.scene.context.rmarkdown) {
case "ioslides_presentation":
if (result.tagName === "SLIDE") return result;
break;
case "slidy_presentation":
if (result.tagName === "DIV" && result.classList.contains("slide"))
return result;
break;
default: return null;
}
result = result.parentElement;
}
return null;
};
/**
* Is this scene visible in the browser?
* @returns { boolean }
*/
rglwidgetClass.prototype.isInBrowserViewport = function() {
var rect = this.canvas.getBoundingClientRect(),
windHeight = (window.innerHeight || document.documentElement.clientHeight),
windWidth = (window.innerWidth || document.documentElement.clientWidth);
if (this.scene.context && this.scene.context.rmarkdown !== null) {
if (this.slide)
return (this.scene.context.rmarkdown === "ioslides_presentation" &&
this.slide.classList.contains("current")) ||
(this.scene.context.rmarkdown === "slidy_presentation" &&
!this.slide.classList.contains("hidden"));
}
return (
rect.top >= -windHeight &&
rect.left >= -windWidth &&
rect.bottom <= 2*windHeight &&
rect.right <= 2*windWidth);
};
/**
* Initialize WebGL
* @returns { Object } the WebGL context
*/
rglwidgetClass.prototype.initGL = function() {
var self = this;
if (this.gl) {
if (!this.drawing && this.gl.isContextLost())
this.restartCanvas();
else
return this.gl;
}
// if (!this.isInBrowserViewport()) return; Return what??? At this point we know this.gl is null.
this.canvas.addEventListener("webglcontextrestored",
this.onContextRestored, false);
this.canvas.addEventListener("webglcontextlost",
this.onContextLost, false);
this.gl = this.canvas.getContext("webgl", this.webGLoptions) ||
this.canvas.getContext("experimental-webgl", this.webGLoptions);
this.index_uint = this.gl.getExtension("OES_element_index_uint");
var save = this.startDrawing();
this.initSphereGL();
Object.keys(this.scene.objects).forEach(function(key){
self.initObj(parseInt(key, 10));
});
this.stopDrawing(save);
return this.gl;
};
/**
* Resize the display to match element
* @param { Object } el - DOM element to match
*/
rglwidgetClass.prototype.resize = function(el) {
this.canvas.width = el.width;
this.canvas.height = el.height;
};
/**
* Draw the whole scene
*/
rglwidgetClass.prototype.drawScene = function() {
var gl = this.gl || this.initGL(),
wasDrawing = this.startDrawing();
if (!wasDrawing) {
if (this.select.state !== "inactive")
this.selectionChanged();
gl.enable(gl.DEPTH_TEST);
gl.depthFunc(gl.LEQUAL);
gl.clearDepth(1.0);
gl.clearColor(1,1,1,1);
gl.depthMask(true); // Must be true before clearing depth buffer
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
this.drawSubscene(this.scene.rootSubscene, true);
this.drawSubscene(this.scene.rootSubscene, false);
}
this.stopDrawing(wasDrawing);
};
/**
* Change the displayed subset
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The subset control data.
*/
rglwidgetClass.prototype.subsetSetter = function(el, control) {
if (typeof control.subscenes === "undefined" ||
control.subscenes === null)
control.subscenes = this.scene.rootSubscene;
var value = Math.round(control.value),
subscenes = [].concat(control.subscenes),
fullset = [].concat(control.fullset),
i, j, entries, subsceneid,
adds = [], deletes = [],
ismissing = function(x) {
return fullset.indexOf(x) < 0;
},
tointeger = function(x) {
return parseInt(x, 10);
};
if (isNaN(value))
value = control.value = 0;
if (control.accumulate)
for (i=0; i <= value; i++)
adds = adds.concat(control.subsets[i]);
else
adds = adds.concat(control.subsets[value]);
deletes = fullset.filter(function(x) { return adds.indexOf(x) < 0; });
for (i = 0; i < subscenes.length; i++) {
subsceneid = subscenes[i];
if (typeof this.getObj(subsceneid) === "undefined")
this.alertOnce("typeof object is undefined");
for (j = 0; j < adds.length; j++)
this.addToSubscene(adds[j], subsceneid);
for (j = 0; j < deletes.length; j++)
this.delFromSubscene(deletes[j], subsceneid);
}
};
/**
* Change the requested property
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The property setter control data.
*/
rglwidgetClass.prototype.propertySetter = function(el, control)  {
var value = control.value,
values = [].concat(control.values),
svals = [].concat(control.param),
direct = values[0] === null,
entries = [].concat(control.entries),
ncol = entries.length,
nrow = values.length/ncol,
properties = this.repeatToLen(control.properties, ncol),
objids = this.repeatToLen(control.objids, ncol),
property, objid = objids[0],
obj = this.getObj(objid),
propvals, i, v1, v2, p, entry, gl, needsBinding,
newprop, newid,
getPropvals = function() {
if (property === "userMatrix")
return obj.par3d.userMatrix.getAsArray();
else if (property === "scale" || property === "FOV" || property === "zoom")
return [].concat(obj.par3d[property]);
else
return [].concat(obj[property]);
};
putPropvals = function(newvals) {
if (newvals.length == 1)
newvals = newvals[0];
if (property === "userMatrix")
obj.par3d.userMatrix.load(newvals);
else if (property === "scale" || property === "FOV" || property === "zoom")
obj.par3d[property] = newvals;
else
obj[property] = newvals;
};
if (direct && typeof value === "undefined")
return;
if (control.interp) {
values = values.slice(0, ncol).concat(values).
concat(values.slice(ncol*(nrow-1), ncol*nrow));
svals = [-Infinity].concat(svals).concat(Infinity);
for (i = 1; i < svals.length; i++) {
if (value <= svals[i]) {
if (svals[i] === Infinity)
p = 1;
else
p = (svals[i] - value)/(svals[i] - svals[i-1]);
break;
}
}
} else if (!direct) {
value = Math.round(value);
}
for (j=0; j<entries.length; j++) {
entry = entries[j];
newprop = properties[j];
newid = objids[j];
if (newprop !== property || newid != objid) {
if (typeof property !== "undefined")
putPropvals(propvals);
property = newprop;
objid = newid;
obj = this.getObj(objid);
propvals = getPropvals();
}
if (control.interp) {
v1 = values[ncol*(i-1) + j];
v2 = values[ncol*i + j];
this.setElement(propvals, entry, p*v1 + (1-p)*v2);
} else if (!direct) {
this.setElement(propvals, entry, values[ncol*value + j]);
} else {
this.setElement(propvals, entry, value[j]);
}
}
putPropvals(propvals);
needsBinding = [];
for (j=0; j < entries.length; j++) {
if (properties[j] === "values" &&
needsBinding.indexOf(objids[j]) === -1) {
needsBinding.push(objids[j]);
}
}
for (j=0; j < needsBinding.length; j++) {
gl = this.gl || this.initGL();
obj = this.getObj(needsBinding[j]);
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
}
};
/**
* Change the requested vertices
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The vertext setter control data.
*/
rglwidgetClass.prototype.vertexSetter = function(el, control)  {
var svals = [].concat(control.param),
j, k, p, a, propvals, stride, ofs, obj, entry,
attrib,
ofss    = {x:"vofs", y:"vofs", z:"vofs",
red:"cofs", green:"cofs", blue:"cofs",
alpha:"cofs", radii:"radofs",
nx:"nofs", ny:"nofs", nz:"nofs",
ox:"oofs", oy:"oofs", oz:"oofs",
ts:"tofs", tt:"tofs"},
pos     = {x:0, y:1, z:2,
red:0, green:1, blue:2,
alpha:3,radii:0,
nx:0, ny:1, nz:2,
ox:0, oy:1, oz:2,
ts:0, tt:1},
values = control.values,
direct = values === null,
ncol,
interp = control.interp,
vertices = [].concat(control.vertices),
attributes = [].concat(control.attributes),
value = control.value, newval, aliases, alias;
ncol = Math.max(vertices.length, attributes.length);
if (!ncol)
return;
vertices = this.repeatToLen(vertices, ncol);
attributes = this.repeatToLen(attributes, ncol);
if (direct)
interp = false;
/* JSON doesn't pass Infinity */
svals[0] = -Infinity;
svals[svals.length - 1] = Infinity;
for (j = 1; j < svals.length; j++) {
if (value <= svals[j]) {
if (interp) {
if (svals[j] === Infinity)
p = 1;
else
p = (svals[j] - value)/(svals[j] - svals[j-1]);
} else {
if (svals[j] - value > value - svals[j-1])
j = j - 1;
}
break;
}
}
obj = this.getObj(control.objid);
// First, make sure color attributes vary in original
if (typeof obj.vOffsets !== "undefined") {
varies = true;
for (k = 0; k < ncol; k++) {
attrib = attributes[k];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[attrib]];
if (ofs < 0) {
switch(attrib) {
case "alpha":
case "red":
case "green":
case "blue":
obj.colors = [obj.colors[0], obj.colors[0]];
break;
}
varies = false;
}
}
}
if (!varies)
this.initObj(control.objid);
}
propvals = obj.values;
aliases = obj.alias;
if (typeof aliases === "undefined")
aliases = [];
for (k=0; k<ncol; k++) {
if (direct) {
newval = value;
} else if (interp) {
newval = p*values[j-1][k] + (1-p)*values[j][k];
} else {
newval = values[j][k];
}       
attrib = attributes[k];
vertex = vertices[k];
alias = aliases[vertex];
if (obj.type === "planes" || obj.type === "clipplanes") {
ofs = ["nx", "ny", "nz", "offset"].indexOf(attrib);
if (ofs >= 0) {
if (ofs < 3) {
if (obj.normals[vertex][ofs] != newval) {  // Assume no aliases here...
obj.normals[vertex][ofs] = newval;
obj.initialized = false;
}
} else {
if (obj.offsets[vertex][0] != newval) {
obj.offsets[vertex][0] = newval;
obj.initialized = false;
}
}
continue;
}
}
// Not a plane setting...
ofs = obj.vOffsets[ofss[attrib]];
if (ofs < 0)
this.alertOnce("Attribute '"+attrib+"' not found in object "+control.objid);
else {
stride = obj.vOffsets.stride;
ofs = ofs + pos[attrib];
entry = vertex*stride + ofs;
propvals[entry] = newval;
if (typeof alias !== "undefined")
for (a = 0; a < alias.length; a++)
propvals[alias[a]*stride + ofs] = newval;
}
}
if (typeof obj.buf !== "undefined") {
var gl = this.gl || this.initGL();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, propvals, gl.STATIC_DRAW);
}
};
/**
* Change the requested vertex properties by age
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The age setter control data.
*/
rglwidgetClass.prototype.ageSetter = function(el, control) {
var objids = [].concat(control.objids),
nobjs = objids.length,
time = control.value,
births = [].concat(control.births),
ages = [].concat(control.ages),
steps = births.length,
j = Array(steps),
p = Array(steps),
i, k, age, j0, propvals, stride, ofs, objid, obj,
attrib, dim, varies, alias, aliases, a, d,
attribs = ["colors", "alpha", "radii", "vertices",
"normals", "origins", "texcoords",
"x", "y", "z",
"red", "green", "blue"],
ofss    = ["cofs", "cofs", "radofs", "vofs",
"nofs", "oofs", "tofs",
"vofs", "vofs", "vofs",
"cofs", "cofs", "cofs"],
dims    = [3,1,1,3,
3,2,2,
1,1,1,
1,1,1],
pos     = [0,3,0,0,
0,0,0,
0,1,2,
0,1,2];
/* Infinity doesn't make it through JSON */
ages[0] = -Infinity;
ages[ages.length-1] = Infinity;
for (i = 0; i < steps; i++) {
if (births[i] !== null) {  // NA in R becomes null
age = time - births[i];
for (j0 = 1; age > ages[j0]; j0++);
if (ages[j0] == Infinity)
p[i] = 1;
else if (ages[j0] > ages[j0-1])
p[i] = (ages[j0] - age)/(ages[j0] - ages[j0-1]);
else
p[i] = 0;
j[i] = j0;
}
}
// First, make sure color attributes vary in original
for (l = 0; l < nobjs; l++) {
objid = objids[l];
obj = this.getObj(objid);
varies = true;
if (typeof obj.vOffsets === "undefined")
continue;
for (k = 0; k < attribs.length; k++) {
attrib = control[attribs[k]];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[k]];
if (ofs < 0) {
switch(attribs[k]) {
case "colors":
case "alpha":
case "red":
case "green":
case "blue":
obj.colors = [obj.colors[0], obj.colors[0]];
break;
}
varies = false;
}
}
}
if (!varies)
this.initObj(objid);
}
for (l = 0; l < nobjs; l++) {
objid = objids[l];
obj = this.getObj(objid);
if (typeof obj.vOffsets === "undefined")
continue;
aliases = obj.alias;
if (typeof aliases === "undefined")
aliases = [];
propvals = obj.values;
stride = obj.vOffsets.stride;
for (k = 0; k < attribs.length; k++) {
attrib = control[attribs[k]];
if (typeof attrib !== "undefined") {
ofs = obj.vOffsets[ofss[k]];
if (ofs >= 0) {
dim = dims[k];
ofs = ofs + pos[k];
for (i = 0; i < steps; i++) {
alias = aliases[i];
if (births[i] !== null) {
for (d=0; d < dim; d++) {
propvals[i*stride + ofs + d] = p[i]*attrib[dim*(j[i]-1) + d] + (1-p[i])*attrib[dim*j[i] + d];
if (typeof alias !== "undefined")
for (a=0; a < alias.length; a++)
propvals[alias[a]*stride + ofs + d] = propvals[i*stride + ofs + d];
}
}
}
} else
this.alertOnce("\'"+attribs[k]+"\' property not found in object "+objid);
}
}
obj.values = propvals;
if (typeof obj.buf !== "undefined") {
gl = this.gl || this.initGL();
gl.bindBuffer(gl.ARRAY_BUFFER, obj.buf);
gl.bufferData(gl.ARRAY_BUFFER, obj.values, gl.STATIC_DRAW);
}
}
};
/**
* Bridge to old style control
* @param { Object } el - Element of the control; not used.
* @param { Object } control - The bridge control data.
*/
rglwidgetClass.prototype.oldBridge = function(el, control) {
var attrname, global = window[control.prefix + "rgl"];
if (global)
for (attrname in global)
this[attrname] = global[attrname];
window[control.prefix + "rgl"] = this;
};
/**
* Set up a player control
* @param { Object } el - The player control element
* @param { Object } control - The player data.
*/
rglwidgetClass.prototype.Player = function(el, control) {
var
self = this,
components = [].concat(control.components),
buttonLabels = [].concat(control.buttonLabels),
Tick = function() { /* "this" will be a timer */
var i,
nominal = this.value,
slider = this.Slider,
labels = this.outputLabels,
output = this.Output,
step;
if (typeof slider !== "undefined" && nominal != slider.value)
slider.value = nominal;
if (typeof output !== "undefined") {
step = Math.round((nominal - output.sliderMin)/output.sliderStep);
if (labels !== null) {
output.innerHTML = labels[step];
} else {
step = step*output.sliderStep + output.sliderMin;
output.innerHTML = step.toPrecision(output.outputPrecision);
}
}
for (i=0; i < this.actions.length; i++) {
this.actions[i].value = nominal;
}
self.applyControls(el, this.actions, false);
self.drawScene();
},
OnSliderInput = function() { /* "this" will be the slider */
this.rgltimer.value = Number(this.value);
this.rgltimer.Tick();
},
addSlider = function(min, max, step, value) {
var slider = document.createElement("input");
slider.type = "range";
slider.min = min;
slider.max = max;
slider.step = step;
slider.value = value;
slider.oninput = OnSliderInput;
slider.sliderActions = control.actions;
slider.sliderScene = this;
slider.className = "rgl-slider";
slider.id = el.id + "-slider";
el.rgltimer.Slider = slider;
slider.rgltimer = el.rgltimer;
el.appendChild(slider);
},
addLabel = function(labels, min, step, precision) {
var output = document.createElement("output");
output.sliderMin = min;
output.sliderStep = step;
output.outputPrecision = precision;
output.className = "rgl-label";
output.id = el.id + "-label";
el.rgltimer.Output = output;
el.rgltimer.outputLabels = labels;
el.appendChild(output);
},
addButton = function(which, label, active) {
var button = document.createElement("input"),
onclicks = {Reverse: function() { this.rgltimer.reverse();},
Play: function() { this.rgltimer.play();
this.value = this.rgltimer.enabled ? this.inactiveValue : this.activeValue; },
Slower: function() { this.rgltimer.slower(); },
Faster: function() { this.rgltimer.faster(); },
Reset: function() { this.rgltimer.reset(); },
Step:  function() { this.rgltimer.step(); }
};
button.rgltimer = el.rgltimer;
button.type = "button";
button.value = label;
button.activeValue = label;
button.inactiveValue = active;
if (which === "Play")
button.rgltimer.PlayButton = button;
button.onclick = onclicks[which];
button.className = "rgl-button";
button.id = el.id + "-" + which;
el.appendChild(button);
};
if (typeof control.reinit !== "undefined" && control.reinit !== null) {
control.actions.reinit = control.reinit;
}
el.rgltimer = new rgltimerClass(Tick, control.start, control.interval, control.stop,
control.step, control.value, control.rate, control.loop, control.actions);
for (var i=0; i < components.length; i++) {
switch(components[i]) {
case "Slider": addSlider(control.start, control.stop,
control.step, control.value);
break;
case "Label": addLabel(control.labels, control.start,
control.step, control.precision);
break;
default:
addButton(components[i], buttonLabels[i], control.pause);
}
}
el.rgltimer.Tick();
};
/**
* Apply all registered controls
* @param { Object } el - DOM element of the control
* @param { Object } x - List of actions to apply
* @param { boolean } [draw=true] - Whether to redraw after applying
*/
rglwidgetClass.prototype.applyControls = function(el, x, draw) {
var self = this, reinit = x.reinit, i, control, type;
for (i = 0; i < x.length; i++) {
control = x[i];
type = control.type;
self[type](el, control);
}
if (typeof reinit !== "undefined" && reinit !== null) {
reinit = [].concat(reinit);
for (i = 0; i < reinit.length; i++)
self.getObj(reinit[i]).initialized = false;
}
if (typeof draw === "undefined" || draw)
self.drawScene();
};
/**
* Handler for scene change
* @param { Object } message - What sort of scene change to do?
*/
rglwidgetClass.prototype.sceneChangeHandler = function(message) {
var self = document.getElementById(message.elementId).rglinstance,
objs = message.objects, mat = message.material,
root = message.rootSubscene,
initSubs = message.initSubscenes,
redraw = message.redrawScene,
skipRedraw = message.skipRedraw,
deletes, subs, allsubs = [], i,j;
if (typeof message.delete !== "undefined") {
deletes = [].concat(message.delete);
if (typeof message.delfromSubscenes !== "undefined")
subs = [].concat(message.delfromSubscenes);
else
subs = [];
for (i = 0; i < deletes.length; i++) {
for (j = 0; j < subs.length; j++) {
self.delFromSubscene(deletes[i], subs[j]);
}
delete self.scene.objects[deletes[i]];
}
}
if (typeof objs !== "undefined") {
Object.keys(objs).forEach(function(key){
key = parseInt(key, 10);
self.scene.objects[key] = objs[key];
self.initObj(key);
var obj = self.getObj(key),
subs = [].concat(obj.inSubscenes), k;
allsubs = allsubs.concat(subs);
for (k = 0; k < subs.length; k++)
self.addToSubscene(key, subs[k]);
});
}
if (typeof mat !== "undefined") {
self.scene.material = mat;
}
if (typeof root !== "undefined") {
self.scene.rootSubscene = root;
}
if (typeof initSubs !== "undefined")
allsubs = allsubs.concat(initSubs);
allsubs = self.unique(allsubs);
for (i = 0; i < allsubs.length; i++) {
self.initSubscene(allsubs[i]);
}
if (typeof skipRedraw !== "undefined") {
root = self.getObj(self.scene.rootSubscene);
root.par3d.skipRedraw = skipRedraw;
}
if (redraw)
self.drawScene();
};
/**
* Set mouse mode for a subscene
* @param { string } mode - name of mode
* @param { number } button - button number (1 to 3)
* @param { number } subscene - subscene id number
* @param { number } stayActive - if truthy, don't clear brush
*/
rglwidgetClass.prototype.setMouseMode = function(mode, button, subscene, stayActive) {
var sub = this.getObj(subscene),
which = ["left", "right", "middle"][button - 1];
if (!stayActive && sub.par3d.mouseMode[which] === "selecting")
this.clearBrush(null);
sub.par3d.mouseMode[which] = mode;
};
/**
* The class of an rgl timer object
* @class
*/
/**
* Construct an rgltimerClass object
* @constructor
* @param { function } Tick - action when timer fires
* @param { number } startTime - nominal start time in seconds
* @param { number } interval - seconds between updates
* @param { number } stopTime - nominal stop time in seconds
* @param { number } stepSize - nominal step size
* @param { number } value - current nominal time
* @param { number } rate - nominal units per second
* @param { string } loop - "none", "cycle" or "oscillate"
* @param { Object } actions - list of actions
*/
rgltimerClass = function(Tick, startTime, interval, stopTime, stepSize, value, rate, loop, actions) {
this.enabled = false;
this.timerId = 0;
/** nominal start time in seconds */
this.startTime = startTime;   
/** current nominal time */      
this.value = value;
/** seconds between updates */                 
this.interval = interval;
/** nominal stop time */           
this.stopTime = stopTime;
/** nominal step size */           
this.stepSize = stepSize;
/** nominal units per second */           
this.rate = rate;
/** "none", "cycle", or "oscillate" */                   
this.loop = loop;
/** real world start time */                   
this.realStart = undefined;
/** multiplier for fast-forward or reverse */         
this.multiplier = 1;                
this.actions = actions;
this.Tick = Tick;
};
/**
* Start playing timer object
*/
rgltimerClass.prototype.play = function() {
if (this.enabled) {
this.enabled = false;
window.clearInterval(this.timerId);
this.timerId = 0;
return;
}
var tick = function(self) {
var now = new Date();
self.value = self.multiplier*self.rate*(now - self.realStart)/1000 + self.startTime;
self.forceToRange();
if (typeof self.Tick !== "undefined") {
self.Tick(self.value);
}
};
this.realStart = new Date() - 1000*(this.value - this.startTime)/this.rate/this.multiplier;
this.timerId = window.setInterval(tick, 1000*this.interval, this);
this.enabled = true;
};
/**
* Force value into legal range
*/
rgltimerClass.prototype.forceToRange = function() {
if (this.value > this.stopTime + this.stepSize/2 || this.value < this.startTime - this.stepSize/2) {
if (!this.loop) {
this.reset();
} else {
var cycle = this.stopTime - this.startTime + this.stepSize,
newval = (this.value - this.startTime) % cycle + this.startTime;
if (newval < this.startTime) {
newval += cycle;
}
this.realStart += (this.value - newval)*1000/this.multiplier/this.rate;
this.value = newval;
}
}
};
/**
* Reset to start values
*/
rgltimerClass.prototype.reset = function() {
this.value = this.startTime;
this.newmultiplier(1);
if (typeof this.Tick !== "undefined") {
this.Tick(this.value);
}
if (this.enabled)
this.play();  /* really pause... */
if (typeof this.PlayButton !== "undefined")
this.PlayButton.value = "Play";
};
/**
* Increase the multiplier to play faster
*/
rgltimerClass.prototype.faster = function() {
this.newmultiplier(Math.SQRT2*this.multiplier);
};
/**
* Decrease the multiplier to play slower
*/
rgltimerClass.prototype.slower = function() {
this.newmultiplier(this.multiplier/Math.SQRT2);
};
/**
* Change sign of multiplier to reverse direction
*/
rgltimerClass.prototype.reverse = function() {
this.newmultiplier(-this.multiplier);
};
/**
* Set multiplier for play speed
* @param { number } newmult - new value
*/
rgltimerClass.prototype.newmultiplier = function(newmult) {
if (newmult != this.multiplier) {
this.realStart += 1000*(this.value - this.startTime)/this.rate*(1/this.multiplier - 1/newmult);
this.multiplier = newmult;
}
};
/**
* Take one step
*/
rgltimerClass.prototype.step = function() {
this.value += this.rate*this.multiplier;
this.forceToRange();
if (typeof this.Tick !== "undefined")
this.Tick(this.value);
};</script>

<script type="text/javascript">
var powergrid_div = document.getElementById("powergrid_div"),
powergrid_rgl = new rglwidgetClass();
powergrid_div.width = 673;
powergrid_div.height = 481;
powergrid_rgl.initialize(powergrid_div,
{"material":{"color":"#000000","alpha":1,"lit":false,"ambient":"#000000","specular":"#FFFFFF","emission":"#000000","shininess":50,"smooth":true,"front":"filled","back":"filled","size":3,"lwd":1,"fog":false,"point_antialias":false,"line_antialias":false,"texture":null,"textype":"rgb","texmipmap":false,"texminfilter":"linear","texmagfilter":"linear","texenvmap":false,"depth_mask":true,"depth_test":"less","isTransparent":false,"polygon_offset":[0,0]},"rootSubscene":6,"objects":{"14":{"id":14,"type":"quads","material":{},"vertices":[[0,0,0],[0.32,0,0],[0.32,0,0.01],[0,0,0.01],[0.32,0,0],[0.32,0.5,0],[0.32,0.5,0.01],[0.32,0,0.01],[0.32,0.5,0],[0,0.5,0],[0,0.5,0.01],[0.32,0.5,0.01],[0,0.5,0],[0,0,0],[0,0,0.01],[0,0.5,0.01],[0,0,0.01],[0.32,0,0.01],[0.32,0.5,0.01],[0,0.5,0.01],[0,0,0],[0.32,0,0],[0.32,0.5,0],[0,0.5,0],[0,0,0.01],[0.32,0,0.01],[0.32,0,0.02],[0,0,0.02],[0.32,0,0.01],[0.32,0.5,0.01],[0.32,0.5,0.02],[0.32,0,0.02],[0.32,0.5,0.01],[0,0.5,0.01],[0,0.5,0.02],[0.32,0.5,0.02],[0,0.5,0.01],[0,0,0.01],[0,0,0.02],[0,0.5,0.02],[0,0,0.02],[0.32,0,0.02],[0.32,0.5,0.02],[0,0.5,0.02],[0,0,0.01],[0.32,0,0.01],[0.32,0.5,0.01],[0,0.5,0.01],[0,0,0.02],[0.32,0,0.02],[0.32,0,0.03],[0,0,0.03],[0.32,0,0.02],[0.32,0.5,0.02],[0.32,0.5,0.03],[0.32,0,0.03],[0.32,0.5,0.02],[0,0.5,0.02],[0,0.5,0.03],[0.32,0.5,0.03],[0,0.5,0.02],[0,0,0.02],[0,0,0.03],[0,0.5,0.03],[0,0,0.03],[0.32,0,0.03],[0.32,0.5,0.03],[0,0.5,0.03],[0,0,0.02],[0.32,0,0.02],[0.32,0.5,0.02],[0,0.5,0.02],[0,0,0.03],[0.32,0,0.03],[0.32,0,0.04],[0,0,0.04],[0.32,0,0.03],[0.32,0.5,0.03],[0.32,0.5,0.04],[0.32,0,0.04],[0.32,0.5,0.03],[0,0.5,0.03],[0,0.5,0.04],[0.32,0.5,0.04],[0,0.5,0.03],[0,0,0.03],[0,0,0.04],[0,0.5,0.04],[0,0,0.04],[0.32,0,0.04],[0.32,0.5,0.04],[0,0.5,0.04],[0,0,0.03],[0.32,0,0.03],[0.32,0.5,0.03],[0,0.5,0.03],[0,0,0.04],[0.32,0,0.04],[0.32,0,0.05],[0,0,0.05],[0.32,0,0.04],[0.32,0.5,0.04],[0.32,0.5,0.05],[0.32,0,0.05],[0.32,0.5,0.04],[0,0.5,0.04],[0,0.5,0.05],[0.32,0.5,0.05],[0,0.5,0.04],[0,0,0.04],[0,0,0.05],[0,0.5,0.05],[0,0,0.05],[0.32,0,0.05],[0.32,0.5,0.05],[0,0.5,0.05],[0,0,0.04],[0.32,0,0.04],[0.32,0.5,0.04],[0,0.5,0.04],[0,0,0.05],[0.32,0,0.05],[0.32,0,0.06],[0,0,0.06],[0.32,0,0.05],[0.32,0.5,0.05],[0.32,0.5,0.06],[0.32,0,0.06],[0.32,0.5,0.05],[0,0.5,0.05],[0,0.5,0.06],[0.32,0.5,0.06],[0,0.5,0.05],[0,0,0.05],[0,0,0.06],[0,0.5,0.06],[0,0,0.06],[0.32,0,0.06],[0.32,0.5,0.06],[0,0.5,0.06],[0,0,0.05],[0.32,0,0.05],[0.32,0.5,0.05],[0,0.5,0.05],[0,0,0.06],[0.32,0,0.06],[0.32,0,0.07],[0,0,0.07],[0.32,0,0.06],[0.32,0.5,0.06],[0.32,0.5,0.07],[0.32,0,0.07],[0.32,0.5,0.06],[0,0.5,0.06],[0,0.5,0.07],[0.32,0.5,0.07],[0,0.5,0.06],[0,0,0.06],[0,0,0.07],[0,0.5,0.07],[0,0,0.07],[0.32,0,0.07],[0.32,0.5,0.07],[0,0.5,0.07],[0,0,0.06],[0.32,0,0.06],[0.32,0.5,0.06],[0,0.5,0.06],[0,0,0.07],[0.32,0,0.07],[0.32,0,0.08],[0,0,0.08],[0.32,0,0.07],[0.32,0.5,0.07],[0.32,0.5,0.08],[0.32,0,0.08],[0.32,0.5,0.07],[0,0.5,0.07],[0,0.5,0.08],[0.32,0.5,0.08],[0,0.5,0.07],[0,0,0.07],[0,0,0.08],[0,0.5,0.08],[0,0,0.08],[0.32,0,0.08],[0.32,0.5,0.08],[0,0.5,0.08],[0,0,0.07],[0.32,0,0.07],[0.32,0.5,0.07],[0,0.5,0.07],[0,0,0.08],[0.32,0,0.08],[0.32,0,0.09],[0,0,0.09],[0.32,0,0.08],[0.32,0.5,0.08],[0.32,0.5,0.09],[0.32,0,0.09],[0.32,0.5,0.08],[0,0.5,0.08],[0,0.5,0.09],[0.32,0.5,0.09],[0,0.5,0.08],[0,0,0.08],[0,0,0.09],[0,0.5,0.09],[0,0,0.09],[0.32,0,0.09],[0.32,0.5,0.09],[0,0.5,0.09],[0,0,0.08],[0.32,0,0.08],[0.32,0.5,0.08],[0,0.5,0.08],[0,0,0.09],[0.32,0,0.09],[0.32,0,0.1],[0,0,0.1],[0.32,0,0.09],[0.32,0.5,0.09],[0.32,0.5,0.1],[0.32,0,0.1],[0.32,0.5,0.09],[0,0.5,0.09],[0,0.5,0.1],[0.32,0.5,0.1],[0,0.5,0.09],[0,0,0.09],[0,0,0.1],[0,0.5,0.1],[0,0,0.1],[0.32,0,0.1],[0.32,0.5,0.1],[0,0.5,0.1],[0,0,0.09],[0.32,0,0.09],[0.32,0.5,0.09],[0,0.5,0.09],[0,0,0.1],[0.32,0,0.1],[0.32,0,0.11],[0,0,0.11],[0.32,0,0.1],[0.32,0.5,0.1],[0.32,0.5,0.11],[0.32,0,0.11],[0.32,0.5,0.1],[0,0.5,0.1],[0,0.5,0.11],[0.32,0.5,0.11],[0,0.5,0.1],[0,0,0.1],[0,0,0.11],[0,0.5,0.11],[0,0,0.11],[0.32,0,0.11],[0.32,0.5,0.11],[0,0.5,0.11],[0,0,0.1],[0.32,0,0.1],[0.32,0.5,0.1],[0,0.5,0.1],[0,0,0.11],[0.32,0,0.11],[0.32,0,0.12],[0,0,0.12],[0.32,0,0.11],[0.32,0.5,0.11],[0.32,0.5,0.12],[0.32,0,0.12],[0.32,0.5,0.11],[0,0.5,0.11],[0,0.5,0.12],[0.32,0.5,0.12],[0,0.5,0.11],[0,0,0.11],[0,0,0.12],[0,0.5,0.12],[0,0,0.12],[0.32,0,0.12],[0.32,0.5,0.12],[0,0.5,0.12],[0,0,0.11],[0.32,0,0.11],[0.32,0.5,0.11],[0,0.5,0.11],[0,0,0.12],[0.32,0,0.12],[0.32,0,0.13],[0,0,0.13],[0.32,0,0.12],[0.32,0.5,0.12],[0.32,0.5,0.13],[0.32,0,0.13],[0.32,0.5,0.12],[0,0.5,0.12],[0,0.5,0.13],[0.32,0.5,0.13],[0,0.5,0.12],[0,0,0.12],[0,0,0.13],[0,0.5,0.13],[0,0,0.13],[0.32,0,0.13],[0.32,0.5,0.13],[0,0.5,0.13],[0,0,0.12],[0.32,0,0.12],[0.32,0.5,0.12],[0,0.5,0.12],[0,0,0.13],[0.32,0,0.13],[0.32,0,0.14],[0,0,0.14],[0.32,0,0.13],[0.32,0.5,0.13],[0.32,0.5,0.14],[0.32,0,0.14],[0.32,0.5,0.13],[0,0.5,0.13],[0,0.5,0.14],[0.32,0.5,0.14],[0,0.5,0.13],[0,0,0.13],[0,0,0.14],[0,0.5,0.14],[0,0,0.14],[0.32,0,0.14],[0.32,0.5,0.14],[0,0.5,0.14],[0,0,0.13],[0.32,0,0.13],[0.32,0.5,0.13],[0,0.5,0.13],[0,0,0.14],[0.32,0,0.14],[0.32,0,0.15],[0,0,0.15],[0.32,0,0.14],[0.32,0.5,0.14],[0.32,0.5,0.15],[0.32,0,0.15],[0.32,0.5,0.14],[0,0.5,0.14],[0,0.5,0.15],[0.32,0.5,0.15],[0,0.5,0.14],[0,0,0.14],[0,0,0.15],[0,0.5,0.15],[0,0,0.15],[0.32,0,0.15],[0.32,0.5,0.15],[0,0.5,0.15],[0,0,0.14],[0.32,0,0.14],[0.32,0.5,0.14],[0,0.5,0.14],[0,0,0.15],[0.32,0,0.15],[0.32,0,0.16],[0,0,0.16],[0.32,0,0.15],[0.32,0.5,0.15],[0.32,0.5,0.16],[0.32,0,0.16],[0.32,0.5,0.15],[0,0.5,0.15],[0,0.5,0.16],[0.32,0.5,0.16],[0,0.5,0.15],[0,0,0.15],[0,0,0.16],[0,0.5,0.16],[0,0,0.16],[0.32,0,0.16],[0.32,0.5,0.16],[0,0.5,0.16],[0,0,0.15],[0.32,0,0.15],[0.32,0.5,0.15],[0,0.5,0.15],[0,0,0.16],[0.32,0,0.16],[0.32,0,0.17],[0,0,0.17],[0.32,0,0.16],[0.32,0.5,0.16],[0.32,0.5,0.17],[0.32,0,0.17],[0.32,0.5,0.16],[0,0.5,0.16],[0,0.5,0.17],[0.32,0.5,0.17],[0,0.5,0.16],[0,0,0.16],[0,0,0.17],[0,0.5,0.17],[0,0,0.17],[0.32,0,0.17],[0.32,0.5,0.17],[0,0.5,0.17],[0,0,0.16],[0.32,0,0.16],[0.32,0.5,0.16],[0,0.5,0.16],[0,0,0.17],[0.32,0,0.17],[0.32,0,0.18],[0,0,0.18],[0.32,0,0.17],[0.32,0.5,0.17],[0.32,0.5,0.18],[0.32,0,0.18],[0.32,0.5,0.17],[0,0.5,0.17],[0,0.5,0.18],[0.32,0.5,0.18],[0,0.5,0.17],[0,0,0.17],[0,0,0.18],[0,0.5,0.18],[0,0,0.18],[0.32,0,0.18],[0.32,0.5,0.18],[0,0.5,0.18],[0,0,0.17],[0.32,0,0.17],[0.32,0.5,0.17],[0,0.5,0.17],[0,0,0.18],[0.32,0,0.18],[0.32,0,0.19],[0,0,0.19],[0.32,0,0.18],[0.32,0.5,0.18],[0.32,0.5,0.19],[0.32,0,0.19],[0.32,0.5,0.18],[0,0.5,0.18],[0,0.5,0.19],[0.32,0.5,0.19],[0,0.5,0.18],[0,0,0.18],[0,0,0.19],[0,0.5,0.19],[0,0,0.19],[0.32,0,0.19],[0.32,0.5,0.19],[0,0.5,0.19],[0,0,0.18],[0.32,0,0.18],[0.32,0.5,0.18],[0,0.5,0.18],[0,0,0.19],[0.32,0,0.19],[0.32,0,0.2],[0,0,0.2],[0.32,0,0.19],[0.32,0.5,0.19],[0.32,0.5,0.2],[0.32,0,0.2],[0.32,0.5,0.19],[0,0.5,0.19],[0,0.5,0.2],[0.32,0.5,0.2],[0,0.5,0.19],[0,0,0.19],[0,0,0.2],[0,0.5,0.2],[0,0,0.2],[0.32,0,0.2],[0.32,0.5,0.2],[0,0.5,0.2],[0,0,0.19],[0.32,0,0.19],[0.32,0.5,0.19],[0,0.5,0.19],[0,0,0.2],[0.32,0,0.2],[0.32,0,0.21],[0,0,0.21],[0.32,0,0.2],[0.32,0.5,0.2],[0.32,0.5,0.21],[0.32,0,0.21],[0.32,0.5,0.2],[0,0.5,0.2],[0,0.5,0.21],[0.32,0.5,0.21],[0,0.5,0.2],[0,0,0.2],[0,0,0.21],[0,0.5,0.21],[0,0,0.21],[0.32,0,0.21],[0.32,0.5,0.21],[0,0.5,0.21],[0,0,0.2],[0.32,0,0.2],[0.32,0.5,0.2],[0,0.5,0.2],[0,0,0.21],[0.32,0,0.21],[0.32,0,0.22],[0,0,0.22],[0.32,0,0.21],[0.32,0.5,0.21],[0.32,0.5,0.22],[0.32,0,0.22],[0.32,0.5,0.21],[0,0.5,0.21],[0,0.5,0.22],[0.32,0.5,0.22],[0,0.5,0.21],[0,0,0.21],[0,0,0.22],[0,0.5,0.22],[0,0,0.22],[0.32,0,0.22],[0.32,0.5,0.22],[0,0.5,0.22],[0,0,0.21],[0.32,0,0.21],[0.32,0.5,0.21],[0,0.5,0.21],[0,0,0.22],[0.32,0,0.22],[0.32,0,0.23],[0,0,0.23],[0.32,0,0.22],[0.32,0.5,0.22],[0.32,0.5,0.23],[0.32,0,0.23],[0.32,0.5,0.22],[0,0.5,0.22],[0,0.5,0.23],[0.32,0.5,0.23],[0,0.5,0.22],[0,0,0.22],[0,0,0.23],[0,0.5,0.23],[0,0,0.23],[0.32,0,0.23],[0.32,0.5,0.23],[0,0.5,0.23],[0,0,0.22],[0.32,0,0.22],[0.32,0.5,0.22],[0,0.5,0.22],[0,0,0.23],[0.32,0,0.23],[0.32,0,0.24],[0,0,0.24],[0.32,0,0.23],[0.32,0.5,0.23],[0.32,0.5,0.24],[0.32,0,0.24],[0.32,0.5,0.23],[0,0.5,0.23],[0,0.5,0.24],[0.32,0.5,0.24],[0,0.5,0.23],[0,0,0.23],[0,0,0.24],[0,0.5,0.24],[0,0,0.24],[0.32,0,0.24],[0.32,0.5,0.24],[0,0.5,0.24],[0,0,0.23],[0.32,0,0.23],[0.32,0.5,0.23],[0,0.5,0.23],[0,0,0.24],[0.32,0,0.24],[0.32,0,0.25],[0,0,0.25],[0.32,0,0.24],[0.32,0.5,0.24],[0.32,0.5,0.25],[0.32,0,0.25],[0.32,0.5,0.24],[0,0.5,0.24],[0,0.5,0.25],[0.32,0.5,0.25],[0,0.5,0.24],[0,0,0.24],[0,0,0.25],[0,0.5,0.25],[0,0,0.25],[0.32,0,0.25],[0.32,0.5,0.25],[0,0.5,0.25],[0,0,0.24],[0.32,0,0.24],[0.32,0.5,0.24],[0,0.5,0.24],[0,0,0.25],[0.32,0,0.25],[0.32,0,0.26],[0,0,0.26],[0.32,0,0.25],[0.32,0.5,0.25],[0.32,0.5,0.26],[0.32,0,0.26],[0.32,0.5,0.25],[0,0.5,0.25],[0,0.5,0.26],[0.32,0.5,0.26],[0,0.5,0.25],[0,0,0.25],[0,0,0.26],[0,0.5,0.26],[0,0,0.26],[0.32,0,0.26],[0.32,0.5,0.26],[0,0.5,0.26],[0,0,0.25],[0.32,0,0.25],[0.32,0.5,0.25],[0,0.5,0.25],[0,0,0.26],[0.32,0,0.26],[0.32,0,0.27],[0,0,0.27],[0.32,0,0.26],[0.32,0.5,0.26],[0.32,0.5,0.27],[0.32,0,0.27],[0.32,0.5,0.26],[0,0.5,0.26],[0,0.5,0.27],[0.32,0.5,0.27],[0,0.5,0.26],[0,0,0.26],[0,0,0.27],[0,0.5,0.27],[0,0,0.27],[0.32,0,0.27],[0.32,0.5,0.27],[0,0.5,0.27],[0,0,0.26],[0.32,0,0.26],[0.32,0.5,0.26],[0,0.5,0.26],[0,0,0.27],[0.32,0,0.27],[0.32,0,0.28],[0,0,0.28],[0.32,0,0.27],[0.32,0.5,0.27],[0.32,0.5,0.28],[0.32,0,0.28],[0.32,0.5,0.27],[0,0.5,0.27],[0,0.5,0.28],[0.32,0.5,0.28],[0,0.5,0.27],[0,0,0.27],[0,0,0.28],[0,0.5,0.28],[0,0,0.28],[0.32,0,0.28],[0.32,0.5,0.28],[0,0.5,0.28],[0,0,0.27],[0.32,0,0.27],[0.32,0.5,0.27],[0,0.5,0.27],[0,0,0.28],[0.32,0,0.28],[0.32,0,0.29],[0,0,0.29],[0.32,0,0.28],[0.32,0.5,0.28],[0.32,0.5,0.29],[0.32,0,0.29],[0.32,0.5,0.28],[0,0.5,0.28],[0,0.5,0.29],[0.32,0.5,0.29],[0,0.5,0.28],[0,0,0.28],[0,0,0.29],[0,0.5,0.29],[0,0,0.29],[0.32,0,0.29],[0.32,0.5,0.29],[0,0.5,0.29],[0,0,0.28],[0.32,0,0.28],[0.32,0.5,0.28],[0,0.5,0.28],[0,0,0.29],[0.32,0,0.29],[0.32,0,0.3],[0,0,0.3],[0.32,0,0.29],[0.32,0.5,0.29],[0.32,0.5,0.3],[0.32,0,0.3],[0.32,0.5,0.29],[0,0.5,0.29],[0,0.5,0.3],[0.32,0.5,0.3],[0,0.5,0.29],[0,0,0.29],[0,0,0.3],[0,0.5,0.3],[0,0,0.3],[0.32,0,0.3],[0.32,0.5,0.3],[0,0.5,0.3],[0,0,0.29],[0.32,0,0.29],[0.32,0.5,0.29],[0,0.5,0.29],[0,0,0.3],[0.32,0,0.3],[0.32,0,0.31],[0,0,0.31],[0.32,0,0.3],[0.32,0.5,0.3],[0.32,0.5,0.31],[0.32,0,0.31],[0.32,0.5,0.3],[0,0.5,0.3],[0,0.5,0.31],[0.32,0.5,0.31],[0,0.5,0.3],[0,0,0.3],[0,0,0.31],[0,0.5,0.31],[0,0,0.31],[0.32,0,0.31],[0.32,0.5,0.31],[0,0.5,0.31],[0,0,0.3],[0.32,0,0.3],[0.32,0.5,0.3],[0,0.5,0.3],[0,0,0.31],[0.32,0,0.31],[0.32,0,0.32],[0,0,0.32],[0.32,0,0.31],[0.32,0.5,0.31],[0.32,0.5,0.32],[0.32,0,0.32],[0.32,0.5,0.31],[0,0.5,0.31],[0,0.5,0.32],[0.32,0.5,0.32],[0,0.5,0.31],[0,0,0.31],[0,0,0.32],[0,0.5,0.32],[0,0,0.32],[0.32,0,0.32],[0.32,0.5,0.32],[0,0.5,0.32],[0,0,0.31],[0.32,0,0.31],[0.32,0.5,0.31],[0,0.5,0.31],[0,0,0.32],[0.32,0,0.32],[0.32,0,0.33],[0,0,0.33],[0.32,0,0.32],[0.32,0.5,0.32],[0.32,0.5,0.33],[0.32,0,0.33],[0.32,0.5,0.32],[0,0.5,0.32],[0,0.5,0.33],[0.32,0.5,0.33],[0,0.5,0.32],[0,0,0.32],[0,0,0.33],[0,0.5,0.33],[0,0,0.33],[0.32,0,0.33],[0.32,0.5,0.33],[0,0.5,0.33],[0,0,0.32],[0.32,0,0.32],[0.32,0.5,0.32],[0,0.5,0.32],[0,0,0.33],[0.32,0,0.33],[0.32,0,0.34],[0,0,0.34],[0.32,0,0.33],[0.32,0.5,0.33],[0.32,0.5,0.34],[0.32,0,0.34],[0.32,0.5,0.33],[0,0.5,0.33],[0,0.5,0.34],[0.32,0.5,0.34],[0,0.5,0.33],[0,0,0.33],[0,0,0.34],[0,0.5,0.34],[0,0,0.34],[0.32,0,0.34],[0.32,0.5,0.34],[0,0.5,0.34],[0,0,0.33],[0.32,0,0.33],[0.32,0.5,0.33],[0,0.5,0.33],[0,0,0.34],[0.32,0,0.34],[0.32,0,0.35],[0,0,0.35],[0.32,0,0.34],[0.32,0.5,0.34],[0.32,0.5,0.35],[0.32,0,0.35],[0.32,0.5,0.34],[0,0.5,0.34],[0,0.5,0.35],[0.32,0.5,0.35],[0,0.5,0.34],[0,0,0.34],[0,0,0.35],[0,0.5,0.35],[0,0,0.35],[0.32,0,0.35],[0.32,0.5,0.35],[0,0.5,0.35],[0,0,0.34],[0.32,0,0.34],[0.32,0.5,0.34],[0,0.5,0.34],[0,0,0.35],[0.32,0,0.35],[0.32,0,0.36],[0,0,0.36],[0.32,0,0.35],[0.32,0.5,0.35],[0.32,0.5,0.36],[0.32,0,0.36],[0.32,0.5,0.35],[0,0.5,0.35],[0,0.5,0.36],[0.32,0.5,0.36],[0,0.5,0.35],[0,0,0.35],[0,0,0.36],[0,0.5,0.36],[0,0,0.36],[0.32,0,0.36],[0.32,0.5,0.36],[0,0.5,0.36],[0,0,0.35],[0.32,0,0.35],[0.32,0.5,0.35],[0,0.5,0.35],[0,0,0.36],[0.32,0,0.36],[0.32,0,0.37],[0,0,0.37],[0.32,0,0.36],[0.32,0.5,0.36],[0.32,0.5,0.37],[0.32,0,0.37],[0.32,0.5,0.36],[0,0.5,0.36],[0,0.5,0.37],[0.32,0.5,0.37],[0,0.5,0.36],[0,0,0.36],[0,0,0.37],[0,0.5,0.37],[0,0,0.37],[0.32,0,0.37],[0.32,0.5,0.37],[0,0.5,0.37],[0,0,0.36],[0.32,0,0.36],[0.32,0.5,0.36],[0,0.5,0.36],[0,0,0.37],[0.32,0,0.37],[0.32,0,0.38],[0,0,0.38],[0.32,0,0.37],[0.32,0.5,0.37],[0.32,0.5,0.38],[0.32,0,0.38],[0.32,0.5,0.37],[0,0.5,0.37],[0,0.5,0.38],[0.32,0.5,0.38],[0,0.5,0.37],[0,0,0.37],[0,0,0.38],[0,0.5,0.38],[0,0,0.38],[0.32,0,0.38],[0.32,0.5,0.38],[0,0.5,0.38],[0,0,0.37],[0.32,0,0.37],[0.32,0.5,0.37],[0,0.5,0.37],[0,0,0.38],[0.32,0,0.38],[0.32,0,0.39],[0,0,0.39],[0.32,0,0.38],[0.32,0.5,0.38],[0.32,0.5,0.39],[0.32,0,0.39],[0.32,0.5,0.38],[0,0.5,0.38],[0,0.5,0.39],[0.32,0.5,0.39],[0,0.5,0.38],[0,0,0.38],[0,0,0.39],[0,0.5,0.39],[0,0,0.39],[0.32,0,0.39],[0.32,0.5,0.39],[0,0.5,0.39],[0,0,0.38],[0.32,0,0.38],[0.32,0.5,0.38],[0,0.5,0.38],[0,0,0.39],[0.32,0,0.39],[0.32,0,0.4],[0,0,0.4],[0.32,0,0.39],[0.32,0.5,0.39],[0.32,0.5,0.4],[0.32,0,0.4],[0.32,0.5,0.39],[0,0.5,0.39],[0,0.5,0.4],[0.32,0.5,0.4],[0,0.5,0.39],[0,0,0.39],[0,0,0.4],[0,0.5,0.4],[0,0,0.4],[0.32,0,0.4],[0.32,0.5,0.4],[0,0.5,0.4],[0,0,0.39],[0.32,0,0.39],[0.32,0.5,0.39],[0,0.5,0.39],[0,0,0.4],[0.32,0,0.4],[0.32,0,0.41],[0,0,0.41],[0.32,0,0.4],[0.32,0.5,0.4],[0.32,0.5,0.41],[0.32,0,0.41],[0.32,0.5,0.4],[0,0.5,0.4],[0,0.5,0.41],[0.32,0.5,0.41],[0,0.5,0.4],[0,0,0.4],[0,0,0.41],[0,0.5,0.41],[0,0,0.41],[0.32,0,0.41],[0.32,0.5,0.41],[0,0.5,0.41],[0,0,0.4],[0.32,0,0.4],[0.32,0.5,0.4],[0,0.5,0.4],[0,0,0.41],[0.32,0,0.41],[0.32,0,0.42],[0,0,0.42],[0.32,0,0.41],[0.32,0.5,0.41],[0.32,0.5,0.42],[0.32,0,0.42],[0.32,0.5,0.41],[0,0.5,0.41],[0,0.5,0.42],[0.32,0.5,0.42],[0,0.5,0.41],[0,0,0.41],[0,0,0.42],[0,0.5,0.42],[0,0,0.42],[0.32,0,0.42],[0.32,0.5,0.42],[0,0.5,0.42],[0,0,0.41],[0.32,0,0.41],[0.32,0.5,0.41],[0,0.5,0.41],[0,0,0.42],[0.32,0,0.42],[0.32,0,0.43],[0,0,0.43],[0.32,0,0.42],[0.32,0.5,0.42],[0.32,0.5,0.43],[0.32,0,0.43],[0.32,0.5,0.42],[0,0.5,0.42],[0,0.5,0.43],[0.32,0.5,0.43],[0,0.5,0.42],[0,0,0.42],[0,0,0.43],[0,0.5,0.43],[0,0,0.43],[0.32,0,0.43],[0.32,0.5,0.43],[0,0.5,0.43],[0,0,0.42],[0.32,0,0.42],[0.32,0.5,0.42],[0,0.5,0.42],[0,0,0.43],[0.32,0,0.43],[0.32,0,0.44],[0,0,0.44],[0.32,0,0.43],[0.32,0.5,0.43],[0.32,0.5,0.44],[0.32,0,0.44],[0.32,0.5,0.43],[0,0.5,0.43],[0,0.5,0.44],[0.32,0.5,0.44],[0,0.5,0.43],[0,0,0.43],[0,0,0.44],[0,0.5,0.44],[0,0,0.44],[0.32,0,0.44],[0.32,0.5,0.44],[0,0.5,0.44],[0,0,0.43],[0.32,0,0.43],[0.32,0.5,0.43],[0,0.5,0.43],[0,0,0.44],[0.32,0,0.44],[0.32,0,0.45],[0,0,0.45],[0.32,0,0.44],[0.32,0.5,0.44],[0.32,0.5,0.45],[0.32,0,0.45],[0.32,0.5,0.44],[0,0.5,0.44],[0,0.5,0.45],[0.32,0.5,0.45],[0,0.5,0.44],[0,0,0.44],[0,0,0.45],[0,0.5,0.45],[0,0,0.45],[0.32,0,0.45],[0.32,0.5,0.45],[0,0.5,0.45],[0,0,0.44],[0.32,0,0.44],[0.32,0.5,0.44],[0,0.5,0.44],[0,0,0.45],[0.32,0,0.45],[0.32,0,0.46],[0,0,0.46],[0.32,0,0.45],[0.32,0.5,0.45],[0.32,0.5,0.46],[0.32,0,0.46],[0.32,0.5,0.45],[0,0.5,0.45],[0,0.5,0.46],[0.32,0.5,0.46],[0,0.5,0.45],[0,0,0.45],[0,0,0.46],[0,0.5,0.46],[0,0,0.46],[0.32,0,0.46],[0.32,0.5,0.46],[0,0.5,0.46],[0,0,0.45],[0.32,0,0.45],[0.32,0.5,0.45],[0,0.5,0.45],[0,0,0.46],[0.32,0,0.46],[0.32,0,0.47],[0,0,0.47],[0.32,0,0.46],[0.32,0.5,0.46],[0.32,0.5,0.47],[0.32,0,0.47],[0.32,0.5,0.46],[0,0.5,0.46],[0,0.5,0.47],[0.32,0.5,0.47],[0,0.5,0.46],[0,0,0.46],[0,0,0.47],[0,0.5,0.47],[0,0,0.47],[0.32,0,0.47],[0.32,0.5,0.47],[0,0.5,0.47],[0,0,0.46],[0.32,0,0.46],[0.32,0.5,0.46],[0,0.5,0.46],[0,0,0.47],[0.32,0,0.47],[0.32,0,0.48],[0,0,0.48],[0.32,0,0.47],[0.32,0.5,0.47],[0.32,0.5,0.48],[0.32,0,0.48],[0.32,0.5,0.47],[0,0.5,0.47],[0,0.5,0.48],[0.32,0.5,0.48],[0,0.5,0.47],[0,0,0.47],[0,0,0.48],[0,0.5,0.48],[0,0,0.48],[0.32,0,0.48],[0.32,0.5,0.48],[0,0.5,0.48],[0,0,0.47],[0.32,0,0.47],[0.32,0.5,0.47],[0,0.5,0.47],[0,0,0.48],[0.32,0,0.48],[0.32,0,0.49],[0,0,0.49],[0.32,0,0.48],[0.32,0.5,0.48],[0.32,0.5,0.49],[0.32,0,0.49],[0.32,0.5,0.48],[0,0.5,0.48],[0,0.5,0.49],[0.32,0.5,0.49],[0,0.5,0.48],[0,0,0.48],[0,0,0.49],[0,0.5,0.49],[0,0,0.49],[0.32,0,0.49],[0.32,0.5,0.49],[0,0.5,0.49],[0,0,0.48],[0.32,0,0.48],[0.32,0.5,0.48],[0,0.5,0.48],[0,0,0.49],[0.32,0,0.49],[0.32,0,0.5],[0,0,0.5],[0.32,0,0.49],[0.32,0.5,0.49],[0.32,0.5,0.5],[0.32,0,0.5],[0.32,0.5,0.49],[0,0.5,0.49],[0,0.5,0.5],[0.32,0.5,0.5],[0,0.5,0.49],[0,0,0.49],[0,0,0.5],[0,0.5,0.5],[0,0,0.5],[0.32,0,0.5],[0.32,0.5,0.5],[0,0.5,0.5],[0,0,0.49],[0.32,0,0.49],[0.32,0.5,0.49],[0,0.5,0.49],[0,0,0.5],[0.32,0,0.5],[0.32,0,0.51],[0,0,0.51],[0.32,0,0.5],[0.32,0.5,0.5],[0.32,0.5,0.51],[0.32,0,0.51],[0.32,0.5,0.5],[0,0.5,0.5],[0,0.5,0.51],[0.32,0.5,0.51],[0,0.5,0.5],[0,0,0.5],[0,0,0.51],[0,0.5,0.51],[0,0,0.51],[0.32,0,0.51],[0.32,0.5,0.51],[0,0.5,0.51],[0,0,0.5],[0.32,0,0.5],[0.32,0.5,0.5],[0,0.5,0.5],[0,0,0.51],[0.32,0,0.51],[0.32,0,0.52],[0,0,0.52],[0.32,0,0.51],[0.32,0.5,0.51],[0.32,0.5,0.52],[0.32,0,0.52],[0.32,0.5,0.51],[0,0.5,0.51],[0,0.5,0.52],[0.32,0.5,0.52],[0,0.5,0.51],[0,0,0.51],[0,0,0.52],[0,0.5,0.52],[0,0,0.52],[0.32,0,0.52],[0.32,0.5,0.52],[0,0.5,0.52],[0,0,0.51],[0.32,0,0.51],[0.32,0.5,0.51],[0,0.5,0.51],[0,0,0.52],[0.32,0,0.52],[0.32,0,0.53],[0,0,0.53],[0.32,0,0.52],[0.32,0.5,0.52],[0.32,0.5,0.53],[0.32,0,0.53],[0.32,0.5,0.52],[0,0.5,0.52],[0,0.5,0.53],[0.32,0.5,0.53],[0,0.5,0.52],[0,0,0.52],[0,0,0.53],[0,0.5,0.53],[0,0,0.53],[0.32,0,0.53],[0.32,0.5,0.53],[0,0.5,0.53],[0,0,0.52],[0.32,0,0.52],[0.32,0.5,0.52],[0,0.5,0.52],[0,0,0.53],[0.32,0,0.53],[0.32,0,0.54],[0,0,0.54],[0.32,0,0.53],[0.32,0.5,0.53],[0.32,0.5,0.54],[0.32,0,0.54],[0.32,0.5,0.53],[0,0.5,0.53],[0,0.5,0.54],[0.32,0.5,0.54],[0,0.5,0.53],[0,0,0.53],[0,0,0.54],[0,0.5,0.54],[0,0,0.54],[0.32,0,0.54],[0.32,0.5,0.54],[0,0.5,0.54],[0,0,0.53],[0.32,0,0.53],[0.32,0.5,0.53],[0,0.5,0.53],[0,0,0.54],[0.32,0,0.54],[0.32,0,0.55],[0,0,0.55],[0.32,0,0.54],[0.32,0.5,0.54],[0.32,0.5,0.55],[0.32,0,0.55],[0.32,0.5,0.54],[0,0.5,0.54],[0,0.5,0.55],[0.32,0.5,0.55],[0,0.5,0.54],[0,0,0.54],[0,0,0.55],[0,0.5,0.55],[0,0,0.55],[0.32,0,0.55],[0.32,0.5,0.55],[0,0.5,0.55],[0,0,0.54],[0.32,0,0.54],[0.32,0.5,0.54],[0,0.5,0.54],[0,0,0.55],[0.32,0,0.55],[0.32,0,0.56],[0,0,0.56],[0.32,0,0.55],[0.32,0.5,0.55],[0.32,0.5,0.56],[0.32,0,0.56],[0.32,0.5,0.55],[0,0.5,0.55],[0,0.5,0.56],[0.32,0.5,0.56],[0,0.5,0.55],[0,0,0.55],[0,0,0.56],[0,0.5,0.56],[0,0,0.56],[0.32,0,0.56],[0.32,0.5,0.56],[0,0.5,0.56],[0,0,0.55],[0.32,0,0.55],[0.32,0.5,0.55],[0,0.5,0.55],[0,0,0.56],[0.32,0,0.56],[0.32,0,0.57],[0,0,0.57],[0.32,0,0.56],[0.32,0.5,0.56],[0.32,0.5,0.57],[0.32,0,0.57],[0.32,0.5,0.56],[0,0.5,0.56],[0,0.5,0.57],[0.32,0.5,0.57],[0,0.5,0.56],[0,0,0.56],[0,0,0.57],[0,0.5,0.57],[0,0,0.57],[0.32,0,0.57],[0.32,0.5,0.57],[0,0.5,0.57],[0,0,0.56],[0.32,0,0.56],[0.32,0.5,0.56],[0,0.5,0.56],[0,0,0.57],[0.32,0,0.57],[0.32,0,0.58],[0,0,0.58],[0.32,0,0.57],[0.32,0.5,0.57],[0.32,0.5,0.58],[0.32,0,0.58],[0.32,0.5,0.57],[0,0.5,0.57],[0,0.5,0.58],[0.32,0.5,0.58],[0,0.5,0.57],[0,0,0.57],[0,0,0.58],[0,0.5,0.58],[0,0,0.58],[0.32,0,0.58],[0.32,0.5,0.58],[0,0.5,0.58],[0,0,0.57],[0.32,0,0.57],[0.32,0.5,0.57],[0,0.5,0.57],[0,0,0.58],[0.32,0,0.58],[0.32,0,0.59],[0,0,0.59],[0.32,0,0.58],[0.32,0.5,0.58],[0.32,0.5,0.59],[0.32,0,0.59],[0.32,0.5,0.58],[0,0.5,0.58],[0,0.5,0.59],[0.32,0.5,0.59],[0,0.5,0.58],[0,0,0.58],[0,0,0.59],[0,0.5,0.59],[0,0,0.59],[0.32,0,0.59],[0.32,0.5,0.59],[0,0.5,0.59],[0,0,0.58],[0.32,0,0.58],[0.32,0.5,0.58],[0,0.5,0.58],[0,0,0.59],[0.32,0,0.59],[0.32,0,0.6],[0,0,0.6],[0.32,0,0.59],[0.32,0.5,0.59],[0.32,0.5,0.6],[0.32,0,0.6],[0.32,0.5,0.59],[0,0.5,0.59],[0,0.5,0.6],[0.32,0.5,0.6],[0,0.5,0.59],[0,0,0.59],[0,0,0.6],[0,0.5,0.6],[0,0,0.6],[0.32,0,0.6],[0.32,0.5,0.6],[0,0.5,0.6],[0,0,0.59],[0.32,0,0.59],[0.32,0.5,0.59],[0,0.5,0.59],[0,0,0.6],[0.32,0,0.6],[0.32,0,0.61],[0,0,0.61],[0.32,0,0.6],[0.32,0.5,0.6],[0.32,0.5,0.61],[0.32,0,0.61],[0.32,0.5,0.6],[0,0.5,0.6],[0,0.5,0.61],[0.32,0.5,0.61],[0,0.5,0.6],[0,0,0.6],[0,0,0.61],[0,0.5,0.61],[0,0,0.61],[0.32,0,0.61],[0.32,0.5,0.61],[0,0.5,0.61],[0,0,0.6],[0.32,0,0.6],[0.32,0.5,0.6],[0,0.5,0.6],[0,0,0.61],[0.32,0,0.61],[0.32,0,0.62],[0,0,0.62],[0.32,0,0.61],[0.32,0.5,0.61],[0.32,0.5,0.62],[0.32,0,0.62],[0.32,0.5,0.61],[0,0.5,0.61],[0,0.5,0.62],[0.32,0.5,0.62],[0,0.5,0.61],[0,0,0.61],[0,0,0.62],[0,0.5,0.62],[0,0,0.62],[0.32,0,0.62],[0.32,0.5,0.62],[0,0.5,0.62],[0,0,0.61],[0.32,0,0.61],[0.32,0.5,0.61],[0,0.5,0.61],[0,0,0.62],[0.32,0,0.62],[0.32,0,0.63],[0,0,0.63],[0.32,0,0.62],[0.32,0.5,0.62],[0.32,0.5,0.63],[0.32,0,0.63],[0.32,0.5,0.62],[0,0.5,0.62],[0,0.5,0.63],[0.32,0.5,0.63],[0,0.5,0.62],[0,0,0.62],[0,0,0.63],[0,0.5,0.63],[0,0,0.63],[0.32,0,0.63],[0.32,0.5,0.63],[0,0.5,0.63],[0,0,0.62],[0.32,0,0.62],[0.32,0.5,0.62],[0,0.5,0.62],[0,0,0.63],[0.32,0,0.63],[0.32,0,0.64],[0,0,0.64],[0.32,0,0.63],[0.32,0.5,0.63],[0.32,0.5,0.64],[0.32,0,0.64],[0.32,0.5,0.63],[0,0.5,0.63],[0,0.5,0.64],[0.32,0.5,0.64],[0,0.5,0.63],[0,0,0.63],[0,0,0.64],[0,0.5,0.64],[0,0,0.64],[0.32,0,0.64],[0.32,0.5,0.64],[0,0.5,0.64],[0,0,0.63],[0.32,0,0.63],[0.32,0.5,0.63],[0,0.5,0.63],[0,0,0.64],[0.32,0,0.64],[0.32,0,0.65],[0,0,0.65],[0.32,0,0.64],[0.32,0.5,0.64],[0.32,0.5,0.65],[0.32,0,0.65],[0.32,0.5,0.64],[0,0.5,0.64],[0,0.5,0.65],[0.32,0.5,0.65],[0,0.5,0.64],[0,0,0.64],[0,0,0.65],[0,0.5,0.65],[0,0,0.65],[0.32,0,0.65],[0.32,0.5,0.65],[0,0.5,0.65],[0,0,0.64],[0.32,0,0.64],[0.32,0.5,0.64],[0,0.5,0.64],[0,0,0.65],[0.32,0,0.65],[0.32,0,0.66],[0,0,0.66],[0.32,0,0.65],[0.32,0.5,0.65],[0.32,0.5,0.66],[0.32,0,0.66],[0.32,0.5,0.65],[0,0.5,0.65],[0,0.5,0.66],[0.32,0.5,0.66],[0,0.5,0.65],[0,0,0.65],[0,0,0.66],[0,0.5,0.66],[0,0,0.66],[0.32,0,0.66],[0.32,0.5,0.66],[0,0.5,0.66],[0,0,0.65],[0.32,0,0.65],[0.32,0.5,0.65],[0,0.5,0.65],[0,0,0.66],[0.32,0,0.66],[0.32,0,0.67],[0,0,0.67],[0.32,0,0.66],[0.32,0.5,0.66],[0.32,0.5,0.67],[0.32,0,0.67],[0.32,0.5,0.66],[0,0.5,0.66],[0,0.5,0.67],[0.32,0.5,0.67],[0,0.5,0.66],[0,0,0.66],[0,0,0.67],[0,0.5,0.67],[0,0,0.67],[0.32,0,0.67],[0.32,0.5,0.67],[0,0.5,0.67],[0,0,0.66],[0.32,0,0.66],[0.32,0.5,0.66],[0,0.5,0.66],[0,0,0.67],[0.32,0,0.67],[0.32,0,0.68],[0,0,0.68],[0.32,0,0.67],[0.32,0.5,0.67],[0.32,0.5,0.68],[0.32,0,0.68],[0.32,0.5,0.67],[0,0.5,0.67],[0,0.5,0.68],[0.32,0.5,0.68],[0,0.5,0.67],[0,0,0.67],[0,0,0.68],[0,0.5,0.68],[0,0,0.68],[0.32,0,0.68],[0.32,0.5,0.68],[0,0.5,0.68],[0,0,0.67],[0.32,0,0.67],[0.32,0.5,0.67],[0,0.5,0.67],[0,0,0.68],[0.32,0,0.68],[0.32,0,0.69],[0,0,0.69],[0.32,0,0.68],[0.32,0.5,0.68],[0.32,0.5,0.69],[0.32,0,0.69],[0.32,0.5,0.68],[0,0.5,0.68],[0,0.5,0.69],[0.32,0.5,0.69],[0,0.5,0.68],[0,0,0.68],[0,0,0.69],[0,0.5,0.69],[0,0,0.69],[0.32,0,0.69],[0.32,0.5,0.69],[0,0.5,0.69],[0,0,0.68],[0.32,0,0.68],[0.32,0.5,0.68],[0,0.5,0.68],[0,0,0.69],[0.32,0,0.69],[0.32,0,0.7],[0,0,0.7],[0.32,0,0.69],[0.32,0.5,0.69],[0.32,0.5,0.7],[0.32,0,0.7],[0.32,0.5,0.69],[0,0.5,0.69],[0,0.5,0.7],[0.32,0.5,0.7],[0,0.5,0.69],[0,0,0.69],[0,0,0.7],[0,0.5,0.7],[0,0,0.7],[0.32,0,0.7],[0.32,0.5,0.7],[0,0.5,0.7],[0,0,0.69],[0.32,0,0.69],[0.32,0.5,0.69],[0,0.5,0.69],[0,0,0.7],[0.32,0,0.7],[0.32,0,0.71],[0,0,0.71],[0.32,0,0.7],[0.32,0.5,0.7],[0.32,0.5,0.71],[0.32,0,0.71],[0.32,0.5,0.7],[0,0.5,0.7],[0,0.5,0.71],[0.32,0.5,0.71],[0,0.5,0.7],[0,0,0.7],[0,0,0.71],[0,0.5,0.71],[0,0,0.71],[0.32,0,0.71],[0.32,0.5,0.71],[0,0.5,0.71],[0,0,0.7],[0.32,0,0.7],[0.32,0.5,0.7],[0,0.5,0.7],[0,0,0.71],[0.32,0,0.71],[0.32,0,0.72],[0,0,0.72],[0.32,0,0.71],[0.32,0.5,0.71],[0.32,0.5,0.72],[0.32,0,0.72],[0.32,0.5,0.71],[0,0.5,0.71],[0,0.5,0.72],[0.32,0.5,0.72],[0,0.5,0.71],[0,0,0.71],[0,0,0.72],[0,0.5,0.72],[0,0,0.72],[0.32,0,0.72],[0.32,0.5,0.72],[0,0.5,0.72],[0,0,0.71],[0.32,0,0.71],[0.32,0.5,0.71],[0,0.5,0.71],[0,0,0.72],[0.32,0,0.72],[0.32,0,0.73],[0,0,0.73],[0.32,0,0.72],[0.32,0.5,0.72],[0.32,0.5,0.73],[0.32,0,0.73],[0.32,0.5,0.72],[0,0.5,0.72],[0,0.5,0.73],[0.32,0.5,0.73],[0,0.5,0.72],[0,0,0.72],[0,0,0.73],[0,0.5,0.73],[0,0,0.73],[0.32,0,0.73],[0.32,0.5,0.73],[0,0.5,0.73],[0,0,0.72],[0.32,0,0.72],[0.32,0.5,0.72],[0,0.5,0.72],[0,0,0.73],[0.32,0,0.73],[0.32,0,0.74],[0,0,0.74],[0.32,0,0.73],[0.32,0.5,0.73],[0.32,0.5,0.74],[0.32,0,0.74],[0.32,0.5,0.73],[0,0.5,0.73],[0,0.5,0.74],[0.32,0.5,0.74],[0,0.5,0.73],[0,0,0.73],[0,0,0.74],[0,0.5,0.74],[0,0,0.74],[0.32,0,0.74],[0.32,0.5,0.74],[0,0.5,0.74],[0,0,0.73],[0.32,0,0.73],[0.32,0.5,0.73],[0,0.5,0.73],[0,0,0.74],[0.32,0,0.74],[0.32,0,0.75],[0,0,0.75],[0.32,0,0.74],[0.32,0.5,0.74],[0.32,0.5,0.75],[0.32,0,0.75],[0.32,0.5,0.74],[0,0.5,0.74],[0,0.5,0.75],[0.32,0.5,0.75],[0,0.5,0.74],[0,0,0.74],[0,0,0.75],[0,0.5,0.75],[0,0,0.75],[0.32,0,0.75],[0.32,0.5,0.75],[0,0.5,0.75],[0,0,0.74],[0.32,0,0.74],[0.32,0.5,0.74],[0,0.5,0.74],[0,0,0.75],[0.32,0,0.75],[0.32,0,0.76],[0,0,0.76],[0.32,0,0.75],[0.32,0.5,0.75],[0.32,0.5,0.76],[0.32,0,0.76],[0.32,0.5,0.75],[0,0.5,0.75],[0,0.5,0.76],[0.32,0.5,0.76],[0,0.5,0.75],[0,0,0.75],[0,0,0.76],[0,0.5,0.76],[0,0,0.76],[0.32,0,0.76],[0.32,0.5,0.76],[0,0.5,0.76],[0,0,0.75],[0.32,0,0.75],[0.32,0.5,0.75],[0,0.5,0.75],[0,0,0.76],[0.32,0,0.76],[0.32,0,0.77],[0,0,0.77],[0.32,0,0.76],[0.32,0.5,0.76],[0.32,0.5,0.77],[0.32,0,0.77],[0.32,0.5,0.76],[0,0.5,0.76],[0,0.5,0.77],[0.32,0.5,0.77],[0,0.5,0.76],[0,0,0.76],[0,0,0.77],[0,0.5,0.77],[0,0,0.77],[0.32,0,0.77],[0.32,0.5,0.77],[0,0.5,0.77],[0,0,0.76],[0.32,0,0.76],[0.32,0.5,0.76],[0,0.5,0.76],[0,0,0.77],[0.32,0,0.77],[0.32,0,0.78],[0,0,0.78],[0.32,0,0.77],[0.32,0.5,0.77],[0.32,0.5,0.78],[0.32,0,0.78],[0.32,0.5,0.77],[0,0.5,0.77],[0,0.5,0.78],[0.32,0.5,0.78],[0,0.5,0.77],[0,0,0.77],[0,0,0.78],[0,0.5,0.78],[0,0,0.78],[0.32,0,0.78],[0.32,0.5,0.78],[0,0.5,0.78],[0,0,0.77],[0.32,0,0.77],[0.32,0.5,0.77],[0,0.5,0.77],[0,0,0.78],[0.32,0,0.78],[0.32,0,0.79],[0,0,0.79],[0.32,0,0.78],[0.32,0.5,0.78],[0.32,0.5,0.79],[0.32,0,0.79],[0.32,0.5,0.78],[0,0.5,0.78],[0,0.5,0.79],[0.32,0.5,0.79],[0,0.5,0.78],[0,0,0.78],[0,0,0.79],[0,0.5,0.79],[0,0,0.79],[0.32,0,0.79],[0.32,0.5,0.79],[0,0.5,0.79],[0,0,0.78],[0.32,0,0.78],[0.32,0.5,0.78],[0,0.5,0.78],[0,0,0.79],[0.32,0,0.79],[0.32,0,0.8],[0,0,0.8],[0.32,0,0.79],[0.32,0.5,0.79],[0.32,0.5,0.8],[0.32,0,0.8],[0.32,0.5,0.79],[0,0.5,0.79],[0,0.5,0.8],[0.32,0.5,0.8],[0,0.5,0.79],[0,0,0.79],[0,0,0.8],[0,0.5,0.8],[0,0,0.8],[0.32,0,0.8],[0.32,0.5,0.8],[0,0.5,0.8],[0,0,0.79],[0.32,0,0.79],[0.32,0.5,0.79],[0,0.5,0.79],[0,0,0.8],[0.32,0,0.8],[0.32,0,0.81],[0,0,0.81],[0.32,0,0.8],[0.32,0.5,0.8],[0.32,0.5,0.81],[0.32,0,0.81],[0.32,0.5,0.8],[0,0.5,0.8],[0,0.5,0.81],[0.32,0.5,0.81],[0,0.5,0.8],[0,0,0.8],[0,0,0.81],[0,0.5,0.81],[0,0,0.81],[0.32,0,0.81],[0.32,0.5,0.81],[0,0.5,0.81],[0,0,0.8],[0.32,0,0.8],[0.32,0.5,0.8],[0,0.5,0.8],[0,0,0.81],[0.32,0,0.81],[0.32,0,0.82],[0,0,0.82],[0.32,0,0.81],[0.32,0.5,0.81],[0.32,0.5,0.82],[0.32,0,0.82],[0.32,0.5,0.81],[0,0.5,0.81],[0,0.5,0.82],[0.32,0.5,0.82],[0,0.5,0.81],[0,0,0.81],[0,0,0.82],[0,0.5,0.82],[0,0,0.82],[0.32,0,0.82],[0.32,0.5,0.82],[0,0.5,0.82],[0,0,0.81],[0.32,0,0.81],[0.32,0.5,0.81],[0,0.5,0.81],[0,0,0.82],[0.32,0,0.82],[0.32,0,0.83],[0,0,0.83],[0.32,0,0.82],[0.32,0.5,0.82],[0.32,0.5,0.83],[0.32,0,0.83],[0.32,0.5,0.82],[0,0.5,0.82],[0,0.5,0.83],[0.32,0.5,0.83],[0,0.5,0.82],[0,0,0.82],[0,0,0.83],[0,0.5,0.83],[0,0,0.83],[0.32,0,0.83],[0.32,0.5,0.83],[0,0.5,0.83],[0,0,0.82],[0.32,0,0.82],[0.32,0.5,0.82],[0,0.5,0.82],[0,0,0.83],[0.32,0,0.83],[0.32,0,0.84],[0,0,0.84],[0.32,0,0.83],[0.32,0.5,0.83],[0.32,0.5,0.84],[0.32,0,0.84],[0.32,0.5,0.83],[0,0.5,0.83],[0,0.5,0.84],[0.32,0.5,0.84],[0,0.5,0.83],[0,0,0.83],[0,0,0.84],[0,0.5,0.84],[0,0,0.84],[0.32,0,0.84],[0.32,0.5,0.84],[0,0.5,0.84],[0,0,0.83],[0.32,0,0.83],[0.32,0.5,0.83],[0,0.5,0.83],[0,0,0.84],[0.32,0,0.84],[0.32,0,0.85],[0,0,0.85],[0.32,0,0.84],[0.32,0.5,0.84],[0.32,0.5,0.85],[0.32,0,0.85],[0.32,0.5,0.84],[0,0.5,0.84],[0,0.5,0.85],[0.32,0.5,0.85],[0,0.5,0.84],[0,0,0.84],[0,0,0.85],[0,0.5,0.85],[0,0,0.85],[0.32,0,0.85],[0.32,0.5,0.85],[0,0.5,0.85],[0,0,0.84],[0.32,0,0.84],[0.32,0.5,0.84],[0,0.5,0.84],[0,0,0.85],[0.32,0,0.85],[0.32,0,0.86],[0,0,0.86],[0.32,0,0.85],[0.32,0.5,0.85],[0.32,0.5,0.86],[0.32,0,0.86],[0.32,0.5,0.85],[0,0.5,0.85],[0,0.5,0.86],[0.32,0.5,0.86],[0,0.5,0.85],[0,0,0.85],[0,0,0.86],[0,0.5,0.86],[0,0,0.86],[0.32,0,0.86],[0.32,0.5,0.86],[0,0.5,0.86],[0,0,0.85],[0.32,0,0.85],[0.32,0.5,0.85],[0,0.5,0.85],[0,0,0.86],[0.32,0,0.86],[0.32,0,0.87],[0,0,0.87],[0.32,0,0.86],[0.32,0.5,0.86],[0.32,0.5,0.87],[0.32,0,0.87],[0.32,0.5,0.86],[0,0.5,0.86],[0,0.5,0.87],[0.32,0.5,0.87],[0,0.5,0.86],[0,0,0.86],[0,0,0.87],[0,0.5,0.87],[0,0,0.87],[0.32,0,0.87],[0.32,0.5,0.87],[0,0.5,0.87],[0,0,0.86],[0.32,0,0.86],[0.32,0.5,0.86],[0,0.5,0.86],[0,0,0.87],[0.32,0,0.87],[0.32,0,0.88],[0,0,0.88],[0.32,0,0.87],[0.32,0.5,0.87],[0.32,0.5,0.88],[0.32,0,0.88],[0.32,0.5,0.87],[0,0.5,0.87],[0,0.5,0.88],[0.32,0.5,0.88],[0,0.5,0.87],[0,0,0.87],[0,0,0.88],[0,0.5,0.88],[0,0,0.88],[0.32,0,0.88],[0.32,0.5,0.88],[0,0.5,0.88],[0,0,0.87],[0.32,0,0.87],[0.32,0.5,0.87],[0,0.5,0.87],[0,0,0.88],[0.32,0,0.88],[0.32,0,0.89],[0,0,0.89],[0.32,0,0.88],[0.32,0.5,0.88],[0.32,0.5,0.89],[0.32,0,0.89],[0.32,0.5,0.88],[0,0.5,0.88],[0,0.5,0.89],[0.32,0.5,0.89],[0,0.5,0.88],[0,0,0.88],[0,0,0.89],[0,0.5,0.89],[0,0,0.89],[0.32,0,0.89],[0.32,0.5,0.89],[0,0.5,0.89],[0,0,0.88],[0.32,0,0.88],[0.32,0.5,0.88],[0,0.5,0.88],[0,0,0.89],[0.32,0,0.89],[0.32,0,0.9],[0,0,0.9],[0.32,0,0.89],[0.32,0.5,0.89],[0.32,0.5,0.9],[0.32,0,0.9],[0.32,0.5,0.89],[0,0.5,0.89],[0,0.5,0.9],[0.32,0.5,0.9],[0,0.5,0.89],[0,0,0.89],[0,0,0.9],[0,0.5,0.9],[0,0,0.9],[0.32,0,0.9],[0.32,0.5,0.9],[0,0.5,0.9],[0,0,0.89],[0.32,0,0.89],[0.32,0.5,0.89],[0,0.5,0.89],[0,0,0.9],[0.32,0,0.9],[0.32,0,0.91],[0,0,0.91],[0.32,0,0.9],[0.32,0.5,0.9],[0.32,0.5,0.91],[0.32,0,0.91],[0.32,0.5,0.9],[0,0.5,0.9],[0,0.5,0.91],[0.32,0.5,0.91],[0,0.5,0.9],[0,0,0.9],[0,0,0.91],[0,0.5,0.91],[0,0,0.91],[0.32,0,0.91],[0.32,0.5,0.91],[0,0.5,0.91],[0,0,0.9],[0.32,0,0.9],[0.32,0.5,0.9],[0,0.5,0.9],[0,0,0.91],[0.32,0,0.91],[0.32,0,0.92],[0,0,0.92],[0.32,0,0.91],[0.32,0.5,0.91],[0.32,0.5,0.92],[0.32,0,0.92],[0.32,0.5,0.91],[0,0.5,0.91],[0,0.5,0.92],[0.32,0.5,0.92],[0,0.5,0.91],[0,0,0.91],[0,0,0.92],[0,0.5,0.92],[0,0,0.92],[0.32,0,0.92],[0.32,0.5,0.92],[0,0.5,0.92],[0,0,0.91],[0.32,0,0.91],[0.32,0.5,0.91],[0,0.5,0.91],[0,0,0.92],[0.32,0,0.92],[0.32,0,0.93],[0,0,0.93],[0.32,0,0.92],[0.32,0.5,0.92],[0.32,0.5,0.93],[0.32,0,0.93],[0.32,0.5,0.92],[0,0.5,0.92],[0,0.5,0.93],[0.32,0.5,0.93],[0,0.5,0.92],[0,0,0.92],[0,0,0.93],[0,0.5,0.93],[0,0,0.93],[0.32,0,0.93],[0.32,0.5,0.93],[0,0.5,0.93],[0,0,0.92],[0.32,0,0.92],[0.32,0.5,0.92],[0,0.5,0.92],[0,0,0.93],[0.32,0,0.93],[0.32,0,0.94],[0,0,0.94],[0.32,0,0.93],[0.32,0.5,0.93],[0.32,0.5,0.94],[0.32,0,0.94],[0.32,0.5,0.93],[0,0.5,0.93],[0,0.5,0.94],[0.32,0.5,0.94],[0,0.5,0.93],[0,0,0.93],[0,0,0.94],[0,0.5,0.94],[0,0,0.94],[0.32,0,0.94],[0.32,0.5,0.94],[0,0.5,0.94],[0,0,0.93],[0.32,0,0.93],[0.32,0.5,0.93],[0,0.5,0.93],[0,0,0.94],[0.32,0,0.94],[0.32,0,0.95],[0,0,0.95],[0.32,0,0.94],[0.32,0.5,0.94],[0.32,0.5,0.95],[0.32,0,0.95],[0.32,0.5,0.94],[0,0.5,0.94],[0,0.5,0.95],[0.32,0.5,0.95],[0,0.5,0.94],[0,0,0.94],[0,0,0.95],[0,0.5,0.95],[0,0,0.95],[0.32,0,0.95],[0.32,0.5,0.95],[0,0.5,0.95],[0,0,0.94],[0.32,0,0.94],[0.32,0.5,0.94],[0,0.5,0.94],[0,0,0.95],[0.32,0,0.95],[0.32,0,0.96],[0,0,0.96],[0.32,0,0.95],[0.32,0.5,0.95],[0.32,0.5,0.96],[0.32,0,0.96],[0.32,0.5,0.95],[0,0.5,0.95],[0,0.5,0.96],[0.32,0.5,0.96],[0,0.5,0.95],[0,0,0.95],[0,0,0.96],[0,0.5,0.96],[0,0,0.96],[0.32,0,0.96],[0.32,0.5,0.96],[0,0.5,0.96],[0,0,0.95],[0.32,0,0.95],[0.32,0.5,0.95],[0,0.5,0.95],[0,0,0.96],[0.32,0,0.96],[0.32,0,0.97],[0,0,0.97],[0.32,0,0.96],[0.32,0.5,0.96],[0.32,0.5,0.97],[0.32,0,0.97],[0.32,0.5,0.96],[0,0.5,0.96],[0,0.5,0.97],[0.32,0.5,0.97],[0,0.5,0.96],[0,0,0.96],[0,0,0.97],[0,0.5,0.97],[0,0,0.97],[0.32,0,0.97],[0.32,0.5,0.97],[0,0.5,0.97],[0,0,0.96],[0.32,0,0.96],[0.32,0.5,0.96],[0,0.5,0.96],[0,0,0.97],[0.32,0,0.97],[0.32,0,0.98],[0,0,0.98],[0.32,0,0.97],[0.32,0.5,0.97],[0.32,0.5,0.98],[0.32,0,0.98],[0.32,0.5,0.97],[0,0.5,0.97],[0,0.5,0.98],[0.32,0.5,0.98],[0,0.5,0.97],[0,0,0.97],[0,0,0.98],[0,0.5,0.98],[0,0,0.98],[0.32,0,0.98],[0.32,0.5,0.98],[0,0.5,0.98],[0,0,0.97],[0.32,0,0.97],[0.32,0.5,0.97],[0,0.5,0.97],[0,0,0.98],[0.32,0,0.98],[0.32,0,0.99],[0,0,0.99],[0.32,0,0.98],[0.32,0.5,0.98],[0.32,0.5,0.99],[0.32,0,0.99],[0.32,0.5,0.98],[0,0.5,0.98],[0,0.5,0.99],[0.32,0.5,0.99],[0,0.5,0.98],[0,0,0.98],[0,0,0.99],[0,0.5,0.99],[0,0,0.99],[0.32,0,0.99],[0.32,0.5,0.99],[0,0.5,0.99],[0,0,0.98],[0.32,0,0.98],[0.32,0.5,0.98],[0,0.5,0.98],[0,0,0.99],[0.32,0,0.99],[0.32,0,1],[0,0,1],[0.32,0,0.99],[0.32,0.5,0.99],[0.32,0.5,1],[0.32,0,1],[0.32,0.5,0.99],[0,0.5,0.99],[0,0.5,1],[0.32,0.5,1],[0,0.5,0.99],[0,0,0.99],[0,0,1],[0,0.5,1],[0,0,1],[0.32,0,1],[0.32,0.5,1],[0,0.5,1],[0,0,0.99],[0.32,0,0.99],[0.32,0.5,0.99],[0,0.5,0.99]],"colors":[[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.01960784,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.03921569,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.05882353,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.07843138,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.09803922,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1176471,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1372549,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1568628,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1764706,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.1960784,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2156863,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.2352941,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.254902,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2745098,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.2941177,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3137255,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3333333,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.3529412,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.372549,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.3921569,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4117647,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4313726,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4509804,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4705882,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.4901961,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.509804,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5294118,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5490196,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5686275,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.5882353,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.6078432,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.627451,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.6470588,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.654902,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6666667,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6784314,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6862745,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.6980392,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7098039,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7215686,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7294118,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7411765,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7529412,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.7647059,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.772549,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7843137,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.7960784,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8039216,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.8156863,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.827451,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8392157,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8470588,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8588235,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8705882,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8823529,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.8901961,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9019608,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9137255,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9215686,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.9333333,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.945098,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9568627,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9647059,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9764706,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[1,0.9882353,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9960784,0.9960784,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9686275,0.9843137,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9372549,0.9686275,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.9058824,0.9568627,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8784314,0.9411765,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8470588,0.9294118,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.8156863,0.9137255,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7843137,0.9019608,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7568628,0.8862745,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.7254902,0.8745098,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6941177,0.8588235,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6627451,0.8470588,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6352941,0.8313726,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.6039216,0.8196079,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.572549,0.8039216,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5450981,0.7921569,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.5137255,0.7764706,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4823529,0.7647059,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4509804,0.7490196,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.4235294,0.7372549,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3921569,0.7215686,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3607843,0.7098039,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3294118,0.6941177,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.3019608,0.682353,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2705882,0.6666667,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2392157,0.654902,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.2117647,0.6392157,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1803922,0.627451,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1490196,0.6117647,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.1176471,0.6,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.09019608,0.5843138,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.05882353,0.572549,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0.02745098,0.5568628,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1],[0,0.5450981,0,1]],"centers":[[0.16,0,0.005],[0.32,0.25,0.005],[0.16,0.5,0.005],[0,0.25,0.005],[0.16,0.25,0.01],[0.16,0.25,0],[0.16,0,0.015],[0.32,0.25,0.015],[0.16,0.5,0.015],[0,0.25,0.015],[0.16,0.25,0.02],[0.16,0.25,0.01],[0.16,0,0.025],[0.32,0.25,0.025],[0.16,0.5,0.025],[0,0.25,0.025],[0.16,0.25,0.03],[0.16,0.25,0.02],[0.16,0,0.035],[0.32,0.25,0.035],[0.16,0.5,0.035],[0,0.25,0.035],[0.16,0.25,0.04],[0.16,0.25,0.03],[0.16,0,0.045],[0.32,0.25,0.045],[0.16,0.5,0.045],[0,0.25,0.045],[0.16,0.25,0.05],[0.16,0.25,0.04],[0.16,0,0.055],[0.32,0.25,0.055],[0.16,0.5,0.055],[0,0.25,0.055],[0.16,0.25,0.06],[0.16,0.25,0.05],[0.16,0,0.065],[0.32,0.25,0.065],[0.16,0.5,0.065],[0,0.25,0.065],[0.16,0.25,0.07],[0.16,0.25,0.06],[0.16,0,0.075],[0.32,0.25,0.075],[0.16,0.5,0.075],[0,0.25,0.075],[0.16,0.25,0.08],[0.16,0.25,0.07],[0.16,0,0.085],[0.32,0.25,0.085],[0.16,0.5,0.085],[0,0.25,0.085],[0.16,0.25,0.09],[0.16,0.25,0.08],[0.16,0,0.095],[0.32,0.25,0.095],[0.16,0.5,0.095],[0,0.25,0.095],[0.16,0.25,0.1],[0.16,0.25,0.09],[0.16,0,0.105],[0.32,0.25,0.105],[0.16,0.5,0.105],[0,0.25,0.105],[0.16,0.25,0.11],[0.16,0.25,0.1],[0.16,0,0.115],[0.32,0.25,0.115],[0.16,0.5,0.115],[0,0.25,0.115],[0.16,0.25,0.12],[0.16,0.25,0.11],[0.16,0,0.125],[0.32,0.25,0.125],[0.16,0.5,0.125],[0,0.25,0.125],[0.16,0.25,0.13],[0.16,0.25,0.12],[0.16,0,0.135],[0.32,0.25,0.135],[0.16,0.5,0.135],[0,0.25,0.135],[0.16,0.25,0.14],[0.16,0.25,0.13],[0.16,0,0.145],[0.32,0.25,0.145],[0.16,0.5,0.145],[0,0.25,0.145],[0.16,0.25,0.15],[0.16,0.25,0.14],[0.16,0,0.155],[0.32,0.25,0.155],[0.16,0.5,0.155],[0,0.25,0.155],[0.16,0.25,0.16],[0.16,0.25,0.15],[0.16,0,0.165],[0.32,0.25,0.165],[0.16,0.5,0.165],[0,0.25,0.165],[0.16,0.25,0.17],[0.16,0.25,0.16],[0.16,0,0.175],[0.32,0.25,0.175],[0.16,0.5,0.175],[0,0.25,0.175],[0.16,0.25,0.18],[0.16,0.25,0.17],[0.16,0,0.185],[0.32,0.25,0.185],[0.16,0.5,0.185],[0,0.25,0.185],[0.16,0.25,0.19],[0.16,0.25,0.18],[0.16,0,0.195],[0.32,0.25,0.195],[0.16,0.5,0.195],[0,0.25,0.195],[0.16,0.25,0.2],[0.16,0.25,0.19],[0.16,0,0.205],[0.32,0.25,0.205],[0.16,0.5,0.205],[0,0.25,0.205],[0.16,0.25,0.21],[0.16,0.25,0.2],[0.16,0,0.215],[0.32,0.25,0.215],[0.16,0.5,0.215],[0,0.25,0.215],[0.16,0.25,0.22],[0.16,0.25,0.21],[0.16,0,0.225],[0.32,0.25,0.225],[0.16,0.5,0.225],[0,0.25,0.225],[0.16,0.25,0.23],[0.16,0.25,0.22],[0.16,0,0.235],[0.32,0.25,0.235],[0.16,0.5,0.235],[0,0.25,0.235],[0.16,0.25,0.24],[0.16,0.25,0.23],[0.16,0,0.245],[0.32,0.25,0.245],[0.16,0.5,0.245],[0,0.25,0.245],[0.16,0.25,0.25],[0.16,0.25,0.24],[0.16,0,0.255],[0.32,0.25,0.255],[0.16,0.5,0.255],[0,0.25,0.255],[0.16,0.25,0.26],[0.16,0.25,0.25],[0.16,0,0.265],[0.32,0.25,0.265],[0.16,0.5,0.265],[0,0.25,0.265],[0.16,0.25,0.27],[0.16,0.25,0.26],[0.16,0,0.275],[0.32,0.25,0.275],[0.16,0.5,0.275],[0,0.25,0.275],[0.16,0.25,0.28],[0.16,0.25,0.27],[0.16,0,0.285],[0.32,0.25,0.285],[0.16,0.5,0.285],[0,0.25,0.285],[0.16,0.25,0.29],[0.16,0.25,0.28],[0.16,0,0.295],[0.32,0.25,0.295],[0.16,0.5,0.295],[0,0.25,0.295],[0.16,0.25,0.3],[0.16,0.25,0.29],[0.16,0,0.305],[0.32,0.25,0.305],[0.16,0.5,0.305],[0,0.25,0.305],[0.16,0.25,0.31],[0.16,0.25,0.3],[0.16,0,0.315],[0.32,0.25,0.315],[0.16,0.5,0.315],[0,0.25,0.315],[0.16,0.25,0.32],[0.16,0.25,0.31],[0.16,0,0.325],[0.32,0.25,0.325],[0.16,0.5,0.325],[0,0.25,0.325],[0.16,0.25,0.33],[0.16,0.25,0.32],[0.16,0,0.335],[0.32,0.25,0.335],[0.16,0.5,0.335],[0,0.25,0.335],[0.16,0.25,0.34],[0.16,0.25,0.33],[0.16,0,0.345],[0.32,0.25,0.345],[0.16,0.5,0.345],[0,0.25,0.345],[0.16,0.25,0.35],[0.16,0.25,0.34],[0.16,0,0.355],[0.32,0.25,0.355],[0.16,0.5,0.355],[0,0.25,0.355],[0.16,0.25,0.36],[0.16,0.25,0.35],[0.16,0,0.365],[0.32,0.25,0.365],[0.16,0.5,0.365],[0,0.25,0.365],[0.16,0.25,0.37],[0.16,0.25,0.36],[0.16,0,0.375],[0.32,0.25,0.375],[0.16,0.5,0.375],[0,0.25,0.375],[0.16,0.25,0.38],[0.16,0.25,0.37],[0.16,0,0.385],[0.32,0.25,0.385],[0.16,0.5,0.385],[0,0.25,0.385],[0.16,0.25,0.39],[0.16,0.25,0.38],[0.16,0,0.395],[0.32,0.25,0.395],[0.16,0.5,0.395],[0,0.25,0.395],[0.16,0.25,0.4],[0.16,0.25,0.39],[0.16,0,0.405],[0.32,0.25,0.405],[0.16,0.5,0.405],[0,0.25,0.405],[0.16,0.25,0.41],[0.16,0.25,0.4],[0.16,0,0.415],[0.32,0.25,0.415],[0.16,0.5,0.415],[0,0.25,0.415],[0.16,0.25,0.42],[0.16,0.25,0.41],[0.16,0,0.425],[0.32,0.25,0.425],[0.16,0.5,0.425],[0,0.25,0.425],[0.16,0.25,0.43],[0.16,0.25,0.42],[0.16,0,0.435],[0.32,0.25,0.435],[0.16,0.5,0.435],[0,0.25,0.435],[0.16,0.25,0.44],[0.16,0.25,0.43],[0.16,0,0.445],[0.32,0.25,0.445],[0.16,0.5,0.445],[0,0.25,0.445],[0.16,0.25,0.45],[0.16,0.25,0.44],[0.16,0,0.455],[0.32,0.25,0.455],[0.16,0.5,0.455],[0,0.25,0.455],[0.16,0.25,0.46],[0.16,0.25,0.45],[0.16,0,0.465],[0.32,0.25,0.465],[0.16,0.5,0.465],[0,0.25,0.465],[0.16,0.25,0.47],[0.16,0.25,0.46],[0.16,0,0.475],[0.32,0.25,0.475],[0.16,0.5,0.475],[0,0.25,0.475],[0.16,0.25,0.48],[0.16,0.25,0.47],[0.16,0,0.485],[0.32,0.25,0.485],[0.16,0.5,0.485],[0,0.25,0.485],[0.16,0.25,0.49],[0.16,0.25,0.48],[0.16,0,0.495],[0.32,0.25,0.495],[0.16,0.5,0.495],[0,0.25,0.495],[0.16,0.25,0.5],[0.16,0.25,0.49],[0.16,0,0.505],[0.32,0.25,0.505],[0.16,0.5,0.505],[0,0.25,0.505],[0.16,0.25,0.51],[0.16,0.25,0.5],[0.16,0,0.515],[0.32,0.25,0.515],[0.16,0.5,0.515],[0,0.25,0.515],[0.16,0.25,0.52],[0.16,0.25,0.51],[0.16,0,0.525],[0.32,0.25,0.525],[0.16,0.5,0.525],[0,0.25,0.525],[0.16,0.25,0.53],[0.16,0.25,0.52],[0.16,0,0.535],[0.32,0.25,0.535],[0.16,0.5,0.535],[0,0.25,0.535],[0.16,0.25,0.54],[0.16,0.25,0.53],[0.16,0,0.545],[0.32,0.25,0.545],[0.16,0.5,0.545],[0,0.25,0.545],[0.16,0.25,0.55],[0.16,0.25,0.54],[0.16,0,0.555],[0.32,0.25,0.555],[0.16,0.5,0.555],[0,0.25,0.555],[0.16,0.25,0.56],[0.16,0.25,0.55],[0.16,0,0.565],[0.32,0.25,0.565],[0.16,0.5,0.565],[0,0.25,0.565],[0.16,0.25,0.57],[0.16,0.25,0.56],[0.16,0,0.575],[0.32,0.25,0.575],[0.16,0.5,0.575],[0,0.25,0.575],[0.16,0.25,0.58],[0.16,0.25,0.57],[0.16,0,0.585],[0.32,0.25,0.585],[0.16,0.5,0.585],[0,0.25,0.585],[0.16,0.25,0.59],[0.16,0.25,0.58],[0.16,0,0.595],[0.32,0.25,0.595],[0.16,0.5,0.595],[0,0.25,0.595],[0.16,0.25,0.6],[0.16,0.25,0.59],[0.16,0,0.605],[0.32,0.25,0.605],[0.16,0.5,0.605],[0,0.25,0.605],[0.16,0.25,0.61],[0.16,0.25,0.6],[0.16,0,0.615],[0.32,0.25,0.615],[0.16,0.5,0.615],[0,0.25,0.615],[0.16,0.25,0.62],[0.16,0.25,0.61],[0.16,0,0.625],[0.32,0.25,0.625],[0.16,0.5,0.625],[0,0.25,0.625],[0.16,0.25,0.63],[0.16,0.25,0.62],[0.16,0,0.635],[0.32,0.25,0.635],[0.16,0.5,0.635],[0,0.25,0.635],[0.16,0.25,0.64],[0.16,0.25,0.63],[0.16,0,0.645],[0.32,0.25,0.645],[0.16,0.5,0.645],[0,0.25,0.645],[0.16,0.25,0.65],[0.16,0.25,0.64],[0.16,0,0.655],[0.32,0.25,0.655],[0.16,0.5,0.655],[0,0.25,0.655],[0.16,0.25,0.66],[0.16,0.25,0.65],[0.16,0,0.665],[0.32,0.25,0.665],[0.16,0.5,0.665],[0,0.25,0.665],[0.16,0.25,0.67],[0.16,0.25,0.66],[0.16,0,0.675],[0.32,0.25,0.675],[0.16,0.5,0.675],[0,0.25,0.675],[0.16,0.25,0.68],[0.16,0.25,0.67],[0.16,0,0.685],[0.32,0.25,0.685],[0.16,0.5,0.685],[0,0.25,0.685],[0.16,0.25,0.69],[0.16,0.25,0.68],[0.16,0,0.695],[0.32,0.25,0.695],[0.16,0.5,0.695],[0,0.25,0.695],[0.16,0.25,0.7],[0.16,0.25,0.69],[0.16,0,0.705],[0.32,0.25,0.705],[0.16,0.5,0.705],[0,0.25,0.705],[0.16,0.25,0.71],[0.16,0.25,0.7],[0.16,0,0.715],[0.32,0.25,0.715],[0.16,0.5,0.715],[0,0.25,0.715],[0.16,0.25,0.72],[0.16,0.25,0.71],[0.16,0,0.725],[0.32,0.25,0.725],[0.16,0.5,0.725],[0,0.25,0.725],[0.16,0.25,0.73],[0.16,0.25,0.72],[0.16,0,0.735],[0.32,0.25,0.735],[0.16,0.5,0.735],[0,0.25,0.735],[0.16,0.25,0.74],[0.16,0.25,0.73],[0.16,0,0.745],[0.32,0.25,0.745],[0.16,0.5,0.745],[0,0.25,0.745],[0.16,0.25,0.75],[0.16,0.25,0.74],[0.16,0,0.755],[0.32,0.25,0.755],[0.16,0.5,0.755],[0,0.25,0.755],[0.16,0.25,0.76],[0.16,0.25,0.75],[0.16,0,0.765],[0.32,0.25,0.765],[0.16,0.5,0.765],[0,0.25,0.765],[0.16,0.25,0.77],[0.16,0.25,0.76],[0.16,0,0.775],[0.32,0.25,0.775],[0.16,0.5,0.775],[0,0.25,0.775],[0.16,0.25,0.78],[0.16,0.25,0.77],[0.16,0,0.785],[0.32,0.25,0.785],[0.16,0.5,0.785],[0,0.25,0.785],[0.16,0.25,0.79],[0.16,0.25,0.78],[0.16,0,0.795],[0.32,0.25,0.795],[0.16,0.5,0.795],[0,0.25,0.795],[0.16,0.25,0.8],[0.16,0.25,0.79],[0.16,0,0.805],[0.32,0.25,0.805],[0.16,0.5,0.805],[0,0.25,0.805],[0.16,0.25,0.81],[0.16,0.25,0.8],[0.16,0,0.815],[0.32,0.25,0.815],[0.16,0.5,0.815],[0,0.25,0.815],[0.16,0.25,0.82],[0.16,0.25,0.81],[0.16,0,0.825],[0.32,0.25,0.825],[0.16,0.5,0.825],[0,0.25,0.825],[0.16,0.25,0.83],[0.16,0.25,0.82],[0.16,0,0.835],[0.32,0.25,0.835],[0.16,0.5,0.835],[0,0.25,0.835],[0.16,0.25,0.84],[0.16,0.25,0.83],[0.16,0,0.845],[0.32,0.25,0.845],[0.16,0.5,0.845],[0,0.25,0.845],[0.16,0.25,0.85],[0.16,0.25,0.84],[0.16,0,0.855],[0.32,0.25,0.855],[0.16,0.5,0.855],[0,0.25,0.855],[0.16,0.25,0.86],[0.16,0.25,0.85],[0.16,0,0.865],[0.32,0.25,0.865],[0.16,0.5,0.865],[0,0.25,0.865],[0.16,0.25,0.87],[0.16,0.25,0.86],[0.16,0,0.875],[0.32,0.25,0.875],[0.16,0.5,0.875],[0,0.25,0.875],[0.16,0.25,0.88],[0.16,0.25,0.87],[0.16,0,0.885],[0.32,0.25,0.885],[0.16,0.5,0.885],[0,0.25,0.885],[0.16,0.25,0.89],[0.16,0.25,0.88],[0.16,0,0.895],[0.32,0.25,0.895],[0.16,0.5,0.895],[0,0.25,0.895],[0.16,0.25,0.9],[0.16,0.25,0.89],[0.16,0,0.905],[0.32,0.25,0.905],[0.16,0.5,0.905],[0,0.25,0.905],[0.16,0.25,0.91],[0.16,0.25,0.9],[0.16,0,0.915],[0.32,0.25,0.915],[0.16,0.5,0.915],[0,0.25,0.915],[0.16,0.25,0.92],[0.16,0.25,0.91],[0.16,0,0.925],[0.32,0.25,0.925],[0.16,0.5,0.925],[0,0.25,0.925],[0.16,0.25,0.93],[0.16,0.25,0.92],[0.16,0,0.935],[0.32,0.25,0.935],[0.16,0.5,0.935],[0,0.25,0.935],[0.16,0.25,0.94],[0.16,0.25,0.93],[0.16,0,0.945],[0.32,0.25,0.945],[0.16,0.5,0.945],[0,0.25,0.945],[0.16,0.25,0.95],[0.16,0.25,0.94],[0.16,0,0.955],[0.32,0.25,0.955],[0.16,0.5,0.955],[0,0.25,0.955],[0.16,0.25,0.96],[0.16,0.25,0.95],[0.16,0,0.965],[0.32,0.25,0.965],[0.16,0.5,0.965],[0,0.25,0.965],[0.16,0.25,0.97],[0.16,0.25,0.96],[0.16,0,0.975],[0.32,0.25,0.975],[0.16,0.5,0.975],[0,0.25,0.975],[0.16,0.25,0.98],[0.16,0.25,0.97],[0.16,0,0.985],[0.32,0.25,0.985],[0.16,0.5,0.985],[0,0.25,0.985],[0.16,0.25,0.99],[0.16,0.25,0.98],[0.16,0,0.995],[0.32,0.25,0.995],[0.16,0.5,0.995],[0,0.25,0.995],[0.16,0.25,1],[0.16,0.25,0.99]],"ignoreExtent":false,"flags":2},"15":{"id":15,"type":"lines","material":{},"vertices":[[0,0,0],[0.32,0,0],[0.32,0,0],[0.32,0,0],[0.32,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,1],[0.32,0,0],[0.32,0,1],[0.32,0,0],[0.32,0,1],[0,0,0],[0,0,1],[0,0,1],[0.32,0,1],[0.32,0,1],[0.32,0,1],[0.32,0,1],[0,0,1],[0,0,1],[0,0,1],[0.32,0,0.1638639],[0.42,0,0.1638639],[0.32,0,0.3728979],[0.42,0,0.3728979],[0.32,0,0.5819319],[0.42,0,0.5819319],[0.32,0,0.790966],[0.42,0,0.790966],[0.32,0,1],[0.42,0,1]],"colors":[[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0.5607843,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1]],"centers":[[0.16,0,0],[0.32,0,0],[0.16,0,0],[0,0,0],[0,0,0.5],[0.32,0,0.5],[0.32,0,0.5],[0,0,0.5],[0.16,0,1],[0.32,0,1],[0.16,0,1],[0,0,1],[0.37,0,0.1638639],[0.37,0,0.3728979],[0.37,0,0.5819319],[0.37,0,0.790966],[0.37,0,1]],"ignoreExtent":false,"flags":64},"16":{"id":16,"type":"text","material":{},"vertices":[[0.44,0,0.1638639],[0.44,0,0.3728979],[0.44,0,0.5819319],[0.44,0,0.790966],[0.44,0,1]],"colors":[[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1]],"texts":[["0.2"],["0.4"],["0.6"],["0.8"],["1"]],"cex":[[1],[1],[1],[1],[1]],"adj":[[0,0.5]],"centers":[[0.44,0,0.1638639],[0.44,0,0.3728979],[0.44,0,0.5819319],[0.44,0,0.790966],[0.44,0,1]],"family":[["sans"],["sans"],["sans"],["sans"],["sans"]],"font":[[1],[1],[1],[1],[1]],"ignoreExtent":false,"flags":2064},"17":{"id":17,"type":"text","material":{},"vertices":[[0.22,-0.08475,-0.1695]],"colors":[[0,0,0,1]],"texts":[[""]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[0.22,-0.08475,-0.1695]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"18":{"id":18,"type":"text","material":{},"vertices":[[-0.07458,0.25,-0.1695]],"colors":[[0,0,0,1]],"texts":[[""]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-0.07458,0.25,-0.1695]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"19":{"id":19,"type":"text","material":{},"vertices":[[-0.07458,-0.08475,0.5]],"colors":[[0,0,0,1]],"texts":[[""]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-0.07458,-0.08475,0.5]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"22":{"id":22,"type":"quads","material":{},"vertices":[[0.575,0.575,0.04321798],[1.425,0.575,0.04321798],[1.425,0.575,0.04321798],[0.575,0.575,0.04321798],[1.575,0.575,0.04321798],[2.425,0.575,0.04321798],[2.425,0.575,0.06596308],[1.575,0.575,0.06596308],[2.575,0.575,0.04321798],[3.425,0.575,0.04321798],[3.425,0.575,0.1021959],[2.575,0.575,0.1021959],[3.575,0.575,0.04321798],[4.425,0.575,0.04321798],[4.425,0.575,0.1075265],[3.575,0.575,0.1075265],[4.575,0.575,0.04321798],[5.425,0.575,0.04321798],[5.425,0.575,0.2316452],[4.575,0.575,0.2316452],[5.575,0.575,0.04321798],[6.425,0.575,0.04321798],[6.425,0.575,0.2258095],[5.575,0.575,0.2258095],[6.575,0.575,0.04321798],[7.425,0.575,0.04321798],[7.425,0.575,0.237572],[6.575,0.575,0.237572],[7.575,0.575,0.04321798],[8.425,0.575,0.04321798],[8.425,0.575,0.3690261],[7.575,0.575,0.3690261],[8.575,0.575,0.04321798],[9.425,0.575,0.04321798],[9.425,0.575,0.3592858],[8.575,0.575,0.3592858],[0.575,1.575,0.04321798],[1.425,1.575,0.04321798],[1.425,1.575,0.05915068],[0.575,1.575,0.05915068],[1.575,1.575,0.04321798],[2.425,1.575,0.04321798],[2.425,1.575,0.1108488],[1.575,1.575,0.1108488],[2.575,1.575,0.04321798],[3.425,1.575,0.04321798],[3.425,1.575,0.2030372],[2.575,1.575,0.2030372],[3.575,1.575,0.04321798],[4.425,1.575,0.04321798],[4.425,1.575,0.217041],[3.575,1.575,0.217041],[4.575,1.575,0.04321798],[5.425,1.575,0.04321798],[5.425,1.575,0.5257641],[4.575,1.575,0.5257641],[5.575,1.575,0.04321798],[6.425,1.575,0.04321798],[6.425,1.575,0.5127993],[5.575,1.575,0.5127993],[6.575,1.575,0.04321798],[7.425,1.575,0.04321798],[7.425,1.575,0.5387281],[6.575,1.575,0.5387281],[7.575,1.575,0.04321798],[8.425,1.575,0.04321798],[8.425,1.575,0.7713429],[7.575,1.575,0.7713429],[8.575,1.575,0.04321798],[9.425,1.575,0.04321798],[9.425,1.575,0.7576808],[8.575,1.575,0.7576808],[0.575,2.575,0.04321798],[1.425,2.575,0.04321798],[1.425,2.575,0.07725442],[0.575,2.575,0.07725442],[1.575,2.575,0.04321798],[2.425,2.575,0.04321798],[2.425,2.575,0.1672267],[1.575,2.575,0.1672267],[2.575,2.575,0.04321798],[3.425,2.575,0.04321798],[3.425,2.575,0.3312058],[2.575,2.575,0.3312058],[3.575,2.575,0.04321798],[4.425,2.575,0.04321798],[4.425,2.575,0.3553697],[3.575,2.575,0.3553697],[4.575,2.575,0.04321798],[5.425,2.575,0.04321798],[5.425,2.575,0.781293],[4.575,2.575,0.781293],[5.575,2.575,0.04321798],[6.425,2.575,0.04321798],[6.425,2.575,0.76804],[5.575,2.575,0.76804],[6.575,2.575,0.04321798],[7.425,2.575,0.04321798],[7.425,2.575,0.7941267],[6.575,2.575,0.7941267],[7.575,2.575,0.04321798],[8.425,2.575,0.04321798],[8.425,2.575,0.9553382],[7.575,2.575,0.9553382],[8.575,2.575,0.04321798],[9.425,2.575,0.04321798],[9.425,2.575,0.9493445],[8.575,2.575,0.9493445],[0.575,3.575,0.04321798],[1.425,3.575,0.04321798],[1.425,3.575,0.1025677],[0.575,3.575,0.1025677],[1.575,3.575,0.04321798],[2.425,3.575,0.04321798],[2.425,3.575,0.2501325],[1.575,3.575,0.2501325],[2.575,3.575,0.04321798],[3.425,3.575,0.04321798],[3.425,3.575,0.5039053],[2.575,3.575,0.5039053],[3.575,3.575,0.04321798],[4.425,3.575,0.04321798],[4.425,3.575,0.537781],[3.575,3.575,0.537781],[4.575,3.575,0.04321798],[5.425,3.575,0.04321798],[5.425,3.575,0.9430385],[4.575,3.575,0.9430385],[5.575,3.575,0.04321798],[6.425,3.575,0.04321798],[6.425,3.575,0.9361931],[5.575,3.575,0.9361931],[6.575,3.575,0.04321798],[7.425,3.575,0.04321798],[7.425,3.575,0.9493053],[6.575,3.575,0.9493053],[7.575,3.575,0.04321798],[8.425,3.575,0.04321798],[8.425,3.575,0.9972544],[7.575,3.575,0.9972544],[8.575,3.575,0.04321798],[9.425,3.575,0.04321798],[9.425,3.575,0.9965154],[8.575,3.575,0.9965154],[0.575,4.575,0.04321798],[1.425,4.575,0.04321798],[1.425,4.575,0.1352801],[0.575,4.575,0.1352801],[1.575,4.575,0.04321798],[2.425,4.575,0.04321798],[2.425,4.575,0.3567418],[1.575,4.575,0.3567418],[2.575,4.575,0.04321798],[3.425,4.575,0.04321798],[3.425,4.575,0.6844499],[2.575,4.575,0.6844499],[3.575,4.575,0.04321798],[4.425,4.575,0.04321798],[4.425,4.575,0.7209224],[3.575,4.575,0.7209224],[4.575,4.575,0.04321798],[5.425,4.575,0.04321798],[5.425,4.575,0.9923994],[4.575,4.575,0.9923994],[5.575,4.575,0.04321798],[6.425,4.575,0.04321798],[6.425,4.575,0.9907784],[5.575,4.575,0.9907784],[6.575,4.575,0.04321798],[7.425,4.575,0.04321798],[7.425,4.575,0.9937668],[6.575,4.575,0.9937668],[7.575,4.575,0.04321798],[8.425,4.575,0.04321798],[8.425,4.575,0.999954],[7.575,4.575,0.999954],[8.575,4.575,0.04321798],[9.425,4.575,0.04321798],[9.425,4.575,0.9999315],[8.575,4.575,0.9999315],[0.575,5.575,0.04321798],[1.425,5.575,0.04321798],[1.425,5.575,0.1997644],[0.575,5.575,0.1997644],[1.575,5.575,0.04321798],[2.425,5.575,0.04321798],[2.425,5.575,0.5454682],[1.575,5.575,0.5454682],[2.575,5.575,0.04321798],[3.425,5.575,0.04321798],[3.425,5.575,0.8883711],[2.575,5.575,0.8883711],[3.575,5.575,0.04321798],[4.425,5.575,0.04321798],[4.425,5.575,0.912325],[3.575,5.575,0.912325],[4.575,5.575,0.04321798],[5.425,5.575,0.04321798],[5.425,5.575,0.9999133],[4.575,5.575,0.9999133],[5.575,5.575,0.04321798],[6.425,5.575,0.04321798],[6.425,5.575,0.9998751],[5.575,5.575,0.9998751],[6.575,5.575,0.04321798],[7.425,5.575,0.04321798],[7.425,5.575,0.9999403],[6.575,5.575,0.9999403],[7.575,5.575,0.04321798],[8.425,5.575,0.04321798],[8.425,5.575,1],[7.575,5.575,1],[8.575,5.575,0.04321798],[9.425,5.575,0.04321798],[9.425,5.575,1],[8.575,5.575,1],[0.575,6.575,0.04321798],[1.425,6.575,0.04321798],[1.425,6.575,0.2790861],[0.575,6.575,0.2790861],[1.575,6.575,0.04321798],[2.425,6.575,0.04321798],[2.425,6.575,0.722779],[1.575,6.575,0.722779],[2.575,6.575,0.04321798],[3.425,6.575,0.04321798],[3.425,6.575,0.9743328],[2.575,6.575,0.9743328],[3.575,6.575,0.04321798],[4.425,6.575,0.04321798],[4.425,6.575,0.9828939],[3.575,6.575,0.9828939],[4.575,6.575,0.04321798],[5.425,6.575,0.04321798],[5.425,6.575,0.9999998],[4.575,6.575,0.9999998],[5.575,6.575,0.04321798],[6.425,6.575,0.04321798],[6.425,6.575,0.9999996],[5.575,6.575,0.9999996],[6.575,6.575,0.04321798],[7.425,6.575,0.04321798],[7.425,6.575,0.9999999],[6.575,6.575,0.9999999],[7.575,6.575,0.04321798],[8.425,6.575,0.04321798],[8.425,6.575,1],[7.575,6.575,1],[8.575,6.575,0.04321798],[9.425,6.575,0.04321798],[9.425,6.575,1],[8.575,6.575,1],[0.575,7.575,0.04321798],[1.425,7.575,0.04321798],[1.425,7.575,0.428676],[0.575,7.575,0.428676],[1.575,7.575,0.04321798],[2.425,7.575,0.04321798],[2.425,7.575,0.9105917],[1.575,7.575,0.9105917],[2.575,7.575,0.04321798],[3.425,7.575,0.04321798],[3.425,7.575,0.9990419],[2.575,7.575,0.9990419],[3.575,7.575,0.04321798],[4.425,7.575,0.04321798],[4.425,7.575,0.9995526],[3.575,7.575,0.9995526],[4.575,7.575,0.04321798],[5.425,7.575,0.04321798],[5.425,7.575,1],[4.575,7.575,1],[5.575,7.575,0.04321798],[6.425,7.575,0.04321798],[6.425,7.575,1],[5.575,7.575,1],[6.575,7.575,0.04321798],[7.425,7.575,0.04321798],[7.425,7.575,1],[6.575,7.575,1],[7.575,7.575,0.04321798],[8.425,7.575,0.04321798],[8.425,7.575,1],[7.575,7.575,1],[8.575,7.575,0.04321798],[9.425,7.575,0.04321798],[9.425,7.575,1],[8.575,7.575,1],[0.575,8.575,0.04321798],[1.425,8.575,0.04321798],[1.425,8.575,0.5885873],[0.575,8.575,0.5885873],[1.575,8.575,0.04321798],[2.425,8.575,0.04321798],[2.425,8.575,0.981856],[1.575,8.575,0.981856],[2.575,8.575,0.04321798],[3.425,8.575,0.04321798],[3.425,8.575,0.9999893],[2.575,8.575,0.9999893],[3.575,8.575,0.04321798],[4.425,8.575,0.04321798],[4.425,8.575,0.9999969],[3.575,8.575,0.9999969],[4.575,8.575,0.04321798],[5.425,8.575,0.04321798],[5.425,8.575,1],[4.575,8.575,1],[5.575,8.575,0.04321798],[6.425,8.575,0.04321798],[6.425,8.575,1],[5.575,8.575,1],[6.575,8.575,0.04321798],[7.425,8.575,0.04321798],[7.425,8.575,1],[6.575,8.575,1],[7.575,8.575,0.04321798],[8.425,8.575,0.04321798],[8.425,8.575,1],[7.575,8.575,1],[8.575,8.575,0.04321798],[9.425,8.575,0.04321798],[9.425,8.575,1],[8.575,8.575,1],[0.575,9.575,0.04321798],[1.425,9.575,0.04321798],[1.425,9.575,0.7135577],[0.575,9.575,0.7135577],[1.575,9.575,0.04321798],[2.425,9.575,0.04321798],[2.425,9.575,0.9967811],[1.575,9.575,0.9967811],[2.575,9.575,0.04321798],[3.425,9.575,0.04321798],[3.425,9.575,0.9999999],[2.575,9.575,0.9999999],[3.575,9.575,0.04321798],[4.425,9.575,0.04321798],[4.425,9.575,1],[3.575,9.575,1],[4.575,9.575,0.04321798],[5.425,9.575,0.04321798],[5.425,9.575,1],[4.575,9.575,1],[5.575,9.575,0.04321798],[6.425,9.575,0.04321798],[6.425,9.575,1],[5.575,9.575,1],[6.575,9.575,0.04321798],[7.425,9.575,0.04321798],[7.425,9.575,1],[6.575,9.575,1],[7.575,9.575,0.04321798],[8.425,9.575,0.04321798],[8.425,9.575,1],[7.575,9.575,1],[8.575,9.575,0.04321798],[9.425,9.575,0.04321798],[9.425,9.575,1],[8.575,9.575,1],[1.425,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,1.425,0.04321798],[1.425,1.425,0.04321798],[2.425,1.425,0.04321798],[1.575,1.425,0.04321798],[1.575,1.425,0.06596308],[2.425,1.425,0.06596308],[3.425,1.425,0.04321798],[2.575,1.425,0.04321798],[2.575,1.425,0.1021959],[3.425,1.425,0.1021959],[4.425,1.425,0.04321798],[3.575,1.425,0.04321798],[3.575,1.425,0.1075265],[4.425,1.425,0.1075265],[5.425,1.425,0.04321798],[4.575,1.425,0.04321798],[4.575,1.425,0.2316452],[5.425,1.425,0.2316452],[6.425,1.425,0.04321798],[5.575,1.425,0.04321798],[5.575,1.425,0.2258095],[6.425,1.425,0.2258095],[7.425,1.425,0.04321798],[6.575,1.425,0.04321798],[6.575,1.425,0.237572],[7.425,1.425,0.237572],[8.425,1.425,0.04321798],[7.575,1.425,0.04321798],[7.575,1.425,0.3690261],[8.425,1.425,0.3690261],[9.425,1.425,0.04321798],[8.575,1.425,0.04321798],[8.575,1.425,0.3592858],[9.425,1.425,0.3592858],[1.425,2.425,0.04321798],[0.575,2.425,0.04321798],[0.575,2.425,0.05915068],[1.425,2.425,0.05915068],[2.425,2.425,0.04321798],[1.575,2.425,0.04321798],[1.575,2.425,0.1108488],[2.425,2.425,0.1108488],[3.425,2.425,0.04321798],[2.575,2.425,0.04321798],[2.575,2.425,0.2030372],[3.425,2.425,0.2030372],[4.425,2.425,0.04321798],[3.575,2.425,0.04321798],[3.575,2.425,0.217041],[4.425,2.425,0.217041],[5.425,2.425,0.04321798],[4.575,2.425,0.04321798],[4.575,2.425,0.5257641],[5.425,2.425,0.5257641],[6.425,2.425,0.04321798],[5.575,2.425,0.04321798],[5.575,2.425,0.5127993],[6.425,2.425,0.5127993],[7.425,2.425,0.04321798],[6.575,2.425,0.04321798],[6.575,2.425,0.5387281],[7.425,2.425,0.5387281],[8.425,2.425,0.04321798],[7.575,2.425,0.04321798],[7.575,2.425,0.7713429],[8.425,2.425,0.7713429],[9.425,2.425,0.04321798],[8.575,2.425,0.04321798],[8.575,2.425,0.7576808],[9.425,2.425,0.7576808],[1.425,3.425,0.04321798],[0.575,3.425,0.04321798],[0.575,3.425,0.07725442],[1.425,3.425,0.07725442],[2.425,3.425,0.04321798],[1.575,3.425,0.04321798],[1.575,3.425,0.1672267],[2.425,3.425,0.1672267],[3.425,3.425,0.04321798],[2.575,3.425,0.04321798],[2.575,3.425,0.3312058],[3.425,3.425,0.3312058],[4.425,3.425,0.04321798],[3.575,3.425,0.04321798],[3.575,3.425,0.3553697],[4.425,3.425,0.3553697],[5.425,3.425,0.04321798],[4.575,3.425,0.04321798],[4.575,3.425,0.781293],[5.425,3.425,0.781293],[6.425,3.425,0.04321798],[5.575,3.425,0.04321798],[5.575,3.425,0.76804],[6.425,3.425,0.76804],[7.425,3.425,0.04321798],[6.575,3.425,0.04321798],[6.575,3.425,0.7941267],[7.425,3.425,0.7941267],[8.425,3.425,0.04321798],[7.575,3.425,0.04321798],[7.575,3.425,0.9553382],[8.425,3.425,0.9553382],[9.425,3.425,0.04321798],[8.575,3.425,0.04321798],[8.575,3.425,0.9493445],[9.425,3.425,0.9493445],[1.425,4.425,0.04321798],[0.575,4.425,0.04321798],[0.575,4.425,0.1025677],[1.425,4.425,0.1025677],[2.425,4.425,0.04321798],[1.575,4.425,0.04321798],[1.575,4.425,0.2501325],[2.425,4.425,0.2501325],[3.425,4.425,0.04321798],[2.575,4.425,0.04321798],[2.575,4.425,0.5039053],[3.425,4.425,0.5039053],[4.425,4.425,0.04321798],[3.575,4.425,0.04321798],[3.575,4.425,0.537781],[4.425,4.425,0.537781],[5.425,4.425,0.04321798],[4.575,4.425,0.04321798],[4.575,4.425,0.9430385],[5.425,4.425,0.9430385],[6.425,4.425,0.04321798],[5.575,4.425,0.04321798],[5.575,4.425,0.9361931],[6.425,4.425,0.9361931],[7.425,4.425,0.04321798],[6.575,4.425,0.04321798],[6.575,4.425,0.9493053],[7.425,4.425,0.9493053],[8.425,4.425,0.04321798],[7.575,4.425,0.04321798],[7.575,4.425,0.9972544],[8.425,4.425,0.9972544],[9.425,4.425,0.04321798],[8.575,4.425,0.04321798],[8.575,4.425,0.9965154],[9.425,4.425,0.9965154],[1.425,5.425,0.04321798],[0.575,5.425,0.04321798],[0.575,5.425,0.1352801],[1.425,5.425,0.1352801],[2.425,5.425,0.04321798],[1.575,5.425,0.04321798],[1.575,5.425,0.3567418],[2.425,5.425,0.3567418],[3.425,5.425,0.04321798],[2.575,5.425,0.04321798],[2.575,5.425,0.6844499],[3.425,5.425,0.6844499],[4.425,5.425,0.04321798],[3.575,5.425,0.04321798],[3.575,5.425,0.7209224],[4.425,5.425,0.7209224],[5.425,5.425,0.04321798],[4.575,5.425,0.04321798],[4.575,5.425,0.9923994],[5.425,5.425,0.9923994],[6.425,5.425,0.04321798],[5.575,5.425,0.04321798],[5.575,5.425,0.9907784],[6.425,5.425,0.9907784],[7.425,5.425,0.04321798],[6.575,5.425,0.04321798],[6.575,5.425,0.9937668],[7.425,5.425,0.9937668],[8.425,5.425,0.04321798],[7.575,5.425,0.04321798],[7.575,5.425,0.999954],[8.425,5.425,0.999954],[9.425,5.425,0.04321798],[8.575,5.425,0.04321798],[8.575,5.425,0.9999315],[9.425,5.425,0.9999315],[1.425,6.425,0.04321798],[0.575,6.425,0.04321798],[0.575,6.425,0.1997644],[1.425,6.425,0.1997644],[2.425,6.425,0.04321798],[1.575,6.425,0.04321798],[1.575,6.425,0.5454682],[2.425,6.425,0.5454682],[3.425,6.425,0.04321798],[2.575,6.425,0.04321798],[2.575,6.425,0.8883711],[3.425,6.425,0.8883711],[4.425,6.425,0.04321798],[3.575,6.425,0.04321798],[3.575,6.425,0.912325],[4.425,6.425,0.912325],[5.425,6.425,0.04321798],[4.575,6.425,0.04321798],[4.575,6.425,0.9999133],[5.425,6.425,0.9999133],[6.425,6.425,0.04321798],[5.575,6.425,0.04321798],[5.575,6.425,0.9998751],[6.425,6.425,0.9998751],[7.425,6.425,0.04321798],[6.575,6.425,0.04321798],[6.575,6.425,0.9999403],[7.425,6.425,0.9999403],[8.425,6.425,0.04321798],[7.575,6.425,0.04321798],[7.575,6.425,1],[8.425,6.425,1],[9.425,6.425,0.04321798],[8.575,6.425,0.04321798],[8.575,6.425,1],[9.425,6.425,1],[1.425,7.425,0.04321798],[0.575,7.425,0.04321798],[0.575,7.425,0.2790861],[1.425,7.425,0.2790861],[2.425,7.425,0.04321798],[1.575,7.425,0.04321798],[1.575,7.425,0.722779],[2.425,7.425,0.722779],[3.425,7.425,0.04321798],[2.575,7.425,0.04321798],[2.575,7.425,0.9743328],[3.425,7.425,0.9743328],[4.425,7.425,0.04321798],[3.575,7.425,0.04321798],[3.575,7.425,0.9828939],[4.425,7.425,0.9828939],[5.425,7.425,0.04321798],[4.575,7.425,0.04321798],[4.575,7.425,0.9999998],[5.425,7.425,0.9999998],[6.425,7.425,0.04321798],[5.575,7.425,0.04321798],[5.575,7.425,0.9999996],[6.425,7.425,0.9999996],[7.425,7.425,0.04321798],[6.575,7.425,0.04321798],[6.575,7.425,0.9999999],[7.425,7.425,0.9999999],[8.425,7.425,0.04321798],[7.575,7.425,0.04321798],[7.575,7.425,1],[8.425,7.425,1],[9.425,7.425,0.04321798],[8.575,7.425,0.04321798],[8.575,7.425,1],[9.425,7.425,1],[1.425,8.425,0.04321798],[0.575,8.425,0.04321798],[0.575,8.425,0.428676],[1.425,8.425,0.428676],[2.425,8.425,0.04321798],[1.575,8.425,0.04321798],[1.575,8.425,0.9105917],[2.425,8.425,0.9105917],[3.425,8.425,0.04321798],[2.575,8.425,0.04321798],[2.575,8.425,0.9990419],[3.425,8.425,0.9990419],[4.425,8.425,0.04321798],[3.575,8.425,0.04321798],[3.575,8.425,0.9995526],[4.425,8.425,0.9995526],[5.425,8.425,0.04321798],[4.575,8.425,0.04321798],[4.575,8.425,1],[5.425,8.425,1],[6.425,8.425,0.04321798],[5.575,8.425,0.04321798],[5.575,8.425,1],[6.425,8.425,1],[7.425,8.425,0.04321798],[6.575,8.425,0.04321798],[6.575,8.425,1],[7.425,8.425,1],[8.425,8.425,0.04321798],[7.575,8.425,0.04321798],[7.575,8.425,1],[8.425,8.425,1],[9.425,8.425,0.04321798],[8.575,8.425,0.04321798],[8.575,8.425,1],[9.425,8.425,1],[1.425,9.425,0.04321798],[0.575,9.425,0.04321798],[0.575,9.425,0.5885873],[1.425,9.425,0.5885873],[2.425,9.425,0.04321798],[1.575,9.425,0.04321798],[1.575,9.425,0.981856],[2.425,9.425,0.981856],[3.425,9.425,0.04321798],[2.575,9.425,0.04321798],[2.575,9.425,0.9999893],[3.425,9.425,0.9999893],[4.425,9.425,0.04321798],[3.575,9.425,0.04321798],[3.575,9.425,0.9999969],[4.425,9.425,0.9999969],[5.425,9.425,0.04321798],[4.575,9.425,0.04321798],[4.575,9.425,1],[5.425,9.425,1],[6.425,9.425,0.04321798],[5.575,9.425,0.04321798],[5.575,9.425,1],[6.425,9.425,1],[7.425,9.425,0.04321798],[6.575,9.425,0.04321798],[6.575,9.425,1],[7.425,9.425,1],[8.425,9.425,0.04321798],[7.575,9.425,0.04321798],[7.575,9.425,1],[8.425,9.425,1],[9.425,9.425,0.04321798],[8.575,9.425,0.04321798],[8.575,9.425,1],[9.425,9.425,1],[1.425,10.425,0.04321798],[0.575,10.425,0.04321798],[0.575,10.425,0.7135577],[1.425,10.425,0.7135577],[2.425,10.425,0.04321798],[1.575,10.425,0.04321798],[1.575,10.425,0.9967811],[2.425,10.425,0.9967811],[3.425,10.425,0.04321798],[2.575,10.425,0.04321798],[2.575,10.425,0.9999999],[3.425,10.425,0.9999999],[4.425,10.425,0.04321798],[3.575,10.425,0.04321798],[3.575,10.425,1],[4.425,10.425,1],[5.425,10.425,0.04321798],[4.575,10.425,0.04321798],[4.575,10.425,1],[5.425,10.425,1],[6.425,10.425,0.04321798],[5.575,10.425,0.04321798],[5.575,10.425,1],[6.425,10.425,1],[7.425,10.425,0.04321798],[6.575,10.425,0.04321798],[6.575,10.425,1],[7.425,10.425,1],[8.425,10.425,0.04321798],[7.575,10.425,0.04321798],[7.575,10.425,1],[8.425,10.425,1],[9.425,10.425,0.04321798],[8.575,10.425,0.04321798],[8.575,10.425,1],[9.425,10.425,1],[0.575,0.575,0.04321798],[0.575,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,0.575,0.04321798],[1.575,0.575,0.04321798],[1.575,1.425,0.04321798],[1.575,1.425,0.06596308],[1.575,0.575,0.06596308],[2.575,0.575,0.04321798],[2.575,1.425,0.04321798],[2.575,1.425,0.1021959],[2.575,0.575,0.1021959],[3.575,0.575,0.04321798],[3.575,1.425,0.04321798],[3.575,1.425,0.1075265],[3.575,0.575,0.1075265],[4.575,0.575,0.04321798],[4.575,1.425,0.04321798],[4.575,1.425,0.2316452],[4.575,0.575,0.2316452],[5.575,0.575,0.04321798],[5.575,1.425,0.04321798],[5.575,1.425,0.2258095],[5.575,0.575,0.2258095],[6.575,0.575,0.04321798],[6.575,1.425,0.04321798],[6.575,1.425,0.237572],[6.575,0.575,0.237572],[7.575,0.575,0.04321798],[7.575,1.425,0.04321798],[7.575,1.425,0.3690261],[7.575,0.575,0.3690261],[8.575,0.575,0.04321798],[8.575,1.425,0.04321798],[8.575,1.425,0.3592858],[8.575,0.575,0.3592858],[0.575,1.575,0.04321798],[0.575,2.425,0.04321798],[0.575,2.425,0.05915068],[0.575,1.575,0.05915068],[1.575,1.575,0.04321798],[1.575,2.425,0.04321798],[1.575,2.425,0.1108488],[1.575,1.575,0.1108488],[2.575,1.575,0.04321798],[2.575,2.425,0.04321798],[2.575,2.425,0.2030372],[2.575,1.575,0.2030372],[3.575,1.575,0.04321798],[3.575,2.425,0.04321798],[3.575,2.425,0.217041],[3.575,1.575,0.217041],[4.575,1.575,0.04321798],[4.575,2.425,0.04321798],[4.575,2.425,0.5257641],[4.575,1.575,0.5257641],[5.575,1.575,0.04321798],[5.575,2.425,0.04321798],[5.575,2.425,0.5127993],[5.575,1.575,0.5127993],[6.575,1.575,0.04321798],[6.575,2.425,0.04321798],[6.575,2.425,0.5387281],[6.575,1.575,0.5387281],[7.575,1.575,0.04321798],[7.575,2.425,0.04321798],[7.575,2.425,0.7713429],[7.575,1.575,0.7713429],[8.575,1.575,0.04321798],[8.575,2.425,0.04321798],[8.575,2.425,0.7576808],[8.575,1.575,0.7576808],[0.575,2.575,0.04321798],[0.575,3.425,0.04321798],[0.575,3.425,0.07725442],[0.575,2.575,0.07725442],[1.575,2.575,0.04321798],[1.575,3.425,0.04321798],[1.575,3.425,0.1672267],[1.575,2.575,0.1672267],[2.575,2.575,0.04321798],[2.575,3.425,0.04321798],[2.575,3.425,0.3312058],[2.575,2.575,0.3312058],[3.575,2.575,0.04321798],[3.575,3.425,0.04321798],[3.575,3.425,0.3553697],[3.575,2.575,0.3553697],[4.575,2.575,0.04321798],[4.575,3.425,0.04321798],[4.575,3.425,0.781293],[4.575,2.575,0.781293],[5.575,2.575,0.04321798],[5.575,3.425,0.04321798],[5.575,3.425,0.76804],[5.575,2.575,0.76804],[6.575,2.575,0.04321798],[6.575,3.425,0.04321798],[6.575,3.425,0.7941267],[6.575,2.575,0.7941267],[7.575,2.575,0.04321798],[7.575,3.425,0.04321798],[7.575,3.425,0.9553382],[7.575,2.575,0.9553382],[8.575,2.575,0.04321798],[8.575,3.425,0.04321798],[8.575,3.425,0.9493445],[8.575,2.575,0.9493445],[0.575,3.575,0.04321798],[0.575,4.425,0.04321798],[0.575,4.425,0.1025677],[0.575,3.575,0.1025677],[1.575,3.575,0.04321798],[1.575,4.425,0.04321798],[1.575,4.425,0.2501325],[1.575,3.575,0.2501325],[2.575,3.575,0.04321798],[2.575,4.425,0.04321798],[2.575,4.425,0.5039053],[2.575,3.575,0.5039053],[3.575,3.575,0.04321798],[3.575,4.425,0.04321798],[3.575,4.425,0.537781],[3.575,3.575,0.537781],[4.575,3.575,0.04321798],[4.575,4.425,0.04321798],[4.575,4.425,0.9430385],[4.575,3.575,0.9430385],[5.575,3.575,0.04321798],[5.575,4.425,0.04321798],[5.575,4.425,0.9361931],[5.575,3.575,0.9361931],[6.575,3.575,0.04321798],[6.575,4.425,0.04321798],[6.575,4.425,0.9493053],[6.575,3.575,0.9493053],[7.575,3.575,0.04321798],[7.575,4.425,0.04321798],[7.575,4.425,0.9972544],[7.575,3.575,0.9972544],[8.575,3.575,0.04321798],[8.575,4.425,0.04321798],[8.575,4.425,0.9965154],[8.575,3.575,0.9965154],[0.575,4.575,0.04321798],[0.575,5.425,0.04321798],[0.575,5.425,0.1352801],[0.575,4.575,0.1352801],[1.575,4.575,0.04321798],[1.575,5.425,0.04321798],[1.575,5.425,0.3567418],[1.575,4.575,0.3567418],[2.575,4.575,0.04321798],[2.575,5.425,0.04321798],[2.575,5.425,0.6844499],[2.575,4.575,0.6844499],[3.575,4.575,0.04321798],[3.575,5.425,0.04321798],[3.575,5.425,0.7209224],[3.575,4.575,0.7209224],[4.575,4.575,0.04321798],[4.575,5.425,0.04321798],[4.575,5.425,0.9923994],[4.575,4.575,0.9923994],[5.575,4.575,0.04321798],[5.575,5.425,0.04321798],[5.575,5.425,0.9907784],[5.575,4.575,0.9907784],[6.575,4.575,0.04321798],[6.575,5.425,0.04321798],[6.575,5.425,0.9937668],[6.575,4.575,0.9937668],[7.575,4.575,0.04321798],[7.575,5.425,0.04321798],[7.575,5.425,0.999954],[7.575,4.575,0.999954],[8.575,4.575,0.04321798],[8.575,5.425,0.04321798],[8.575,5.425,0.9999315],[8.575,4.575,0.9999315],[0.575,5.575,0.04321798],[0.575,6.425,0.04321798],[0.575,6.425,0.1997644],[0.575,5.575,0.1997644],[1.575,5.575,0.04321798],[1.575,6.425,0.04321798],[1.575,6.425,0.5454682],[1.575,5.575,0.5454682],[2.575,5.575,0.04321798],[2.575,6.425,0.04321798],[2.575,6.425,0.8883711],[2.575,5.575,0.8883711],[3.575,5.575,0.04321798],[3.575,6.425,0.04321798],[3.575,6.425,0.912325],[3.575,5.575,0.912325],[4.575,5.575,0.04321798],[4.575,6.425,0.04321798],[4.575,6.425,0.9999133],[4.575,5.575,0.9999133],[5.575,5.575,0.04321798],[5.575,6.425,0.04321798],[5.575,6.425,0.9998751],[5.575,5.575,0.9998751],[6.575,5.575,0.04321798],[6.575,6.425,0.04321798],[6.575,6.425,0.9999403],[6.575,5.575,0.9999403],[7.575,5.575,0.04321798],[7.575,6.425,0.04321798],[7.575,6.425,1],[7.575,5.575,1],[8.575,5.575,0.04321798],[8.575,6.425,0.04321798],[8.575,6.425,1],[8.575,5.575,1],[0.575,6.575,0.04321798],[0.575,7.425,0.04321798],[0.575,7.425,0.2790861],[0.575,6.575,0.2790861],[1.575,6.575,0.04321798],[1.575,7.425,0.04321798],[1.575,7.425,0.722779],[1.575,6.575,0.722779],[2.575,6.575,0.04321798],[2.575,7.425,0.04321798],[2.575,7.425,0.9743328],[2.575,6.575,0.9743328],[3.575,6.575,0.04321798],[3.575,7.425,0.04321798],[3.575,7.425,0.9828939],[3.575,6.575,0.9828939],[4.575,6.575,0.04321798],[4.575,7.425,0.04321798],[4.575,7.425,0.9999998],[4.575,6.575,0.9999998],[5.575,6.575,0.04321798],[5.575,7.425,0.04321798],[5.575,7.425,0.9999996],[5.575,6.575,0.9999996],[6.575,6.575,0.04321798],[6.575,7.425,0.04321798],[6.575,7.425,0.9999999],[6.575,6.575,0.9999999],[7.575,6.575,0.04321798],[7.575,7.425,0.04321798],[7.575,7.425,1],[7.575,6.575,1],[8.575,6.575,0.04321798],[8.575,7.425,0.04321798],[8.575,7.425,1],[8.575,6.575,1],[0.575,7.575,0.04321798],[0.575,8.425,0.04321798],[0.575,8.425,0.428676],[0.575,7.575,0.428676],[1.575,7.575,0.04321798],[1.575,8.425,0.04321798],[1.575,8.425,0.9105917],[1.575,7.575,0.9105917],[2.575,7.575,0.04321798],[2.575,8.425,0.04321798],[2.575,8.425,0.9990419],[2.575,7.575,0.9990419],[3.575,7.575,0.04321798],[3.575,8.425,0.04321798],[3.575,8.425,0.9995526],[3.575,7.575,0.9995526],[4.575,7.575,0.04321798],[4.575,8.425,0.04321798],[4.575,8.425,1],[4.575,7.575,1],[5.575,7.575,0.04321798],[5.575,8.425,0.04321798],[5.575,8.425,1],[5.575,7.575,1],[6.575,7.575,0.04321798],[6.575,8.425,0.04321798],[6.575,8.425,1],[6.575,7.575,1],[7.575,7.575,0.04321798],[7.575,8.425,0.04321798],[7.575,8.425,1],[7.575,7.575,1],[8.575,7.575,0.04321798],[8.575,8.425,0.04321798],[8.575,8.425,1],[8.575,7.575,1],[0.575,8.575,0.04321798],[0.575,9.425,0.04321798],[0.575,9.425,0.5885873],[0.575,8.575,0.5885873],[1.575,8.575,0.04321798],[1.575,9.425,0.04321798],[1.575,9.425,0.981856],[1.575,8.575,0.981856],[2.575,8.575,0.04321798],[2.575,9.425,0.04321798],[2.575,9.425,0.9999893],[2.575,8.575,0.9999893],[3.575,8.575,0.04321798],[3.575,9.425,0.04321798],[3.575,9.425,0.9999969],[3.575,8.575,0.9999969],[4.575,8.575,0.04321798],[4.575,9.425,0.04321798],[4.575,9.425,1],[4.575,8.575,1],[5.575,8.575,0.04321798],[5.575,9.425,0.04321798],[5.575,9.425,1],[5.575,8.575,1],[6.575,8.575,0.04321798],[6.575,9.425,0.04321798],[6.575,9.425,1],[6.575,8.575,1],[7.575,8.575,0.04321798],[7.575,9.425,0.04321798],[7.575,9.425,1],[7.575,8.575,1],[8.575,8.575,0.04321798],[8.575,9.425,0.04321798],[8.575,9.425,1],[8.575,8.575,1],[0.575,9.575,0.04321798],[0.575,10.425,0.04321798],[0.575,10.425,0.7135577],[0.575,9.575,0.7135577],[1.575,9.575,0.04321798],[1.575,10.425,0.04321798],[1.575,10.425,0.9967811],[1.575,9.575,0.9967811],[2.575,9.575,0.04321798],[2.575,10.425,0.04321798],[2.575,10.425,0.9999999],[2.575,9.575,0.9999999],[3.575,9.575,0.04321798],[3.575,10.425,0.04321798],[3.575,10.425,1],[3.575,9.575,1],[4.575,9.575,0.04321798],[4.575,10.425,0.04321798],[4.575,10.425,1],[4.575,9.575,1],[5.575,9.575,0.04321798],[5.575,10.425,0.04321798],[5.575,10.425,1],[5.575,9.575,1],[6.575,9.575,0.04321798],[6.575,10.425,0.04321798],[6.575,10.425,1],[6.575,9.575,1],[7.575,9.575,0.04321798],[7.575,10.425,0.04321798],[7.575,10.425,1],[7.575,9.575,1],[8.575,9.575,0.04321798],[8.575,10.425,0.04321798],[8.575,10.425,1],[8.575,9.575,1],[1.425,0.575,0.04321798],[1.425,1.425,0.04321798],[1.425,1.425,0.04321798],[1.425,0.575,0.04321798],[2.425,0.575,0.04321798],[2.425,1.425,0.04321798],[2.425,1.425,0.06596308],[2.425,0.575,0.06596308],[3.425,0.575,0.04321798],[3.425,1.425,0.04321798],[3.425,1.425,0.1021959],[3.425,0.575,0.1021959],[4.425,0.575,0.04321798],[4.425,1.425,0.04321798],[4.425,1.425,0.1075265],[4.425,0.575,0.1075265],[5.425,0.575,0.04321798],[5.425,1.425,0.04321798],[5.425,1.425,0.2316452],[5.425,0.575,0.2316452],[6.425,0.575,0.04321798],[6.425,1.425,0.04321798],[6.425,1.425,0.2258095],[6.425,0.575,0.2258095],[7.425,0.575,0.04321798],[7.425,1.425,0.04321798],[7.425,1.425,0.237572],[7.425,0.575,0.237572],[8.425,0.575,0.04321798],[8.425,1.425,0.04321798],[8.425,1.425,0.3690261],[8.425,0.575,0.3690261],[9.425,0.575,0.04321798],[9.425,1.425,0.04321798],[9.425,1.425,0.3592858],[9.425,0.575,0.3592858],[1.425,1.575,0.04321798],[1.425,2.425,0.04321798],[1.425,2.425,0.05915068],[1.425,1.575,0.05915068],[2.425,1.575,0.04321798],[2.425,2.425,0.04321798],[2.425,2.425,0.1108488],[2.425,1.575,0.1108488],[3.425,1.575,0.04321798],[3.425,2.425,0.04321798],[3.425,2.425,0.2030372],[3.425,1.575,0.2030372],[4.425,1.575,0.04321798],[4.425,2.425,0.04321798],[4.425,2.425,0.217041],[4.425,1.575,0.217041],[5.425,1.575,0.04321798],[5.425,2.425,0.04321798],[5.425,2.425,0.5257641],[5.425,1.575,0.5257641],[6.425,1.575,0.04321798],[6.425,2.425,0.04321798],[6.425,2.425,0.5127993],[6.425,1.575,0.5127993],[7.425,1.575,0.04321798],[7.425,2.425,0.04321798],[7.425,2.425,0.5387281],[7.425,1.575,0.5387281],[8.425,1.575,0.04321798],[8.425,2.425,0.04321798],[8.425,2.425,0.7713429],[8.425,1.575,0.7713429],[9.425,1.575,0.04321798],[9.425,2.425,0.04321798],[9.425,2.425,0.7576808],[9.425,1.575,0.7576808],[1.425,2.575,0.04321798],[1.425,3.425,0.04321798],[1.425,3.425,0.07725442],[1.425,2.575,0.07725442],[2.425,2.575,0.04321798],[2.425,3.425,0.04321798],[2.425,3.425,0.1672267],[2.425,2.575,0.1672267],[3.425,2.575,0.04321798],[3.425,3.425,0.04321798],[3.425,3.425,0.3312058],[3.425,2.575,0.3312058],[4.425,2.575,0.04321798],[4.425,3.425,0.04321798],[4.425,3.425,0.3553697],[4.425,2.575,0.3553697],[5.425,2.575,0.04321798],[5.425,3.425,0.04321798],[5.425,3.425,0.781293],[5.425,2.575,0.781293],[6.425,2.575,0.04321798],[6.425,3.425,0.04321798],[6.425,3.425,0.76804],[6.425,2.575,0.76804],[7.425,2.575,0.04321798],[7.425,3.425,0.04321798],[7.425,3.425,0.7941267],[7.425,2.575,0.7941267],[8.425,2.575,0.04321798],[8.425,3.425,0.04321798],[8.425,3.425,0.9553382],[8.425,2.575,0.9553382],[9.425,2.575,0.04321798],[9.425,3.425,0.04321798],[9.425,3.425,0.9493445],[9.425,2.575,0.9493445],[1.425,3.575,0.04321798],[1.425,4.425,0.04321798],[1.425,4.425,0.1025677],[1.425,3.575,0.1025677],[2.425,3.575,0.04321798],[2.425,4.425,0.04321798],[2.425,4.425,0.2501325],[2.425,3.575,0.2501325],[3.425,3.575,0.04321798],[3.425,4.425,0.04321798],[3.425,4.425,0.5039053],[3.425,3.575,0.5039053],[4.425,3.575,0.04321798],[4.425,4.425,0.04321798],[4.425,4.425,0.537781],[4.425,3.575,0.537781],[5.425,3.575,0.04321798],[5.425,4.425,0.04321798],[5.425,4.425,0.9430385],[5.425,3.575,0.9430385],[6.425,3.575,0.04321798],[6.425,4.425,0.04321798],[6.425,4.425,0.9361931],[6.425,3.575,0.9361931],[7.425,3.575,0.04321798],[7.425,4.425,0.04321798],[7.425,4.425,0.9493053],[7.425,3.575,0.9493053],[8.425,3.575,0.04321798],[8.425,4.425,0.04321798],[8.425,4.425,0.9972544],[8.425,3.575,0.9972544],[9.425,3.575,0.04321798],[9.425,4.425,0.04321798],[9.425,4.425,0.9965154],[9.425,3.575,0.9965154],[1.425,4.575,0.04321798],[1.425,5.425,0.04321798],[1.425,5.425,0.1352801],[1.425,4.575,0.1352801],[2.425,4.575,0.04321798],[2.425,5.425,0.04321798],[2.425,5.425,0.3567418],[2.425,4.575,0.3567418],[3.425,4.575,0.04321798],[3.425,5.425,0.04321798],[3.425,5.425,0.6844499],[3.425,4.575,0.6844499],[4.425,4.575,0.04321798],[4.425,5.425,0.04321798],[4.425,5.425,0.7209224],[4.425,4.575,0.7209224],[5.425,4.575,0.04321798],[5.425,5.425,0.04321798],[5.425,5.425,0.9923994],[5.425,4.575,0.9923994],[6.425,4.575,0.04321798],[6.425,5.425,0.04321798],[6.425,5.425,0.9907784],[6.425,4.575,0.9907784],[7.425,4.575,0.04321798],[7.425,5.425,0.04321798],[7.425,5.425,0.9937668],[7.425,4.575,0.9937668],[8.425,4.575,0.04321798],[8.425,5.425,0.04321798],[8.425,5.425,0.999954],[8.425,4.575,0.999954],[9.425,4.575,0.04321798],[9.425,5.425,0.04321798],[9.425,5.425,0.9999315],[9.425,4.575,0.9999315],[1.425,5.575,0.04321798],[1.425,6.425,0.04321798],[1.425,6.425,0.1997644],[1.425,5.575,0.1997644],[2.425,5.575,0.04321798],[2.425,6.425,0.04321798],[2.425,6.425,0.5454682],[2.425,5.575,0.5454682],[3.425,5.575,0.04321798],[3.425,6.425,0.04321798],[3.425,6.425,0.8883711],[3.425,5.575,0.8883711],[4.425,5.575,0.04321798],[4.425,6.425,0.04321798],[4.425,6.425,0.912325],[4.425,5.575,0.912325],[5.425,5.575,0.04321798],[5.425,6.425,0.04321798],[5.425,6.425,0.9999133],[5.425,5.575,0.9999133],[6.425,5.575,0.04321798],[6.425,6.425,0.04321798],[6.425,6.425,0.9998751],[6.425,5.575,0.9998751],[7.425,5.575,0.04321798],[7.425,6.425,0.04321798],[7.425,6.425,0.9999403],[7.425,5.575,0.9999403],[8.425,5.575,0.04321798],[8.425,6.425,0.04321798],[8.425,6.425,1],[8.425,5.575,1],[9.425,5.575,0.04321798],[9.425,6.425,0.04321798],[9.425,6.425,1],[9.425,5.575,1],[1.425,6.575,0.04321798],[1.425,7.425,0.04321798],[1.425,7.425,0.2790861],[1.425,6.575,0.2790861],[2.425,6.575,0.04321798],[2.425,7.425,0.04321798],[2.425,7.425,0.722779],[2.425,6.575,0.722779],[3.425,6.575,0.04321798],[3.425,7.425,0.04321798],[3.425,7.425,0.9743328],[3.425,6.575,0.9743328],[4.425,6.575,0.04321798],[4.425,7.425,0.04321798],[4.425,7.425,0.9828939],[4.425,6.575,0.9828939],[5.425,6.575,0.04321798],[5.425,7.425,0.04321798],[5.425,7.425,0.9999998],[5.425,6.575,0.9999998],[6.425,6.575,0.04321798],[6.425,7.425,0.04321798],[6.425,7.425,0.9999996],[6.425,6.575,0.9999996],[7.425,6.575,0.04321798],[7.425,7.425,0.04321798],[7.425,7.425,0.9999999],[7.425,6.575,0.9999999],[8.425,6.575,0.04321798],[8.425,7.425,0.04321798],[8.425,7.425,1],[8.425,6.575,1],[9.425,6.575,0.04321798],[9.425,7.425,0.04321798],[9.425,7.425,1],[9.425,6.575,1],[1.425,7.575,0.04321798],[1.425,8.425,0.04321798],[1.425,8.425,0.428676],[1.425,7.575,0.428676],[2.425,7.575,0.04321798],[2.425,8.425,0.04321798],[2.425,8.425,0.9105917],[2.425,7.575,0.9105917],[3.425,7.575,0.04321798],[3.425,8.425,0.04321798],[3.425,8.425,0.9990419],[3.425,7.575,0.9990419],[4.425,7.575,0.04321798],[4.425,8.425,0.04321798],[4.425,8.425,0.9995526],[4.425,7.575,0.9995526],[5.425,7.575,0.04321798],[5.425,8.425,0.04321798],[5.425,8.425,1],[5.425,7.575,1],[6.425,7.575,0.04321798],[6.425,8.425,0.04321798],[6.425,8.425,1],[6.425,7.575,1],[7.425,7.575,0.04321798],[7.425,8.425,0.04321798],[7.425,8.425,1],[7.425,7.575,1],[8.425,7.575,0.04321798],[8.425,8.425,0.04321798],[8.425,8.425,1],[8.425,7.575,1],[9.425,7.575,0.04321798],[9.425,8.425,0.04321798],[9.425,8.425,1],[9.425,7.575,1],[1.425,8.575,0.04321798],[1.425,9.425,0.04321798],[1.425,9.425,0.5885873],[1.425,8.575,0.5885873],[2.425,8.575,0.04321798],[2.425,9.425,0.04321798],[2.425,9.425,0.981856],[2.425,8.575,0.981856],[3.425,8.575,0.04321798],[3.425,9.425,0.04321798],[3.425,9.425,0.9999893],[3.425,8.575,0.9999893],[4.425,8.575,0.04321798],[4.425,9.425,0.04321798],[4.425,9.425,0.9999969],[4.425,8.575,0.9999969],[5.425,8.575,0.04321798],[5.425,9.425,0.04321798],[5.425,9.425,1],[5.425,8.575,1],[6.425,8.575,0.04321798],[6.425,9.425,0.04321798],[6.425,9.425,1],[6.425,8.575,1],[7.425,8.575,0.04321798],[7.425,9.425,0.04321798],[7.425,9.425,1],[7.425,8.575,1],[8.425,8.575,0.04321798],[8.425,9.425,0.04321798],[8.425,9.425,1],[8.425,8.575,1],[9.425,8.575,0.04321798],[9.425,9.425,0.04321798],[9.425,9.425,1],[9.425,8.575,1],[1.425,9.575,0.04321798],[1.425,10.425,0.04321798],[1.425,10.425,0.7135577],[1.425,9.575,0.7135577],[2.425,9.575,0.04321798],[2.425,10.425,0.04321798],[2.425,10.425,0.9967811],[2.425,9.575,0.9967811],[3.425,9.575,0.04321798],[3.425,10.425,0.04321798],[3.425,10.425,0.9999999],[3.425,9.575,0.9999999],[4.425,9.575,0.04321798],[4.425,10.425,0.04321798],[4.425,10.425,1],[4.425,9.575,1],[5.425,9.575,0.04321798],[5.425,10.425,0.04321798],[5.425,10.425,1],[5.425,9.575,1],[6.425,9.575,0.04321798],[6.425,10.425,0.04321798],[6.425,10.425,1],[6.425,9.575,1],[7.425,9.575,0.04321798],[7.425,10.425,0.04321798],[7.425,10.425,1],[7.425,9.575,1],[8.425,9.575,0.04321798],[8.425,10.425,0.04321798],[8.425,10.425,1],[8.425,9.575,1],[9.425,9.575,0.04321798],[9.425,10.425,0.04321798],[9.425,10.425,1],[9.425,9.575,1],[0.575,0.575,0.04321798],[1.425,0.575,0.04321798],[1.425,1.425,0.04321798],[0.575,1.425,0.04321798],[1.575,0.575,0.06596308],[2.425,0.575,0.06596308],[2.425,1.425,0.06596308],[1.575,1.425,0.06596308],[2.575,0.575,0.1021959],[3.425,0.575,0.1021959],[3.425,1.425,0.1021959],[2.575,1.425,0.1021959],[3.575,0.575,0.1075265],[4.425,0.575,0.1075265],[4.425,1.425,0.1075265],[3.575,1.425,0.1075265],[4.575,0.575,0.2316452],[5.425,0.575,0.2316452],[5.425,1.425,0.2316452],[4.575,1.425,0.2316452],[5.575,0.575,0.2258095],[6.425,0.575,0.2258095],[6.425,1.425,0.2258095],[5.575,1.425,0.2258095],[6.575,0.575,0.237572],[7.425,0.575,0.237572],[7.425,1.425,0.237572],[6.575,1.425,0.237572],[7.575,0.575,0.3690261],[8.425,0.575,0.3690261],[8.425,1.425,0.3690261],[7.575,1.425,0.3690261],[8.575,0.575,0.3592858],[9.425,0.575,0.3592858],[9.425,1.425,0.3592858],[8.575,1.425,0.3592858],[0.575,1.575,0.05915068],[1.425,1.575,0.05915068],[1.425,2.425,0.05915068],[0.575,2.425,0.05915068],[1.575,1.575,0.1108488],[2.425,1.575,0.1108488],[2.425,2.425,0.1108488],[1.575,2.425,0.1108488],[2.575,1.575,0.2030372],[3.425,1.575,0.2030372],[3.425,2.425,0.2030372],[2.575,2.425,0.2030372],[3.575,1.575,0.217041],[4.425,1.575,0.217041],[4.425,2.425,0.217041],[3.575,2.425,0.217041],[4.575,1.575,0.5257641],[5.425,1.575,0.5257641],[5.425,2.425,0.5257641],[4.575,2.425,0.5257641],[5.575,1.575,0.5127993],[6.425,1.575,0.5127993],[6.425,2.425,0.5127993],[5.575,2.425,0.5127993],[6.575,1.575,0.5387281],[7.425,1.575,0.5387281],[7.425,2.425,0.5387281],[6.575,2.425,0.5387281],[7.575,1.575,0.7713429],[8.425,1.575,0.7713429],[8.425,2.425,0.7713429],[7.575,2.425,0.7713429],[8.575,1.575,0.7576808],[9.425,1.575,0.7576808],[9.425,2.425,0.7576808],[8.575,2.425,0.7576808],[0.575,2.575,0.07725442],[1.425,2.575,0.07725442],[1.425,3.425,0.07725442],[0.575,3.425,0.07725442],[1.575,2.575,0.1672267],[2.425,2.575,0.1672267],[2.425,3.425,0.1672267],[1.575,3.425,0.1672267],[2.575,2.575,0.3312058],[3.425,2.575,0.3312058],[3.425,3.425,0.3312058],[2.575,3.425,0.3312058],[3.575,2.575,0.3553697],[4.425,2.575,0.3553697],[4.425,3.425,0.3553697],[3.575,3.425,0.3553697],[4.575,2.575,0.781293],[5.425,2.575,0.781293],[5.425,3.425,0.781293],[4.575,3.425,0.781293],[5.575,2.575,0.76804],[6.425,2.575,0.76804],[6.425,3.425,0.76804],[5.575,3.425,0.76804],[6.575,2.575,0.7941267],[7.425,2.575,0.7941267],[7.425,3.425,0.7941267],[6.575,3.425,0.7941267],[7.575,2.575,0.9553382],[8.425,2.575,0.9553382],[8.425,3.425,0.9553382],[7.575,3.425,0.9553382],[8.575,2.575,0.9493445],[9.425,2.575,0.9493445],[9.425,3.425,0.9493445],[8.575,3.425,0.9493445],[0.575,3.575,0.1025677],[1.425,3.575,0.1025677],[1.425,4.425,0.1025677],[0.575,4.425,0.1025677],[1.575,3.575,0.2501325],[2.425,3.575,0.2501325],[2.425,4.425,0.2501325],[1.575,4.425,0.2501325],[2.575,3.575,0.5039053],[3.425,3.575,0.5039053],[3.425,4.425,0.5039053],[2.575,4.425,0.5039053],[3.575,3.575,0.537781],[4.425,3.575,0.537781],[4.425,4.425,0.537781],[3.575,4.425,0.537781],[4.575,3.575,0.9430385],[5.425,3.575,0.9430385],[5.425,4.425,0.9430385],[4.575,4.425,0.9430385],[5.575,3.575,0.9361931],[6.425,3.575,0.9361931],[6.425,4.425,0.9361931],[5.575,4.425,0.9361931],[6.575,3.575,0.9493053],[7.425,3.575,0.9493053],[7.425,4.425,0.9493053],[6.575,4.425,0.9493053],[7.575,3.575,0.9972544],[8.425,3.575,0.9972544],[8.425,4.425,0.9972544],[7.575,4.425,0.9972544],[8.575,3.575,0.9965154],[9.425,3.575,0.9965154],[9.425,4.425,0.9965154],[8.575,4.425,0.9965154],[0.575,4.575,0.1352801],[1.425,4.575,0.1352801],[1.425,5.425,0.1352801],[0.575,5.425,0.1352801],[1.575,4.575,0.3567418],[2.425,4.575,0.3567418],[2.425,5.425,0.3567418],[1.575,5.425,0.3567418],[2.575,4.575,0.6844499],[3.425,4.575,0.6844499],[3.425,5.425,0.6844499],[2.575,5.425,0.6844499],[3.575,4.575,0.7209224],[4.425,4.575,0.7209224],[4.425,5.425,0.7209224],[3.575,5.425,0.7209224],[4.575,4.575,0.9923994],[5.425,4.575,0.9923994],[5.425,5.425,0.9923994],[4.575,5.425,0.9923994],[5.575,4.575,0.9907784],[6.425,4.575,0.9907784],[6.425,5.425,0.9907784],[5.575,5.425,0.9907784],[6.575,4.575,0.9937668],[7.425,4.575,0.9937668],[7.425,5.425,0.9937668],[6.575,5.425,0.9937668],[7.575,4.575,0.999954],[8.425,4.575,0.999954],[8.425,5.425,0.999954],[7.575,5.425,0.999954],[8.575,4.575,0.9999315],[9.425,4.575,0.9999315],[9.425,5.425,0.9999315],[8.575,5.425,0.9999315],[0.575,5.575,0.1997644],[1.425,5.575,0.1997644],[1.425,6.425,0.1997644],[0.575,6.425,0.1997644],[1.575,5.575,0.5454682],[2.425,5.575,0.5454682],[2.425,6.425,0.5454682],[1.575,6.425,0.5454682],[2.575,5.575,0.8883711],[3.425,5.575,0.8883711],[3.425,6.425,0.8883711],[2.575,6.425,0.8883711],[3.575,5.575,0.912325],[4.425,5.575,0.912325],[4.425,6.425,0.912325],[3.575,6.425,0.912325],[4.575,5.575,0.9999133],[5.425,5.575,0.9999133],[5.425,6.425,0.9999133],[4.575,6.425,0.9999133],[5.575,5.575,0.9998751],[6.425,5.575,0.9998751],[6.425,6.425,0.9998751],[5.575,6.425,0.9998751],[6.575,5.575,0.9999403],[7.425,5.575,0.9999403],[7.425,6.425,0.9999403],[6.575,6.425,0.9999403],[7.575,5.575,1],[8.425,5.575,1],[8.425,6.425,1],[7.575,6.425,1],[8.575,5.575,1],[9.425,5.575,1],[9.425,6.425,1],[8.575,6.425,1],[0.575,6.575,0.2790861],[1.425,6.575,0.2790861],[1.425,7.425,0.2790861],[0.575,7.425,0.2790861],[1.575,6.575,0.722779],[2.425,6.575,0.722779],[2.425,7.425,0.722779],[1.575,7.425,0.722779],[2.575,6.575,0.9743328],[3.425,6.575,0.9743328],[3.425,7.425,0.9743328],[2.575,7.425,0.9743328],[3.575,6.575,0.9828939],[4.425,6.575,0.9828939],[4.425,7.425,0.9828939],[3.575,7.425,0.9828939],[4.575,6.575,0.9999998],[5.425,6.575,0.9999998],[5.425,7.425,0.9999998],[4.575,7.425,0.9999998],[5.575,6.575,0.9999996],[6.425,6.575,0.9999996],[6.425,7.425,0.9999996],[5.575,7.425,0.9999996],[6.575,6.575,0.9999999],[7.425,6.575,0.9999999],[7.425,7.425,0.9999999],[6.575,7.425,0.9999999],[7.575,6.575,1],[8.425,6.575,1],[8.425,7.425,1],[7.575,7.425,1],[8.575,6.575,1],[9.425,6.575,1],[9.425,7.425,1],[8.575,7.425,1],[0.575,7.575,0.428676],[1.425,7.575,0.428676],[1.425,8.425,0.428676],[0.575,8.425,0.428676],[1.575,7.575,0.9105917],[2.425,7.575,0.9105917],[2.425,8.425,0.9105917],[1.575,8.425,0.9105917],[2.575,7.575,0.9990419],[3.425,7.575,0.9990419],[3.425,8.425,0.9990419],[2.575,8.425,0.9990419],[3.575,7.575,0.9995526],[4.425,7.575,0.9995526],[4.425,8.425,0.9995526],[3.575,8.425,0.9995526],[4.575,7.575,1],[5.425,7.575,1],[5.425,8.425,1],[4.575,8.425,1],[5.575,7.575,1],[6.425,7.575,1],[6.425,8.425,1],[5.575,8.425,1],[6.575,7.575,1],[7.425,7.575,1],[7.425,8.425,1],[6.575,8.425,1],[7.575,7.575,1],[8.425,7.575,1],[8.425,8.425,1],[7.575,8.425,1],[8.575,7.575,1],[9.425,7.575,1],[9.425,8.425,1],[8.575,8.425,1],[0.575,8.575,0.5885873],[1.425,8.575,0.5885873],[1.425,9.425,0.5885873],[0.575,9.425,0.5885873],[1.575,8.575,0.981856],[2.425,8.575,0.981856],[2.425,9.425,0.981856],[1.575,9.425,0.981856],[2.575,8.575,0.9999893],[3.425,8.575,0.9999893],[3.425,9.425,0.9999893],[2.575,9.425,0.9999893],[3.575,8.575,0.9999969],[4.425,8.575,0.9999969],[4.425,9.425,0.9999969],[3.575,9.425,0.9999969],[4.575,8.575,1],[5.425,8.575,1],[5.425,9.425,1],[4.575,9.425,1],[5.575,8.575,1],[6.425,8.575,1],[6.425,9.425,1],[5.575,9.425,1],[6.575,8.575,1],[7.425,8.575,1],[7.425,9.425,1],[6.575,9.425,1],[7.575,8.575,1],[8.425,8.575,1],[8.425,9.425,1],[7.575,9.425,1],[8.575,8.575,1],[9.425,8.575,1],[9.425,9.425,1],[8.575,9.425,1],[0.575,9.575,0.7135577],[1.425,9.575,0.7135577],[1.425,10.425,0.7135577],[0.575,10.425,0.7135577],[1.575,9.575,0.9967811],[2.425,9.575,0.9967811],[2.425,10.425,0.9967811],[1.575,10.425,0.9967811],[2.575,9.575,0.9999999],[3.425,9.575,0.9999999],[3.425,10.425,0.9999999],[2.575,10.425,0.9999999],[3.575,9.575,1],[4.425,9.575,1],[4.425,10.425,1],[3.575,10.425,1],[4.575,9.575,1],[5.425,9.575,1],[5.425,10.425,1],[4.575,10.425,1],[5.575,9.575,1],[6.425,9.575,1],[6.425,10.425,1],[5.575,10.425,1],[6.575,9.575,1],[7.425,9.575,1],[7.425,10.425,1],[6.575,10.425,1],[7.575,9.575,1],[8.425,9.575,1],[8.425,10.425,1],[7.575,10.425,1],[8.575,9.575,1],[9.425,9.575,1],[9.425,10.425,1],[8.575,10.425,1]],"colors":[[0.2156863,0,0,1],[0.2156863,0,0,1],[0.2156863,0,0,1],[0.2156863,0,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.1568628,0.1882353,0,1],[0.1568628,0.1882353,0,1],[0.1568628,0.1882353,0,1],[0.1568628,0.1882353,0,1],[0.1686275,0.1921569,0,1],[0.1686275,0.1921569,0,1],[0.1686275,0.1921569,0,1],[0.1686275,0.1921569,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.1490196,0.1843137,0,1],[0.1490196,0.1843137,0,1],[0.1490196,0.1843137,0,1],[0.1490196,0.1843137,0,1],[0.1647059,0.1921569,0,1],[0.1647059,0.1921569,0,1],[0.1647059,0.1921569,0,1],[0.1647059,0.1921569,0,1],[0.1411765,0.1803922,0,1],[0.1411765,0.1803922,0,1],[0.1411765,0.1803922,0,1],[0.1411765,0.1803922,0,1],[0.03137255,0.1333333,0,1],[0.03137255,0.1333333,0,1],[0.03137255,0.1333333,0,1],[0.03137255,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.04705882,0.1372549,0,1],[0.04705882,0.1372549,0,1],[0.04705882,0.1372549,0,1],[0.04705882,0.1372549,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.2156863,0,1],[0.2156863,0.2156863,0,1],[0.2156863,0.2156863,0,1],[0.2156863,0.2156863,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.07843138,0.1529412,0,1],[0.07843138,0.1529412,0,1],[0.07843138,0.1529412,0,1],[0.07843138,0.1529412,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.1019608,0,1],[0.2156863,0.1019608,0,1],[0.2156863,0.1019608,0,1],[0.2156863,0.1019608,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.01960784,0.1254902,0,1],[0.01960784,0.1254902,0,1],[0.01960784,0.1254902,0,1],[0.01960784,0.1254902,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.1529412,0,1],[0.2156863,0.1529412,0,1],[0.2156863,0.1529412,0,1],[0.2156863,0.1529412,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0.2156863,0.1921569,0,1],[0.2156863,0.1921569,0,1],[0.2156863,0.1921569,0,1],[0.2156863,0.1921569,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0.1960784,0.2039216,0,1],[0.1960784,0.2039216,0,1],[0.1960784,0.2039216,0,1],[0.1960784,0.2039216,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0.8823529,0,0,1],[0.8823529,0,0,1],[0.8823529,0,0,1],[0.8823529,0,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.6392157,0.7686275,0,1],[0.6392157,0.7686275,0,1],[0.6392157,0.7686275,0,1],[0.6392157,0.7686275,0,1],[0.6901961,0.7960784,0,1],[0.6901961,0.7960784,0,1],[0.6901961,0.7960784,0,1],[0.6901961,0.7960784,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.6117647,0.7568628,0,1],[0.6117647,0.7568628,0,1],[0.6117647,0.7568628,0,1],[0.6117647,0.7568628,0,1],[0.6666667,0.7803922,0,1],[0.6666667,0.7803922,0,1],[0.6666667,0.7803922,0,1],[0.6666667,0.7803922,0,1],[0.5843138,0.7450981,0,1],[0.5843138,0.7450981,0,1],[0.5843138,0.7450981,0,1],[0.5843138,0.7450981,0,1],[0.1294118,0.5372549,0,1],[0.1294118,0.5372549,0,1],[0.1294118,0.5372549,0,1],[0.1294118,0.5372549,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1882353,0.5647059,0,1],[0.1882353,0.5647059,0,1],[0.1882353,0.5647059,0,1],[0.1882353,0.5647059,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8784314,0.8784314,0,1],[0.8784314,0.8784314,0,1],[0.8784314,0.8784314,0,1],[0.8784314,0.8784314,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.3176471,0.6235294,0,1],[0.3176471,0.6235294,0,1],[0.3176471,0.6235294,0,1],[0.3176471,0.6235294,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.4156863,0,1],[0.8823529,0.4156863,0,1],[0.8823529,0.4156863,0,1],[0.8823529,0.4156863,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.07843138,0.5137255,0,1],[0.07843138,0.5137255,0,1],[0.07843138,0.5137255,0,1],[0.07843138,0.5137255,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.6235294,0,1],[0.8823529,0.6235294,0,1],[0.8823529,0.6235294,0,1],[0.8823529,0.6235294,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0.8823529,0.7843137,0,1],[0.8823529,0.7843137,0,1],[0.8823529,0.7843137,0,1],[0.8823529,0.7843137,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0.8,0.8431373,0,1],[0.8,0.8431373,0,1],[0.8,0.8431373,0,1],[0.8,0.8431373,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0.8823529,0,0,1],[0.8823529,0,0,1],[0.8823529,0,0,1],[0.8823529,0,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.03529412,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3294118,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.3098039,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.345098,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5686275,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.01568628,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.2941177,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7176471,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7098039,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.6392157,0.7686275,0,1],[0.6392157,0.7686275,0,1],[0.6392157,0.7686275,0,1],[0.6392157,0.7686275,0,1],[0.6901961,0.7960784,0,1],[0.6901961,0.7960784,0,1],[0.6901961,0.7960784,0,1],[0.6901961,0.7960784,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.05098039,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.2078431,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5019608,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.6117647,0.7568628,0,1],[0.6117647,0.7568628,0,1],[0.6117647,0.7568628,0,1],[0.6117647,0.7568628,0,1],[0.6666667,0.7803922,0,1],[0.6666667,0.7803922,0,1],[0.6666667,0.7803922,0,1],[0.6666667,0.7803922,0,1],[0.5843138,0.7450981,0,1],[0.5843138,0.7450981,0,1],[0.5843138,0.7450981,0,1],[0.5843138,0.7450981,0,1],[0.1294118,0.5372549,0,1],[0.1294118,0.5372549,0,1],[0.1294118,0.5372549,0,1],[0.1294118,0.5372549,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.1019608,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.3647059,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7019608,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1882353,0.5647059,0,1],[0.1882353,0.5647059,0,1],[0.1882353,0.5647059,0,1],[0.1882353,0.5647059,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.1607843,0.5529412,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.1568628,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8823529,0.5529412,0,1],[0.8784314,0.8784314,0,1],[0.8784314,0.8784314,0,1],[0.8784314,0.8784314,0,1],[0.8784314,0.8784314,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.2745098,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.8823529,0.7411765,0,1],[0.3176471,0.6235294,0,1],[0.3176471,0.6235294,0,1],[0.3176471,0.6235294,0,1],[0.3176471,0.6235294,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.4156863,0,1],[0.8823529,0.4156863,0,1],[0.8823529,0.4156863,0,1],[0.8823529,0.4156863,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.772549,0.827451,0,1],[0.07843138,0.5137255,0,1],[0.07843138,0.5137255,0,1],[0.07843138,0.5137255,0,1],[0.07843138,0.5137255,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.8823529,0.6235294,0,1],[0.8823529,0.6235294,0,1],[0.8823529,0.6235294,0,1],[0.8823529,0.6235294,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.2666667,0.6,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0.8823529,0.7843137,0,1],[0.8823529,0.7843137,0,1],[0.8823529,0.7843137,0,1],[0.8823529,0.7843137,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.05098039,0.5058824,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0.8,0.8431373,0,1],[0.8,0.8431373,0,1],[0.8,0.8431373,0,1],[0.8,0.8431373,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0.02352941,0.4901961,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0,0.4784314,0,1],[0.2156863,0,0,1],[0.2156863,0,0,1],[0.2156863,0,0,1],[0.2156863,0,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.007843138,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07843138,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.07450981,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.08627451,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1372549,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.003921569,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.07058824,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.1764706,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.1568628,0.1882353,0,1],[0.1568628,0.1882353,0,1],[0.1568628,0.1882353,0,1],[0.1568628,0.1882353,0,1],[0.1686275,0.1921569,0,1],[0.1686275,0.1921569,0,1],[0.1686275,0.1921569,0,1],[0.1686275,0.1921569,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.01176471,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.05098039,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1215686,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.1490196,0.1843137,0,1],[0.1490196,0.1843137,0,1],[0.1490196,0.1843137,0,1],[0.1490196,0.1843137,0,1],[0.1647059,0.1921569,0,1],[0.1647059,0.1921569,0,1],[0.1647059,0.1921569,0,1],[0.1647059,0.1921569,0,1],[0.1411765,0.1803922,0,1],[0.1411765,0.1803922,0,1],[0.1411765,0.1803922,0,1],[0.1411765,0.1803922,0,1],[0.03137255,0.1333333,0,1],[0.03137255,0.1333333,0,1],[0.03137255,0.1333333,0,1],[0.03137255,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.02352941,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.09019608,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.172549,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.04705882,0.1372549,0,1],[0.04705882,0.1372549,0,1],[0.04705882,0.1372549,0,1],[0.04705882,0.1372549,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.03921569,0.1333333,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.03921569,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.1333333,0,1],[0.2156863,0.2156863,0,1],[0.2156863,0.2156863,0,1],[0.2156863,0.2156863,0,1],[0.2156863,0.2156863,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.06666667,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.2156863,0.1803922,0,1],[0.07843138,0.1529412,0,1],[0.07843138,0.1529412,0,1],[0.07843138,0.1529412,0,1],[0.07843138,0.1529412,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.1019608,0,1],[0.2156863,0.1019608,0,1],[0.2156863,0.1019608,0,1],[0.2156863,0.1019608,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.1882353,0.2039216,0,1],[0.01960784,0.1254902,0,1],[0.01960784,0.1254902,0,1],[0.01960784,0.1254902,0,1],[0.01960784,0.1254902,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.2156863,0.1529412,0,1],[0.2156863,0.1529412,0,1],[0.2156863,0.1529412,0,1],[0.2156863,0.1529412,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.06666667,0.145098,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0.2156863,0.1921569,0,1],[0.2156863,0.1921569,0,1],[0.2156863,0.1921569,0,1],[0.2156863,0.1921569,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.01176471,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0.1960784,0.2039216,0,1],[0.1960784,0.2039216,0,1],[0.1960784,0.2039216,0,1],[0.1960784,0.2039216,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0.007843138,0.1215686,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0,0.1176471,0,1],[0.572549,0,0,1],[0.572549,0,0,1],[0.572549,0,0,1],[0.572549,0,0,1],[0.572549,0.02352941,0,1],[0.572549,0.02352941,0,1],[0.572549,0.02352941,0,1],[0.572549,0.02352941,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.2156863,0,1],[0.572549,0.2156863,0,1],[0.572549,0.2156863,0,1],[0.572549,0.2156863,0,1],[0.572549,0.2039216,0,1],[0.572549,0.2039216,0,1],[0.572549,0.2039216,0,1],[0.572549,0.2039216,0,1],[0.572549,0.2235294,0,1],[0.572549,0.2235294,0,1],[0.572549,0.2235294,0,1],[0.572549,0.2235294,0,1],[0.572549,0.372549,0,1],[0.572549,0.372549,0,1],[0.572549,0.372549,0,1],[0.572549,0.372549,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.01176471,0,1],[0.572549,0.01176471,0,1],[0.572549,0.01176471,0,1],[0.572549,0.01176471,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.1803922,0,1],[0.572549,0.1803922,0,1],[0.572549,0.1803922,0,1],[0.572549,0.1803922,0,1],[0.572549,0.1921569,0,1],[0.572549,0.1921569,0,1],[0.572549,0.1921569,0,1],[0.572549,0.1921569,0,1],[0.572549,0.4666667,0,1],[0.572549,0.4666667,0,1],[0.572549,0.4666667,0,1],[0.572549,0.4666667,0,1],[0.572549,0.4627451,0,1],[0.572549,0.4627451,0,1],[0.572549,0.4627451,0,1],[0.572549,0.4627451,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.4156863,0.5019608,0,1],[0.4156863,0.5019608,0,1],[0.4156863,0.5019608,0,1],[0.4156863,0.5019608,0,1],[0.4509804,0.5176471,0,1],[0.4509804,0.5176471,0,1],[0.4509804,0.5176471,0,1],[0.4509804,0.5176471,0,1],[0.572549,0.03529412,0,1],[0.572549,0.03529412,0,1],[0.572549,0.03529412,0,1],[0.572549,0.03529412,0,1],[0.572549,0.1333333,0,1],[0.572549,0.1333333,0,1],[0.572549,0.1333333,0,1],[0.572549,0.1333333,0,1],[0.572549,0.3254902,0,1],[0.572549,0.3254902,0,1],[0.572549,0.3254902,0,1],[0.572549,0.3254902,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.4,0.4941176,0,1],[0.4,0.4941176,0,1],[0.4,0.4941176,0,1],[0.4,0.4941176,0,1],[0.4352941,0.509804,0,1],[0.4352941,0.509804,0,1],[0.4352941,0.509804,0,1],[0.4352941,0.509804,0,1],[0.3803922,0.4862745,0,1],[0.3803922,0.4862745,0,1],[0.3803922,0.4862745,0,1],[0.3803922,0.4862745,0,1],[0.08627451,0.3529412,0,1],[0.08627451,0.3529412,0,1],[0.08627451,0.3529412,0,1],[0.08627451,0.3529412,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.06666667,0,1],[0.572549,0.2352941,0,1],[0.572549,0.2352941,0,1],[0.572549,0.2352941,0,1],[0.572549,0.2352941,0,1],[0.572549,0.4588235,0,1],[0.572549,0.4588235,0,1],[0.572549,0.4588235,0,1],[0.572549,0.4588235,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.1215686,0.3686275,0,1],[0.1215686,0.3686275,0,1],[0.1215686,0.3686275,0,1],[0.1215686,0.3686275,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.1019608,0.3607843,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.572549,0.1019608,0,1],[0.572549,0.1019608,0,1],[0.572549,0.1019608,0,1],[0.572549,0.1019608,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.3607843,0,1],[0.572549,0.572549,0,1],[0.572549,0.572549,0,1],[0.572549,0.572549,0,1],[0.572549,0.572549,0,1],[0.5058824,0.5411765,0,1],[0.5058824,0.5411765,0,1],[0.5058824,0.5411765,0,1],[0.5058824,0.5411765,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.572549,0.1803922,0,1],[0.572549,0.1803922,0,1],[0.572549,0.1803922,0,1],[0.572549,0.1803922,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.572549,0.4823529,0,1],[0.2078431,0.4078431,0,1],[0.2078431,0.4078431,0,1],[0.2078431,0.4078431,0,1],[0.2078431,0.4078431,0,1],[0.172549,0.3921569,0,1],[0.172549,0.3921569,0,1],[0.172549,0.3921569,0,1],[0.172549,0.3921569,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.572549,0.2705882,0,1],[0.572549,0.2705882,0,1],[0.572549,0.2705882,0,1],[0.572549,0.2705882,0,1],[0.5058824,0.5411765,0,1],[0.5058824,0.5411765,0,1],[0.5058824,0.5411765,0,1],[0.5058824,0.5411765,0,1],[0.05098039,0.3372549,0,1],[0.05098039,0.3372549,0,1],[0.05098039,0.3372549,0,1],[0.05098039,0.3372549,0,1],[0.03529412,0.3294118,0,1],[0.03529412,0.3294118,0,1],[0.03529412,0.3294118,0,1],[0.03529412,0.3294118,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.572549,0.4078431,0,1],[0.572549,0.4078431,0,1],[0.572549,0.4078431,0,1],[0.572549,0.4078431,0,1],[0.172549,0.3921569,0,1],[0.172549,0.3921569,0,1],[0.172549,0.3921569,0,1],[0.172549,0.3921569,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0.572549,0.509804,0,1],[0.572549,0.509804,0,1],[0.572549,0.509804,0,1],[0.572549,0.509804,0,1],[0.03529412,0.3294118,0,1],[0.03529412,0.3294118,0,1],[0.03529412,0.3294118,0,1],[0.03529412,0.3294118,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0.5215687,0.5490196,0,1],[0.5215687,0.5490196,0,1],[0.5215687,0.5490196,0,1],[0.5215687,0.5490196,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0.01568628,0.3215686,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1],[0,0.3137255,0,1]],"centers":[[1,0.575,0.04321798],[2,0.575,0.05459053],[3,0.575,0.07270695],[4,0.575,0.07537223],[5,0.575,0.1374316],[6,0.575,0.1345138],[7,0.575,0.140395],[8,0.575,0.206122],[9,0.575,0.2012519],[1,1.575,0.05118433],[2,1.575,0.0770334],[3,1.575,0.1231276],[4,1.575,0.1301295],[5,1.575,0.2844911],[6,1.575,0.2780086],[7,1.575,0.2909731],[8,1.575,0.4072804],[9,1.575,0.4004494],[1,2.575,0.0602362],[2,2.575,0.1052224],[3,2.575,0.1872119],[4,2.575,0.1992938],[5,2.575,0.4122555],[6,2.575,0.405629],[7,2.575,0.4186724],[8,2.575,0.4992781],[9,2.575,0.4962812],[1,3.575,0.07289286],[2,3.575,0.1466753],[3,3.575,0.2735617],[4,3.575,0.2904995],[5,3.575,0.4931282],[6,3.575,0.4897055],[7,3.575,0.4962616],[8,3.575,0.5202362],[9,3.575,0.5198667],[1,4.575,0.08924903],[2,4.575,0.1999799],[3,4.575,0.363834],[4,4.575,0.3820702],[5,4.575,0.5178087],[6,4.575,0.5169982],[7,4.575,0.5184923],[8,4.575,0.5215859],[9,4.575,0.5215747],[1,5.575,0.1214912],[2,5.575,0.2943431],[3,5.575,0.4657945],[4,5.575,0.4777715],[5,5.575,0.5215657],[6,5.575,0.5215465],[7,5.575,0.5215791],[8,5.575,0.5216089],[9,5.575,0.5216089],[1,6.575,0.161152],[2,6.575,0.3829985],[3,6.575,0.5087754],[4,6.575,0.5130559],[5,6.575,0.5216089],[6,6.575,0.5216088],[7,6.575,0.5216089],[8,6.575,0.5216089],[9,6.575,0.5216089],[1,7.575,0.235947],[2,7.575,0.4769049],[3,7.575,0.5211299],[4,7.575,0.5213853],[5,7.575,0.5216089],[6,7.575,0.5216089],[7,7.575,0.5216089],[8,7.575,0.5216089],[9,7.575,0.5216089],[1,8.575,0.3159027],[2,8.575,0.512537],[3,8.575,0.5216036],[4,8.575,0.5216074],[5,8.575,0.5216089],[6,8.575,0.5216089],[7,8.575,0.5216089],[8,8.575,0.5216089],[9,8.575,0.5216089],[1,9.575,0.3783879],[2,9.575,0.5199995],[3,9.575,0.5216089],[4,9.575,0.5216089],[5,9.575,0.5216089],[6,9.575,0.5216089],[7,9.575,0.5216089],[8,9.575,0.5216089],[9,9.575,0.5216089],[1,1.425,0.04321798],[2,1.425,0.05459053],[3,1.425,0.07270695],[4,1.425,0.07537223],[5,1.425,0.1374316],[6,1.425,0.1345138],[7,1.425,0.140395],[8,1.425,0.206122],[9,1.425,0.2012519],[1,2.425,0.05118433],[2,2.425,0.0770334],[3,2.425,0.1231276],[4,2.425,0.1301295],[5,2.425,0.2844911],[6,2.425,0.2780086],[7,2.425,0.2909731],[8,2.425,0.4072804],[9,2.425,0.4004494],[1,3.425,0.0602362],[2,3.425,0.1052224],[3,3.425,0.1872119],[4,3.425,0.1992938],[5,3.425,0.4122555],[6,3.425,0.405629],[7,3.425,0.4186724],[8,3.425,0.4992781],[9,3.425,0.4962812],[1,4.425,0.07289286],[2,4.425,0.1466753],[3,4.425,0.2735617],[4,4.425,0.2904995],[5,4.425,0.4931282],[6,4.425,0.4897055],[7,4.425,0.4962616],[8,4.425,0.5202362],[9,4.425,0.5198667],[1,5.425,0.08924903],[2,5.425,0.1999799],[3,5.425,0.363834],[4,5.425,0.3820702],[5,5.425,0.5178087],[6,5.425,0.5169982],[7,5.425,0.5184923],[8,5.425,0.5215859],[9,5.425,0.5215747],[1,6.425,0.1214912],[2,6.425,0.2943431],[3,6.425,0.4657945],[4,6.425,0.4777715],[5,6.425,0.5215657],[6,6.425,0.5215465],[7,6.425,0.5215791],[8,6.425,0.5216089],[9,6.425,0.5216089],[1,7.425,0.161152],[2,7.425,0.3829985],[3,7.425,0.5087754],[4,7.425,0.5130559],[5,7.425,0.5216089],[6,7.425,0.5216088],[7,7.425,0.5216089],[8,7.425,0.5216089],[9,7.425,0.5216089],[1,8.425,0.235947],[2,8.425,0.4769049],[3,8.425,0.5211299],[4,8.425,0.5213853],[5,8.425,0.5216089],[6,8.425,0.5216089],[7,8.425,0.5216089],[8,8.425,0.5216089],[9,8.425,0.5216089],[1,9.425,0.3159027],[2,9.425,0.512537],[3,9.425,0.5216036],[4,9.425,0.5216074],[5,9.425,0.5216089],[6,9.425,0.5216089],[7,9.425,0.5216089],[8,9.425,0.5216089],[9,9.425,0.5216089],[1,10.425,0.3783879],[2,10.425,0.5199995],[3,10.425,0.5216089],[4,10.425,0.5216089],[5,10.425,0.5216089],[6,10.425,0.5216089],[7,10.425,0.5216089],[8,10.425,0.5216089],[9,10.425,0.5216089],[0.575,1,0.04321798],[1.575,1,0.05459053],[2.575,1,0.07270695],[3.575,1,0.07537223],[4.575,1,0.1374316],[5.575,1,0.1345138],[6.575,1,0.140395],[7.575,1,0.206122],[8.575,1,0.2012519],[0.575,2,0.05118433],[1.575,2,0.0770334],[2.575,2,0.1231276],[3.575,2,0.1301295],[4.575,2,0.2844911],[5.575,2,0.2780086],[6.575,2,0.2909731],[7.575,2,0.4072804],[8.575,2,0.4004494],[0.575,3,0.0602362],[1.575,3,0.1052224],[2.575,3,0.1872119],[3.575,3,0.1992938],[4.575,3,0.4122555],[5.575,3,0.405629],[6.575,3,0.4186724],[7.575,3,0.4992781],[8.575,3,0.4962812],[0.575,4,0.07289286],[1.575,4,0.1466753],[2.575,4,0.2735617],[3.575,4,0.2904995],[4.575,4,0.4931282],[5.575,4,0.4897055],[6.575,4,0.4962616],[7.575,4,0.5202362],[8.575,4,0.5198667],[0.575,5,0.08924903],[1.575,5,0.1999799],[2.575,5,0.363834],[3.575,5,0.3820702],[4.575,5,0.5178087],[5.575,5,0.5169982],[6.575,5,0.5184923],[7.575,5,0.5215859],[8.575,5,0.5215747],[0.575,6,0.1214912],[1.575,6,0.2943431],[2.575,6,0.4657945],[3.575,6,0.4777715],[4.575,6,0.5215657],[5.575,6,0.5215465],[6.575,6,0.5215791],[7.575,6,0.5216089],[8.575,6,0.5216089],[0.575,7,0.161152],[1.575,7,0.3829985],[2.575,7,0.5087754],[3.575,7,0.5130559],[4.575,7,0.5216089],[5.575,7,0.5216088],[6.575,7,0.5216089],[7.575,7,0.5216089],[8.575,7,0.5216089],[0.575,8,0.235947],[1.575,8,0.4769049],[2.575,8,0.5211299],[3.575,8,0.5213853],[4.575,8,0.5216089],[5.575,8,0.5216089],[6.575,8,0.5216089],[7.575,8,0.5216089],[8.575,8,0.5216089],[0.575,9,0.3159027],[1.575,9,0.512537],[2.575,9,0.5216036],[3.575,9,0.5216074],[4.575,9,0.5216089],[5.575,9,0.5216089],[6.575,9,0.5216089],[7.575,9,0.5216089],[8.575,9,0.5216089],[0.575,10,0.3783879],[1.575,10,0.5199995],[2.575,10,0.5216089],[3.575,10,0.5216089],[4.575,10,0.5216089],[5.575,10,0.5216089],[6.575,10,0.5216089],[7.575,10,0.5216089],[8.575,10,0.5216089],[1.425,1,0.04321798],[2.425,1,0.05459053],[3.425,1,0.07270695],[4.425,1,0.07537223],[5.425,1,0.1374316],[6.425,1,0.1345138],[7.425,1,0.140395],[8.425,1,0.206122],[9.425,1,0.2012519],[1.425,2,0.05118433],[2.425,2,0.0770334],[3.425,2,0.1231276],[4.425,2,0.1301295],[5.425,2,0.2844911],[6.425,2,0.2780086],[7.425,2,0.2909731],[8.425,2,0.4072804],[9.425,2,0.4004494],[1.425,3,0.0602362],[2.425,3,0.1052224],[3.425,3,0.1872119],[4.425,3,0.1992938],[5.425,3,0.4122555],[6.425,3,0.405629],[7.425,3,0.4186724],[8.425,3,0.4992781],[9.425,3,0.4962812],[1.425,4,0.07289286],[2.425,4,0.1466753],[3.425,4,0.2735617],[4.425,4,0.2904995],[5.425,4,0.4931282],[6.425,4,0.4897055],[7.425,4,0.4962616],[8.425,4,0.5202362],[9.425,4,0.5198667],[1.425,5,0.08924903],[2.425,5,0.1999799],[3.425,5,0.363834],[4.425,5,0.3820702],[5.425,5,0.5178087],[6.425,5,0.5169982],[7.425,5,0.5184923],[8.425,5,0.5215859],[9.425,5,0.5215747],[1.425,6,0.1214912],[2.425,6,0.2943431],[3.425,6,0.4657945],[4.425,6,0.4777715],[5.425,6,0.5215657],[6.425,6,0.5215465],[7.425,6,0.5215791],[8.425,6,0.5216089],[9.425,6,0.5216089],[1.425,7,0.161152],[2.425,7,0.3829985],[3.425,7,0.5087754],[4.425,7,0.5130559],[5.425,7,0.5216089],[6.425,7,0.5216088],[7.425,7,0.5216089],[8.425,7,0.5216089],[9.425,7,0.5216089],[1.425,8,0.235947],[2.425,8,0.4769049],[3.425,8,0.5211299],[4.425,8,0.5213853],[5.425,8,0.5216089],[6.425,8,0.5216089],[7.425,8,0.5216089],[8.425,8,0.5216089],[9.425,8,0.5216089],[1.425,9,0.3159027],[2.425,9,0.512537],[3.425,9,0.5216036],[4.425,9,0.5216074],[5.425,9,0.5216089],[6.425,9,0.5216089],[7.425,9,0.5216089],[8.425,9,0.5216089],[9.425,9,0.5216089],[1.425,10,0.3783879],[2.425,10,0.5199995],[3.425,10,0.5216089],[4.425,10,0.5216089],[5.425,10,0.5216089],[6.425,10,0.5216089],[7.425,10,0.5216089],[8.425,10,0.5216089],[9.425,10,0.5216089],[1,0.9999999,0.04321798],[2,0.9999999,0.06596308],[3,0.9999999,0.1021959],[4,0.9999999,0.1075265],[5,0.9999999,0.2316452],[6,0.9999999,0.2258095],[7,0.9999999,0.237572],[8,0.9999999,0.3690261],[9,0.9999999,0.3592858],[1,2,0.05915068],[2,2,0.1108488],[3,2,0.2030372],[4,2,0.217041],[5,2,0.5257641],[6,2,0.5127993],[7,2,0.5387281],[8,2,0.7713429],[9,2,0.7576808],[1,3,0.07725442],[2,3,0.1672267],[3,3,0.3312058],[4,3,0.3553697],[5,3,0.781293],[6,3,0.76804],[7,3,0.7941267],[8,3,0.9553382],[9,3,0.9493445],[1,4,0.1025677],[2,4,0.2501325],[3,4,0.5039053],[4,4,0.537781],[5,4,0.9430385],[6,4,0.9361931],[7,4,0.9493053],[8,4,0.9972544],[9,4,0.9965154],[1,5,0.1352801],[2,5,0.3567418],[3,5,0.6844499],[4,5,0.7209224],[5,5,0.9923994],[6,5,0.9907784],[7,5,0.9937668],[8,5,0.999954],[9,5,0.9999315],[1,6,0.1997644],[2,6,0.5454682],[3,6,0.8883711],[4,6,0.912325],[5,6,0.9999133],[6,6,0.9998751],[7,6,0.9999403],[8,6,1],[9,6,1],[1,7,0.2790861],[2,7,0.722779],[3,7,0.9743328],[4,7,0.9828939],[5,7,0.9999998],[6,7,0.9999996],[7,7,0.9999999],[8,7,1],[9,7,1],[1,8,0.428676],[2,8,0.9105917],[3,8,0.9990419],[4,8,0.9995526],[5,8,1],[6,8,1],[7,8,1],[8,8,1],[9,8,1],[1,9,0.5885873],[2,9,0.981856],[3,9,0.9999893],[4,9,0.9999969],[5,9,1],[6,9,1],[7,9,1],[8,9,1],[9,9,1],[1,10,0.7135577],[2,10,0.9967811],[3,10,0.9999999],[4,10,1],[5,10,1],[6,10,1],[7,10,1],[8,10,1],[9,10,1]],"ignoreExtent":false,"flags":2},"23":{"id":23,"type":"linestrip","material":{},"vertices":[[0.575,0.575,0.04321798],[1.425,0.575,0.04321798],[1.425,0.575,0.04321798],[0.575,0.575,0.04321798],[0.575,0.575,0.04321798],["NaN","NaN","NaN"],[1.575,0.575,0.04321798],[2.425,0.575,0.04321798],[2.425,0.575,0.06596308],[1.575,0.575,0.06596308],[1.575,0.575,0.04321798],["NaN","NaN","NaN"],[2.575,0.575,0.04321798],[3.425,0.575,0.04321798],[3.425,0.575,0.1021959],[2.575,0.575,0.1021959],[2.575,0.575,0.04321798],["NaN","NaN","NaN"],[3.575,0.575,0.04321798],[4.425,0.575,0.04321798],[4.425,0.575,0.1075265],[3.575,0.575,0.1075265],[3.575,0.575,0.04321798],["NaN","NaN","NaN"],[4.575,0.575,0.04321798],[5.425,0.575,0.04321798],[5.425,0.575,0.2316452],[4.575,0.575,0.2316452],[4.575,0.575,0.04321798],["NaN","NaN","NaN"],[5.575,0.575,0.04321798],[6.425,0.575,0.04321798],[6.425,0.575,0.2258095],[5.575,0.575,0.2258095],[5.575,0.575,0.04321798],["NaN","NaN","NaN"],[6.575,0.575,0.04321798],[7.425,0.575,0.04321798],[7.425,0.575,0.237572],[6.575,0.575,0.237572],[6.575,0.575,0.04321798],["NaN","NaN","NaN"],[7.575,0.575,0.04321798],[8.425,0.575,0.04321798],[8.425,0.575,0.3690261],[7.575,0.575,0.3690261],[7.575,0.575,0.04321798],["NaN","NaN","NaN"],[8.575,0.575,0.04321798],[9.425,0.575,0.04321798],[9.425,0.575,0.3592858],[8.575,0.575,0.3592858],[8.575,0.575,0.04321798],["NaN","NaN","NaN"],[0.575,1.575,0.04321798],[1.425,1.575,0.04321798],[1.425,1.575,0.05915068],[0.575,1.575,0.05915068],[0.575,1.575,0.04321798],["NaN","NaN","NaN"],[1.575,1.575,0.04321798],[2.425,1.575,0.04321798],[2.425,1.575,0.1108488],[1.575,1.575,0.1108488],[1.575,1.575,0.04321798],["NaN","NaN","NaN"],[2.575,1.575,0.04321798],[3.425,1.575,0.04321798],[3.425,1.575,0.2030372],[2.575,1.575,0.2030372],[2.575,1.575,0.04321798],["NaN","NaN","NaN"],[3.575,1.575,0.04321798],[4.425,1.575,0.04321798],[4.425,1.575,0.217041],[3.575,1.575,0.217041],[3.575,1.575,0.04321798],["NaN","NaN","NaN"],[4.575,1.575,0.04321798],[5.425,1.575,0.04321798],[5.425,1.575,0.5257641],[4.575,1.575,0.5257641],[4.575,1.575,0.04321798],["NaN","NaN","NaN"],[5.575,1.575,0.04321798],[6.425,1.575,0.04321798],[6.425,1.575,0.5127993],[5.575,1.575,0.5127993],[5.575,1.575,0.04321798],["NaN","NaN","NaN"],[6.575,1.575,0.04321798],[7.425,1.575,0.04321798],[7.425,1.575,0.5387281],[6.575,1.575,0.5387281],[6.575,1.575,0.04321798],["NaN","NaN","NaN"],[7.575,1.575,0.04321798],[8.425,1.575,0.04321798],[8.425,1.575,0.7713429],[7.575,1.575,0.7713429],[7.575,1.575,0.04321798],["NaN","NaN","NaN"],[8.575,1.575,0.04321798],[9.425,1.575,0.04321798],[9.425,1.575,0.7576808],[8.575,1.575,0.7576808],[8.575,1.575,0.04321798],["NaN","NaN","NaN"],[0.575,2.575,0.04321798],[1.425,2.575,0.04321798],[1.425,2.575,0.07725442],[0.575,2.575,0.07725442],[0.575,2.575,0.04321798],["NaN","NaN","NaN"],[1.575,2.575,0.04321798],[2.425,2.575,0.04321798],[2.425,2.575,0.1672267],[1.575,2.575,0.1672267],[1.575,2.575,0.04321798],["NaN","NaN","NaN"],[2.575,2.575,0.04321798],[3.425,2.575,0.04321798],[3.425,2.575,0.3312058],[2.575,2.575,0.3312058],[2.575,2.575,0.04321798],["NaN","NaN","NaN"],[3.575,2.575,0.04321798],[4.425,2.575,0.04321798],[4.425,2.575,0.3553697],[3.575,2.575,0.3553697],[3.575,2.575,0.04321798],["NaN","NaN","NaN"],[4.575,2.575,0.04321798],[5.425,2.575,0.04321798],[5.425,2.575,0.781293],[4.575,2.575,0.781293],[4.575,2.575,0.04321798],["NaN","NaN","NaN"],[5.575,2.575,0.04321798],[6.425,2.575,0.04321798],[6.425,2.575,0.76804],[5.575,2.575,0.76804],[5.575,2.575,0.04321798],["NaN","NaN","NaN"],[6.575,2.575,0.04321798],[7.425,2.575,0.04321798],[7.425,2.575,0.7941267],[6.575,2.575,0.7941267],[6.575,2.575,0.04321798],["NaN","NaN","NaN"],[7.575,2.575,0.04321798],[8.425,2.575,0.04321798],[8.425,2.575,0.9553382],[7.575,2.575,0.9553382],[7.575,2.575,0.04321798],["NaN","NaN","NaN"],[8.575,2.575,0.04321798],[9.425,2.575,0.04321798],[9.425,2.575,0.9493445],[8.575,2.575,0.9493445],[8.575,2.575,0.04321798],["NaN","NaN","NaN"],[0.575,3.575,0.04321798],[1.425,3.575,0.04321798],[1.425,3.575,0.1025677],[0.575,3.575,0.1025677],[0.575,3.575,0.04321798],["NaN","NaN","NaN"],[1.575,3.575,0.04321798],[2.425,3.575,0.04321798],[2.425,3.575,0.2501325],[1.575,3.575,0.2501325],[1.575,3.575,0.04321798],["NaN","NaN","NaN"],[2.575,3.575,0.04321798],[3.425,3.575,0.04321798],[3.425,3.575,0.5039053],[2.575,3.575,0.5039053],[2.575,3.575,0.04321798],["NaN","NaN","NaN"],[3.575,3.575,0.04321798],[4.425,3.575,0.04321798],[4.425,3.575,0.537781],[3.575,3.575,0.537781],[3.575,3.575,0.04321798],["NaN","NaN","NaN"],[4.575,3.575,0.04321798],[5.425,3.575,0.04321798],[5.425,3.575,0.9430385],[4.575,3.575,0.9430385],[4.575,3.575,0.04321798],["NaN","NaN","NaN"],[5.575,3.575,0.04321798],[6.425,3.575,0.04321798],[6.425,3.575,0.9361931],[5.575,3.575,0.9361931],[5.575,3.575,0.04321798],["NaN","NaN","NaN"],[6.575,3.575,0.04321798],[7.425,3.575,0.04321798],[7.425,3.575,0.9493053],[6.575,3.575,0.9493053],[6.575,3.575,0.04321798],["NaN","NaN","NaN"],[7.575,3.575,0.04321798],[8.425,3.575,0.04321798],[8.425,3.575,0.9972544],[7.575,3.575,0.9972544],[7.575,3.575,0.04321798],["NaN","NaN","NaN"],[8.575,3.575,0.04321798],[9.425,3.575,0.04321798],[9.425,3.575,0.9965154],[8.575,3.575,0.9965154],[8.575,3.575,0.04321798],["NaN","NaN","NaN"],[0.575,4.575,0.04321798],[1.425,4.575,0.04321798],[1.425,4.575,0.1352801],[0.575,4.575,0.1352801],[0.575,4.575,0.04321798],["NaN","NaN","NaN"],[1.575,4.575,0.04321798],[2.425,4.575,0.04321798],[2.425,4.575,0.3567418],[1.575,4.575,0.3567418],[1.575,4.575,0.04321798],["NaN","NaN","NaN"],[2.575,4.575,0.04321798],[3.425,4.575,0.04321798],[3.425,4.575,0.6844499],[2.575,4.575,0.6844499],[2.575,4.575,0.04321798],["NaN","NaN","NaN"],[3.575,4.575,0.04321798],[4.425,4.575,0.04321798],[4.425,4.575,0.7209224],[3.575,4.575,0.7209224],[3.575,4.575,0.04321798],["NaN","NaN","NaN"],[4.575,4.575,0.04321798],[5.425,4.575,0.04321798],[5.425,4.575,0.9923994],[4.575,4.575,0.9923994],[4.575,4.575,0.04321798],["NaN","NaN","NaN"],[5.575,4.575,0.04321798],[6.425,4.575,0.04321798],[6.425,4.575,0.9907784],[5.575,4.575,0.9907784],[5.575,4.575,0.04321798],["NaN","NaN","NaN"],[6.575,4.575,0.04321798],[7.425,4.575,0.04321798],[7.425,4.575,0.9937668],[6.575,4.575,0.9937668],[6.575,4.575,0.04321798],["NaN","NaN","NaN"],[7.575,4.575,0.04321798],[8.425,4.575,0.04321798],[8.425,4.575,0.999954],[7.575,4.575,0.999954],[7.575,4.575,0.04321798],["NaN","NaN","NaN"],[8.575,4.575,0.04321798],[9.425,4.575,0.04321798],[9.425,4.575,0.9999315],[8.575,4.575,0.9999315],[8.575,4.575,0.04321798],["NaN","NaN","NaN"],[0.575,5.575,0.04321798],[1.425,5.575,0.04321798],[1.425,5.575,0.1997644],[0.575,5.575,0.1997644],[0.575,5.575,0.04321798],["NaN","NaN","NaN"],[1.575,5.575,0.04321798],[2.425,5.575,0.04321798],[2.425,5.575,0.5454682],[1.575,5.575,0.5454682],[1.575,5.575,0.04321798],["NaN","NaN","NaN"],[2.575,5.575,0.04321798],[3.425,5.575,0.04321798],[3.425,5.575,0.8883711],[2.575,5.575,0.8883711],[2.575,5.575,0.04321798],["NaN","NaN","NaN"],[3.575,5.575,0.04321798],[4.425,5.575,0.04321798],[4.425,5.575,0.912325],[3.575,5.575,0.912325],[3.575,5.575,0.04321798],["NaN","NaN","NaN"],[4.575,5.575,0.04321798],[5.425,5.575,0.04321798],[5.425,5.575,0.9999133],[4.575,5.575,0.9999133],[4.575,5.575,0.04321798],["NaN","NaN","NaN"],[5.575,5.575,0.04321798],[6.425,5.575,0.04321798],[6.425,5.575,0.9998751],[5.575,5.575,0.9998751],[5.575,5.575,0.04321798],["NaN","NaN","NaN"],[6.575,5.575,0.04321798],[7.425,5.575,0.04321798],[7.425,5.575,0.9999403],[6.575,5.575,0.9999403],[6.575,5.575,0.04321798],["NaN","NaN","NaN"],[7.575,5.575,0.04321798],[8.425,5.575,0.04321798],[8.425,5.575,1],[7.575,5.575,1],[7.575,5.575,0.04321798],["NaN","NaN","NaN"],[8.575,5.575,0.04321798],[9.425,5.575,0.04321798],[9.425,5.575,1],[8.575,5.575,1],[8.575,5.575,0.04321798],["NaN","NaN","NaN"],[0.575,6.575,0.04321798],[1.425,6.575,0.04321798],[1.425,6.575,0.2790861],[0.575,6.575,0.2790861],[0.575,6.575,0.04321798],["NaN","NaN","NaN"],[1.575,6.575,0.04321798],[2.425,6.575,0.04321798],[2.425,6.575,0.722779],[1.575,6.575,0.722779],[1.575,6.575,0.04321798],["NaN","NaN","NaN"],[2.575,6.575,0.04321798],[3.425,6.575,0.04321798],[3.425,6.575,0.9743328],[2.575,6.575,0.9743328],[2.575,6.575,0.04321798],["NaN","NaN","NaN"],[3.575,6.575,0.04321798],[4.425,6.575,0.04321798],[4.425,6.575,0.9828939],[3.575,6.575,0.9828939],[3.575,6.575,0.04321798],["NaN","NaN","NaN"],[4.575,6.575,0.04321798],[5.425,6.575,0.04321798],[5.425,6.575,0.9999998],[4.575,6.575,0.9999998],[4.575,6.575,0.04321798],["NaN","NaN","NaN"],[5.575,6.575,0.04321798],[6.425,6.575,0.04321798],[6.425,6.575,0.9999996],[5.575,6.575,0.9999996],[5.575,6.575,0.04321798],["NaN","NaN","NaN"],[6.575,6.575,0.04321798],[7.425,6.575,0.04321798],[7.425,6.575,0.9999999],[6.575,6.575,0.9999999],[6.575,6.575,0.04321798],["NaN","NaN","NaN"],[7.575,6.575,0.04321798],[8.425,6.575,0.04321798],[8.425,6.575,1],[7.575,6.575,1],[7.575,6.575,0.04321798],["NaN","NaN","NaN"],[8.575,6.575,0.04321798],[9.425,6.575,0.04321798],[9.425,6.575,1],[8.575,6.575,1],[8.575,6.575,0.04321798],["NaN","NaN","NaN"],[0.575,7.575,0.04321798],[1.425,7.575,0.04321798],[1.425,7.575,0.428676],[0.575,7.575,0.428676],[0.575,7.575,0.04321798],["NaN","NaN","NaN"],[1.575,7.575,0.04321798],[2.425,7.575,0.04321798],[2.425,7.575,0.9105917],[1.575,7.575,0.9105917],[1.575,7.575,0.04321798],["NaN","NaN","NaN"],[2.575,7.575,0.04321798],[3.425,7.575,0.04321798],[3.425,7.575,0.9990419],[2.575,7.575,0.9990419],[2.575,7.575,0.04321798],["NaN","NaN","NaN"],[3.575,7.575,0.04321798],[4.425,7.575,0.04321798],[4.425,7.575,0.9995526],[3.575,7.575,0.9995526],[3.575,7.575,0.04321798],["NaN","NaN","NaN"],[4.575,7.575,0.04321798],[5.425,7.575,0.04321798],[5.425,7.575,1],[4.575,7.575,1],[4.575,7.575,0.04321798],["NaN","NaN","NaN"],[5.575,7.575,0.04321798],[6.425,7.575,0.04321798],[6.425,7.575,1],[5.575,7.575,1],[5.575,7.575,0.04321798],["NaN","NaN","NaN"],[6.575,7.575,0.04321798],[7.425,7.575,0.04321798],[7.425,7.575,1],[6.575,7.575,1],[6.575,7.575,0.04321798],["NaN","NaN","NaN"],[7.575,7.575,0.04321798],[8.425,7.575,0.04321798],[8.425,7.575,1],[7.575,7.575,1],[7.575,7.575,0.04321798],["NaN","NaN","NaN"],[8.575,7.575,0.04321798],[9.425,7.575,0.04321798],[9.425,7.575,1],[8.575,7.575,1],[8.575,7.575,0.04321798],["NaN","NaN","NaN"],[0.575,8.575,0.04321798],[1.425,8.575,0.04321798],[1.425,8.575,0.5885873],[0.575,8.575,0.5885873],[0.575,8.575,0.04321798],["NaN","NaN","NaN"],[1.575,8.575,0.04321798],[2.425,8.575,0.04321798],[2.425,8.575,0.981856],[1.575,8.575,0.981856],[1.575,8.575,0.04321798],["NaN","NaN","NaN"],[2.575,8.575,0.04321798],[3.425,8.575,0.04321798],[3.425,8.575,0.9999893],[2.575,8.575,0.9999893],[2.575,8.575,0.04321798],["NaN","NaN","NaN"],[3.575,8.575,0.04321798],[4.425,8.575,0.04321798],[4.425,8.575,0.9999969],[3.575,8.575,0.9999969],[3.575,8.575,0.04321798],["NaN","NaN","NaN"],[4.575,8.575,0.04321798],[5.425,8.575,0.04321798],[5.425,8.575,1],[4.575,8.575,1],[4.575,8.575,0.04321798],["NaN","NaN","NaN"],[5.575,8.575,0.04321798],[6.425,8.575,0.04321798],[6.425,8.575,1],[5.575,8.575,1],[5.575,8.575,0.04321798],["NaN","NaN","NaN"],[6.575,8.575,0.04321798],[7.425,8.575,0.04321798],[7.425,8.575,1],[6.575,8.575,1],[6.575,8.575,0.04321798],["NaN","NaN","NaN"],[7.575,8.575,0.04321798],[8.425,8.575,0.04321798],[8.425,8.575,1],[7.575,8.575,1],[7.575,8.575,0.04321798],["NaN","NaN","NaN"],[8.575,8.575,0.04321798],[9.425,8.575,0.04321798],[9.425,8.575,1],[8.575,8.575,1],[8.575,8.575,0.04321798],["NaN","NaN","NaN"],[0.575,9.575,0.04321798],[1.425,9.575,0.04321798],[1.425,9.575,0.7135577],[0.575,9.575,0.7135577],[0.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.575,9.575,0.04321798],[2.425,9.575,0.04321798],[2.425,9.575,0.9967811],[1.575,9.575,0.9967811],[1.575,9.575,0.04321798],["NaN","NaN","NaN"],[2.575,9.575,0.04321798],[3.425,9.575,0.04321798],[3.425,9.575,0.9999999],[2.575,9.575,0.9999999],[2.575,9.575,0.04321798],["NaN","NaN","NaN"],[3.575,9.575,0.04321798],[4.425,9.575,0.04321798],[4.425,9.575,1],[3.575,9.575,1],[3.575,9.575,0.04321798],["NaN","NaN","NaN"],[4.575,9.575,0.04321798],[5.425,9.575,0.04321798],[5.425,9.575,1],[4.575,9.575,1],[4.575,9.575,0.04321798],["NaN","NaN","NaN"],[5.575,9.575,0.04321798],[6.425,9.575,0.04321798],[6.425,9.575,1],[5.575,9.575,1],[5.575,9.575,0.04321798],["NaN","NaN","NaN"],[6.575,9.575,0.04321798],[7.425,9.575,0.04321798],[7.425,9.575,1],[6.575,9.575,1],[6.575,9.575,0.04321798],["NaN","NaN","NaN"],[7.575,9.575,0.04321798],[8.425,9.575,0.04321798],[8.425,9.575,1],[7.575,9.575,1],[7.575,9.575,0.04321798],["NaN","NaN","NaN"],[8.575,9.575,0.04321798],[9.425,9.575,0.04321798],[9.425,9.575,1],[8.575,9.575,1],[8.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.425,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,1.425,0.04321798],[1.425,1.425,0.04321798],[1.425,1.425,0.04321798],["NaN","NaN","NaN"],[2.425,1.425,0.04321798],[1.575,1.425,0.04321798],[1.575,1.425,0.06596308],[2.425,1.425,0.06596308],[2.425,1.425,0.04321798],["NaN","NaN","NaN"],[3.425,1.425,0.04321798],[2.575,1.425,0.04321798],[2.575,1.425,0.1021959],[3.425,1.425,0.1021959],[3.425,1.425,0.04321798],["NaN","NaN","NaN"],[4.425,1.425,0.04321798],[3.575,1.425,0.04321798],[3.575,1.425,0.1075265],[4.425,1.425,0.1075265],[4.425,1.425,0.04321798],["NaN","NaN","NaN"],[5.425,1.425,0.04321798],[4.575,1.425,0.04321798],[4.575,1.425,0.2316452],[5.425,1.425,0.2316452],[5.425,1.425,0.04321798],["NaN","NaN","NaN"],[6.425,1.425,0.04321798],[5.575,1.425,0.04321798],[5.575,1.425,0.2258095],[6.425,1.425,0.2258095],[6.425,1.425,0.04321798],["NaN","NaN","NaN"],[7.425,1.425,0.04321798],[6.575,1.425,0.04321798],[6.575,1.425,0.237572],[7.425,1.425,0.237572],[7.425,1.425,0.04321798],["NaN","NaN","NaN"],[8.425,1.425,0.04321798],[7.575,1.425,0.04321798],[7.575,1.425,0.3690261],[8.425,1.425,0.3690261],[8.425,1.425,0.04321798],["NaN","NaN","NaN"],[9.425,1.425,0.04321798],[8.575,1.425,0.04321798],[8.575,1.425,0.3592858],[9.425,1.425,0.3592858],[9.425,1.425,0.04321798],["NaN","NaN","NaN"],[1.425,2.425,0.04321798],[0.575,2.425,0.04321798],[0.575,2.425,0.05915068],[1.425,2.425,0.05915068],[1.425,2.425,0.04321798],["NaN","NaN","NaN"],[2.425,2.425,0.04321798],[1.575,2.425,0.04321798],[1.575,2.425,0.1108488],[2.425,2.425,0.1108488],[2.425,2.425,0.04321798],["NaN","NaN","NaN"],[3.425,2.425,0.04321798],[2.575,2.425,0.04321798],[2.575,2.425,0.2030372],[3.425,2.425,0.2030372],[3.425,2.425,0.04321798],["NaN","NaN","NaN"],[4.425,2.425,0.04321798],[3.575,2.425,0.04321798],[3.575,2.425,0.217041],[4.425,2.425,0.217041],[4.425,2.425,0.04321798],["NaN","NaN","NaN"],[5.425,2.425,0.04321798],[4.575,2.425,0.04321798],[4.575,2.425,0.5257641],[5.425,2.425,0.5257641],[5.425,2.425,0.04321798],["NaN","NaN","NaN"],[6.425,2.425,0.04321798],[5.575,2.425,0.04321798],[5.575,2.425,0.5127993],[6.425,2.425,0.5127993],[6.425,2.425,0.04321798],["NaN","NaN","NaN"],[7.425,2.425,0.04321798],[6.575,2.425,0.04321798],[6.575,2.425,0.5387281],[7.425,2.425,0.5387281],[7.425,2.425,0.04321798],["NaN","NaN","NaN"],[8.425,2.425,0.04321798],[7.575,2.425,0.04321798],[7.575,2.425,0.7713429],[8.425,2.425,0.7713429],[8.425,2.425,0.04321798],["NaN","NaN","NaN"],[9.425,2.425,0.04321798],[8.575,2.425,0.04321798],[8.575,2.425,0.7576808],[9.425,2.425,0.7576808],[9.425,2.425,0.04321798],["NaN","NaN","NaN"],[1.425,3.425,0.04321798],[0.575,3.425,0.04321798],[0.575,3.425,0.07725442],[1.425,3.425,0.07725442],[1.425,3.425,0.04321798],["NaN","NaN","NaN"],[2.425,3.425,0.04321798],[1.575,3.425,0.04321798],[1.575,3.425,0.1672267],[2.425,3.425,0.1672267],[2.425,3.425,0.04321798],["NaN","NaN","NaN"],[3.425,3.425,0.04321798],[2.575,3.425,0.04321798],[2.575,3.425,0.3312058],[3.425,3.425,0.3312058],[3.425,3.425,0.04321798],["NaN","NaN","NaN"],[4.425,3.425,0.04321798],[3.575,3.425,0.04321798],[3.575,3.425,0.3553697],[4.425,3.425,0.3553697],[4.425,3.425,0.04321798],["NaN","NaN","NaN"],[5.425,3.425,0.04321798],[4.575,3.425,0.04321798],[4.575,3.425,0.781293],[5.425,3.425,0.781293],[5.425,3.425,0.04321798],["NaN","NaN","NaN"],[6.425,3.425,0.04321798],[5.575,3.425,0.04321798],[5.575,3.425,0.76804],[6.425,3.425,0.76804],[6.425,3.425,0.04321798],["NaN","NaN","NaN"],[7.425,3.425,0.04321798],[6.575,3.425,0.04321798],[6.575,3.425,0.7941267],[7.425,3.425,0.7941267],[7.425,3.425,0.04321798],["NaN","NaN","NaN"],[8.425,3.425,0.04321798],[7.575,3.425,0.04321798],[7.575,3.425,0.9553382],[8.425,3.425,0.9553382],[8.425,3.425,0.04321798],["NaN","NaN","NaN"],[9.425,3.425,0.04321798],[8.575,3.425,0.04321798],[8.575,3.425,0.9493445],[9.425,3.425,0.9493445],[9.425,3.425,0.04321798],["NaN","NaN","NaN"],[1.425,4.425,0.04321798],[0.575,4.425,0.04321798],[0.575,4.425,0.1025677],[1.425,4.425,0.1025677],[1.425,4.425,0.04321798],["NaN","NaN","NaN"],[2.425,4.425,0.04321798],[1.575,4.425,0.04321798],[1.575,4.425,0.2501325],[2.425,4.425,0.2501325],[2.425,4.425,0.04321798],["NaN","NaN","NaN"],[3.425,4.425,0.04321798],[2.575,4.425,0.04321798],[2.575,4.425,0.5039053],[3.425,4.425,0.5039053],[3.425,4.425,0.04321798],["NaN","NaN","NaN"],[4.425,4.425,0.04321798],[3.575,4.425,0.04321798],[3.575,4.425,0.537781],[4.425,4.425,0.537781],[4.425,4.425,0.04321798],["NaN","NaN","NaN"],[5.425,4.425,0.04321798],[4.575,4.425,0.04321798],[4.575,4.425,0.9430385],[5.425,4.425,0.9430385],[5.425,4.425,0.04321798],["NaN","NaN","NaN"],[6.425,4.425,0.04321798],[5.575,4.425,0.04321798],[5.575,4.425,0.9361931],[6.425,4.425,0.9361931],[6.425,4.425,0.04321798],["NaN","NaN","NaN"],[7.425,4.425,0.04321798],[6.575,4.425,0.04321798],[6.575,4.425,0.9493053],[7.425,4.425,0.9493053],[7.425,4.425,0.04321798],["NaN","NaN","NaN"],[8.425,4.425,0.04321798],[7.575,4.425,0.04321798],[7.575,4.425,0.9972544],[8.425,4.425,0.9972544],[8.425,4.425,0.04321798],["NaN","NaN","NaN"],[9.425,4.425,0.04321798],[8.575,4.425,0.04321798],[8.575,4.425,0.9965154],[9.425,4.425,0.9965154],[9.425,4.425,0.04321798],["NaN","NaN","NaN"],[1.425,5.425,0.04321798],[0.575,5.425,0.04321798],[0.575,5.425,0.1352801],[1.425,5.425,0.1352801],[1.425,5.425,0.04321798],["NaN","NaN","NaN"],[2.425,5.425,0.04321798],[1.575,5.425,0.04321798],[1.575,5.425,0.3567418],[2.425,5.425,0.3567418],[2.425,5.425,0.04321798],["NaN","NaN","NaN"],[3.425,5.425,0.04321798],[2.575,5.425,0.04321798],[2.575,5.425,0.6844499],[3.425,5.425,0.6844499],[3.425,5.425,0.04321798],["NaN","NaN","NaN"],[4.425,5.425,0.04321798],[3.575,5.425,0.04321798],[3.575,5.425,0.7209224],[4.425,5.425,0.7209224],[4.425,5.425,0.04321798],["NaN","NaN","NaN"],[5.425,5.425,0.04321798],[4.575,5.425,0.04321798],[4.575,5.425,0.9923994],[5.425,5.425,0.9923994],[5.425,5.425,0.04321798],["NaN","NaN","NaN"],[6.425,5.425,0.04321798],[5.575,5.425,0.04321798],[5.575,5.425,0.9907784],[6.425,5.425,0.9907784],[6.425,5.425,0.04321798],["NaN","NaN","NaN"],[7.425,5.425,0.04321798],[6.575,5.425,0.04321798],[6.575,5.425,0.9937668],[7.425,5.425,0.9937668],[7.425,5.425,0.04321798],["NaN","NaN","NaN"],[8.425,5.425,0.04321798],[7.575,5.425,0.04321798],[7.575,5.425,0.999954],[8.425,5.425,0.999954],[8.425,5.425,0.04321798],["NaN","NaN","NaN"],[9.425,5.425,0.04321798],[8.575,5.425,0.04321798],[8.575,5.425,0.9999315],[9.425,5.425,0.9999315],[9.425,5.425,0.04321798],["NaN","NaN","NaN"],[1.425,6.425,0.04321798],[0.575,6.425,0.04321798],[0.575,6.425,0.1997644],[1.425,6.425,0.1997644],[1.425,6.425,0.04321798],["NaN","NaN","NaN"],[2.425,6.425,0.04321798],[1.575,6.425,0.04321798],[1.575,6.425,0.5454682],[2.425,6.425,0.5454682],[2.425,6.425,0.04321798],["NaN","NaN","NaN"],[3.425,6.425,0.04321798],[2.575,6.425,0.04321798],[2.575,6.425,0.8883711],[3.425,6.425,0.8883711],[3.425,6.425,0.04321798],["NaN","NaN","NaN"],[4.425,6.425,0.04321798],[3.575,6.425,0.04321798],[3.575,6.425,0.912325],[4.425,6.425,0.912325],[4.425,6.425,0.04321798],["NaN","NaN","NaN"],[5.425,6.425,0.04321798],[4.575,6.425,0.04321798],[4.575,6.425,0.9999133],[5.425,6.425,0.9999133],[5.425,6.425,0.04321798],["NaN","NaN","NaN"],[6.425,6.425,0.04321798],[5.575,6.425,0.04321798],[5.575,6.425,0.9998751],[6.425,6.425,0.9998751],[6.425,6.425,0.04321798],["NaN","NaN","NaN"],[7.425,6.425,0.04321798],[6.575,6.425,0.04321798],[6.575,6.425,0.9999403],[7.425,6.425,0.9999403],[7.425,6.425,0.04321798],["NaN","NaN","NaN"],[8.425,6.425,0.04321798],[7.575,6.425,0.04321798],[7.575,6.425,1],[8.425,6.425,1],[8.425,6.425,0.04321798],["NaN","NaN","NaN"],[9.425,6.425,0.04321798],[8.575,6.425,0.04321798],[8.575,6.425,1],[9.425,6.425,1],[9.425,6.425,0.04321798],["NaN","NaN","NaN"],[1.425,7.425,0.04321798],[0.575,7.425,0.04321798],[0.575,7.425,0.2790861],[1.425,7.425,0.2790861],[1.425,7.425,0.04321798],["NaN","NaN","NaN"],[2.425,7.425,0.04321798],[1.575,7.425,0.04321798],[1.575,7.425,0.722779],[2.425,7.425,0.722779],[2.425,7.425,0.04321798],["NaN","NaN","NaN"],[3.425,7.425,0.04321798],[2.575,7.425,0.04321798],[2.575,7.425,0.9743328],[3.425,7.425,0.9743328],[3.425,7.425,0.04321798],["NaN","NaN","NaN"],[4.425,7.425,0.04321798],[3.575,7.425,0.04321798],[3.575,7.425,0.9828939],[4.425,7.425,0.9828939],[4.425,7.425,0.04321798],["NaN","NaN","NaN"],[5.425,7.425,0.04321798],[4.575,7.425,0.04321798],[4.575,7.425,0.9999998],[5.425,7.425,0.9999998],[5.425,7.425,0.04321798],["NaN","NaN","NaN"],[6.425,7.425,0.04321798],[5.575,7.425,0.04321798],[5.575,7.425,0.9999996],[6.425,7.425,0.9999996],[6.425,7.425,0.04321798],["NaN","NaN","NaN"],[7.425,7.425,0.04321798],[6.575,7.425,0.04321798],[6.575,7.425,0.9999999],[7.425,7.425,0.9999999],[7.425,7.425,0.04321798],["NaN","NaN","NaN"],[8.425,7.425,0.04321798],[7.575,7.425,0.04321798],[7.575,7.425,1],[8.425,7.425,1],[8.425,7.425,0.04321798],["NaN","NaN","NaN"],[9.425,7.425,0.04321798],[8.575,7.425,0.04321798],[8.575,7.425,1],[9.425,7.425,1],[9.425,7.425,0.04321798],["NaN","NaN","NaN"],[1.425,8.425,0.04321798],[0.575,8.425,0.04321798],[0.575,8.425,0.428676],[1.425,8.425,0.428676],[1.425,8.425,0.04321798],["NaN","NaN","NaN"],[2.425,8.425,0.04321798],[1.575,8.425,0.04321798],[1.575,8.425,0.9105917],[2.425,8.425,0.9105917],[2.425,8.425,0.04321798],["NaN","NaN","NaN"],[3.425,8.425,0.04321798],[2.575,8.425,0.04321798],[2.575,8.425,0.9990419],[3.425,8.425,0.9990419],[3.425,8.425,0.04321798],["NaN","NaN","NaN"],[4.425,8.425,0.04321798],[3.575,8.425,0.04321798],[3.575,8.425,0.9995526],[4.425,8.425,0.9995526],[4.425,8.425,0.04321798],["NaN","NaN","NaN"],[5.425,8.425,0.04321798],[4.575,8.425,0.04321798],[4.575,8.425,1],[5.425,8.425,1],[5.425,8.425,0.04321798],["NaN","NaN","NaN"],[6.425,8.425,0.04321798],[5.575,8.425,0.04321798],[5.575,8.425,1],[6.425,8.425,1],[6.425,8.425,0.04321798],["NaN","NaN","NaN"],[7.425,8.425,0.04321798],[6.575,8.425,0.04321798],[6.575,8.425,1],[7.425,8.425,1],[7.425,8.425,0.04321798],["NaN","NaN","NaN"],[8.425,8.425,0.04321798],[7.575,8.425,0.04321798],[7.575,8.425,1],[8.425,8.425,1],[8.425,8.425,0.04321798],["NaN","NaN","NaN"],[9.425,8.425,0.04321798],[8.575,8.425,0.04321798],[8.575,8.425,1],[9.425,8.425,1],[9.425,8.425,0.04321798],["NaN","NaN","NaN"],[1.425,9.425,0.04321798],[0.575,9.425,0.04321798],[0.575,9.425,0.5885873],[1.425,9.425,0.5885873],[1.425,9.425,0.04321798],["NaN","NaN","NaN"],[2.425,9.425,0.04321798],[1.575,9.425,0.04321798],[1.575,9.425,0.981856],[2.425,9.425,0.981856],[2.425,9.425,0.04321798],["NaN","NaN","NaN"],[3.425,9.425,0.04321798],[2.575,9.425,0.04321798],[2.575,9.425,0.9999893],[3.425,9.425,0.9999893],[3.425,9.425,0.04321798],["NaN","NaN","NaN"],[4.425,9.425,0.04321798],[3.575,9.425,0.04321798],[3.575,9.425,0.9999969],[4.425,9.425,0.9999969],[4.425,9.425,0.04321798],["NaN","NaN","NaN"],[5.425,9.425,0.04321798],[4.575,9.425,0.04321798],[4.575,9.425,1],[5.425,9.425,1],[5.425,9.425,0.04321798],["NaN","NaN","NaN"],[6.425,9.425,0.04321798],[5.575,9.425,0.04321798],[5.575,9.425,1],[6.425,9.425,1],[6.425,9.425,0.04321798],["NaN","NaN","NaN"],[7.425,9.425,0.04321798],[6.575,9.425,0.04321798],[6.575,9.425,1],[7.425,9.425,1],[7.425,9.425,0.04321798],["NaN","NaN","NaN"],[8.425,9.425,0.04321798],[7.575,9.425,0.04321798],[7.575,9.425,1],[8.425,9.425,1],[8.425,9.425,0.04321798],["NaN","NaN","NaN"],[9.425,9.425,0.04321798],[8.575,9.425,0.04321798],[8.575,9.425,1],[9.425,9.425,1],[9.425,9.425,0.04321798],["NaN","NaN","NaN"],[1.425,10.425,0.04321798],[0.575,10.425,0.04321798],[0.575,10.425,0.7135577],[1.425,10.425,0.7135577],[1.425,10.425,0.04321798],["NaN","NaN","NaN"],[2.425,10.425,0.04321798],[1.575,10.425,0.04321798],[1.575,10.425,0.9967811],[2.425,10.425,0.9967811],[2.425,10.425,0.04321798],["NaN","NaN","NaN"],[3.425,10.425,0.04321798],[2.575,10.425,0.04321798],[2.575,10.425,0.9999999],[3.425,10.425,0.9999999],[3.425,10.425,0.04321798],["NaN","NaN","NaN"],[4.425,10.425,0.04321798],[3.575,10.425,0.04321798],[3.575,10.425,1],[4.425,10.425,1],[4.425,10.425,0.04321798],["NaN","NaN","NaN"],[5.425,10.425,0.04321798],[4.575,10.425,0.04321798],[4.575,10.425,1],[5.425,10.425,1],[5.425,10.425,0.04321798],["NaN","NaN","NaN"],[6.425,10.425,0.04321798],[5.575,10.425,0.04321798],[5.575,10.425,1],[6.425,10.425,1],[6.425,10.425,0.04321798],["NaN","NaN","NaN"],[7.425,10.425,0.04321798],[6.575,10.425,0.04321798],[6.575,10.425,1],[7.425,10.425,1],[7.425,10.425,0.04321798],["NaN","NaN","NaN"],[8.425,10.425,0.04321798],[7.575,10.425,0.04321798],[7.575,10.425,1],[8.425,10.425,1],[8.425,10.425,0.04321798],["NaN","NaN","NaN"],[9.425,10.425,0.04321798],[8.575,10.425,0.04321798],[8.575,10.425,1],[9.425,10.425,1],[9.425,10.425,0.04321798],["NaN","NaN","NaN"],[0.575,0.575,0.04321798],[0.575,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,0.575,0.04321798],[0.575,0.575,0.04321798],["NaN","NaN","NaN"],[1.575,0.575,0.04321798],[1.575,1.425,0.04321798],[1.575,1.425,0.06596308],[1.575,0.575,0.06596308],[1.575,0.575,0.04321798],["NaN","NaN","NaN"],[2.575,0.575,0.04321798],[2.575,1.425,0.04321798],[2.575,1.425,0.1021959],[2.575,0.575,0.1021959],[2.575,0.575,0.04321798],["NaN","NaN","NaN"],[3.575,0.575,0.04321798],[3.575,1.425,0.04321798],[3.575,1.425,0.1075265],[3.575,0.575,0.1075265],[3.575,0.575,0.04321798],["NaN","NaN","NaN"],[4.575,0.575,0.04321798],[4.575,1.425,0.04321798],[4.575,1.425,0.2316452],[4.575,0.575,0.2316452],[4.575,0.575,0.04321798],["NaN","NaN","NaN"],[5.575,0.575,0.04321798],[5.575,1.425,0.04321798],[5.575,1.425,0.2258095],[5.575,0.575,0.2258095],[5.575,0.575,0.04321798],["NaN","NaN","NaN"],[6.575,0.575,0.04321798],[6.575,1.425,0.04321798],[6.575,1.425,0.237572],[6.575,0.575,0.237572],[6.575,0.575,0.04321798],["NaN","NaN","NaN"],[7.575,0.575,0.04321798],[7.575,1.425,0.04321798],[7.575,1.425,0.3690261],[7.575,0.575,0.3690261],[7.575,0.575,0.04321798],["NaN","NaN","NaN"],[8.575,0.575,0.04321798],[8.575,1.425,0.04321798],[8.575,1.425,0.3592858],[8.575,0.575,0.3592858],[8.575,0.575,0.04321798],["NaN","NaN","NaN"],[0.575,1.575,0.04321798],[0.575,2.425,0.04321798],[0.575,2.425,0.05915068],[0.575,1.575,0.05915068],[0.575,1.575,0.04321798],["NaN","NaN","NaN"],[1.575,1.575,0.04321798],[1.575,2.425,0.04321798],[1.575,2.425,0.1108488],[1.575,1.575,0.1108488],[1.575,1.575,0.04321798],["NaN","NaN","NaN"],[2.575,1.575,0.04321798],[2.575,2.425,0.04321798],[2.575,2.425,0.2030372],[2.575,1.575,0.2030372],[2.575,1.575,0.04321798],["NaN","NaN","NaN"],[3.575,1.575,0.04321798],[3.575,2.425,0.04321798],[3.575,2.425,0.217041],[3.575,1.575,0.217041],[3.575,1.575,0.04321798],["NaN","NaN","NaN"],[4.575,1.575,0.04321798],[4.575,2.425,0.04321798],[4.575,2.425,0.5257641],[4.575,1.575,0.5257641],[4.575,1.575,0.04321798],["NaN","NaN","NaN"],[5.575,1.575,0.04321798],[5.575,2.425,0.04321798],[5.575,2.425,0.5127993],[5.575,1.575,0.5127993],[5.575,1.575,0.04321798],["NaN","NaN","NaN"],[6.575,1.575,0.04321798],[6.575,2.425,0.04321798],[6.575,2.425,0.5387281],[6.575,1.575,0.5387281],[6.575,1.575,0.04321798],["NaN","NaN","NaN"],[7.575,1.575,0.04321798],[7.575,2.425,0.04321798],[7.575,2.425,0.7713429],[7.575,1.575,0.7713429],[7.575,1.575,0.04321798],["NaN","NaN","NaN"],[8.575,1.575,0.04321798],[8.575,2.425,0.04321798],[8.575,2.425,0.7576808],[8.575,1.575,0.7576808],[8.575,1.575,0.04321798],["NaN","NaN","NaN"],[0.575,2.575,0.04321798],[0.575,3.425,0.04321798],[0.575,3.425,0.07725442],[0.575,2.575,0.07725442],[0.575,2.575,0.04321798],["NaN","NaN","NaN"],[1.575,2.575,0.04321798],[1.575,3.425,0.04321798],[1.575,3.425,0.1672267],[1.575,2.575,0.1672267],[1.575,2.575,0.04321798],["NaN","NaN","NaN"],[2.575,2.575,0.04321798],[2.575,3.425,0.04321798],[2.575,3.425,0.3312058],[2.575,2.575,0.3312058],[2.575,2.575,0.04321798],["NaN","NaN","NaN"],[3.575,2.575,0.04321798],[3.575,3.425,0.04321798],[3.575,3.425,0.3553697],[3.575,2.575,0.3553697],[3.575,2.575,0.04321798],["NaN","NaN","NaN"],[4.575,2.575,0.04321798],[4.575,3.425,0.04321798],[4.575,3.425,0.781293],[4.575,2.575,0.781293],[4.575,2.575,0.04321798],["NaN","NaN","NaN"],[5.575,2.575,0.04321798],[5.575,3.425,0.04321798],[5.575,3.425,0.76804],[5.575,2.575,0.76804],[5.575,2.575,0.04321798],["NaN","NaN","NaN"],[6.575,2.575,0.04321798],[6.575,3.425,0.04321798],[6.575,3.425,0.7941267],[6.575,2.575,0.7941267],[6.575,2.575,0.04321798],["NaN","NaN","NaN"],[7.575,2.575,0.04321798],[7.575,3.425,0.04321798],[7.575,3.425,0.9553382],[7.575,2.575,0.9553382],[7.575,2.575,0.04321798],["NaN","NaN","NaN"],[8.575,2.575,0.04321798],[8.575,3.425,0.04321798],[8.575,3.425,0.9493445],[8.575,2.575,0.9493445],[8.575,2.575,0.04321798],["NaN","NaN","NaN"],[0.575,3.575,0.04321798],[0.575,4.425,0.04321798],[0.575,4.425,0.1025677],[0.575,3.575,0.1025677],[0.575,3.575,0.04321798],["NaN","NaN","NaN"],[1.575,3.575,0.04321798],[1.575,4.425,0.04321798],[1.575,4.425,0.2501325],[1.575,3.575,0.2501325],[1.575,3.575,0.04321798],["NaN","NaN","NaN"],[2.575,3.575,0.04321798],[2.575,4.425,0.04321798],[2.575,4.425,0.5039053],[2.575,3.575,0.5039053],[2.575,3.575,0.04321798],["NaN","NaN","NaN"],[3.575,3.575,0.04321798],[3.575,4.425,0.04321798],[3.575,4.425,0.537781],[3.575,3.575,0.537781],[3.575,3.575,0.04321798],["NaN","NaN","NaN"],[4.575,3.575,0.04321798],[4.575,4.425,0.04321798],[4.575,4.425,0.9430385],[4.575,3.575,0.9430385],[4.575,3.575,0.04321798],["NaN","NaN","NaN"],[5.575,3.575,0.04321798],[5.575,4.425,0.04321798],[5.575,4.425,0.9361931],[5.575,3.575,0.9361931],[5.575,3.575,0.04321798],["NaN","NaN","NaN"],[6.575,3.575,0.04321798],[6.575,4.425,0.04321798],[6.575,4.425,0.9493053],[6.575,3.575,0.9493053],[6.575,3.575,0.04321798],["NaN","NaN","NaN"],[7.575,3.575,0.04321798],[7.575,4.425,0.04321798],[7.575,4.425,0.9972544],[7.575,3.575,0.9972544],[7.575,3.575,0.04321798],["NaN","NaN","NaN"],[8.575,3.575,0.04321798],[8.575,4.425,0.04321798],[8.575,4.425,0.9965154],[8.575,3.575,0.9965154],[8.575,3.575,0.04321798],["NaN","NaN","NaN"],[0.575,4.575,0.04321798],[0.575,5.425,0.04321798],[0.575,5.425,0.1352801],[0.575,4.575,0.1352801],[0.575,4.575,0.04321798],["NaN","NaN","NaN"],[1.575,4.575,0.04321798],[1.575,5.425,0.04321798],[1.575,5.425,0.3567418],[1.575,4.575,0.3567418],[1.575,4.575,0.04321798],["NaN","NaN","NaN"],[2.575,4.575,0.04321798],[2.575,5.425,0.04321798],[2.575,5.425,0.6844499],[2.575,4.575,0.6844499],[2.575,4.575,0.04321798],["NaN","NaN","NaN"],[3.575,4.575,0.04321798],[3.575,5.425,0.04321798],[3.575,5.425,0.7209224],[3.575,4.575,0.7209224],[3.575,4.575,0.04321798],["NaN","NaN","NaN"],[4.575,4.575,0.04321798],[4.575,5.425,0.04321798],[4.575,5.425,0.9923994],[4.575,4.575,0.9923994],[4.575,4.575,0.04321798],["NaN","NaN","NaN"],[5.575,4.575,0.04321798],[5.575,5.425,0.04321798],[5.575,5.425,0.9907784],[5.575,4.575,0.9907784],[5.575,4.575,0.04321798],["NaN","NaN","NaN"],[6.575,4.575,0.04321798],[6.575,5.425,0.04321798],[6.575,5.425,0.9937668],[6.575,4.575,0.9937668],[6.575,4.575,0.04321798],["NaN","NaN","NaN"],[7.575,4.575,0.04321798],[7.575,5.425,0.04321798],[7.575,5.425,0.999954],[7.575,4.575,0.999954],[7.575,4.575,0.04321798],["NaN","NaN","NaN"],[8.575,4.575,0.04321798],[8.575,5.425,0.04321798],[8.575,5.425,0.9999315],[8.575,4.575,0.9999315],[8.575,4.575,0.04321798],["NaN","NaN","NaN"],[0.575,5.575,0.04321798],[0.575,6.425,0.04321798],[0.575,6.425,0.1997644],[0.575,5.575,0.1997644],[0.575,5.575,0.04321798],["NaN","NaN","NaN"],[1.575,5.575,0.04321798],[1.575,6.425,0.04321798],[1.575,6.425,0.5454682],[1.575,5.575,0.5454682],[1.575,5.575,0.04321798],["NaN","NaN","NaN"],[2.575,5.575,0.04321798],[2.575,6.425,0.04321798],[2.575,6.425,0.8883711],[2.575,5.575,0.8883711],[2.575,5.575,0.04321798],["NaN","NaN","NaN"],[3.575,5.575,0.04321798],[3.575,6.425,0.04321798],[3.575,6.425,0.912325],[3.575,5.575,0.912325],[3.575,5.575,0.04321798],["NaN","NaN","NaN"],[4.575,5.575,0.04321798],[4.575,6.425,0.04321798],[4.575,6.425,0.9999133],[4.575,5.575,0.9999133],[4.575,5.575,0.04321798],["NaN","NaN","NaN"],[5.575,5.575,0.04321798],[5.575,6.425,0.04321798],[5.575,6.425,0.9998751],[5.575,5.575,0.9998751],[5.575,5.575,0.04321798],["NaN","NaN","NaN"],[6.575,5.575,0.04321798],[6.575,6.425,0.04321798],[6.575,6.425,0.9999403],[6.575,5.575,0.9999403],[6.575,5.575,0.04321798],["NaN","NaN","NaN"],[7.575,5.575,0.04321798],[7.575,6.425,0.04321798],[7.575,6.425,1],[7.575,5.575,1],[7.575,5.575,0.04321798],["NaN","NaN","NaN"],[8.575,5.575,0.04321798],[8.575,6.425,0.04321798],[8.575,6.425,1],[8.575,5.575,1],[8.575,5.575,0.04321798],["NaN","NaN","NaN"],[0.575,6.575,0.04321798],[0.575,7.425,0.04321798],[0.575,7.425,0.2790861],[0.575,6.575,0.2790861],[0.575,6.575,0.04321798],["NaN","NaN","NaN"],[1.575,6.575,0.04321798],[1.575,7.425,0.04321798],[1.575,7.425,0.722779],[1.575,6.575,0.722779],[1.575,6.575,0.04321798],["NaN","NaN","NaN"],[2.575,6.575,0.04321798],[2.575,7.425,0.04321798],[2.575,7.425,0.9743328],[2.575,6.575,0.9743328],[2.575,6.575,0.04321798],["NaN","NaN","NaN"],[3.575,6.575,0.04321798],[3.575,7.425,0.04321798],[3.575,7.425,0.9828939],[3.575,6.575,0.9828939],[3.575,6.575,0.04321798],["NaN","NaN","NaN"],[4.575,6.575,0.04321798],[4.575,7.425,0.04321798],[4.575,7.425,0.9999998],[4.575,6.575,0.9999998],[4.575,6.575,0.04321798],["NaN","NaN","NaN"],[5.575,6.575,0.04321798],[5.575,7.425,0.04321798],[5.575,7.425,0.9999996],[5.575,6.575,0.9999996],[5.575,6.575,0.04321798],["NaN","NaN","NaN"],[6.575,6.575,0.04321798],[6.575,7.425,0.04321798],[6.575,7.425,0.9999999],[6.575,6.575,0.9999999],[6.575,6.575,0.04321798],["NaN","NaN","NaN"],[7.575,6.575,0.04321798],[7.575,7.425,0.04321798],[7.575,7.425,1],[7.575,6.575,1],[7.575,6.575,0.04321798],["NaN","NaN","NaN"],[8.575,6.575,0.04321798],[8.575,7.425,0.04321798],[8.575,7.425,1],[8.575,6.575,1],[8.575,6.575,0.04321798],["NaN","NaN","NaN"],[0.575,7.575,0.04321798],[0.575,8.425,0.04321798],[0.575,8.425,0.428676],[0.575,7.575,0.428676],[0.575,7.575,0.04321798],["NaN","NaN","NaN"],[1.575,7.575,0.04321798],[1.575,8.425,0.04321798],[1.575,8.425,0.9105917],[1.575,7.575,0.9105917],[1.575,7.575,0.04321798],["NaN","NaN","NaN"],[2.575,7.575,0.04321798],[2.575,8.425,0.04321798],[2.575,8.425,0.9990419],[2.575,7.575,0.9990419],[2.575,7.575,0.04321798],["NaN","NaN","NaN"],[3.575,7.575,0.04321798],[3.575,8.425,0.04321798],[3.575,8.425,0.9995526],[3.575,7.575,0.9995526],[3.575,7.575,0.04321798],["NaN","NaN","NaN"],[4.575,7.575,0.04321798],[4.575,8.425,0.04321798],[4.575,8.425,1],[4.575,7.575,1],[4.575,7.575,0.04321798],["NaN","NaN","NaN"],[5.575,7.575,0.04321798],[5.575,8.425,0.04321798],[5.575,8.425,1],[5.575,7.575,1],[5.575,7.575,0.04321798],["NaN","NaN","NaN"],[6.575,7.575,0.04321798],[6.575,8.425,0.04321798],[6.575,8.425,1],[6.575,7.575,1],[6.575,7.575,0.04321798],["NaN","NaN","NaN"],[7.575,7.575,0.04321798],[7.575,8.425,0.04321798],[7.575,8.425,1],[7.575,7.575,1],[7.575,7.575,0.04321798],["NaN","NaN","NaN"],[8.575,7.575,0.04321798],[8.575,8.425,0.04321798],[8.575,8.425,1],[8.575,7.575,1],[8.575,7.575,0.04321798],["NaN","NaN","NaN"],[0.575,8.575,0.04321798],[0.575,9.425,0.04321798],[0.575,9.425,0.5885873],[0.575,8.575,0.5885873],[0.575,8.575,0.04321798],["NaN","NaN","NaN"],[1.575,8.575,0.04321798],[1.575,9.425,0.04321798],[1.575,9.425,0.981856],[1.575,8.575,0.981856],[1.575,8.575,0.04321798],["NaN","NaN","NaN"],[2.575,8.575,0.04321798],[2.575,9.425,0.04321798],[2.575,9.425,0.9999893],[2.575,8.575,0.9999893],[2.575,8.575,0.04321798],["NaN","NaN","NaN"],[3.575,8.575,0.04321798],[3.575,9.425,0.04321798],[3.575,9.425,0.9999969],[3.575,8.575,0.9999969],[3.575,8.575,0.04321798],["NaN","NaN","NaN"],[4.575,8.575,0.04321798],[4.575,9.425,0.04321798],[4.575,9.425,1],[4.575,8.575,1],[4.575,8.575,0.04321798],["NaN","NaN","NaN"],[5.575,8.575,0.04321798],[5.575,9.425,0.04321798],[5.575,9.425,1],[5.575,8.575,1],[5.575,8.575,0.04321798],["NaN","NaN","NaN"],[6.575,8.575,0.04321798],[6.575,9.425,0.04321798],[6.575,9.425,1],[6.575,8.575,1],[6.575,8.575,0.04321798],["NaN","NaN","NaN"],[7.575,8.575,0.04321798],[7.575,9.425,0.04321798],[7.575,9.425,1],[7.575,8.575,1],[7.575,8.575,0.04321798],["NaN","NaN","NaN"],[8.575,8.575,0.04321798],[8.575,9.425,0.04321798],[8.575,9.425,1],[8.575,8.575,1],[8.575,8.575,0.04321798],["NaN","NaN","NaN"],[0.575,9.575,0.04321798],[0.575,10.425,0.04321798],[0.575,10.425,0.7135577],[0.575,9.575,0.7135577],[0.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.575,9.575,0.04321798],[1.575,10.425,0.04321798],[1.575,10.425,0.9967811],[1.575,9.575,0.9967811],[1.575,9.575,0.04321798],["NaN","NaN","NaN"],[2.575,9.575,0.04321798],[2.575,10.425,0.04321798],[2.575,10.425,0.9999999],[2.575,9.575,0.9999999],[2.575,9.575,0.04321798],["NaN","NaN","NaN"],[3.575,9.575,0.04321798],[3.575,10.425,0.04321798],[3.575,10.425,1],[3.575,9.575,1],[3.575,9.575,0.04321798],["NaN","NaN","NaN"],[4.575,9.575,0.04321798],[4.575,10.425,0.04321798],[4.575,10.425,1],[4.575,9.575,1],[4.575,9.575,0.04321798],["NaN","NaN","NaN"],[5.575,9.575,0.04321798],[5.575,10.425,0.04321798],[5.575,10.425,1],[5.575,9.575,1],[5.575,9.575,0.04321798],["NaN","NaN","NaN"],[6.575,9.575,0.04321798],[6.575,10.425,0.04321798],[6.575,10.425,1],[6.575,9.575,1],[6.575,9.575,0.04321798],["NaN","NaN","NaN"],[7.575,9.575,0.04321798],[7.575,10.425,0.04321798],[7.575,10.425,1],[7.575,9.575,1],[7.575,9.575,0.04321798],["NaN","NaN","NaN"],[8.575,9.575,0.04321798],[8.575,10.425,0.04321798],[8.575,10.425,1],[8.575,9.575,1],[8.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.425,0.575,0.04321798],[1.425,1.425,0.04321798],[1.425,1.425,0.04321798],[1.425,0.575,0.04321798],[1.425,0.575,0.04321798],["NaN","NaN","NaN"],[2.425,0.575,0.04321798],[2.425,1.425,0.04321798],[2.425,1.425,0.06596308],[2.425,0.575,0.06596308],[2.425,0.575,0.04321798],["NaN","NaN","NaN"],[3.425,0.575,0.04321798],[3.425,1.425,0.04321798],[3.425,1.425,0.1021959],[3.425,0.575,0.1021959],[3.425,0.575,0.04321798],["NaN","NaN","NaN"],[4.425,0.575,0.04321798],[4.425,1.425,0.04321798],[4.425,1.425,0.1075265],[4.425,0.575,0.1075265],[4.425,0.575,0.04321798],["NaN","NaN","NaN"],[5.425,0.575,0.04321798],[5.425,1.425,0.04321798],[5.425,1.425,0.2316452],[5.425,0.575,0.2316452],[5.425,0.575,0.04321798],["NaN","NaN","NaN"],[6.425,0.575,0.04321798],[6.425,1.425,0.04321798],[6.425,1.425,0.2258095],[6.425,0.575,0.2258095],[6.425,0.575,0.04321798],["NaN","NaN","NaN"],[7.425,0.575,0.04321798],[7.425,1.425,0.04321798],[7.425,1.425,0.237572],[7.425,0.575,0.237572],[7.425,0.575,0.04321798],["NaN","NaN","NaN"],[8.425,0.575,0.04321798],[8.425,1.425,0.04321798],[8.425,1.425,0.3690261],[8.425,0.575,0.3690261],[8.425,0.575,0.04321798],["NaN","NaN","NaN"],[9.425,0.575,0.04321798],[9.425,1.425,0.04321798],[9.425,1.425,0.3592858],[9.425,0.575,0.3592858],[9.425,0.575,0.04321798],["NaN","NaN","NaN"],[1.425,1.575,0.04321798],[1.425,2.425,0.04321798],[1.425,2.425,0.05915068],[1.425,1.575,0.05915068],[1.425,1.575,0.04321798],["NaN","NaN","NaN"],[2.425,1.575,0.04321798],[2.425,2.425,0.04321798],[2.425,2.425,0.1108488],[2.425,1.575,0.1108488],[2.425,1.575,0.04321798],["NaN","NaN","NaN"],[3.425,1.575,0.04321798],[3.425,2.425,0.04321798],[3.425,2.425,0.2030372],[3.425,1.575,0.2030372],[3.425,1.575,0.04321798],["NaN","NaN","NaN"],[4.425,1.575,0.04321798],[4.425,2.425,0.04321798],[4.425,2.425,0.217041],[4.425,1.575,0.217041],[4.425,1.575,0.04321798],["NaN","NaN","NaN"],[5.425,1.575,0.04321798],[5.425,2.425,0.04321798],[5.425,2.425,0.5257641],[5.425,1.575,0.5257641],[5.425,1.575,0.04321798],["NaN","NaN","NaN"],[6.425,1.575,0.04321798],[6.425,2.425,0.04321798],[6.425,2.425,0.5127993],[6.425,1.575,0.5127993],[6.425,1.575,0.04321798],["NaN","NaN","NaN"],[7.425,1.575,0.04321798],[7.425,2.425,0.04321798],[7.425,2.425,0.5387281],[7.425,1.575,0.5387281],[7.425,1.575,0.04321798],["NaN","NaN","NaN"],[8.425,1.575,0.04321798],[8.425,2.425,0.04321798],[8.425,2.425,0.7713429],[8.425,1.575,0.7713429],[8.425,1.575,0.04321798],["NaN","NaN","NaN"],[9.425,1.575,0.04321798],[9.425,2.425,0.04321798],[9.425,2.425,0.7576808],[9.425,1.575,0.7576808],[9.425,1.575,0.04321798],["NaN","NaN","NaN"],[1.425,2.575,0.04321798],[1.425,3.425,0.04321798],[1.425,3.425,0.07725442],[1.425,2.575,0.07725442],[1.425,2.575,0.04321798],["NaN","NaN","NaN"],[2.425,2.575,0.04321798],[2.425,3.425,0.04321798],[2.425,3.425,0.1672267],[2.425,2.575,0.1672267],[2.425,2.575,0.04321798],["NaN","NaN","NaN"],[3.425,2.575,0.04321798],[3.425,3.425,0.04321798],[3.425,3.425,0.3312058],[3.425,2.575,0.3312058],[3.425,2.575,0.04321798],["NaN","NaN","NaN"],[4.425,2.575,0.04321798],[4.425,3.425,0.04321798],[4.425,3.425,0.3553697],[4.425,2.575,0.3553697],[4.425,2.575,0.04321798],["NaN","NaN","NaN"],[5.425,2.575,0.04321798],[5.425,3.425,0.04321798],[5.425,3.425,0.781293],[5.425,2.575,0.781293],[5.425,2.575,0.04321798],["NaN","NaN","NaN"],[6.425,2.575,0.04321798],[6.425,3.425,0.04321798],[6.425,3.425,0.76804],[6.425,2.575,0.76804],[6.425,2.575,0.04321798],["NaN","NaN","NaN"],[7.425,2.575,0.04321798],[7.425,3.425,0.04321798],[7.425,3.425,0.7941267],[7.425,2.575,0.7941267],[7.425,2.575,0.04321798],["NaN","NaN","NaN"],[8.425,2.575,0.04321798],[8.425,3.425,0.04321798],[8.425,3.425,0.9553382],[8.425,2.575,0.9553382],[8.425,2.575,0.04321798],["NaN","NaN","NaN"],[9.425,2.575,0.04321798],[9.425,3.425,0.04321798],[9.425,3.425,0.9493445],[9.425,2.575,0.9493445],[9.425,2.575,0.04321798],["NaN","NaN","NaN"],[1.425,3.575,0.04321798],[1.425,4.425,0.04321798],[1.425,4.425,0.1025677],[1.425,3.575,0.1025677],[1.425,3.575,0.04321798],["NaN","NaN","NaN"],[2.425,3.575,0.04321798],[2.425,4.425,0.04321798],[2.425,4.425,0.2501325],[2.425,3.575,0.2501325],[2.425,3.575,0.04321798],["NaN","NaN","NaN"],[3.425,3.575,0.04321798],[3.425,4.425,0.04321798],[3.425,4.425,0.5039053],[3.425,3.575,0.5039053],[3.425,3.575,0.04321798],["NaN","NaN","NaN"],[4.425,3.575,0.04321798],[4.425,4.425,0.04321798],[4.425,4.425,0.537781],[4.425,3.575,0.537781],[4.425,3.575,0.04321798],["NaN","NaN","NaN"],[5.425,3.575,0.04321798],[5.425,4.425,0.04321798],[5.425,4.425,0.9430385],[5.425,3.575,0.9430385],[5.425,3.575,0.04321798],["NaN","NaN","NaN"],[6.425,3.575,0.04321798],[6.425,4.425,0.04321798],[6.425,4.425,0.9361931],[6.425,3.575,0.9361931],[6.425,3.575,0.04321798],["NaN","NaN","NaN"],[7.425,3.575,0.04321798],[7.425,4.425,0.04321798],[7.425,4.425,0.9493053],[7.425,3.575,0.9493053],[7.425,3.575,0.04321798],["NaN","NaN","NaN"],[8.425,3.575,0.04321798],[8.425,4.425,0.04321798],[8.425,4.425,0.9972544],[8.425,3.575,0.9972544],[8.425,3.575,0.04321798],["NaN","NaN","NaN"],[9.425,3.575,0.04321798],[9.425,4.425,0.04321798],[9.425,4.425,0.9965154],[9.425,3.575,0.9965154],[9.425,3.575,0.04321798],["NaN","NaN","NaN"],[1.425,4.575,0.04321798],[1.425,5.425,0.04321798],[1.425,5.425,0.1352801],[1.425,4.575,0.1352801],[1.425,4.575,0.04321798],["NaN","NaN","NaN"],[2.425,4.575,0.04321798],[2.425,5.425,0.04321798],[2.425,5.425,0.3567418],[2.425,4.575,0.3567418],[2.425,4.575,0.04321798],["NaN","NaN","NaN"],[3.425,4.575,0.04321798],[3.425,5.425,0.04321798],[3.425,5.425,0.6844499],[3.425,4.575,0.6844499],[3.425,4.575,0.04321798],["NaN","NaN","NaN"],[4.425,4.575,0.04321798],[4.425,5.425,0.04321798],[4.425,5.425,0.7209224],[4.425,4.575,0.7209224],[4.425,4.575,0.04321798],["NaN","NaN","NaN"],[5.425,4.575,0.04321798],[5.425,5.425,0.04321798],[5.425,5.425,0.9923994],[5.425,4.575,0.9923994],[5.425,4.575,0.04321798],["NaN","NaN","NaN"],[6.425,4.575,0.04321798],[6.425,5.425,0.04321798],[6.425,5.425,0.9907784],[6.425,4.575,0.9907784],[6.425,4.575,0.04321798],["NaN","NaN","NaN"],[7.425,4.575,0.04321798],[7.425,5.425,0.04321798],[7.425,5.425,0.9937668],[7.425,4.575,0.9937668],[7.425,4.575,0.04321798],["NaN","NaN","NaN"],[8.425,4.575,0.04321798],[8.425,5.425,0.04321798],[8.425,5.425,0.999954],[8.425,4.575,0.999954],[8.425,4.575,0.04321798],["NaN","NaN","NaN"],[9.425,4.575,0.04321798],[9.425,5.425,0.04321798],[9.425,5.425,0.9999315],[9.425,4.575,0.9999315],[9.425,4.575,0.04321798],["NaN","NaN","NaN"],[1.425,5.575,0.04321798],[1.425,6.425,0.04321798],[1.425,6.425,0.1997644],[1.425,5.575,0.1997644],[1.425,5.575,0.04321798],["NaN","NaN","NaN"],[2.425,5.575,0.04321798],[2.425,6.425,0.04321798],[2.425,6.425,0.5454682],[2.425,5.575,0.5454682],[2.425,5.575,0.04321798],["NaN","NaN","NaN"],[3.425,5.575,0.04321798],[3.425,6.425,0.04321798],[3.425,6.425,0.8883711],[3.425,5.575,0.8883711],[3.425,5.575,0.04321798],["NaN","NaN","NaN"],[4.425,5.575,0.04321798],[4.425,6.425,0.04321798],[4.425,6.425,0.912325],[4.425,5.575,0.912325],[4.425,5.575,0.04321798],["NaN","NaN","NaN"],[5.425,5.575,0.04321798],[5.425,6.425,0.04321798],[5.425,6.425,0.9999133],[5.425,5.575,0.9999133],[5.425,5.575,0.04321798],["NaN","NaN","NaN"],[6.425,5.575,0.04321798],[6.425,6.425,0.04321798],[6.425,6.425,0.9998751],[6.425,5.575,0.9998751],[6.425,5.575,0.04321798],["NaN","NaN","NaN"],[7.425,5.575,0.04321798],[7.425,6.425,0.04321798],[7.425,6.425,0.9999403],[7.425,5.575,0.9999403],[7.425,5.575,0.04321798],["NaN","NaN","NaN"],[8.425,5.575,0.04321798],[8.425,6.425,0.04321798],[8.425,6.425,1],[8.425,5.575,1],[8.425,5.575,0.04321798],["NaN","NaN","NaN"],[9.425,5.575,0.04321798],[9.425,6.425,0.04321798],[9.425,6.425,1],[9.425,5.575,1],[9.425,5.575,0.04321798],["NaN","NaN","NaN"],[1.425,6.575,0.04321798],[1.425,7.425,0.04321798],[1.425,7.425,0.2790861],[1.425,6.575,0.2790861],[1.425,6.575,0.04321798],["NaN","NaN","NaN"],[2.425,6.575,0.04321798],[2.425,7.425,0.04321798],[2.425,7.425,0.722779],[2.425,6.575,0.722779],[2.425,6.575,0.04321798],["NaN","NaN","NaN"],[3.425,6.575,0.04321798],[3.425,7.425,0.04321798],[3.425,7.425,0.9743328],[3.425,6.575,0.9743328],[3.425,6.575,0.04321798],["NaN","NaN","NaN"],[4.425,6.575,0.04321798],[4.425,7.425,0.04321798],[4.425,7.425,0.9828939],[4.425,6.575,0.9828939],[4.425,6.575,0.04321798],["NaN","NaN","NaN"],[5.425,6.575,0.04321798],[5.425,7.425,0.04321798],[5.425,7.425,0.9999998],[5.425,6.575,0.9999998],[5.425,6.575,0.04321798],["NaN","NaN","NaN"],[6.425,6.575,0.04321798],[6.425,7.425,0.04321798],[6.425,7.425,0.9999996],[6.425,6.575,0.9999996],[6.425,6.575,0.04321798],["NaN","NaN","NaN"],[7.425,6.575,0.04321798],[7.425,7.425,0.04321798],[7.425,7.425,0.9999999],[7.425,6.575,0.9999999],[7.425,6.575,0.04321798],["NaN","NaN","NaN"],[8.425,6.575,0.04321798],[8.425,7.425,0.04321798],[8.425,7.425,1],[8.425,6.575,1],[8.425,6.575,0.04321798],["NaN","NaN","NaN"],[9.425,6.575,0.04321798],[9.425,7.425,0.04321798],[9.425,7.425,1],[9.425,6.575,1],[9.425,6.575,0.04321798],["NaN","NaN","NaN"],[1.425,7.575,0.04321798],[1.425,8.425,0.04321798],[1.425,8.425,0.428676],[1.425,7.575,0.428676],[1.425,7.575,0.04321798],["NaN","NaN","NaN"],[2.425,7.575,0.04321798],[2.425,8.425,0.04321798],[2.425,8.425,0.9105917],[2.425,7.575,0.9105917],[2.425,7.575,0.04321798],["NaN","NaN","NaN"],[3.425,7.575,0.04321798],[3.425,8.425,0.04321798],[3.425,8.425,0.9990419],[3.425,7.575,0.9990419],[3.425,7.575,0.04321798],["NaN","NaN","NaN"],[4.425,7.575,0.04321798],[4.425,8.425,0.04321798],[4.425,8.425,0.9995526],[4.425,7.575,0.9995526],[4.425,7.575,0.04321798],["NaN","NaN","NaN"],[5.425,7.575,0.04321798],[5.425,8.425,0.04321798],[5.425,8.425,1],[5.425,7.575,1],[5.425,7.575,0.04321798],["NaN","NaN","NaN"],[6.425,7.575,0.04321798],[6.425,8.425,0.04321798],[6.425,8.425,1],[6.425,7.575,1],[6.425,7.575,0.04321798],["NaN","NaN","NaN"],[7.425,7.575,0.04321798],[7.425,8.425,0.04321798],[7.425,8.425,1],[7.425,7.575,1],[7.425,7.575,0.04321798],["NaN","NaN","NaN"],[8.425,7.575,0.04321798],[8.425,8.425,0.04321798],[8.425,8.425,1],[8.425,7.575,1],[8.425,7.575,0.04321798],["NaN","NaN","NaN"],[9.425,7.575,0.04321798],[9.425,8.425,0.04321798],[9.425,8.425,1],[9.425,7.575,1],[9.425,7.575,0.04321798],["NaN","NaN","NaN"],[1.425,8.575,0.04321798],[1.425,9.425,0.04321798],[1.425,9.425,0.5885873],[1.425,8.575,0.5885873],[1.425,8.575,0.04321798],["NaN","NaN","NaN"],[2.425,8.575,0.04321798],[2.425,9.425,0.04321798],[2.425,9.425,0.981856],[2.425,8.575,0.981856],[2.425,8.575,0.04321798],["NaN","NaN","NaN"],[3.425,8.575,0.04321798],[3.425,9.425,0.04321798],[3.425,9.425,0.9999893],[3.425,8.575,0.9999893],[3.425,8.575,0.04321798],["NaN","NaN","NaN"],[4.425,8.575,0.04321798],[4.425,9.425,0.04321798],[4.425,9.425,0.9999969],[4.425,8.575,0.9999969],[4.425,8.575,0.04321798],["NaN","NaN","NaN"],[5.425,8.575,0.04321798],[5.425,9.425,0.04321798],[5.425,9.425,1],[5.425,8.575,1],[5.425,8.575,0.04321798],["NaN","NaN","NaN"],[6.425,8.575,0.04321798],[6.425,9.425,0.04321798],[6.425,9.425,1],[6.425,8.575,1],[6.425,8.575,0.04321798],["NaN","NaN","NaN"],[7.425,8.575,0.04321798],[7.425,9.425,0.04321798],[7.425,9.425,1],[7.425,8.575,1],[7.425,8.575,0.04321798],["NaN","NaN","NaN"],[8.425,8.575,0.04321798],[8.425,9.425,0.04321798],[8.425,9.425,1],[8.425,8.575,1],[8.425,8.575,0.04321798],["NaN","NaN","NaN"],[9.425,8.575,0.04321798],[9.425,9.425,0.04321798],[9.425,9.425,1],[9.425,8.575,1],[9.425,8.575,0.04321798],["NaN","NaN","NaN"],[1.425,9.575,0.04321798],[1.425,10.425,0.04321798],[1.425,10.425,0.7135577],[1.425,9.575,0.7135577],[1.425,9.575,0.04321798],["NaN","NaN","NaN"],[2.425,9.575,0.04321798],[2.425,10.425,0.04321798],[2.425,10.425,0.9967811],[2.425,9.575,0.9967811],[2.425,9.575,0.04321798],["NaN","NaN","NaN"],[3.425,9.575,0.04321798],[3.425,10.425,0.04321798],[3.425,10.425,0.9999999],[3.425,9.575,0.9999999],[3.425,9.575,0.04321798],["NaN","NaN","NaN"],[4.425,9.575,0.04321798],[4.425,10.425,0.04321798],[4.425,10.425,1],[4.425,9.575,1],[4.425,9.575,0.04321798],["NaN","NaN","NaN"],[5.425,9.575,0.04321798],[5.425,10.425,0.04321798],[5.425,10.425,1],[5.425,9.575,1],[5.425,9.575,0.04321798],["NaN","NaN","NaN"],[6.425,9.575,0.04321798],[6.425,10.425,0.04321798],[6.425,10.425,1],[6.425,9.575,1],[6.425,9.575,0.04321798],["NaN","NaN","NaN"],[7.425,9.575,0.04321798],[7.425,10.425,0.04321798],[7.425,10.425,1],[7.425,9.575,1],[7.425,9.575,0.04321798],["NaN","NaN","NaN"],[8.425,9.575,0.04321798],[8.425,10.425,0.04321798],[8.425,10.425,1],[8.425,9.575,1],[8.425,9.575,0.04321798],["NaN","NaN","NaN"],[9.425,9.575,0.04321798],[9.425,10.425,0.04321798],[9.425,10.425,1],[9.425,9.575,1],[9.425,9.575,0.04321798],["NaN","NaN","NaN"],[0.575,0.575,0.04321798],[1.425,0.575,0.04321798],[1.425,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,0.575,0.04321798],["NaN","NaN","NaN"],[1.575,0.575,0.06596308],[2.425,0.575,0.06596308],[2.425,1.425,0.06596308],[1.575,1.425,0.06596308],[1.575,0.575,0.06596308],["NaN","NaN","NaN"],[2.575,0.575,0.1021959],[3.425,0.575,0.1021959],[3.425,1.425,0.1021959],[2.575,1.425,0.1021959],[2.575,0.575,0.1021959],["NaN","NaN","NaN"],[3.575,0.575,0.1075265],[4.425,0.575,0.1075265],[4.425,1.425,0.1075265],[3.575,1.425,0.1075265],[3.575,0.575,0.1075265],["NaN","NaN","NaN"],[4.575,0.575,0.2316452],[5.425,0.575,0.2316452],[5.425,1.425,0.2316452],[4.575,1.425,0.2316452],[4.575,0.575,0.2316452],["NaN","NaN","NaN"],[5.575,0.575,0.2258095],[6.425,0.575,0.2258095],[6.425,1.425,0.2258095],[5.575,1.425,0.2258095],[5.575,0.575,0.2258095],["NaN","NaN","NaN"],[6.575,0.575,0.237572],[7.425,0.575,0.237572],[7.425,1.425,0.237572],[6.575,1.425,0.237572],[6.575,0.575,0.237572],["NaN","NaN","NaN"],[7.575,0.575,0.3690261],[8.425,0.575,0.3690261],[8.425,1.425,0.3690261],[7.575,1.425,0.3690261],[7.575,0.575,0.3690261],["NaN","NaN","NaN"],[8.575,0.575,0.3592858],[9.425,0.575,0.3592858],[9.425,1.425,0.3592858],[8.575,1.425,0.3592858],[8.575,0.575,0.3592858],["NaN","NaN","NaN"],[0.575,1.575,0.05915068],[1.425,1.575,0.05915068],[1.425,2.425,0.05915068],[0.575,2.425,0.05915068],[0.575,1.575,0.05915068],["NaN","NaN","NaN"],[1.575,1.575,0.1108488],[2.425,1.575,0.1108488],[2.425,2.425,0.1108488],[1.575,2.425,0.1108488],[1.575,1.575,0.1108488],["NaN","NaN","NaN"],[2.575,1.575,0.2030372],[3.425,1.575,0.2030372],[3.425,2.425,0.2030372],[2.575,2.425,0.2030372],[2.575,1.575,0.2030372],["NaN","NaN","NaN"],[3.575,1.575,0.217041],[4.425,1.575,0.217041],[4.425,2.425,0.217041],[3.575,2.425,0.217041],[3.575,1.575,0.217041],["NaN","NaN","NaN"],[4.575,1.575,0.5257641],[5.425,1.575,0.5257641],[5.425,2.425,0.5257641],[4.575,2.425,0.5257641],[4.575,1.575,0.5257641],["NaN","NaN","NaN"],[5.575,1.575,0.5127993],[6.425,1.575,0.5127993],[6.425,2.425,0.5127993],[5.575,2.425,0.5127993],[5.575,1.575,0.5127993],["NaN","NaN","NaN"],[6.575,1.575,0.5387281],[7.425,1.575,0.5387281],[7.425,2.425,0.5387281],[6.575,2.425,0.5387281],[6.575,1.575,0.5387281],["NaN","NaN","NaN"],[7.575,1.575,0.7713429],[8.425,1.575,0.7713429],[8.425,2.425,0.7713429],[7.575,2.425,0.7713429],[7.575,1.575,0.7713429],["NaN","NaN","NaN"],[8.575,1.575,0.7576808],[9.425,1.575,0.7576808],[9.425,2.425,0.7576808],[8.575,2.425,0.7576808],[8.575,1.575,0.7576808],["NaN","NaN","NaN"],[0.575,2.575,0.07725442],[1.425,2.575,0.07725442],[1.425,3.425,0.07725442],[0.575,3.425,0.07725442],[0.575,2.575,0.07725442],["NaN","NaN","NaN"],[1.575,2.575,0.1672267],[2.425,2.575,0.1672267],[2.425,3.425,0.1672267],[1.575,3.425,0.1672267],[1.575,2.575,0.1672267],["NaN","NaN","NaN"],[2.575,2.575,0.3312058],[3.425,2.575,0.3312058],[3.425,3.425,0.3312058],[2.575,3.425,0.3312058],[2.575,2.575,0.3312058],["NaN","NaN","NaN"],[3.575,2.575,0.3553697],[4.425,2.575,0.3553697],[4.425,3.425,0.3553697],[3.575,3.425,0.3553697],[3.575,2.575,0.3553697],["NaN","NaN","NaN"],[4.575,2.575,0.781293],[5.425,2.575,0.781293],[5.425,3.425,0.781293],[4.575,3.425,0.781293],[4.575,2.575,0.781293],["NaN","NaN","NaN"],[5.575,2.575,0.76804],[6.425,2.575,0.76804],[6.425,3.425,0.76804],[5.575,3.425,0.76804],[5.575,2.575,0.76804],["NaN","NaN","NaN"],[6.575,2.575,0.7941267],[7.425,2.575,0.7941267],[7.425,3.425,0.7941267],[6.575,3.425,0.7941267],[6.575,2.575,0.7941267],["NaN","NaN","NaN"],[7.575,2.575,0.9553382],[8.425,2.575,0.9553382],[8.425,3.425,0.9553382],[7.575,3.425,0.9553382],[7.575,2.575,0.9553382],["NaN","NaN","NaN"],[8.575,2.575,0.9493445],[9.425,2.575,0.9493445],[9.425,3.425,0.9493445],[8.575,3.425,0.9493445],[8.575,2.575,0.9493445],["NaN","NaN","NaN"],[0.575,3.575,0.1025677],[1.425,3.575,0.1025677],[1.425,4.425,0.1025677],[0.575,4.425,0.1025677],[0.575,3.575,0.1025677],["NaN","NaN","NaN"],[1.575,3.575,0.2501325],[2.425,3.575,0.2501325],[2.425,4.425,0.2501325],[1.575,4.425,0.2501325],[1.575,3.575,0.2501325],["NaN","NaN","NaN"],[2.575,3.575,0.5039053],[3.425,3.575,0.5039053],[3.425,4.425,0.5039053],[2.575,4.425,0.5039053],[2.575,3.575,0.5039053],["NaN","NaN","NaN"],[3.575,3.575,0.537781],[4.425,3.575,0.537781],[4.425,4.425,0.537781],[3.575,4.425,0.537781],[3.575,3.575,0.537781],["NaN","NaN","NaN"],[4.575,3.575,0.9430385],[5.425,3.575,0.9430385],[5.425,4.425,0.9430385],[4.575,4.425,0.9430385],[4.575,3.575,0.9430385],["NaN","NaN","NaN"],[5.575,3.575,0.9361931],[6.425,3.575,0.9361931],[6.425,4.425,0.9361931],[5.575,4.425,0.9361931],[5.575,3.575,0.9361931],["NaN","NaN","NaN"],[6.575,3.575,0.9493053],[7.425,3.575,0.9493053],[7.425,4.425,0.9493053],[6.575,4.425,0.9493053],[6.575,3.575,0.9493053],["NaN","NaN","NaN"],[7.575,3.575,0.9972544],[8.425,3.575,0.9972544],[8.425,4.425,0.9972544],[7.575,4.425,0.9972544],[7.575,3.575,0.9972544],["NaN","NaN","NaN"],[8.575,3.575,0.9965154],[9.425,3.575,0.9965154],[9.425,4.425,0.9965154],[8.575,4.425,0.9965154],[8.575,3.575,0.9965154],["NaN","NaN","NaN"],[0.575,4.575,0.1352801],[1.425,4.575,0.1352801],[1.425,5.425,0.1352801],[0.575,5.425,0.1352801],[0.575,4.575,0.1352801],["NaN","NaN","NaN"],[1.575,4.575,0.3567418],[2.425,4.575,0.3567418],[2.425,5.425,0.3567418],[1.575,5.425,0.3567418],[1.575,4.575,0.3567418],["NaN","NaN","NaN"],[2.575,4.575,0.6844499],[3.425,4.575,0.6844499],[3.425,5.425,0.6844499],[2.575,5.425,0.6844499],[2.575,4.575,0.6844499],["NaN","NaN","NaN"],[3.575,4.575,0.7209224],[4.425,4.575,0.7209224],[4.425,5.425,0.7209224],[3.575,5.425,0.7209224],[3.575,4.575,0.7209224],["NaN","NaN","NaN"],[4.575,4.575,0.9923994],[5.425,4.575,0.9923994],[5.425,5.425,0.9923994],[4.575,5.425,0.9923994],[4.575,4.575,0.9923994],["NaN","NaN","NaN"],[5.575,4.575,0.9907784],[6.425,4.575,0.9907784],[6.425,5.425,0.9907784],[5.575,5.425,0.9907784],[5.575,4.575,0.9907784],["NaN","NaN","NaN"],[6.575,4.575,0.9937668],[7.425,4.575,0.9937668],[7.425,5.425,0.9937668],[6.575,5.425,0.9937668],[6.575,4.575,0.9937668],["NaN","NaN","NaN"],[7.575,4.575,0.999954],[8.425,4.575,0.999954],[8.425,5.425,0.999954],[7.575,5.425,0.999954],[7.575,4.575,0.999954],["NaN","NaN","NaN"],[8.575,4.575,0.9999315],[9.425,4.575,0.9999315],[9.425,5.425,0.9999315],[8.575,5.425,0.9999315],[8.575,4.575,0.9999315],["NaN","NaN","NaN"],[0.575,5.575,0.1997644],[1.425,5.575,0.1997644],[1.425,6.425,0.1997644],[0.575,6.425,0.1997644],[0.575,5.575,0.1997644],["NaN","NaN","NaN"],[1.575,5.575,0.5454682],[2.425,5.575,0.5454682],[2.425,6.425,0.5454682],[1.575,6.425,0.5454682],[1.575,5.575,0.5454682],["NaN","NaN","NaN"],[2.575,5.575,0.8883711],[3.425,5.575,0.8883711],[3.425,6.425,0.8883711],[2.575,6.425,0.8883711],[2.575,5.575,0.8883711],["NaN","NaN","NaN"],[3.575,5.575,0.912325],[4.425,5.575,0.912325],[4.425,6.425,0.912325],[3.575,6.425,0.912325],[3.575,5.575,0.912325],["NaN","NaN","NaN"],[4.575,5.575,0.9999133],[5.425,5.575,0.9999133],[5.425,6.425,0.9999133],[4.575,6.425,0.9999133],[4.575,5.575,0.9999133],["NaN","NaN","NaN"],[5.575,5.575,0.9998751],[6.425,5.575,0.9998751],[6.425,6.425,0.9998751],[5.575,6.425,0.9998751],[5.575,5.575,0.9998751],["NaN","NaN","NaN"],[6.575,5.575,0.9999403],[7.425,5.575,0.9999403],[7.425,6.425,0.9999403],[6.575,6.425,0.9999403],[6.575,5.575,0.9999403],["NaN","NaN","NaN"],[7.575,5.575,1],[8.425,5.575,1],[8.425,6.425,1],[7.575,6.425,1],[7.575,5.575,1],["NaN","NaN","NaN"],[8.575,5.575,1],[9.425,5.575,1],[9.425,6.425,1],[8.575,6.425,1],[8.575,5.575,1],["NaN","NaN","NaN"],[0.575,6.575,0.2790861],[1.425,6.575,0.2790861],[1.425,7.425,0.2790861],[0.575,7.425,0.2790861],[0.575,6.575,0.2790861],["NaN","NaN","NaN"],[1.575,6.575,0.722779],[2.425,6.575,0.722779],[2.425,7.425,0.722779],[1.575,7.425,0.722779],[1.575,6.575,0.722779],["NaN","NaN","NaN"],[2.575,6.575,0.9743328],[3.425,6.575,0.9743328],[3.425,7.425,0.9743328],[2.575,7.425,0.9743328],[2.575,6.575,0.9743328],["NaN","NaN","NaN"],[3.575,6.575,0.9828939],[4.425,6.575,0.9828939],[4.425,7.425,0.9828939],[3.575,7.425,0.9828939],[3.575,6.575,0.9828939],["NaN","NaN","NaN"],[4.575,6.575,0.9999998],[5.425,6.575,0.9999998],[5.425,7.425,0.9999998],[4.575,7.425,0.9999998],[4.575,6.575,0.9999998],["NaN","NaN","NaN"],[5.575,6.575,0.9999996],[6.425,6.575,0.9999996],[6.425,7.425,0.9999996],[5.575,7.425,0.9999996],[5.575,6.575,0.9999996],["NaN","NaN","NaN"],[6.575,6.575,0.9999999],[7.425,6.575,0.9999999],[7.425,7.425,0.9999999],[6.575,7.425,0.9999999],[6.575,6.575,0.9999999],["NaN","NaN","NaN"],[7.575,6.575,1],[8.425,6.575,1],[8.425,7.425,1],[7.575,7.425,1],[7.575,6.575,1],["NaN","NaN","NaN"],[8.575,6.575,1],[9.425,6.575,1],[9.425,7.425,1],[8.575,7.425,1],[8.575,6.575,1],["NaN","NaN","NaN"],[0.575,7.575,0.428676],[1.425,7.575,0.428676],[1.425,8.425,0.428676],[0.575,8.425,0.428676],[0.575,7.575,0.428676],["NaN","NaN","NaN"],[1.575,7.575,0.9105917],[2.425,7.575,0.9105917],[2.425,8.425,0.9105917],[1.575,8.425,0.9105917],[1.575,7.575,0.9105917],["NaN","NaN","NaN"],[2.575,7.575,0.9990419],[3.425,7.575,0.9990419],[3.425,8.425,0.9990419],[2.575,8.425,0.9990419],[2.575,7.575,0.9990419],["NaN","NaN","NaN"],[3.575,7.575,0.9995526],[4.425,7.575,0.9995526],[4.425,8.425,0.9995526],[3.575,8.425,0.9995526],[3.575,7.575,0.9995526],["NaN","NaN","NaN"],[4.575,7.575,1],[5.425,7.575,1],[5.425,8.425,1],[4.575,8.425,1],[4.575,7.575,1],["NaN","NaN","NaN"],[5.575,7.575,1],[6.425,7.575,1],[6.425,8.425,1],[5.575,8.425,1],[5.575,7.575,1],["NaN","NaN","NaN"],[6.575,7.575,1],[7.425,7.575,1],[7.425,8.425,1],[6.575,8.425,1],[6.575,7.575,1],["NaN","NaN","NaN"],[7.575,7.575,1],[8.425,7.575,1],[8.425,8.425,1],[7.575,8.425,1],[7.575,7.575,1],["NaN","NaN","NaN"],[8.575,7.575,1],[9.425,7.575,1],[9.425,8.425,1],[8.575,8.425,1],[8.575,7.575,1],["NaN","NaN","NaN"],[0.575,8.575,0.5885873],[1.425,8.575,0.5885873],[1.425,9.425,0.5885873],[0.575,9.425,0.5885873],[0.575,8.575,0.5885873],["NaN","NaN","NaN"],[1.575,8.575,0.981856],[2.425,8.575,0.981856],[2.425,9.425,0.981856],[1.575,9.425,0.981856],[1.575,8.575,0.981856],["NaN","NaN","NaN"],[2.575,8.575,0.9999893],[3.425,8.575,0.9999893],[3.425,9.425,0.9999893],[2.575,9.425,0.9999893],[2.575,8.575,0.9999893],["NaN","NaN","NaN"],[3.575,8.575,0.9999969],[4.425,8.575,0.9999969],[4.425,9.425,0.9999969],[3.575,9.425,0.9999969],[3.575,8.575,0.9999969],["NaN","NaN","NaN"],[4.575,8.575,1],[5.425,8.575,1],[5.425,9.425,1],[4.575,9.425,1],[4.575,8.575,1],["NaN","NaN","NaN"],[5.575,8.575,1],[6.425,8.575,1],[6.425,9.425,1],[5.575,9.425,1],[5.575,8.575,1],["NaN","NaN","NaN"],[6.575,8.575,1],[7.425,8.575,1],[7.425,9.425,1],[6.575,9.425,1],[6.575,8.575,1],["NaN","NaN","NaN"],[7.575,8.575,1],[8.425,8.575,1],[8.425,9.425,1],[7.575,9.425,1],[7.575,8.575,1],["NaN","NaN","NaN"],[8.575,8.575,1],[9.425,8.575,1],[9.425,9.425,1],[8.575,9.425,1],[8.575,8.575,1],["NaN","NaN","NaN"],[0.575,9.575,0.7135577],[1.425,9.575,0.7135577],[1.425,10.425,0.7135577],[0.575,10.425,0.7135577],[0.575,9.575,0.7135577],["NaN","NaN","NaN"],[1.575,9.575,0.9967811],[2.425,9.575,0.9967811],[2.425,10.425,0.9967811],[1.575,10.425,0.9967811],[1.575,9.575,0.9967811],["NaN","NaN","NaN"],[2.575,9.575,0.9999999],[3.425,9.575,0.9999999],[3.425,10.425,0.9999999],[2.575,10.425,0.9999999],[2.575,9.575,0.9999999],["NaN","NaN","NaN"],[3.575,9.575,1],[4.425,9.575,1],[4.425,10.425,1],[3.575,10.425,1],[3.575,9.575,1],["NaN","NaN","NaN"],[4.575,9.575,1],[5.425,9.575,1],[5.425,10.425,1],[4.575,10.425,1],[4.575,9.575,1],["NaN","NaN","NaN"],[5.575,9.575,1],[6.425,9.575,1],[6.425,10.425,1],[5.575,10.425,1],[5.575,9.575,1],["NaN","NaN","NaN"],[6.575,9.575,1],[7.425,9.575,1],[7.425,10.425,1],[6.575,10.425,1],[6.575,9.575,1],["NaN","NaN","NaN"],[7.575,9.575,1],[8.425,9.575,1],[8.425,10.425,1],[7.575,10.425,1],[7.575,9.575,1],["NaN","NaN","NaN"],[8.575,9.575,1],[9.425,9.575,1],[9.425,10.425,1],[8.575,10.425,1],[8.575,9.575,1],["NaN","NaN","NaN"]],"colors":[[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1]],"centers":[[0.575,0.575,0.04321798],[1.425,0.575,0.04321798],[1.425,0.575,0.04321798],[0.575,0.575,0.04321798],[0.575,0.575,0.04321798],["NaN","NaN","NaN"],[1.575,0.575,0.04321798],[2.425,0.575,0.04321798],[2.425,0.575,0.06596308],[1.575,0.575,0.06596308],[1.575,0.575,0.04321798],["NaN","NaN","NaN"],[2.575,0.575,0.04321798],[3.425,0.575,0.04321798],[3.425,0.575,0.1021959],[2.575,0.575,0.1021959],[2.575,0.575,0.04321798],["NaN","NaN","NaN"],[3.575,0.575,0.04321798],[4.425,0.575,0.04321798],[4.425,0.575,0.1075265],[3.575,0.575,0.1075265],[3.575,0.575,0.04321798],["NaN","NaN","NaN"],[4.575,0.575,0.04321798],[5.425,0.575,0.04321798],[5.425,0.575,0.2316452],[4.575,0.575,0.2316452],[4.575,0.575,0.04321798],["NaN","NaN","NaN"],[5.575,0.575,0.04321798],[6.425,0.575,0.04321798],[6.425,0.575,0.2258095],[5.575,0.575,0.2258095],[5.575,0.575,0.04321798],["NaN","NaN","NaN"],[6.575,0.575,0.04321798],[7.425,0.575,0.04321798],[7.425,0.575,0.237572],[6.575,0.575,0.237572],[6.575,0.575,0.04321798],["NaN","NaN","NaN"],[7.575,0.575,0.04321798],[8.425,0.575,0.04321798],[8.425,0.575,0.3690261],[7.575,0.575,0.3690261],[7.575,0.575,0.04321798],["NaN","NaN","NaN"],[8.575,0.575,0.04321798],[9.425,0.575,0.04321798],[9.425,0.575,0.3592858],[8.575,0.575,0.3592858],[8.575,0.575,0.04321798],["NaN","NaN","NaN"],[0.575,1.575,0.04321798],[1.425,1.575,0.04321798],[1.425,1.575,0.05915068],[0.575,1.575,0.05915068],[0.575,1.575,0.04321798],["NaN","NaN","NaN"],[1.575,1.575,0.04321798],[2.425,1.575,0.04321798],[2.425,1.575,0.1108488],[1.575,1.575,0.1108488],[1.575,1.575,0.04321798],["NaN","NaN","NaN"],[2.575,1.575,0.04321798],[3.425,1.575,0.04321798],[3.425,1.575,0.2030372],[2.575,1.575,0.2030372],[2.575,1.575,0.04321798],["NaN","NaN","NaN"],[3.575,1.575,0.04321798],[4.425,1.575,0.04321798],[4.425,1.575,0.217041],[3.575,1.575,0.217041],[3.575,1.575,0.04321798],["NaN","NaN","NaN"],[4.575,1.575,0.04321798],[5.425,1.575,0.04321798],[5.425,1.575,0.5257641],[4.575,1.575,0.5257641],[4.575,1.575,0.04321798],["NaN","NaN","NaN"],[5.575,1.575,0.04321798],[6.425,1.575,0.04321798],[6.425,1.575,0.5127993],[5.575,1.575,0.5127993],[5.575,1.575,0.04321798],["NaN","NaN","NaN"],[6.575,1.575,0.04321798],[7.425,1.575,0.04321798],[7.425,1.575,0.5387281],[6.575,1.575,0.5387281],[6.575,1.575,0.04321798],["NaN","NaN","NaN"],[7.575,1.575,0.04321798],[8.425,1.575,0.04321798],[8.425,1.575,0.7713429],[7.575,1.575,0.7713429],[7.575,1.575,0.04321798],["NaN","NaN","NaN"],[8.575,1.575,0.04321798],[9.425,1.575,0.04321798],[9.425,1.575,0.7576808],[8.575,1.575,0.7576808],[8.575,1.575,0.04321798],["NaN","NaN","NaN"],[0.575,2.575,0.04321798],[1.425,2.575,0.04321798],[1.425,2.575,0.07725442],[0.575,2.575,0.07725442],[0.575,2.575,0.04321798],["NaN","NaN","NaN"],[1.575,2.575,0.04321798],[2.425,2.575,0.04321798],[2.425,2.575,0.1672267],[1.575,2.575,0.1672267],[1.575,2.575,0.04321798],["NaN","NaN","NaN"],[2.575,2.575,0.04321798],[3.425,2.575,0.04321798],[3.425,2.575,0.3312058],[2.575,2.575,0.3312058],[2.575,2.575,0.04321798],["NaN","NaN","NaN"],[3.575,2.575,0.04321798],[4.425,2.575,0.04321798],[4.425,2.575,0.3553697],[3.575,2.575,0.3553697],[3.575,2.575,0.04321798],["NaN","NaN","NaN"],[4.575,2.575,0.04321798],[5.425,2.575,0.04321798],[5.425,2.575,0.781293],[4.575,2.575,0.781293],[4.575,2.575,0.04321798],["NaN","NaN","NaN"],[5.575,2.575,0.04321798],[6.425,2.575,0.04321798],[6.425,2.575,0.76804],[5.575,2.575,0.76804],[5.575,2.575,0.04321798],["NaN","NaN","NaN"],[6.575,2.575,0.04321798],[7.425,2.575,0.04321798],[7.425,2.575,0.7941267],[6.575,2.575,0.7941267],[6.575,2.575,0.04321798],["NaN","NaN","NaN"],[7.575,2.575,0.04321798],[8.425,2.575,0.04321798],[8.425,2.575,0.9553382],[7.575,2.575,0.9553382],[7.575,2.575,0.04321798],["NaN","NaN","NaN"],[8.575,2.575,0.04321798],[9.425,2.575,0.04321798],[9.425,2.575,0.9493445],[8.575,2.575,0.9493445],[8.575,2.575,0.04321798],["NaN","NaN","NaN"],[0.575,3.575,0.04321798],[1.425,3.575,0.04321798],[1.425,3.575,0.1025677],[0.575,3.575,0.1025677],[0.575,3.575,0.04321798],["NaN","NaN","NaN"],[1.575,3.575,0.04321798],[2.425,3.575,0.04321798],[2.425,3.575,0.2501325],[1.575,3.575,0.2501325],[1.575,3.575,0.04321798],["NaN","NaN","NaN"],[2.575,3.575,0.04321798],[3.425,3.575,0.04321798],[3.425,3.575,0.5039053],[2.575,3.575,0.5039053],[2.575,3.575,0.04321798],["NaN","NaN","NaN"],[3.575,3.575,0.04321798],[4.425,3.575,0.04321798],[4.425,3.575,0.537781],[3.575,3.575,0.537781],[3.575,3.575,0.04321798],["NaN","NaN","NaN"],[4.575,3.575,0.04321798],[5.425,3.575,0.04321798],[5.425,3.575,0.9430385],[4.575,3.575,0.9430385],[4.575,3.575,0.04321798],["NaN","NaN","NaN"],[5.575,3.575,0.04321798],[6.425,3.575,0.04321798],[6.425,3.575,0.9361931],[5.575,3.575,0.9361931],[5.575,3.575,0.04321798],["NaN","NaN","NaN"],[6.575,3.575,0.04321798],[7.425,3.575,0.04321798],[7.425,3.575,0.9493053],[6.575,3.575,0.9493053],[6.575,3.575,0.04321798],["NaN","NaN","NaN"],[7.575,3.575,0.04321798],[8.425,3.575,0.04321798],[8.425,3.575,0.9972544],[7.575,3.575,0.9972544],[7.575,3.575,0.04321798],["NaN","NaN","NaN"],[8.575,3.575,0.04321798],[9.425,3.575,0.04321798],[9.425,3.575,0.9965154],[8.575,3.575,0.9965154],[8.575,3.575,0.04321798],["NaN","NaN","NaN"],[0.575,4.575,0.04321798],[1.425,4.575,0.04321798],[1.425,4.575,0.1352801],[0.575,4.575,0.1352801],[0.575,4.575,0.04321798],["NaN","NaN","NaN"],[1.575,4.575,0.04321798],[2.425,4.575,0.04321798],[2.425,4.575,0.3567418],[1.575,4.575,0.3567418],[1.575,4.575,0.04321798],["NaN","NaN","NaN"],[2.575,4.575,0.04321798],[3.425,4.575,0.04321798],[3.425,4.575,0.6844499],[2.575,4.575,0.6844499],[2.575,4.575,0.04321798],["NaN","NaN","NaN"],[3.575,4.575,0.04321798],[4.425,4.575,0.04321798],[4.425,4.575,0.7209224],[3.575,4.575,0.7209224],[3.575,4.575,0.04321798],["NaN","NaN","NaN"],[4.575,4.575,0.04321798],[5.425,4.575,0.04321798],[5.425,4.575,0.9923994],[4.575,4.575,0.9923994],[4.575,4.575,0.04321798],["NaN","NaN","NaN"],[5.575,4.575,0.04321798],[6.425,4.575,0.04321798],[6.425,4.575,0.9907784],[5.575,4.575,0.9907784],[5.575,4.575,0.04321798],["NaN","NaN","NaN"],[6.575,4.575,0.04321798],[7.425,4.575,0.04321798],[7.425,4.575,0.9937668],[6.575,4.575,0.9937668],[6.575,4.575,0.04321798],["NaN","NaN","NaN"],[7.575,4.575,0.04321798],[8.425,4.575,0.04321798],[8.425,4.575,0.999954],[7.575,4.575,0.999954],[7.575,4.575,0.04321798],["NaN","NaN","NaN"],[8.575,4.575,0.04321798],[9.425,4.575,0.04321798],[9.425,4.575,0.9999315],[8.575,4.575,0.9999315],[8.575,4.575,0.04321798],["NaN","NaN","NaN"],[0.575,5.575,0.04321798],[1.425,5.575,0.04321798],[1.425,5.575,0.1997644],[0.575,5.575,0.1997644],[0.575,5.575,0.04321798],["NaN","NaN","NaN"],[1.575,5.575,0.04321798],[2.425,5.575,0.04321798],[2.425,5.575,0.5454682],[1.575,5.575,0.5454682],[1.575,5.575,0.04321798],["NaN","NaN","NaN"],[2.575,5.575,0.04321798],[3.425,5.575,0.04321798],[3.425,5.575,0.8883711],[2.575,5.575,0.8883711],[2.575,5.575,0.04321798],["NaN","NaN","NaN"],[3.575,5.575,0.04321798],[4.425,5.575,0.04321798],[4.425,5.575,0.912325],[3.575,5.575,0.912325],[3.575,5.575,0.04321798],["NaN","NaN","NaN"],[4.575,5.575,0.04321798],[5.425,5.575,0.04321798],[5.425,5.575,0.9999133],[4.575,5.575,0.9999133],[4.575,5.575,0.04321798],["NaN","NaN","NaN"],[5.575,5.575,0.04321798],[6.425,5.575,0.04321798],[6.425,5.575,0.9998751],[5.575,5.575,0.9998751],[5.575,5.575,0.04321798],["NaN","NaN","NaN"],[6.575,5.575,0.04321798],[7.425,5.575,0.04321798],[7.425,5.575,0.9999403],[6.575,5.575,0.9999403],[6.575,5.575,0.04321798],["NaN","NaN","NaN"],[7.575,5.575,0.04321798],[8.425,5.575,0.04321798],[8.425,5.575,1],[7.575,5.575,1],[7.575,5.575,0.04321798],["NaN","NaN","NaN"],[8.575,5.575,0.04321798],[9.425,5.575,0.04321798],[9.425,5.575,1],[8.575,5.575,1],[8.575,5.575,0.04321798],["NaN","NaN","NaN"],[0.575,6.575,0.04321798],[1.425,6.575,0.04321798],[1.425,6.575,0.2790861],[0.575,6.575,0.2790861],[0.575,6.575,0.04321798],["NaN","NaN","NaN"],[1.575,6.575,0.04321798],[2.425,6.575,0.04321798],[2.425,6.575,0.722779],[1.575,6.575,0.722779],[1.575,6.575,0.04321798],["NaN","NaN","NaN"],[2.575,6.575,0.04321798],[3.425,6.575,0.04321798],[3.425,6.575,0.9743328],[2.575,6.575,0.9743328],[2.575,6.575,0.04321798],["NaN","NaN","NaN"],[3.575,6.575,0.04321798],[4.425,6.575,0.04321798],[4.425,6.575,0.9828939],[3.575,6.575,0.9828939],[3.575,6.575,0.04321798],["NaN","NaN","NaN"],[4.575,6.575,0.04321798],[5.425,6.575,0.04321798],[5.425,6.575,0.9999998],[4.575,6.575,0.9999998],[4.575,6.575,0.04321798],["NaN","NaN","NaN"],[5.575,6.575,0.04321798],[6.425,6.575,0.04321798],[6.425,6.575,0.9999996],[5.575,6.575,0.9999996],[5.575,6.575,0.04321798],["NaN","NaN","NaN"],[6.575,6.575,0.04321798],[7.425,6.575,0.04321798],[7.425,6.575,0.9999999],[6.575,6.575,0.9999999],[6.575,6.575,0.04321798],["NaN","NaN","NaN"],[7.575,6.575,0.04321798],[8.425,6.575,0.04321798],[8.425,6.575,1],[7.575,6.575,1],[7.575,6.575,0.04321798],["NaN","NaN","NaN"],[8.575,6.575,0.04321798],[9.425,6.575,0.04321798],[9.425,6.575,1],[8.575,6.575,1],[8.575,6.575,0.04321798],["NaN","NaN","NaN"],[0.575,7.575,0.04321798],[1.425,7.575,0.04321798],[1.425,7.575,0.428676],[0.575,7.575,0.428676],[0.575,7.575,0.04321798],["NaN","NaN","NaN"],[1.575,7.575,0.04321798],[2.425,7.575,0.04321798],[2.425,7.575,0.9105917],[1.575,7.575,0.9105917],[1.575,7.575,0.04321798],["NaN","NaN","NaN"],[2.575,7.575,0.04321798],[3.425,7.575,0.04321798],[3.425,7.575,0.9990419],[2.575,7.575,0.9990419],[2.575,7.575,0.04321798],["NaN","NaN","NaN"],[3.575,7.575,0.04321798],[4.425,7.575,0.04321798],[4.425,7.575,0.9995526],[3.575,7.575,0.9995526],[3.575,7.575,0.04321798],["NaN","NaN","NaN"],[4.575,7.575,0.04321798],[5.425,7.575,0.04321798],[5.425,7.575,1],[4.575,7.575,1],[4.575,7.575,0.04321798],["NaN","NaN","NaN"],[5.575,7.575,0.04321798],[6.425,7.575,0.04321798],[6.425,7.575,1],[5.575,7.575,1],[5.575,7.575,0.04321798],["NaN","NaN","NaN"],[6.575,7.575,0.04321798],[7.425,7.575,0.04321798],[7.425,7.575,1],[6.575,7.575,1],[6.575,7.575,0.04321798],["NaN","NaN","NaN"],[7.575,7.575,0.04321798],[8.425,7.575,0.04321798],[8.425,7.575,1],[7.575,7.575,1],[7.575,7.575,0.04321798],["NaN","NaN","NaN"],[8.575,7.575,0.04321798],[9.425,7.575,0.04321798],[9.425,7.575,1],[8.575,7.575,1],[8.575,7.575,0.04321798],["NaN","NaN","NaN"],[0.575,8.575,0.04321798],[1.425,8.575,0.04321798],[1.425,8.575,0.5885873],[0.575,8.575,0.5885873],[0.575,8.575,0.04321798],["NaN","NaN","NaN"],[1.575,8.575,0.04321798],[2.425,8.575,0.04321798],[2.425,8.575,0.981856],[1.575,8.575,0.981856],[1.575,8.575,0.04321798],["NaN","NaN","NaN"],[2.575,8.575,0.04321798],[3.425,8.575,0.04321798],[3.425,8.575,0.9999893],[2.575,8.575,0.9999893],[2.575,8.575,0.04321798],["NaN","NaN","NaN"],[3.575,8.575,0.04321798],[4.425,8.575,0.04321798],[4.425,8.575,0.9999969],[3.575,8.575,0.9999969],[3.575,8.575,0.04321798],["NaN","NaN","NaN"],[4.575,8.575,0.04321798],[5.425,8.575,0.04321798],[5.425,8.575,1],[4.575,8.575,1],[4.575,8.575,0.04321798],["NaN","NaN","NaN"],[5.575,8.575,0.04321798],[6.425,8.575,0.04321798],[6.425,8.575,1],[5.575,8.575,1],[5.575,8.575,0.04321798],["NaN","NaN","NaN"],[6.575,8.575,0.04321798],[7.425,8.575,0.04321798],[7.425,8.575,1],[6.575,8.575,1],[6.575,8.575,0.04321798],["NaN","NaN","NaN"],[7.575,8.575,0.04321798],[8.425,8.575,0.04321798],[8.425,8.575,1],[7.575,8.575,1],[7.575,8.575,0.04321798],["NaN","NaN","NaN"],[8.575,8.575,0.04321798],[9.425,8.575,0.04321798],[9.425,8.575,1],[8.575,8.575,1],[8.575,8.575,0.04321798],["NaN","NaN","NaN"],[0.575,9.575,0.04321798],[1.425,9.575,0.04321798],[1.425,9.575,0.7135577],[0.575,9.575,0.7135577],[0.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.575,9.575,0.04321798],[2.425,9.575,0.04321798],[2.425,9.575,0.9967811],[1.575,9.575,0.9967811],[1.575,9.575,0.04321798],["NaN","NaN","NaN"],[2.575,9.575,0.04321798],[3.425,9.575,0.04321798],[3.425,9.575,0.9999999],[2.575,9.575,0.9999999],[2.575,9.575,0.04321798],["NaN","NaN","NaN"],[3.575,9.575,0.04321798],[4.425,9.575,0.04321798],[4.425,9.575,1],[3.575,9.575,1],[3.575,9.575,0.04321798],["NaN","NaN","NaN"],[4.575,9.575,0.04321798],[5.425,9.575,0.04321798],[5.425,9.575,1],[4.575,9.575,1],[4.575,9.575,0.04321798],["NaN","NaN","NaN"],[5.575,9.575,0.04321798],[6.425,9.575,0.04321798],[6.425,9.575,1],[5.575,9.575,1],[5.575,9.575,0.04321798],["NaN","NaN","NaN"],[6.575,9.575,0.04321798],[7.425,9.575,0.04321798],[7.425,9.575,1],[6.575,9.575,1],[6.575,9.575,0.04321798],["NaN","NaN","NaN"],[7.575,9.575,0.04321798],[8.425,9.575,0.04321798],[8.425,9.575,1],[7.575,9.575,1],[7.575,9.575,0.04321798],["NaN","NaN","NaN"],[8.575,9.575,0.04321798],[9.425,9.575,0.04321798],[9.425,9.575,1],[8.575,9.575,1],[8.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.425,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,1.425,0.04321798],[1.425,1.425,0.04321798],[1.425,1.425,0.04321798],["NaN","NaN","NaN"],[2.425,1.425,0.04321798],[1.575,1.425,0.04321798],[1.575,1.425,0.06596308],[2.425,1.425,0.06596308],[2.425,1.425,0.04321798],["NaN","NaN","NaN"],[3.425,1.425,0.04321798],[2.575,1.425,0.04321798],[2.575,1.425,0.1021959],[3.425,1.425,0.1021959],[3.425,1.425,0.04321798],["NaN","NaN","NaN"],[4.425,1.425,0.04321798],[3.575,1.425,0.04321798],[3.575,1.425,0.1075265],[4.425,1.425,0.1075265],[4.425,1.425,0.04321798],["NaN","NaN","NaN"],[5.425,1.425,0.04321798],[4.575,1.425,0.04321798],[4.575,1.425,0.2316452],[5.425,1.425,0.2316452],[5.425,1.425,0.04321798],["NaN","NaN","NaN"],[6.425,1.425,0.04321798],[5.575,1.425,0.04321798],[5.575,1.425,0.2258095],[6.425,1.425,0.2258095],[6.425,1.425,0.04321798],["NaN","NaN","NaN"],[7.425,1.425,0.04321798],[6.575,1.425,0.04321798],[6.575,1.425,0.237572],[7.425,1.425,0.237572],[7.425,1.425,0.04321798],["NaN","NaN","NaN"],[8.425,1.425,0.04321798],[7.575,1.425,0.04321798],[7.575,1.425,0.3690261],[8.425,1.425,0.3690261],[8.425,1.425,0.04321798],["NaN","NaN","NaN"],[9.425,1.425,0.04321798],[8.575,1.425,0.04321798],[8.575,1.425,0.3592858],[9.425,1.425,0.3592858],[9.425,1.425,0.04321798],["NaN","NaN","NaN"],[1.425,2.425,0.04321798],[0.575,2.425,0.04321798],[0.575,2.425,0.05915068],[1.425,2.425,0.05915068],[1.425,2.425,0.04321798],["NaN","NaN","NaN"],[2.425,2.425,0.04321798],[1.575,2.425,0.04321798],[1.575,2.425,0.1108488],[2.425,2.425,0.1108488],[2.425,2.425,0.04321798],["NaN","NaN","NaN"],[3.425,2.425,0.04321798],[2.575,2.425,0.04321798],[2.575,2.425,0.2030372],[3.425,2.425,0.2030372],[3.425,2.425,0.04321798],["NaN","NaN","NaN"],[4.425,2.425,0.04321798],[3.575,2.425,0.04321798],[3.575,2.425,0.217041],[4.425,2.425,0.217041],[4.425,2.425,0.04321798],["NaN","NaN","NaN"],[5.425,2.425,0.04321798],[4.575,2.425,0.04321798],[4.575,2.425,0.5257641],[5.425,2.425,0.5257641],[5.425,2.425,0.04321798],["NaN","NaN","NaN"],[6.425,2.425,0.04321798],[5.575,2.425,0.04321798],[5.575,2.425,0.5127993],[6.425,2.425,0.5127993],[6.425,2.425,0.04321798],["NaN","NaN","NaN"],[7.425,2.425,0.04321798],[6.575,2.425,0.04321798],[6.575,2.425,0.5387281],[7.425,2.425,0.5387281],[7.425,2.425,0.04321798],["NaN","NaN","NaN"],[8.425,2.425,0.04321798],[7.575,2.425,0.04321798],[7.575,2.425,0.7713429],[8.425,2.425,0.7713429],[8.425,2.425,0.04321798],["NaN","NaN","NaN"],[9.425,2.425,0.04321798],[8.575,2.425,0.04321798],[8.575,2.425,0.7576808],[9.425,2.425,0.7576808],[9.425,2.425,0.04321798],["NaN","NaN","NaN"],[1.425,3.425,0.04321798],[0.575,3.425,0.04321798],[0.575,3.425,0.07725442],[1.425,3.425,0.07725442],[1.425,3.425,0.04321798],["NaN","NaN","NaN"],[2.425,3.425,0.04321798],[1.575,3.425,0.04321798],[1.575,3.425,0.1672267],[2.425,3.425,0.1672267],[2.425,3.425,0.04321798],["NaN","NaN","NaN"],[3.425,3.425,0.04321798],[2.575,3.425,0.04321798],[2.575,3.425,0.3312058],[3.425,3.425,0.3312058],[3.425,3.425,0.04321798],["NaN","NaN","NaN"],[4.425,3.425,0.04321798],[3.575,3.425,0.04321798],[3.575,3.425,0.3553697],[4.425,3.425,0.3553697],[4.425,3.425,0.04321798],["NaN","NaN","NaN"],[5.425,3.425,0.04321798],[4.575,3.425,0.04321798],[4.575,3.425,0.781293],[5.425,3.425,0.781293],[5.425,3.425,0.04321798],["NaN","NaN","NaN"],[6.425,3.425,0.04321798],[5.575,3.425,0.04321798],[5.575,3.425,0.76804],[6.425,3.425,0.76804],[6.425,3.425,0.04321798],["NaN","NaN","NaN"],[7.425,3.425,0.04321798],[6.575,3.425,0.04321798],[6.575,3.425,0.7941267],[7.425,3.425,0.7941267],[7.425,3.425,0.04321798],["NaN","NaN","NaN"],[8.425,3.425,0.04321798],[7.575,3.425,0.04321798],[7.575,3.425,0.9553382],[8.425,3.425,0.9553382],[8.425,3.425,0.04321798],["NaN","NaN","NaN"],[9.425,3.425,0.04321798],[8.575,3.425,0.04321798],[8.575,3.425,0.9493445],[9.425,3.425,0.9493445],[9.425,3.425,0.04321798],["NaN","NaN","NaN"],[1.425,4.425,0.04321798],[0.575,4.425,0.04321798],[0.575,4.425,0.1025677],[1.425,4.425,0.1025677],[1.425,4.425,0.04321798],["NaN","NaN","NaN"],[2.425,4.425,0.04321798],[1.575,4.425,0.04321798],[1.575,4.425,0.2501325],[2.425,4.425,0.2501325],[2.425,4.425,0.04321798],["NaN","NaN","NaN"],[3.425,4.425,0.04321798],[2.575,4.425,0.04321798],[2.575,4.425,0.5039053],[3.425,4.425,0.5039053],[3.425,4.425,0.04321798],["NaN","NaN","NaN"],[4.425,4.425,0.04321798],[3.575,4.425,0.04321798],[3.575,4.425,0.537781],[4.425,4.425,0.537781],[4.425,4.425,0.04321798],["NaN","NaN","NaN"],[5.425,4.425,0.04321798],[4.575,4.425,0.04321798],[4.575,4.425,0.9430385],[5.425,4.425,0.9430385],[5.425,4.425,0.04321798],["NaN","NaN","NaN"],[6.425,4.425,0.04321798],[5.575,4.425,0.04321798],[5.575,4.425,0.9361931],[6.425,4.425,0.9361931],[6.425,4.425,0.04321798],["NaN","NaN","NaN"],[7.425,4.425,0.04321798],[6.575,4.425,0.04321798],[6.575,4.425,0.9493053],[7.425,4.425,0.9493053],[7.425,4.425,0.04321798],["NaN","NaN","NaN"],[8.425,4.425,0.04321798],[7.575,4.425,0.04321798],[7.575,4.425,0.9972544],[8.425,4.425,0.9972544],[8.425,4.425,0.04321798],["NaN","NaN","NaN"],[9.425,4.425,0.04321798],[8.575,4.425,0.04321798],[8.575,4.425,0.9965154],[9.425,4.425,0.9965154],[9.425,4.425,0.04321798],["NaN","NaN","NaN"],[1.425,5.425,0.04321798],[0.575,5.425,0.04321798],[0.575,5.425,0.1352801],[1.425,5.425,0.1352801],[1.425,5.425,0.04321798],["NaN","NaN","NaN"],[2.425,5.425,0.04321798],[1.575,5.425,0.04321798],[1.575,5.425,0.3567418],[2.425,5.425,0.3567418],[2.425,5.425,0.04321798],["NaN","NaN","NaN"],[3.425,5.425,0.04321798],[2.575,5.425,0.04321798],[2.575,5.425,0.6844499],[3.425,5.425,0.6844499],[3.425,5.425,0.04321798],["NaN","NaN","NaN"],[4.425,5.425,0.04321798],[3.575,5.425,0.04321798],[3.575,5.425,0.7209224],[4.425,5.425,0.7209224],[4.425,5.425,0.04321798],["NaN","NaN","NaN"],[5.425,5.425,0.04321798],[4.575,5.425,0.04321798],[4.575,5.425,0.9923994],[5.425,5.425,0.9923994],[5.425,5.425,0.04321798],["NaN","NaN","NaN"],[6.425,5.425,0.04321798],[5.575,5.425,0.04321798],[5.575,5.425,0.9907784],[6.425,5.425,0.9907784],[6.425,5.425,0.04321798],["NaN","NaN","NaN"],[7.425,5.425,0.04321798],[6.575,5.425,0.04321798],[6.575,5.425,0.9937668],[7.425,5.425,0.9937668],[7.425,5.425,0.04321798],["NaN","NaN","NaN"],[8.425,5.425,0.04321798],[7.575,5.425,0.04321798],[7.575,5.425,0.999954],[8.425,5.425,0.999954],[8.425,5.425,0.04321798],["NaN","NaN","NaN"],[9.425,5.425,0.04321798],[8.575,5.425,0.04321798],[8.575,5.425,0.9999315],[9.425,5.425,0.9999315],[9.425,5.425,0.04321798],["NaN","NaN","NaN"],[1.425,6.425,0.04321798],[0.575,6.425,0.04321798],[0.575,6.425,0.1997644],[1.425,6.425,0.1997644],[1.425,6.425,0.04321798],["NaN","NaN","NaN"],[2.425,6.425,0.04321798],[1.575,6.425,0.04321798],[1.575,6.425,0.5454682],[2.425,6.425,0.5454682],[2.425,6.425,0.04321798],["NaN","NaN","NaN"],[3.425,6.425,0.04321798],[2.575,6.425,0.04321798],[2.575,6.425,0.8883711],[3.425,6.425,0.8883711],[3.425,6.425,0.04321798],["NaN","NaN","NaN"],[4.425,6.425,0.04321798],[3.575,6.425,0.04321798],[3.575,6.425,0.912325],[4.425,6.425,0.912325],[4.425,6.425,0.04321798],["NaN","NaN","NaN"],[5.425,6.425,0.04321798],[4.575,6.425,0.04321798],[4.575,6.425,0.9999133],[5.425,6.425,0.9999133],[5.425,6.425,0.04321798],["NaN","NaN","NaN"],[6.425,6.425,0.04321798],[5.575,6.425,0.04321798],[5.575,6.425,0.9998751],[6.425,6.425,0.9998751],[6.425,6.425,0.04321798],["NaN","NaN","NaN"],[7.425,6.425,0.04321798],[6.575,6.425,0.04321798],[6.575,6.425,0.9999403],[7.425,6.425,0.9999403],[7.425,6.425,0.04321798],["NaN","NaN","NaN"],[8.425,6.425,0.04321798],[7.575,6.425,0.04321798],[7.575,6.425,1],[8.425,6.425,1],[8.425,6.425,0.04321798],["NaN","NaN","NaN"],[9.425,6.425,0.04321798],[8.575,6.425,0.04321798],[8.575,6.425,1],[9.425,6.425,1],[9.425,6.425,0.04321798],["NaN","NaN","NaN"],[1.425,7.425,0.04321798],[0.575,7.425,0.04321798],[0.575,7.425,0.2790861],[1.425,7.425,0.2790861],[1.425,7.425,0.04321798],["NaN","NaN","NaN"],[2.425,7.425,0.04321798],[1.575,7.425,0.04321798],[1.575,7.425,0.722779],[2.425,7.425,0.722779],[2.425,7.425,0.04321798],["NaN","NaN","NaN"],[3.425,7.425,0.04321798],[2.575,7.425,0.04321798],[2.575,7.425,0.9743328],[3.425,7.425,0.9743328],[3.425,7.425,0.04321798],["NaN","NaN","NaN"],[4.425,7.425,0.04321798],[3.575,7.425,0.04321798],[3.575,7.425,0.9828939],[4.425,7.425,0.9828939],[4.425,7.425,0.04321798],["NaN","NaN","NaN"],[5.425,7.425,0.04321798],[4.575,7.425,0.04321798],[4.575,7.425,0.9999998],[5.425,7.425,0.9999998],[5.425,7.425,0.04321798],["NaN","NaN","NaN"],[6.425,7.425,0.04321798],[5.575,7.425,0.04321798],[5.575,7.425,0.9999996],[6.425,7.425,0.9999996],[6.425,7.425,0.04321798],["NaN","NaN","NaN"],[7.425,7.425,0.04321798],[6.575,7.425,0.04321798],[6.575,7.425,0.9999999],[7.425,7.425,0.9999999],[7.425,7.425,0.04321798],["NaN","NaN","NaN"],[8.425,7.425,0.04321798],[7.575,7.425,0.04321798],[7.575,7.425,1],[8.425,7.425,1],[8.425,7.425,0.04321798],["NaN","NaN","NaN"],[9.425,7.425,0.04321798],[8.575,7.425,0.04321798],[8.575,7.425,1],[9.425,7.425,1],[9.425,7.425,0.04321798],["NaN","NaN","NaN"],[1.425,8.425,0.04321798],[0.575,8.425,0.04321798],[0.575,8.425,0.428676],[1.425,8.425,0.428676],[1.425,8.425,0.04321798],["NaN","NaN","NaN"],[2.425,8.425,0.04321798],[1.575,8.425,0.04321798],[1.575,8.425,0.9105917],[2.425,8.425,0.9105917],[2.425,8.425,0.04321798],["NaN","NaN","NaN"],[3.425,8.425,0.04321798],[2.575,8.425,0.04321798],[2.575,8.425,0.9990419],[3.425,8.425,0.9990419],[3.425,8.425,0.04321798],["NaN","NaN","NaN"],[4.425,8.425,0.04321798],[3.575,8.425,0.04321798],[3.575,8.425,0.9995526],[4.425,8.425,0.9995526],[4.425,8.425,0.04321798],["NaN","NaN","NaN"],[5.425,8.425,0.04321798],[4.575,8.425,0.04321798],[4.575,8.425,1],[5.425,8.425,1],[5.425,8.425,0.04321798],["NaN","NaN","NaN"],[6.425,8.425,0.04321798],[5.575,8.425,0.04321798],[5.575,8.425,1],[6.425,8.425,1],[6.425,8.425,0.04321798],["NaN","NaN","NaN"],[7.425,8.425,0.04321798],[6.575,8.425,0.04321798],[6.575,8.425,1],[7.425,8.425,1],[7.425,8.425,0.04321798],["NaN","NaN","NaN"],[8.425,8.425,0.04321798],[7.575,8.425,0.04321798],[7.575,8.425,1],[8.425,8.425,1],[8.425,8.425,0.04321798],["NaN","NaN","NaN"],[9.425,8.425,0.04321798],[8.575,8.425,0.04321798],[8.575,8.425,1],[9.425,8.425,1],[9.425,8.425,0.04321798],["NaN","NaN","NaN"],[1.425,9.425,0.04321798],[0.575,9.425,0.04321798],[0.575,9.425,0.5885873],[1.425,9.425,0.5885873],[1.425,9.425,0.04321798],["NaN","NaN","NaN"],[2.425,9.425,0.04321798],[1.575,9.425,0.04321798],[1.575,9.425,0.981856],[2.425,9.425,0.981856],[2.425,9.425,0.04321798],["NaN","NaN","NaN"],[3.425,9.425,0.04321798],[2.575,9.425,0.04321798],[2.575,9.425,0.9999893],[3.425,9.425,0.9999893],[3.425,9.425,0.04321798],["NaN","NaN","NaN"],[4.425,9.425,0.04321798],[3.575,9.425,0.04321798],[3.575,9.425,0.9999969],[4.425,9.425,0.9999969],[4.425,9.425,0.04321798],["NaN","NaN","NaN"],[5.425,9.425,0.04321798],[4.575,9.425,0.04321798],[4.575,9.425,1],[5.425,9.425,1],[5.425,9.425,0.04321798],["NaN","NaN","NaN"],[6.425,9.425,0.04321798],[5.575,9.425,0.04321798],[5.575,9.425,1],[6.425,9.425,1],[6.425,9.425,0.04321798],["NaN","NaN","NaN"],[7.425,9.425,0.04321798],[6.575,9.425,0.04321798],[6.575,9.425,1],[7.425,9.425,1],[7.425,9.425,0.04321798],["NaN","NaN","NaN"],[8.425,9.425,0.04321798],[7.575,9.425,0.04321798],[7.575,9.425,1],[8.425,9.425,1],[8.425,9.425,0.04321798],["NaN","NaN","NaN"],[9.425,9.425,0.04321798],[8.575,9.425,0.04321798],[8.575,9.425,1],[9.425,9.425,1],[9.425,9.425,0.04321798],["NaN","NaN","NaN"],[1.425,10.425,0.04321798],[0.575,10.425,0.04321798],[0.575,10.425,0.7135577],[1.425,10.425,0.7135577],[1.425,10.425,0.04321798],["NaN","NaN","NaN"],[2.425,10.425,0.04321798],[1.575,10.425,0.04321798],[1.575,10.425,0.9967811],[2.425,10.425,0.9967811],[2.425,10.425,0.04321798],["NaN","NaN","NaN"],[3.425,10.425,0.04321798],[2.575,10.425,0.04321798],[2.575,10.425,0.9999999],[3.425,10.425,0.9999999],[3.425,10.425,0.04321798],["NaN","NaN","NaN"],[4.425,10.425,0.04321798],[3.575,10.425,0.04321798],[3.575,10.425,1],[4.425,10.425,1],[4.425,10.425,0.04321798],["NaN","NaN","NaN"],[5.425,10.425,0.04321798],[4.575,10.425,0.04321798],[4.575,10.425,1],[5.425,10.425,1],[5.425,10.425,0.04321798],["NaN","NaN","NaN"],[6.425,10.425,0.04321798],[5.575,10.425,0.04321798],[5.575,10.425,1],[6.425,10.425,1],[6.425,10.425,0.04321798],["NaN","NaN","NaN"],[7.425,10.425,0.04321798],[6.575,10.425,0.04321798],[6.575,10.425,1],[7.425,10.425,1],[7.425,10.425,0.04321798],["NaN","NaN","NaN"],[8.425,10.425,0.04321798],[7.575,10.425,0.04321798],[7.575,10.425,1],[8.425,10.425,1],[8.425,10.425,0.04321798],["NaN","NaN","NaN"],[9.425,10.425,0.04321798],[8.575,10.425,0.04321798],[8.575,10.425,1],[9.425,10.425,1],[9.425,10.425,0.04321798],["NaN","NaN","NaN"],[0.575,0.575,0.04321798],[0.575,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,0.575,0.04321798],[0.575,0.575,0.04321798],["NaN","NaN","NaN"],[1.575,0.575,0.04321798],[1.575,1.425,0.04321798],[1.575,1.425,0.06596308],[1.575,0.575,0.06596308],[1.575,0.575,0.04321798],["NaN","NaN","NaN"],[2.575,0.575,0.04321798],[2.575,1.425,0.04321798],[2.575,1.425,0.1021959],[2.575,0.575,0.1021959],[2.575,0.575,0.04321798],["NaN","NaN","NaN"],[3.575,0.575,0.04321798],[3.575,1.425,0.04321798],[3.575,1.425,0.1075265],[3.575,0.575,0.1075265],[3.575,0.575,0.04321798],["NaN","NaN","NaN"],[4.575,0.575,0.04321798],[4.575,1.425,0.04321798],[4.575,1.425,0.2316452],[4.575,0.575,0.2316452],[4.575,0.575,0.04321798],["NaN","NaN","NaN"],[5.575,0.575,0.04321798],[5.575,1.425,0.04321798],[5.575,1.425,0.2258095],[5.575,0.575,0.2258095],[5.575,0.575,0.04321798],["NaN","NaN","NaN"],[6.575,0.575,0.04321798],[6.575,1.425,0.04321798],[6.575,1.425,0.237572],[6.575,0.575,0.237572],[6.575,0.575,0.04321798],["NaN","NaN","NaN"],[7.575,0.575,0.04321798],[7.575,1.425,0.04321798],[7.575,1.425,0.3690261],[7.575,0.575,0.3690261],[7.575,0.575,0.04321798],["NaN","NaN","NaN"],[8.575,0.575,0.04321798],[8.575,1.425,0.04321798],[8.575,1.425,0.3592858],[8.575,0.575,0.3592858],[8.575,0.575,0.04321798],["NaN","NaN","NaN"],[0.575,1.575,0.04321798],[0.575,2.425,0.04321798],[0.575,2.425,0.05915068],[0.575,1.575,0.05915068],[0.575,1.575,0.04321798],["NaN","NaN","NaN"],[1.575,1.575,0.04321798],[1.575,2.425,0.04321798],[1.575,2.425,0.1108488],[1.575,1.575,0.1108488],[1.575,1.575,0.04321798],["NaN","NaN","NaN"],[2.575,1.575,0.04321798],[2.575,2.425,0.04321798],[2.575,2.425,0.2030372],[2.575,1.575,0.2030372],[2.575,1.575,0.04321798],["NaN","NaN","NaN"],[3.575,1.575,0.04321798],[3.575,2.425,0.04321798],[3.575,2.425,0.217041],[3.575,1.575,0.217041],[3.575,1.575,0.04321798],["NaN","NaN","NaN"],[4.575,1.575,0.04321798],[4.575,2.425,0.04321798],[4.575,2.425,0.5257641],[4.575,1.575,0.5257641],[4.575,1.575,0.04321798],["NaN","NaN","NaN"],[5.575,1.575,0.04321798],[5.575,2.425,0.04321798],[5.575,2.425,0.5127993],[5.575,1.575,0.5127993],[5.575,1.575,0.04321798],["NaN","NaN","NaN"],[6.575,1.575,0.04321798],[6.575,2.425,0.04321798],[6.575,2.425,0.5387281],[6.575,1.575,0.5387281],[6.575,1.575,0.04321798],["NaN","NaN","NaN"],[7.575,1.575,0.04321798],[7.575,2.425,0.04321798],[7.575,2.425,0.7713429],[7.575,1.575,0.7713429],[7.575,1.575,0.04321798],["NaN","NaN","NaN"],[8.575,1.575,0.04321798],[8.575,2.425,0.04321798],[8.575,2.425,0.7576808],[8.575,1.575,0.7576808],[8.575,1.575,0.04321798],["NaN","NaN","NaN"],[0.575,2.575,0.04321798],[0.575,3.425,0.04321798],[0.575,3.425,0.07725442],[0.575,2.575,0.07725442],[0.575,2.575,0.04321798],["NaN","NaN","NaN"],[1.575,2.575,0.04321798],[1.575,3.425,0.04321798],[1.575,3.425,0.1672267],[1.575,2.575,0.1672267],[1.575,2.575,0.04321798],["NaN","NaN","NaN"],[2.575,2.575,0.04321798],[2.575,3.425,0.04321798],[2.575,3.425,0.3312058],[2.575,2.575,0.3312058],[2.575,2.575,0.04321798],["NaN","NaN","NaN"],[3.575,2.575,0.04321798],[3.575,3.425,0.04321798],[3.575,3.425,0.3553697],[3.575,2.575,0.3553697],[3.575,2.575,0.04321798],["NaN","NaN","NaN"],[4.575,2.575,0.04321798],[4.575,3.425,0.04321798],[4.575,3.425,0.781293],[4.575,2.575,0.781293],[4.575,2.575,0.04321798],["NaN","NaN","NaN"],[5.575,2.575,0.04321798],[5.575,3.425,0.04321798],[5.575,3.425,0.76804],[5.575,2.575,0.76804],[5.575,2.575,0.04321798],["NaN","NaN","NaN"],[6.575,2.575,0.04321798],[6.575,3.425,0.04321798],[6.575,3.425,0.7941267],[6.575,2.575,0.7941267],[6.575,2.575,0.04321798],["NaN","NaN","NaN"],[7.575,2.575,0.04321798],[7.575,3.425,0.04321798],[7.575,3.425,0.9553382],[7.575,2.575,0.9553382],[7.575,2.575,0.04321798],["NaN","NaN","NaN"],[8.575,2.575,0.04321798],[8.575,3.425,0.04321798],[8.575,3.425,0.9493445],[8.575,2.575,0.9493445],[8.575,2.575,0.04321798],["NaN","NaN","NaN"],[0.575,3.575,0.04321798],[0.575,4.425,0.04321798],[0.575,4.425,0.1025677],[0.575,3.575,0.1025677],[0.575,3.575,0.04321798],["NaN","NaN","NaN"],[1.575,3.575,0.04321798],[1.575,4.425,0.04321798],[1.575,4.425,0.2501325],[1.575,3.575,0.2501325],[1.575,3.575,0.04321798],["NaN","NaN","NaN"],[2.575,3.575,0.04321798],[2.575,4.425,0.04321798],[2.575,4.425,0.5039053],[2.575,3.575,0.5039053],[2.575,3.575,0.04321798],["NaN","NaN","NaN"],[3.575,3.575,0.04321798],[3.575,4.425,0.04321798],[3.575,4.425,0.537781],[3.575,3.575,0.537781],[3.575,3.575,0.04321798],["NaN","NaN","NaN"],[4.575,3.575,0.04321798],[4.575,4.425,0.04321798],[4.575,4.425,0.9430385],[4.575,3.575,0.9430385],[4.575,3.575,0.04321798],["NaN","NaN","NaN"],[5.575,3.575,0.04321798],[5.575,4.425,0.04321798],[5.575,4.425,0.9361931],[5.575,3.575,0.9361931],[5.575,3.575,0.04321798],["NaN","NaN","NaN"],[6.575,3.575,0.04321798],[6.575,4.425,0.04321798],[6.575,4.425,0.9493053],[6.575,3.575,0.9493053],[6.575,3.575,0.04321798],["NaN","NaN","NaN"],[7.575,3.575,0.04321798],[7.575,4.425,0.04321798],[7.575,4.425,0.9972544],[7.575,3.575,0.9972544],[7.575,3.575,0.04321798],["NaN","NaN","NaN"],[8.575,3.575,0.04321798],[8.575,4.425,0.04321798],[8.575,4.425,0.9965154],[8.575,3.575,0.9965154],[8.575,3.575,0.04321798],["NaN","NaN","NaN"],[0.575,4.575,0.04321798],[0.575,5.425,0.04321798],[0.575,5.425,0.1352801],[0.575,4.575,0.1352801],[0.575,4.575,0.04321798],["NaN","NaN","NaN"],[1.575,4.575,0.04321798],[1.575,5.425,0.04321798],[1.575,5.425,0.3567418],[1.575,4.575,0.3567418],[1.575,4.575,0.04321798],["NaN","NaN","NaN"],[2.575,4.575,0.04321798],[2.575,5.425,0.04321798],[2.575,5.425,0.6844499],[2.575,4.575,0.6844499],[2.575,4.575,0.04321798],["NaN","NaN","NaN"],[3.575,4.575,0.04321798],[3.575,5.425,0.04321798],[3.575,5.425,0.7209224],[3.575,4.575,0.7209224],[3.575,4.575,0.04321798],["NaN","NaN","NaN"],[4.575,4.575,0.04321798],[4.575,5.425,0.04321798],[4.575,5.425,0.9923994],[4.575,4.575,0.9923994],[4.575,4.575,0.04321798],["NaN","NaN","NaN"],[5.575,4.575,0.04321798],[5.575,5.425,0.04321798],[5.575,5.425,0.9907784],[5.575,4.575,0.9907784],[5.575,4.575,0.04321798],["NaN","NaN","NaN"],[6.575,4.575,0.04321798],[6.575,5.425,0.04321798],[6.575,5.425,0.9937668],[6.575,4.575,0.9937668],[6.575,4.575,0.04321798],["NaN","NaN","NaN"],[7.575,4.575,0.04321798],[7.575,5.425,0.04321798],[7.575,5.425,0.999954],[7.575,4.575,0.999954],[7.575,4.575,0.04321798],["NaN","NaN","NaN"],[8.575,4.575,0.04321798],[8.575,5.425,0.04321798],[8.575,5.425,0.9999315],[8.575,4.575,0.9999315],[8.575,4.575,0.04321798],["NaN","NaN","NaN"],[0.575,5.575,0.04321798],[0.575,6.425,0.04321798],[0.575,6.425,0.1997644],[0.575,5.575,0.1997644],[0.575,5.575,0.04321798],["NaN","NaN","NaN"],[1.575,5.575,0.04321798],[1.575,6.425,0.04321798],[1.575,6.425,0.5454682],[1.575,5.575,0.5454682],[1.575,5.575,0.04321798],["NaN","NaN","NaN"],[2.575,5.575,0.04321798],[2.575,6.425,0.04321798],[2.575,6.425,0.8883711],[2.575,5.575,0.8883711],[2.575,5.575,0.04321798],["NaN","NaN","NaN"],[3.575,5.575,0.04321798],[3.575,6.425,0.04321798],[3.575,6.425,0.912325],[3.575,5.575,0.912325],[3.575,5.575,0.04321798],["NaN","NaN","NaN"],[4.575,5.575,0.04321798],[4.575,6.425,0.04321798],[4.575,6.425,0.9999133],[4.575,5.575,0.9999133],[4.575,5.575,0.04321798],["NaN","NaN","NaN"],[5.575,5.575,0.04321798],[5.575,6.425,0.04321798],[5.575,6.425,0.9998751],[5.575,5.575,0.9998751],[5.575,5.575,0.04321798],["NaN","NaN","NaN"],[6.575,5.575,0.04321798],[6.575,6.425,0.04321798],[6.575,6.425,0.9999403],[6.575,5.575,0.9999403],[6.575,5.575,0.04321798],["NaN","NaN","NaN"],[7.575,5.575,0.04321798],[7.575,6.425,0.04321798],[7.575,6.425,1],[7.575,5.575,1],[7.575,5.575,0.04321798],["NaN","NaN","NaN"],[8.575,5.575,0.04321798],[8.575,6.425,0.04321798],[8.575,6.425,1],[8.575,5.575,1],[8.575,5.575,0.04321798],["NaN","NaN","NaN"],[0.575,6.575,0.04321798],[0.575,7.425,0.04321798],[0.575,7.425,0.2790861],[0.575,6.575,0.2790861],[0.575,6.575,0.04321798],["NaN","NaN","NaN"],[1.575,6.575,0.04321798],[1.575,7.425,0.04321798],[1.575,7.425,0.722779],[1.575,6.575,0.722779],[1.575,6.575,0.04321798],["NaN","NaN","NaN"],[2.575,6.575,0.04321798],[2.575,7.425,0.04321798],[2.575,7.425,0.9743328],[2.575,6.575,0.9743328],[2.575,6.575,0.04321798],["NaN","NaN","NaN"],[3.575,6.575,0.04321798],[3.575,7.425,0.04321798],[3.575,7.425,0.9828939],[3.575,6.575,0.9828939],[3.575,6.575,0.04321798],["NaN","NaN","NaN"],[4.575,6.575,0.04321798],[4.575,7.425,0.04321798],[4.575,7.425,0.9999998],[4.575,6.575,0.9999998],[4.575,6.575,0.04321798],["NaN","NaN","NaN"],[5.575,6.575,0.04321798],[5.575,7.425,0.04321798],[5.575,7.425,0.9999996],[5.575,6.575,0.9999996],[5.575,6.575,0.04321798],["NaN","NaN","NaN"],[6.575,6.575,0.04321798],[6.575,7.425,0.04321798],[6.575,7.425,0.9999999],[6.575,6.575,0.9999999],[6.575,6.575,0.04321798],["NaN","NaN","NaN"],[7.575,6.575,0.04321798],[7.575,7.425,0.04321798],[7.575,7.425,1],[7.575,6.575,1],[7.575,6.575,0.04321798],["NaN","NaN","NaN"],[8.575,6.575,0.04321798],[8.575,7.425,0.04321798],[8.575,7.425,1],[8.575,6.575,1],[8.575,6.575,0.04321798],["NaN","NaN","NaN"],[0.575,7.575,0.04321798],[0.575,8.425,0.04321798],[0.575,8.425,0.428676],[0.575,7.575,0.428676],[0.575,7.575,0.04321798],["NaN","NaN","NaN"],[1.575,7.575,0.04321798],[1.575,8.425,0.04321798],[1.575,8.425,0.9105917],[1.575,7.575,0.9105917],[1.575,7.575,0.04321798],["NaN","NaN","NaN"],[2.575,7.575,0.04321798],[2.575,8.425,0.04321798],[2.575,8.425,0.9990419],[2.575,7.575,0.9990419],[2.575,7.575,0.04321798],["NaN","NaN","NaN"],[3.575,7.575,0.04321798],[3.575,8.425,0.04321798],[3.575,8.425,0.9995526],[3.575,7.575,0.9995526],[3.575,7.575,0.04321798],["NaN","NaN","NaN"],[4.575,7.575,0.04321798],[4.575,8.425,0.04321798],[4.575,8.425,1],[4.575,7.575,1],[4.575,7.575,0.04321798],["NaN","NaN","NaN"],[5.575,7.575,0.04321798],[5.575,8.425,0.04321798],[5.575,8.425,1],[5.575,7.575,1],[5.575,7.575,0.04321798],["NaN","NaN","NaN"],[6.575,7.575,0.04321798],[6.575,8.425,0.04321798],[6.575,8.425,1],[6.575,7.575,1],[6.575,7.575,0.04321798],["NaN","NaN","NaN"],[7.575,7.575,0.04321798],[7.575,8.425,0.04321798],[7.575,8.425,1],[7.575,7.575,1],[7.575,7.575,0.04321798],["NaN","NaN","NaN"],[8.575,7.575,0.04321798],[8.575,8.425,0.04321798],[8.575,8.425,1],[8.575,7.575,1],[8.575,7.575,0.04321798],["NaN","NaN","NaN"],[0.575,8.575,0.04321798],[0.575,9.425,0.04321798],[0.575,9.425,0.5885873],[0.575,8.575,0.5885873],[0.575,8.575,0.04321798],["NaN","NaN","NaN"],[1.575,8.575,0.04321798],[1.575,9.425,0.04321798],[1.575,9.425,0.981856],[1.575,8.575,0.981856],[1.575,8.575,0.04321798],["NaN","NaN","NaN"],[2.575,8.575,0.04321798],[2.575,9.425,0.04321798],[2.575,9.425,0.9999893],[2.575,8.575,0.9999893],[2.575,8.575,0.04321798],["NaN","NaN","NaN"],[3.575,8.575,0.04321798],[3.575,9.425,0.04321798],[3.575,9.425,0.9999969],[3.575,8.575,0.9999969],[3.575,8.575,0.04321798],["NaN","NaN","NaN"],[4.575,8.575,0.04321798],[4.575,9.425,0.04321798],[4.575,9.425,1],[4.575,8.575,1],[4.575,8.575,0.04321798],["NaN","NaN","NaN"],[5.575,8.575,0.04321798],[5.575,9.425,0.04321798],[5.575,9.425,1],[5.575,8.575,1],[5.575,8.575,0.04321798],["NaN","NaN","NaN"],[6.575,8.575,0.04321798],[6.575,9.425,0.04321798],[6.575,9.425,1],[6.575,8.575,1],[6.575,8.575,0.04321798],["NaN","NaN","NaN"],[7.575,8.575,0.04321798],[7.575,9.425,0.04321798],[7.575,9.425,1],[7.575,8.575,1],[7.575,8.575,0.04321798],["NaN","NaN","NaN"],[8.575,8.575,0.04321798],[8.575,9.425,0.04321798],[8.575,9.425,1],[8.575,8.575,1],[8.575,8.575,0.04321798],["NaN","NaN","NaN"],[0.575,9.575,0.04321798],[0.575,10.425,0.04321798],[0.575,10.425,0.7135577],[0.575,9.575,0.7135577],[0.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.575,9.575,0.04321798],[1.575,10.425,0.04321798],[1.575,10.425,0.9967811],[1.575,9.575,0.9967811],[1.575,9.575,0.04321798],["NaN","NaN","NaN"],[2.575,9.575,0.04321798],[2.575,10.425,0.04321798],[2.575,10.425,0.9999999],[2.575,9.575,0.9999999],[2.575,9.575,0.04321798],["NaN","NaN","NaN"],[3.575,9.575,0.04321798],[3.575,10.425,0.04321798],[3.575,10.425,1],[3.575,9.575,1],[3.575,9.575,0.04321798],["NaN","NaN","NaN"],[4.575,9.575,0.04321798],[4.575,10.425,0.04321798],[4.575,10.425,1],[4.575,9.575,1],[4.575,9.575,0.04321798],["NaN","NaN","NaN"],[5.575,9.575,0.04321798],[5.575,10.425,0.04321798],[5.575,10.425,1],[5.575,9.575,1],[5.575,9.575,0.04321798],["NaN","NaN","NaN"],[6.575,9.575,0.04321798],[6.575,10.425,0.04321798],[6.575,10.425,1],[6.575,9.575,1],[6.575,9.575,0.04321798],["NaN","NaN","NaN"],[7.575,9.575,0.04321798],[7.575,10.425,0.04321798],[7.575,10.425,1],[7.575,9.575,1],[7.575,9.575,0.04321798],["NaN","NaN","NaN"],[8.575,9.575,0.04321798],[8.575,10.425,0.04321798],[8.575,10.425,1],[8.575,9.575,1],[8.575,9.575,0.04321798],["NaN","NaN","NaN"],[1.425,0.575,0.04321798],[1.425,1.425,0.04321798],[1.425,1.425,0.04321798],[1.425,0.575,0.04321798],[1.425,0.575,0.04321798],["NaN","NaN","NaN"],[2.425,0.575,0.04321798],[2.425,1.425,0.04321798],[2.425,1.425,0.06596308],[2.425,0.575,0.06596308],[2.425,0.575,0.04321798],["NaN","NaN","NaN"],[3.425,0.575,0.04321798],[3.425,1.425,0.04321798],[3.425,1.425,0.1021959],[3.425,0.575,0.1021959],[3.425,0.575,0.04321798],["NaN","NaN","NaN"],[4.425,0.575,0.04321798],[4.425,1.425,0.04321798],[4.425,1.425,0.1075265],[4.425,0.575,0.1075265],[4.425,0.575,0.04321798],["NaN","NaN","NaN"],[5.425,0.575,0.04321798],[5.425,1.425,0.04321798],[5.425,1.425,0.2316452],[5.425,0.575,0.2316452],[5.425,0.575,0.04321798],["NaN","NaN","NaN"],[6.425,0.575,0.04321798],[6.425,1.425,0.04321798],[6.425,1.425,0.2258095],[6.425,0.575,0.2258095],[6.425,0.575,0.04321798],["NaN","NaN","NaN"],[7.425,0.575,0.04321798],[7.425,1.425,0.04321798],[7.425,1.425,0.237572],[7.425,0.575,0.237572],[7.425,0.575,0.04321798],["NaN","NaN","NaN"],[8.425,0.575,0.04321798],[8.425,1.425,0.04321798],[8.425,1.425,0.3690261],[8.425,0.575,0.3690261],[8.425,0.575,0.04321798],["NaN","NaN","NaN"],[9.425,0.575,0.04321798],[9.425,1.425,0.04321798],[9.425,1.425,0.3592858],[9.425,0.575,0.3592858],[9.425,0.575,0.04321798],["NaN","NaN","NaN"],[1.425,1.575,0.04321798],[1.425,2.425,0.04321798],[1.425,2.425,0.05915068],[1.425,1.575,0.05915068],[1.425,1.575,0.04321798],["NaN","NaN","NaN"],[2.425,1.575,0.04321798],[2.425,2.425,0.04321798],[2.425,2.425,0.1108488],[2.425,1.575,0.1108488],[2.425,1.575,0.04321798],["NaN","NaN","NaN"],[3.425,1.575,0.04321798],[3.425,2.425,0.04321798],[3.425,2.425,0.2030372],[3.425,1.575,0.2030372],[3.425,1.575,0.04321798],["NaN","NaN","NaN"],[4.425,1.575,0.04321798],[4.425,2.425,0.04321798],[4.425,2.425,0.217041],[4.425,1.575,0.217041],[4.425,1.575,0.04321798],["NaN","NaN","NaN"],[5.425,1.575,0.04321798],[5.425,2.425,0.04321798],[5.425,2.425,0.5257641],[5.425,1.575,0.5257641],[5.425,1.575,0.04321798],["NaN","NaN","NaN"],[6.425,1.575,0.04321798],[6.425,2.425,0.04321798],[6.425,2.425,0.5127993],[6.425,1.575,0.5127993],[6.425,1.575,0.04321798],["NaN","NaN","NaN"],[7.425,1.575,0.04321798],[7.425,2.425,0.04321798],[7.425,2.425,0.5387281],[7.425,1.575,0.5387281],[7.425,1.575,0.04321798],["NaN","NaN","NaN"],[8.425,1.575,0.04321798],[8.425,2.425,0.04321798],[8.425,2.425,0.7713429],[8.425,1.575,0.7713429],[8.425,1.575,0.04321798],["NaN","NaN","NaN"],[9.425,1.575,0.04321798],[9.425,2.425,0.04321798],[9.425,2.425,0.7576808],[9.425,1.575,0.7576808],[9.425,1.575,0.04321798],["NaN","NaN","NaN"],[1.425,2.575,0.04321798],[1.425,3.425,0.04321798],[1.425,3.425,0.07725442],[1.425,2.575,0.07725442],[1.425,2.575,0.04321798],["NaN","NaN","NaN"],[2.425,2.575,0.04321798],[2.425,3.425,0.04321798],[2.425,3.425,0.1672267],[2.425,2.575,0.1672267],[2.425,2.575,0.04321798],["NaN","NaN","NaN"],[3.425,2.575,0.04321798],[3.425,3.425,0.04321798],[3.425,3.425,0.3312058],[3.425,2.575,0.3312058],[3.425,2.575,0.04321798],["NaN","NaN","NaN"],[4.425,2.575,0.04321798],[4.425,3.425,0.04321798],[4.425,3.425,0.3553697],[4.425,2.575,0.3553697],[4.425,2.575,0.04321798],["NaN","NaN","NaN"],[5.425,2.575,0.04321798],[5.425,3.425,0.04321798],[5.425,3.425,0.781293],[5.425,2.575,0.781293],[5.425,2.575,0.04321798],["NaN","NaN","NaN"],[6.425,2.575,0.04321798],[6.425,3.425,0.04321798],[6.425,3.425,0.76804],[6.425,2.575,0.76804],[6.425,2.575,0.04321798],["NaN","NaN","NaN"],[7.425,2.575,0.04321798],[7.425,3.425,0.04321798],[7.425,3.425,0.7941267],[7.425,2.575,0.7941267],[7.425,2.575,0.04321798],["NaN","NaN","NaN"],[8.425,2.575,0.04321798],[8.425,3.425,0.04321798],[8.425,3.425,0.9553382],[8.425,2.575,0.9553382],[8.425,2.575,0.04321798],["NaN","NaN","NaN"],[9.425,2.575,0.04321798],[9.425,3.425,0.04321798],[9.425,3.425,0.9493445],[9.425,2.575,0.9493445],[9.425,2.575,0.04321798],["NaN","NaN","NaN"],[1.425,3.575,0.04321798],[1.425,4.425,0.04321798],[1.425,4.425,0.1025677],[1.425,3.575,0.1025677],[1.425,3.575,0.04321798],["NaN","NaN","NaN"],[2.425,3.575,0.04321798],[2.425,4.425,0.04321798],[2.425,4.425,0.2501325],[2.425,3.575,0.2501325],[2.425,3.575,0.04321798],["NaN","NaN","NaN"],[3.425,3.575,0.04321798],[3.425,4.425,0.04321798],[3.425,4.425,0.5039053],[3.425,3.575,0.5039053],[3.425,3.575,0.04321798],["NaN","NaN","NaN"],[4.425,3.575,0.04321798],[4.425,4.425,0.04321798],[4.425,4.425,0.537781],[4.425,3.575,0.537781],[4.425,3.575,0.04321798],["NaN","NaN","NaN"],[5.425,3.575,0.04321798],[5.425,4.425,0.04321798],[5.425,4.425,0.9430385],[5.425,3.575,0.9430385],[5.425,3.575,0.04321798],["NaN","NaN","NaN"],[6.425,3.575,0.04321798],[6.425,4.425,0.04321798],[6.425,4.425,0.9361931],[6.425,3.575,0.9361931],[6.425,3.575,0.04321798],["NaN","NaN","NaN"],[7.425,3.575,0.04321798],[7.425,4.425,0.04321798],[7.425,4.425,0.9493053],[7.425,3.575,0.9493053],[7.425,3.575,0.04321798],["NaN","NaN","NaN"],[8.425,3.575,0.04321798],[8.425,4.425,0.04321798],[8.425,4.425,0.9972544],[8.425,3.575,0.9972544],[8.425,3.575,0.04321798],["NaN","NaN","NaN"],[9.425,3.575,0.04321798],[9.425,4.425,0.04321798],[9.425,4.425,0.9965154],[9.425,3.575,0.9965154],[9.425,3.575,0.04321798],["NaN","NaN","NaN"],[1.425,4.575,0.04321798],[1.425,5.425,0.04321798],[1.425,5.425,0.1352801],[1.425,4.575,0.1352801],[1.425,4.575,0.04321798],["NaN","NaN","NaN"],[2.425,4.575,0.04321798],[2.425,5.425,0.04321798],[2.425,5.425,0.3567418],[2.425,4.575,0.3567418],[2.425,4.575,0.04321798],["NaN","NaN","NaN"],[3.425,4.575,0.04321798],[3.425,5.425,0.04321798],[3.425,5.425,0.6844499],[3.425,4.575,0.6844499],[3.425,4.575,0.04321798],["NaN","NaN","NaN"],[4.425,4.575,0.04321798],[4.425,5.425,0.04321798],[4.425,5.425,0.7209224],[4.425,4.575,0.7209224],[4.425,4.575,0.04321798],["NaN","NaN","NaN"],[5.425,4.575,0.04321798],[5.425,5.425,0.04321798],[5.425,5.425,0.9923994],[5.425,4.575,0.9923994],[5.425,4.575,0.04321798],["NaN","NaN","NaN"],[6.425,4.575,0.04321798],[6.425,5.425,0.04321798],[6.425,5.425,0.9907784],[6.425,4.575,0.9907784],[6.425,4.575,0.04321798],["NaN","NaN","NaN"],[7.425,4.575,0.04321798],[7.425,5.425,0.04321798],[7.425,5.425,0.9937668],[7.425,4.575,0.9937668],[7.425,4.575,0.04321798],["NaN","NaN","NaN"],[8.425,4.575,0.04321798],[8.425,5.425,0.04321798],[8.425,5.425,0.999954],[8.425,4.575,0.999954],[8.425,4.575,0.04321798],["NaN","NaN","NaN"],[9.425,4.575,0.04321798],[9.425,5.425,0.04321798],[9.425,5.425,0.9999315],[9.425,4.575,0.9999315],[9.425,4.575,0.04321798],["NaN","NaN","NaN"],[1.425,5.575,0.04321798],[1.425,6.425,0.04321798],[1.425,6.425,0.1997644],[1.425,5.575,0.1997644],[1.425,5.575,0.04321798],["NaN","NaN","NaN"],[2.425,5.575,0.04321798],[2.425,6.425,0.04321798],[2.425,6.425,0.5454682],[2.425,5.575,0.5454682],[2.425,5.575,0.04321798],["NaN","NaN","NaN"],[3.425,5.575,0.04321798],[3.425,6.425,0.04321798],[3.425,6.425,0.8883711],[3.425,5.575,0.8883711],[3.425,5.575,0.04321798],["NaN","NaN","NaN"],[4.425,5.575,0.04321798],[4.425,6.425,0.04321798],[4.425,6.425,0.912325],[4.425,5.575,0.912325],[4.425,5.575,0.04321798],["NaN","NaN","NaN"],[5.425,5.575,0.04321798],[5.425,6.425,0.04321798],[5.425,6.425,0.9999133],[5.425,5.575,0.9999133],[5.425,5.575,0.04321798],["NaN","NaN","NaN"],[6.425,5.575,0.04321798],[6.425,6.425,0.04321798],[6.425,6.425,0.9998751],[6.425,5.575,0.9998751],[6.425,5.575,0.04321798],["NaN","NaN","NaN"],[7.425,5.575,0.04321798],[7.425,6.425,0.04321798],[7.425,6.425,0.9999403],[7.425,5.575,0.9999403],[7.425,5.575,0.04321798],["NaN","NaN","NaN"],[8.425,5.575,0.04321798],[8.425,6.425,0.04321798],[8.425,6.425,1],[8.425,5.575,1],[8.425,5.575,0.04321798],["NaN","NaN","NaN"],[9.425,5.575,0.04321798],[9.425,6.425,0.04321798],[9.425,6.425,1],[9.425,5.575,1],[9.425,5.575,0.04321798],["NaN","NaN","NaN"],[1.425,6.575,0.04321798],[1.425,7.425,0.04321798],[1.425,7.425,0.2790861],[1.425,6.575,0.2790861],[1.425,6.575,0.04321798],["NaN","NaN","NaN"],[2.425,6.575,0.04321798],[2.425,7.425,0.04321798],[2.425,7.425,0.722779],[2.425,6.575,0.722779],[2.425,6.575,0.04321798],["NaN","NaN","NaN"],[3.425,6.575,0.04321798],[3.425,7.425,0.04321798],[3.425,7.425,0.9743328],[3.425,6.575,0.9743328],[3.425,6.575,0.04321798],["NaN","NaN","NaN"],[4.425,6.575,0.04321798],[4.425,7.425,0.04321798],[4.425,7.425,0.9828939],[4.425,6.575,0.9828939],[4.425,6.575,0.04321798],["NaN","NaN","NaN"],[5.425,6.575,0.04321798],[5.425,7.425,0.04321798],[5.425,7.425,0.9999998],[5.425,6.575,0.9999998],[5.425,6.575,0.04321798],["NaN","NaN","NaN"],[6.425,6.575,0.04321798],[6.425,7.425,0.04321798],[6.425,7.425,0.9999996],[6.425,6.575,0.9999996],[6.425,6.575,0.04321798],["NaN","NaN","NaN"],[7.425,6.575,0.04321798],[7.425,7.425,0.04321798],[7.425,7.425,0.9999999],[7.425,6.575,0.9999999],[7.425,6.575,0.04321798],["NaN","NaN","NaN"],[8.425,6.575,0.04321798],[8.425,7.425,0.04321798],[8.425,7.425,1],[8.425,6.575,1],[8.425,6.575,0.04321798],["NaN","NaN","NaN"],[9.425,6.575,0.04321798],[9.425,7.425,0.04321798],[9.425,7.425,1],[9.425,6.575,1],[9.425,6.575,0.04321798],["NaN","NaN","NaN"],[1.425,7.575,0.04321798],[1.425,8.425,0.04321798],[1.425,8.425,0.428676],[1.425,7.575,0.428676],[1.425,7.575,0.04321798],["NaN","NaN","NaN"],[2.425,7.575,0.04321798],[2.425,8.425,0.04321798],[2.425,8.425,0.9105917],[2.425,7.575,0.9105917],[2.425,7.575,0.04321798],["NaN","NaN","NaN"],[3.425,7.575,0.04321798],[3.425,8.425,0.04321798],[3.425,8.425,0.9990419],[3.425,7.575,0.9990419],[3.425,7.575,0.04321798],["NaN","NaN","NaN"],[4.425,7.575,0.04321798],[4.425,8.425,0.04321798],[4.425,8.425,0.9995526],[4.425,7.575,0.9995526],[4.425,7.575,0.04321798],["NaN","NaN","NaN"],[5.425,7.575,0.04321798],[5.425,8.425,0.04321798],[5.425,8.425,1],[5.425,7.575,1],[5.425,7.575,0.04321798],["NaN","NaN","NaN"],[6.425,7.575,0.04321798],[6.425,8.425,0.04321798],[6.425,8.425,1],[6.425,7.575,1],[6.425,7.575,0.04321798],["NaN","NaN","NaN"],[7.425,7.575,0.04321798],[7.425,8.425,0.04321798],[7.425,8.425,1],[7.425,7.575,1],[7.425,7.575,0.04321798],["NaN","NaN","NaN"],[8.425,7.575,0.04321798],[8.425,8.425,0.04321798],[8.425,8.425,1],[8.425,7.575,1],[8.425,7.575,0.04321798],["NaN","NaN","NaN"],[9.425,7.575,0.04321798],[9.425,8.425,0.04321798],[9.425,8.425,1],[9.425,7.575,1],[9.425,7.575,0.04321798],["NaN","NaN","NaN"],[1.425,8.575,0.04321798],[1.425,9.425,0.04321798],[1.425,9.425,0.5885873],[1.425,8.575,0.5885873],[1.425,8.575,0.04321798],["NaN","NaN","NaN"],[2.425,8.575,0.04321798],[2.425,9.425,0.04321798],[2.425,9.425,0.981856],[2.425,8.575,0.981856],[2.425,8.575,0.04321798],["NaN","NaN","NaN"],[3.425,8.575,0.04321798],[3.425,9.425,0.04321798],[3.425,9.425,0.9999893],[3.425,8.575,0.9999893],[3.425,8.575,0.04321798],["NaN","NaN","NaN"],[4.425,8.575,0.04321798],[4.425,9.425,0.04321798],[4.425,9.425,0.9999969],[4.425,8.575,0.9999969],[4.425,8.575,0.04321798],["NaN","NaN","NaN"],[5.425,8.575,0.04321798],[5.425,9.425,0.04321798],[5.425,9.425,1],[5.425,8.575,1],[5.425,8.575,0.04321798],["NaN","NaN","NaN"],[6.425,8.575,0.04321798],[6.425,9.425,0.04321798],[6.425,9.425,1],[6.425,8.575,1],[6.425,8.575,0.04321798],["NaN","NaN","NaN"],[7.425,8.575,0.04321798],[7.425,9.425,0.04321798],[7.425,9.425,1],[7.425,8.575,1],[7.425,8.575,0.04321798],["NaN","NaN","NaN"],[8.425,8.575,0.04321798],[8.425,9.425,0.04321798],[8.425,9.425,1],[8.425,8.575,1],[8.425,8.575,0.04321798],["NaN","NaN","NaN"],[9.425,8.575,0.04321798],[9.425,9.425,0.04321798],[9.425,9.425,1],[9.425,8.575,1],[9.425,8.575,0.04321798],["NaN","NaN","NaN"],[1.425,9.575,0.04321798],[1.425,10.425,0.04321798],[1.425,10.425,0.7135577],[1.425,9.575,0.7135577],[1.425,9.575,0.04321798],["NaN","NaN","NaN"],[2.425,9.575,0.04321798],[2.425,10.425,0.04321798],[2.425,10.425,0.9967811],[2.425,9.575,0.9967811],[2.425,9.575,0.04321798],["NaN","NaN","NaN"],[3.425,9.575,0.04321798],[3.425,10.425,0.04321798],[3.425,10.425,0.9999999],[3.425,9.575,0.9999999],[3.425,9.575,0.04321798],["NaN","NaN","NaN"],[4.425,9.575,0.04321798],[4.425,10.425,0.04321798],[4.425,10.425,1],[4.425,9.575,1],[4.425,9.575,0.04321798],["NaN","NaN","NaN"],[5.425,9.575,0.04321798],[5.425,10.425,0.04321798],[5.425,10.425,1],[5.425,9.575,1],[5.425,9.575,0.04321798],["NaN","NaN","NaN"],[6.425,9.575,0.04321798],[6.425,10.425,0.04321798],[6.425,10.425,1],[6.425,9.575,1],[6.425,9.575,0.04321798],["NaN","NaN","NaN"],[7.425,9.575,0.04321798],[7.425,10.425,0.04321798],[7.425,10.425,1],[7.425,9.575,1],[7.425,9.575,0.04321798],["NaN","NaN","NaN"],[8.425,9.575,0.04321798],[8.425,10.425,0.04321798],[8.425,10.425,1],[8.425,9.575,1],[8.425,9.575,0.04321798],["NaN","NaN","NaN"],[9.425,9.575,0.04321798],[9.425,10.425,0.04321798],[9.425,10.425,1],[9.425,9.575,1],[9.425,9.575,0.04321798],["NaN","NaN","NaN"],[0.575,0.575,0.04321798],[1.425,0.575,0.04321798],[1.425,1.425,0.04321798],[0.575,1.425,0.04321798],[0.575,0.575,0.04321798],["NaN","NaN","NaN"],[1.575,0.575,0.06596308],[2.425,0.575,0.06596308],[2.425,1.425,0.06596308],[1.575,1.425,0.06596308],[1.575,0.575,0.06596308],["NaN","NaN","NaN"],[2.575,0.575,0.1021959],[3.425,0.575,0.1021959],[3.425,1.425,0.1021959],[2.575,1.425,0.1021959],[2.575,0.575,0.1021959],["NaN","NaN","NaN"],[3.575,0.575,0.1075265],[4.425,0.575,0.1075265],[4.425,1.425,0.1075265],[3.575,1.425,0.1075265],[3.575,0.575,0.1075265],["NaN","NaN","NaN"],[4.575,0.575,0.2316452],[5.425,0.575,0.2316452],[5.425,1.425,0.2316452],[4.575,1.425,0.2316452],[4.575,0.575,0.2316452],["NaN","NaN","NaN"],[5.575,0.575,0.2258095],[6.425,0.575,0.2258095],[6.425,1.425,0.2258095],[5.575,1.425,0.2258095],[5.575,0.575,0.2258095],["NaN","NaN","NaN"],[6.575,0.575,0.237572],[7.425,0.575,0.237572],[7.425,1.425,0.237572],[6.575,1.425,0.237572],[6.575,0.575,0.237572],["NaN","NaN","NaN"],[7.575,0.575,0.3690261],[8.425,0.575,0.3690261],[8.425,1.425,0.3690261],[7.575,1.425,0.3690261],[7.575,0.575,0.3690261],["NaN","NaN","NaN"],[8.575,0.575,0.3592858],[9.425,0.575,0.3592858],[9.425,1.425,0.3592858],[8.575,1.425,0.3592858],[8.575,0.575,0.3592858],["NaN","NaN","NaN"],[0.575,1.575,0.05915068],[1.425,1.575,0.05915068],[1.425,2.425,0.05915068],[0.575,2.425,0.05915068],[0.575,1.575,0.05915068],["NaN","NaN","NaN"],[1.575,1.575,0.1108488],[2.425,1.575,0.1108488],[2.425,2.425,0.1108488],[1.575,2.425,0.1108488],[1.575,1.575,0.1108488],["NaN","NaN","NaN"],[2.575,1.575,0.2030372],[3.425,1.575,0.2030372],[3.425,2.425,0.2030372],[2.575,2.425,0.2030372],[2.575,1.575,0.2030372],["NaN","NaN","NaN"],[3.575,1.575,0.217041],[4.425,1.575,0.217041],[4.425,2.425,0.217041],[3.575,2.425,0.217041],[3.575,1.575,0.217041],["NaN","NaN","NaN"],[4.575,1.575,0.5257641],[5.425,1.575,0.5257641],[5.425,2.425,0.5257641],[4.575,2.425,0.5257641],[4.575,1.575,0.5257641],["NaN","NaN","NaN"],[5.575,1.575,0.5127993],[6.425,1.575,0.5127993],[6.425,2.425,0.5127993],[5.575,2.425,0.5127993],[5.575,1.575,0.5127993],["NaN","NaN","NaN"],[6.575,1.575,0.5387281],[7.425,1.575,0.5387281],[7.425,2.425,0.5387281],[6.575,2.425,0.5387281],[6.575,1.575,0.5387281],["NaN","NaN","NaN"],[7.575,1.575,0.7713429],[8.425,1.575,0.7713429],[8.425,2.425,0.7713429],[7.575,2.425,0.7713429],[7.575,1.575,0.7713429],["NaN","NaN","NaN"],[8.575,1.575,0.7576808],[9.425,1.575,0.7576808],[9.425,2.425,0.7576808],[8.575,2.425,0.7576808],[8.575,1.575,0.7576808],["NaN","NaN","NaN"],[0.575,2.575,0.07725442],[1.425,2.575,0.07725442],[1.425,3.425,0.07725442],[0.575,3.425,0.07725442],[0.575,2.575,0.07725442],["NaN","NaN","NaN"],[1.575,2.575,0.1672267],[2.425,2.575,0.1672267],[2.425,3.425,0.1672267],[1.575,3.425,0.1672267],[1.575,2.575,0.1672267],["NaN","NaN","NaN"],[2.575,2.575,0.3312058],[3.425,2.575,0.3312058],[3.425,3.425,0.3312058],[2.575,3.425,0.3312058],[2.575,2.575,0.3312058],["NaN","NaN","NaN"],[3.575,2.575,0.3553697],[4.425,2.575,0.3553697],[4.425,3.425,0.3553697],[3.575,3.425,0.3553697],[3.575,2.575,0.3553697],["NaN","NaN","NaN"],[4.575,2.575,0.781293],[5.425,2.575,0.781293],[5.425,3.425,0.781293],[4.575,3.425,0.781293],[4.575,2.575,0.781293],["NaN","NaN","NaN"],[5.575,2.575,0.76804],[6.425,2.575,0.76804],[6.425,3.425,0.76804],[5.575,3.425,0.76804],[5.575,2.575,0.76804],["NaN","NaN","NaN"],[6.575,2.575,0.7941267],[7.425,2.575,0.7941267],[7.425,3.425,0.7941267],[6.575,3.425,0.7941267],[6.575,2.575,0.7941267],["NaN","NaN","NaN"],[7.575,2.575,0.9553382],[8.425,2.575,0.9553382],[8.425,3.425,0.9553382],[7.575,3.425,0.9553382],[7.575,2.575,0.9553382],["NaN","NaN","NaN"],[8.575,2.575,0.9493445],[9.425,2.575,0.9493445],[9.425,3.425,0.9493445],[8.575,3.425,0.9493445],[8.575,2.575,0.9493445],["NaN","NaN","NaN"],[0.575,3.575,0.1025677],[1.425,3.575,0.1025677],[1.425,4.425,0.1025677],[0.575,4.425,0.1025677],[0.575,3.575,0.1025677],["NaN","NaN","NaN"],[1.575,3.575,0.2501325],[2.425,3.575,0.2501325],[2.425,4.425,0.2501325],[1.575,4.425,0.2501325],[1.575,3.575,0.2501325],["NaN","NaN","NaN"],[2.575,3.575,0.5039053],[3.425,3.575,0.5039053],[3.425,4.425,0.5039053],[2.575,4.425,0.5039053],[2.575,3.575,0.5039053],["NaN","NaN","NaN"],[3.575,3.575,0.537781],[4.425,3.575,0.537781],[4.425,4.425,0.537781],[3.575,4.425,0.537781],[3.575,3.575,0.537781],["NaN","NaN","NaN"],[4.575,3.575,0.9430385],[5.425,3.575,0.9430385],[5.425,4.425,0.9430385],[4.575,4.425,0.9430385],[4.575,3.575,0.9430385],["NaN","NaN","NaN"],[5.575,3.575,0.9361931],[6.425,3.575,0.9361931],[6.425,4.425,0.9361931],[5.575,4.425,0.9361931],[5.575,3.575,0.9361931],["NaN","NaN","NaN"],[6.575,3.575,0.9493053],[7.425,3.575,0.9493053],[7.425,4.425,0.9493053],[6.575,4.425,0.9493053],[6.575,3.575,0.9493053],["NaN","NaN","NaN"],[7.575,3.575,0.9972544],[8.425,3.575,0.9972544],[8.425,4.425,0.9972544],[7.575,4.425,0.9972544],[7.575,3.575,0.9972544],["NaN","NaN","NaN"],[8.575,3.575,0.9965154],[9.425,3.575,0.9965154],[9.425,4.425,0.9965154],[8.575,4.425,0.9965154],[8.575,3.575,0.9965154],["NaN","NaN","NaN"],[0.575,4.575,0.1352801],[1.425,4.575,0.1352801],[1.425,5.425,0.1352801],[0.575,5.425,0.1352801],[0.575,4.575,0.1352801],["NaN","NaN","NaN"],[1.575,4.575,0.3567418],[2.425,4.575,0.3567418],[2.425,5.425,0.3567418],[1.575,5.425,0.3567418],[1.575,4.575,0.3567418],["NaN","NaN","NaN"],[2.575,4.575,0.6844499],[3.425,4.575,0.6844499],[3.425,5.425,0.6844499],[2.575,5.425,0.6844499],[2.575,4.575,0.6844499],["NaN","NaN","NaN"],[3.575,4.575,0.7209224],[4.425,4.575,0.7209224],[4.425,5.425,0.7209224],[3.575,5.425,0.7209224],[3.575,4.575,0.7209224],["NaN","NaN","NaN"],[4.575,4.575,0.9923994],[5.425,4.575,0.9923994],[5.425,5.425,0.9923994],[4.575,5.425,0.9923994],[4.575,4.575,0.9923994],["NaN","NaN","NaN"],[5.575,4.575,0.9907784],[6.425,4.575,0.9907784],[6.425,5.425,0.9907784],[5.575,5.425,0.9907784],[5.575,4.575,0.9907784],["NaN","NaN","NaN"],[6.575,4.575,0.9937668],[7.425,4.575,0.9937668],[7.425,5.425,0.9937668],[6.575,5.425,0.9937668],[6.575,4.575,0.9937668],["NaN","NaN","NaN"],[7.575,4.575,0.999954],[8.425,4.575,0.999954],[8.425,5.425,0.999954],[7.575,5.425,0.999954],[7.575,4.575,0.999954],["NaN","NaN","NaN"],[8.575,4.575,0.9999315],[9.425,4.575,0.9999315],[9.425,5.425,0.9999315],[8.575,5.425,0.9999315],[8.575,4.575,0.9999315],["NaN","NaN","NaN"],[0.575,5.575,0.1997644],[1.425,5.575,0.1997644],[1.425,6.425,0.1997644],[0.575,6.425,0.1997644],[0.575,5.575,0.1997644],["NaN","NaN","NaN"],[1.575,5.575,0.5454682],[2.425,5.575,0.5454682],[2.425,6.425,0.5454682],[1.575,6.425,0.5454682],[1.575,5.575,0.5454682],["NaN","NaN","NaN"],[2.575,5.575,0.8883711],[3.425,5.575,0.8883711],[3.425,6.425,0.8883711],[2.575,6.425,0.8883711],[2.575,5.575,0.8883711],["NaN","NaN","NaN"],[3.575,5.575,0.912325],[4.425,5.575,0.912325],[4.425,6.425,0.912325],[3.575,6.425,0.912325],[3.575,5.575,0.912325],["NaN","NaN","NaN"],[4.575,5.575,0.9999133],[5.425,5.575,0.9999133],[5.425,6.425,0.9999133],[4.575,6.425,0.9999133],[4.575,5.575,0.9999133],["NaN","NaN","NaN"],[5.575,5.575,0.9998751],[6.425,5.575,0.9998751],[6.425,6.425,0.9998751],[5.575,6.425,0.9998751],[5.575,5.575,0.9998751],["NaN","NaN","NaN"],[6.575,5.575,0.9999403],[7.425,5.575,0.9999403],[7.425,6.425,0.9999403],[6.575,6.425,0.9999403],[6.575,5.575,0.9999403],["NaN","NaN","NaN"],[7.575,5.575,1],[8.425,5.575,1],[8.425,6.425,1],[7.575,6.425,1],[7.575,5.575,1],["NaN","NaN","NaN"],[8.575,5.575,1],[9.425,5.575,1],[9.425,6.425,1],[8.575,6.425,1],[8.575,5.575,1],["NaN","NaN","NaN"],[0.575,6.575,0.2790861],[1.425,6.575,0.2790861],[1.425,7.425,0.2790861],[0.575,7.425,0.2790861],[0.575,6.575,0.2790861],["NaN","NaN","NaN"],[1.575,6.575,0.722779],[2.425,6.575,0.722779],[2.425,7.425,0.722779],[1.575,7.425,0.722779],[1.575,6.575,0.722779],["NaN","NaN","NaN"],[2.575,6.575,0.9743328],[3.425,6.575,0.9743328],[3.425,7.425,0.9743328],[2.575,7.425,0.9743328],[2.575,6.575,0.9743328],["NaN","NaN","NaN"],[3.575,6.575,0.9828939],[4.425,6.575,0.9828939],[4.425,7.425,0.9828939],[3.575,7.425,0.9828939],[3.575,6.575,0.9828939],["NaN","NaN","NaN"],[4.575,6.575,0.9999998],[5.425,6.575,0.9999998],[5.425,7.425,0.9999998],[4.575,7.425,0.9999998],[4.575,6.575,0.9999998],["NaN","NaN","NaN"],[5.575,6.575,0.9999996],[6.425,6.575,0.9999996],[6.425,7.425,0.9999996],[5.575,7.425,0.9999996],[5.575,6.575,0.9999996],["NaN","NaN","NaN"],[6.575,6.575,0.9999999],[7.425,6.575,0.9999999],[7.425,7.425,0.9999999],[6.575,7.425,0.9999999],[6.575,6.575,0.9999999],["NaN","NaN","NaN"],[7.575,6.575,1],[8.425,6.575,1],[8.425,7.425,1],[7.575,7.425,1],[7.575,6.575,1],["NaN","NaN","NaN"],[8.575,6.575,1],[9.425,6.575,1],[9.425,7.425,1],[8.575,7.425,1],[8.575,6.575,1],["NaN","NaN","NaN"],[0.575,7.575,0.428676],[1.425,7.575,0.428676],[1.425,8.425,0.428676],[0.575,8.425,0.428676],[0.575,7.575,0.428676],["NaN","NaN","NaN"],[1.575,7.575,0.9105917],[2.425,7.575,0.9105917],[2.425,8.425,0.9105917],[1.575,8.425,0.9105917],[1.575,7.575,0.9105917],["NaN","NaN","NaN"],[2.575,7.575,0.9990419],[3.425,7.575,0.9990419],[3.425,8.425,0.9990419],[2.575,8.425,0.9990419],[2.575,7.575,0.9990419],["NaN","NaN","NaN"],[3.575,7.575,0.9995526],[4.425,7.575,0.9995526],[4.425,8.425,0.9995526],[3.575,8.425,0.9995526],[3.575,7.575,0.9995526],["NaN","NaN","NaN"],[4.575,7.575,1],[5.425,7.575,1],[5.425,8.425,1],[4.575,8.425,1],[4.575,7.575,1],["NaN","NaN","NaN"],[5.575,7.575,1],[6.425,7.575,1],[6.425,8.425,1],[5.575,8.425,1],[5.575,7.575,1],["NaN","NaN","NaN"],[6.575,7.575,1],[7.425,7.575,1],[7.425,8.425,1],[6.575,8.425,1],[6.575,7.575,1],["NaN","NaN","NaN"],[7.575,7.575,1],[8.425,7.575,1],[8.425,8.425,1],[7.575,8.425,1],[7.575,7.575,1],["NaN","NaN","NaN"],[8.575,7.575,1],[9.425,7.575,1],[9.425,8.425,1],[8.575,8.425,1],[8.575,7.575,1],["NaN","NaN","NaN"],[0.575,8.575,0.5885873],[1.425,8.575,0.5885873],[1.425,9.425,0.5885873],[0.575,9.425,0.5885873],[0.575,8.575,0.5885873],["NaN","NaN","NaN"],[1.575,8.575,0.981856],[2.425,8.575,0.981856],[2.425,9.425,0.981856],[1.575,9.425,0.981856],[1.575,8.575,0.981856],["NaN","NaN","NaN"],[2.575,8.575,0.9999893],[3.425,8.575,0.9999893],[3.425,9.425,0.9999893],[2.575,9.425,0.9999893],[2.575,8.575,0.9999893],["NaN","NaN","NaN"],[3.575,8.575,0.9999969],[4.425,8.575,0.9999969],[4.425,9.425,0.9999969],[3.575,9.425,0.9999969],[3.575,8.575,0.9999969],["NaN","NaN","NaN"],[4.575,8.575,1],[5.425,8.575,1],[5.425,9.425,1],[4.575,9.425,1],[4.575,8.575,1],["NaN","NaN","NaN"],[5.575,8.575,1],[6.425,8.575,1],[6.425,9.425,1],[5.575,9.425,1],[5.575,8.575,1],["NaN","NaN","NaN"],[6.575,8.575,1],[7.425,8.575,1],[7.425,9.425,1],[6.575,9.425,1],[6.575,8.575,1],["NaN","NaN","NaN"],[7.575,8.575,1],[8.425,8.575,1],[8.425,9.425,1],[7.575,9.425,1],[7.575,8.575,1],["NaN","NaN","NaN"],[8.575,8.575,1],[9.425,8.575,1],[9.425,9.425,1],[8.575,9.425,1],[8.575,8.575,1],["NaN","NaN","NaN"],[0.575,9.575,0.7135577],[1.425,9.575,0.7135577],[1.425,10.425,0.7135577],[0.575,10.425,0.7135577],[0.575,9.575,0.7135577],["NaN","NaN","NaN"],[1.575,9.575,0.9967811],[2.425,9.575,0.9967811],[2.425,10.425,0.9967811],[1.575,10.425,0.9967811],[1.575,9.575,0.9967811],["NaN","NaN","NaN"],[2.575,9.575,0.9999999],[3.425,9.575,0.9999999],[3.425,10.425,0.9999999],[2.575,10.425,0.9999999],[2.575,9.575,0.9999999],["NaN","NaN","NaN"],[3.575,9.575,1],[4.425,9.575,1],[4.425,10.425,1],[3.575,10.425,1],[3.575,9.575,1],["NaN","NaN","NaN"],[4.575,9.575,1],[5.425,9.575,1],[5.425,10.425,1],[4.575,10.425,1],[4.575,9.575,1],["NaN","NaN","NaN"],[5.575,9.575,1],[6.425,9.575,1],[6.425,10.425,1],[5.575,10.425,1],[5.575,9.575,1],["NaN","NaN","NaN"],[6.575,9.575,1],[7.425,9.575,1],[7.425,10.425,1],[6.575,10.425,1],[6.575,9.575,1],["NaN","NaN","NaN"],[7.575,9.575,1],[8.425,9.575,1],[8.425,10.425,1],[7.575,10.425,1],[7.575,9.575,1],["NaN","NaN","NaN"],[8.575,9.575,1],[9.425,9.575,1],[9.425,10.425,1],[8.575,10.425,1],[8.575,9.575,1],["NaN","NaN","NaN"]],"ignoreExtent":false,"flags":64},"24":{"id":24,"type":"text","material":{},"vertices":[[1,-1,0],[2,-1,0],[3,-1,0],[4,-1,0],[5,-1,0],[6,-1,0],[7,-1,0],[8,-1,0],[9,-1,0],[0.5,1,0],[0.5,2,0],[0.5,3,0],[0.5,4,0],[0.5,5,0],[0.5,6,0],[0.5,7,0],[0.5,8,0],[0.5,9,0],[0.5,10,0]],"colors":[[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1]],"texts":[["1213-Ep-9-KODE"],["8_9-DiHETrE"],["SUM-TriHOME"],["8-HETE"],["5_6-DiHETrE"],["9-HEPE"],["1415-EpETrE"],["14_15-DiHETE"],["9-HETE"],["3"],["6"],["10"],["16"],["24"],["40"],["60"],["100"],["150"],["200"]],"cex":[[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]],"adj":[[1,0.5]],"centers":[[1,-1,0],[2,-1,0],[3,-1,0],[4,-1,0],[5,-1,0],[6,-1,0],[7,-1,0],[8,-1,0],[9,-1,0],[0.5,1,0],[0.5,2,0],[0.5,3,0],[0.5,4,0],[0.5,5,0],[0.5,6,0],[0.5,7,0],[0.5,8,0],[0.5,9,0],[0.5,10,0]],"family":[["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"],["sans"]],"font":[[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1]],"ignoreExtent":false,"flags":2064},"25":{"id":25,"type":"text","material":{},"vertices":[[4.9625,12.36154,1.1695]],"colors":[[0,0,0,1]],"texts":[["Power Array (Interactive 3D Model)"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[4.9625,12.36154,1.1695]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"26":{"id":26,"type":"text","material":{},"vertices":[[4.9625,-2.936538,-0.1695]],"colors":[[0,0,0,1]],"texts":[[""]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[4.9625,-2.936538,-0.1695]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"27":{"id":27,"type":"text","material":{},"vertices":[[-1.012788,4.7125,-0.1695]],"colors":[[0,0,0,1]],"texts":[[""]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-1.012788,4.7125,-0.1695]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"28":{"id":28,"type":"text","material":{},"vertices":[[-1.012788,-2.936538,0.5]],"colors":[[0,0,0,1]],"texts":[[""]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-1.012788,-2.936538,0.5]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"29":{"id":29,"type":"lines","material":{},"vertices":[[0.366125,-1.171375,0],[0.366125,-1.171375,1],[0.366125,-1.171375,0],[0.1363062,-1.465569,0],[0.366125,-1.171375,0.2],[0.1363062,-1.465569,0.2],[0.366125,-1.171375,0.4],[0.1363062,-1.465569,0.4],[0.366125,-1.171375,0.6],[0.1363062,-1.465569,0.6],[0.366125,-1.171375,0.8],[0.1363062,-1.465569,0.8],[0.366125,-1.171375,1],[0.1363062,-1.465569,1]],"colors":[[0,0,0,1]],"centers":[[0.366125,-1.171375,0.5],[0.2512156,-1.318472,0],[0.2512156,-1.318472,0.2],[0.2512156,-1.318472,0.4],[0.2512156,-1.318472,0.6],[0.2512156,-1.318472,0.8],[0.2512156,-1.318472,1]],"ignoreExtent":true,"flags":64},"30":{"id":30,"type":"text","material":{},"vertices":[[-0.3233313,-2.053956,0],[-0.3233313,-2.053956,0.2],[-0.3233313,-2.053956,0.4],[-0.3233313,-2.053956,0.6],[-0.3233313,-2.053956,0.8],[-0.3233313,-2.053956,1]],"colors":[[0,0,0,1]],"texts":[["0.0"],["0.2"],["0.4"],["0.6"],["0.8"],["1.0"]],"cex":[[1]],"adj":[[0.5,0.5]],"centers":[[-0.3233313,-2.053956,0],[-0.3233313,-2.053956,0.2],[-0.3233313,-2.053956,0.4],[-0.3233313,-2.053956,0.6],[-0.3233313,-2.053956,0.8],[-0.3233313,-2.053956,1]],"family":[["sans"]],"font":[[1]],"ignoreExtent":true,"flags":2064},"10":{"id":10,"type":"light","vertices":[[0,0,1]],"colors":[[1,1,1,1],[1,1,1,1],[1,1,1,1]],"viewpoint":true,"finite":false},"9":{"id":9,"type":"background","material":{"lit":true,"fog":true},"colors":[[0.2980392,0.2980392,0.2980392,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"11":{"id":11,"type":"background","material":{"back":"lines"},"colors":[[1,1,1,1]],"centers":[[0,0,0]],"sphere":false,"fogtype":"none","flags":0},"12":{"id":12,"type":"subscene","par3d":{"antialias":8,"FOV":0,"ignoreExtent":false,"listeners":12,"mouseMode":{"left":"trackball","right":"zoom","middle":"user","wheel":"pull"},"observer":[0,0,20.17666],"modelMatrix":[[4,0,0,-0.88],[0,-1.748456e-07,20,-10],[0,-4,-8.742278e-07,-19.17666],[0,0,0,1]],"projMatrix":[["NaN",0,0,"NaN"],["NaN",0.1982489,0,"NaN"],["NaN",0,-0.09912442,"NaN"],["NaN",0,0,"NaN"]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,-4.371139e-08,1,0],[0,-1,-4.371139e-08,0],[0,0,0,1]],"userProjection":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],"scale":[4,4,20],"viewport":{"x":0,"y":0,"width":0,"height":0},"zoom":0.5,"bbox":[0,0.44,0,0.5,0,1],"windowRect":[100,100,772,580],"family":"sans","font":1,"cex":1,"useFreeType":false,"fontname":"TT Arial","maxClipPlanes":8,"glVersion":4.6,"activeSubscene":0},"embeddings":{"viewport":"replace","projection":"replace","model":"replace","mouse":"inherit"},"objects":[14,15,16,17,18,19,10],"parent":6,"subscenes":[],"flags":2642},"13":{"id":13,"type":"subscene","par3d":{"antialias":8,"FOV":30,"ignoreExtent":false,"listeners":13,"mouseMode":{"left":"trackball","right":"zoom","middle":"user","wheel":"pull"},"observer":[0,0,7.104477],"modelMatrix":[[0.2222222,0,0,-1.102778],[0,0.06840403,1.964277,-1.304493],[0,-0.1879385,0.7149385,-6.576286],[0,0,0,1]],"projMatrix":[["NaN",0,"NaN",0],["NaN",3.732051,"NaN",0],["NaN",0,"NaN",-25.61082],["NaN",0,"NaN",0]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,0.3420201,0.9396926,0],[0,-0.9396926,0.3420201,0],[0,0,0,1]],"userProjection":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],"scale":[0.2222222,0.2,2.09034],"viewport":{"x":0,"y":0,"width":0,"height":0},"zoom":1,"bbox":[0.5,9.425,-1,10.425,0,1],"windowRect":[100,100,772,580],"family":"sans","font":1,"cex":1,"useFreeType":false,"fontname":"TT Arial","maxClipPlanes":8,"glVersion":4.6,"activeSubscene":0},"embeddings":{"viewport":"replace","projection":"replace","model":"replace","mouse":"inherit"},"objects":[22,23,24,25,26,27,28,29,30,10],"parent":6,"subscenes":[],"flags":2642},"6":{"id":6,"type":"subscene","par3d":{"antialias":8,"FOV":30,"ignoreExtent":false,"listeners":6,"mouseMode":{"left":"trackball","right":"zoom","middle":"user","wheel":"pull"},"observer":[0,0,3.863703],"modelMatrix":[[1,0,0,0],[0,0.3420202,0.9396926,0],[0,-0.9396926,0.3420202,-3.863703],[0,0,0,1]],"projMatrix":[[2.665751,0,0,0],[0,3.732051,0,0],[0,0,-3.863703,-13.9282],[0,0,-1,0]],"skipRedraw":false,"userMatrix":[[1,0,0,0],[0,0.3420201,0.9396926,0],[0,-0.9396926,0.3420201,0],[0,0,0,1]],"userProjection":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],"scale":[1,1,1],"viewport":{"x":0,"y":0,"width":1,"height":1},"zoom":1,"bbox":[3.402823e+38,-3.402823e+38,3.402823e+38,-3.402823e+38,3.402823e+38,-3.402823e+38],"windowRect":[100,100,772,580],"family":"sans","font":1,"cex":1,"useFreeType":false,"fontname":"TT Arial","maxClipPlanes":8,"glVersion":4.6,"activeSubscene":0},"embeddings":{"viewport":"replace","projection":"replace","model":"replace","mouse":"replace"},"objects":[11,10,12,13],"subscenes":[12,13],"flags":2642}},"snapshot":"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAqAAAAHgCAIAAAD17khjAAAAHXRFWHRTb2Z0d2FyZQBSL1JHTCBwYWNrYWdlL2xpYnBuZ7GveO8AAAecSURBVHic7dVBCQAwDMDA+jfdmhgMwp2C/DILAOTM7wAA4D2DB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4AggweAIIMHgCCDB4CgA6dpGKtfMTX3AAAAAElFTkSuQmCC","width":673,"height":481,"sphereVerts":{"vb":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1],[0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1950903,-0.3826834,-0.5555702,-0.7071068,-0.8314696,-0.9238795,-0.9807853,-1,-0.9807853,-0.9238795,-0.8314696,-0.7071068,-0.5555702,-0.3826834,-0.1950903,-0,-0,-0.18024,-0.3535534,-0.51328,-0.6532815,-0.7681778,-0.8535534,-0.9061274,-0.9238795,-0.9061274,-0.8535534,-0.7681778,-0.6532815,-0.51328,-0.3535534,-0.18024,-0,-0,-0.1379497,-0.2705981,-0.3928475,-0.5,-0.5879378,-0.6532815,-0.6935199,-0.7071068,-0.6935199,-0.6532815,-0.5879378,-0.5,-0.3928475,-0.2705981,-0.1379497,-0,-0,-0.07465783,-0.1464466,-0.2126075,-0.2705981,-0.3181896,-0.3535534,-0.3753303,-0.3826834,-0.3753303,-0.3535534,-0.3181896,-0.2705981,-0.2126075,-0.1464466,-0.07465783,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.07465783,0.1464466,0.2126075,0.2705981,0.3181896,0.3535534,0.3753303,0.3826834,0.3753303,0.3535534,0.3181896,0.2705981,0.2126075,0.1464466,0.07465783,0,0,0.1379497,0.2705981,0.3928475,0.5,0.5879378,0.6532815,0.6935199,0.7071068,0.6935199,0.6532815,0.5879378,0.5,0.3928475,0.2705981,0.1379497,0,0,0.18024,0.3535534,0.51328,0.6532815,0.7681778,0.8535534,0.9061274,0.9238795,0.9061274,0.8535534,0.7681778,0.6532815,0.51328,0.3535534,0.18024,0,0,0.1950903,0.3826834,0.5555702,0.7071068,0.8314696,0.9238795,0.9807853,1,0.9807853,0.9238795,0.8314696,0.7071068,0.5555702,0.3826834,0.1950903,0]],"it":[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270],[17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288],[18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271]],"material":[],"normals":null,"texcoords":[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.1875,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.3125,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.4375,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.5625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.625,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.6875,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.75,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.8125,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.875,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,0.9375,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1,0,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,0.625,0.6875,0.75,0.8125,0.875,0.9375,1]],"meshColor":"vertices"},"context":{"shiny":false,"rmarkdown":"md_document"},"crosstalk":{"key":[],"group":[],"id":[],"options":[]}});
powergrid_rgl.prefix = "powergrid_";
</script>
<p id="powergrid_debug">
You must enable Javascript to view this page properly.
</p>
<script>powergrid_rgl.start();</script>

------------------------------------------------------------------------

<br>

I used the following R libraries to program rMetabolomics:
==========================================================

ggplot2  
reshape2  
plyr  
scales  
gridExtra  
grid  
gtable  
pwr  
plot3D  
rgl  
plot3Drgl

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19042)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] plot3Drgl_1.0.1      plot3D_1.3           pwr_1.3-0           
    ##  [4] gtable_0.3.0         gridExtra_2.3        scales_1.1.1        
    ##  [7] plyr_1.8.6           reshape2_1.4.4       ggplot2_3.3.2       
    ## [10] rgl_0.100.54         knitr_1.29           RevoUtils_11.0.2    
    ## [13] RevoUtilsMath_11.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0        xfun_0.15               purrr_0.3.4            
    ##  [4] colorspace_1.4-1        vctrs_0.3.2             generics_0.0.2         
    ##  [7] miniUI_0.1.1.1          htmltools_0.5.0         yaml_2.2.1             
    ## [10] rlang_0.4.7             manipulateWidget_0.10.1 later_1.1.0.1          
    ## [13] pillar_1.4.6            glue_1.4.1              withr_2.2.0            
    ## [16] lifecycle_0.2.0         stringr_1.4.0           munsell_0.5.0          
    ## [19] htmlwidgets_1.5.1       evaluate_0.14           labeling_0.3           
    ## [22] misc3d_0.8-4            fastmap_1.0.1           httpuv_1.5.4           
    ## [25] crosstalk_1.1.0.1       Rcpp_1.0.5              xtable_1.8-4           
    ## [28] promises_1.1.1          webshot_0.5.2           jsonlite_1.7.0         
    ## [31] farver_2.0.3            mime_0.9                digest_0.6.25          
    ## [34] stringi_1.4.6           dplyr_1.0.0             shiny_1.5.0            
    ## [37] tools_4.0.2             magrittr_1.5            tibble_3.0.3           
    ## [40] crayon_1.3.4            pkgconfig_2.0.3         ellipsis_0.3.1         
    ## [43] rmarkdown_2.3           R6_2.4.1                compiler_4.0.2
