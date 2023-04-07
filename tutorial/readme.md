# GateMeClass Tutorial

## INSTALLATION:

Before using GateMeClass, is strongly reconmended to run across the following steps.

### Step 1 "R-Libraries":

```
library(stringi)
library(data.table)
library(moments)
library(mclust)
library(stringr)
library(scales)
library(caret)
library(batchelor)
library(recipes)
library(dplyr)

```

### Step 2 "Source":

To use the functions of the package, it's mandatory to execute the line:

```
source("GateMeClass.R")

```

The file called "GateMeClass.R" has to be present in the Working Directory.

## EXECUTION:

### Classification using a manually defined marker table:

The first thing to do is to assign the dataset of interest to a variable (se) and to extract both the expression matrix (m) and the labels (lab) from it.

Both the dataset and the marker table have to be in the working directory.

```
se <- readRDS("Levine32.Rds")    
m <- se@assays@data$exprs    
lab <- se$labels
```

Then, to use the function of GateMeClass, assign to the "testing_set" variable the matrix (m) and to the variable "gate" the marker table.
```
testing_set<-m
gate<-read_xls("Levine32_marker_table.xls")
```

```
res <- GateMeClass_annotate(testing_set,
                            marker_table = gate,
                            reject_option = F,
                            GMM_parameterization = "V",
                            RSS = T,
                            k = 20,				
                            sampling = 0.1,
                            verbose = T,
                            seed = 1)
```
The output of this function will be the "res" object. This object contains labels, marker-table and cell-signatures. To access to one of theese elements use "$".

