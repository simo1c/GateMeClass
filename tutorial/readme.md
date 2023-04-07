# GateMeClass Tutorial

## INSTALLATION:

Before using GateMeClass, is strongly reconmended to run across the following steps.

### Step 1 "R-Libraries":

```
install.packages("stringi")
install.packages("data.table")
install.packages("moments")
install.packages("mclust")
install.packages("stringr")
install.packages("scales")
install.packages("recipes")
install.packages("dplyr")
install.packages("caret", dependencies = c("Depends", "Suggests"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("batchelor")

```

### Step 2 "Source":

To use the functions of the package, it's mandatory to execute the line:

```
source("GateMeClass.R")

```

The file called "GateMeClass.R" has to be present in the Working Directory.

## EXECUTION:

### Classification using a manually defined marker table:

The annotation function of GateMeClass is able to attonate the cells of a dataset of interest using a manual gated marker table.
The first thing to do is to assign the dataset of interest to a variable (se) and to extract both the expression matrix (m) and the labels (lab) from it.
Both the dataset and the manually gated marker table have to be in the working directory.

```
se <- readRDS("Levine32.Rds")    
m <- se@assays@data$exprs    
lab <- se$labels
```
 
To use the annotation function of GateMeClass, you have to set the input variables. Assign to the "testing_set" variable the matrix (m) and to the variable "gate" the manually gated marker table.

```
testing_set<-m
gate<-read_xls("Levine32_marker_table.xls")
```

After this, you can run the annotation function. [WITH DEFAULT PARAMETERS...]


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


The output of this function will be saved in "res". This object contains different informations: 1) labels; 2)marker-table; 3) cell-signatures. To access one of these elements you can use the dollar sign "$".


### Classification extracting the marker table from a reference annotated dataset:

[Training Function Annotation function]


### All-In-One

[PARAMETERS]

