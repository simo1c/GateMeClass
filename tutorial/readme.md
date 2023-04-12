# GateMeClass Tutorial

## DESCRIPTION:

The concept behind GateMeClass was to build a user friendly tool (in R programming language), able to fully automate both the definition of a marker table and the cell annotation. GateMeClass stands for “Gate-Mining-and-Classification” and it's a tool applicable to cytometry data.

GateMeClass is composed by two main modules, the training module, and the annotation module. Both modules can be used individually or in sequence according to user needs.

GateMeClass is able to learn a marker table from an external annotated dataset and is able to apply the cell labels to a dataset of interest.

GateMeClass gives you also the option to manually define a marker table without an annotated dataset and to apply it directly on the dataset of interest to perform cells classification.

Here we present an easy tutorial to guide you through its applicability on real biological datasets. 


## INSTALLATION:

Before using GateMeClass, it is strongly reconmended to follow these steps.

### Step 1 "R-Libraries":

This is a list of all the packages that are fundamental to be able to run GateMeClass.

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

In order to run GateMeClass, you have to import the file called *GateMeClass.R*, which has to be present in the Working Directory.
It contains the architecture of the tool, its functions with all the default parameters and all the control procedures.
To import all the functions of the tool, it's mandatory to execute the line:

```
source("GateMeClass.R")

```


## EXECUTION:

### Classification using a manually defined marker table:

The annotation function of GateMeClass is able to annotate the cells of a dataset of interest using a manual gated marker table.
The first thing to do is to assign the dataset of interest to a variable (*se*) and to extract both the expression matrix (*m*) and the labels (*lab*) from it.
Both the dataset and the manually gated marker table have to be in the working directory.

```
se <- readRDS("Levine32.Rds")    
m <- se@assays@data$exprs    
lab <- se$labels
```
 
To use the annotation function of GateMeClass, you have to set correctly the input variables.
Assign to the *testing_set* variable the matrix of expression (*m*) and to the variable *gate* the manually created marker table.

```
testing_set<-m
gate<-read_xls("Levine32_marker_table.xls")
```

After this, you can run the annotation function.


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


Parameters description:
```
testing_set         : Matrix of expression
marker_table        : Manually curated table of markers
reject_option       : default "T"
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: default "V" (Variance) or "E" (Equal)
RSS                 : RSS (Ranked-Set-Sampling) observations, used to have a better resolution: default "T"
k                   : k-NN (k-Nearest-Neighbour), used to refine the uncertain labels to the most similar already annotated
sampling            : default "0.1"
verbose             : default "T"
seed                : default = 1
```

The output of this function will be saved in a variable (*res*). 
This object is relatively complex, because it contains different informations.


1) labels           (class: character)
2) marker_table     (class: tbl_df ; tbl ; data.frame)
3) cell_signatures  (class: data.frame)


These different fields are easily accessible with the "$" character:





### Classification extracting the marker table from a reference annotated dataset:

[Training Function Annotation function]

Estrazione tabella  (funzione di training)
Richiamo tabella estratta (annotation)

### All-In-One


Se l’utente preferisce fare tutto in una volta --> Definizione parametri per farla one-shot (annotation)

[PARAMETERS]


