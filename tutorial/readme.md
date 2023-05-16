# GateMeClass Tutorial

## DESCRIPTION

The idea behind GateMeClass was to create a user-friendly tool, implemented in the R programming language, that automates the process of defining a marker table and cell annotation. GateMeClass, short for "Gate-Mining-and-Classification," is specifically designed for cytometry data.

GateMeClass consists of two main modules: the training module and the annotation module. These modules can be utilized independently or in a sequence depending on the user's requirements.

GateMeClass has the capability to learn a marker table from an externally annotated dataset and apply the corresponding cell labels to a different dataset of interest.

Additionally, GateMeClass provides the option to manually define a marker table without relying on an annotated dataset. This marker table can then be directly applied to the dataset of interest for cell classification.

To demonstrate its practicality with real biological datasets, we have prepared an easy-to-follow tutorial.

## INSTALLATION

It is highly recommended to follow these steps prior to utilizing GateMeClass.

### Step 1 "R-libraries"

This is a list of all the R packages that are fundamental to be able to run GateMeClass.

```
install.packages("stringi")
install.packages("data.table")
install.packages("moments")
install.packages("mclust")
install.packages("stringr")
install.packages("scales")
install.packages("recipes")
install.packages("dplyr")
install.packages("RcppAlgos")
install.packages("caret", dependencies = c("Depends", "Suggests"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("batchelor")

```


### Step 2 "Source"

In order to run GateMeClass, you have to import the file called *GateMeClass.R*, which has to be present in the Working Directory.

It contains the architecture of the tool, its functions with default parameters and all the control procedures.

To import all the functions of the tool, it's mandatory to execute the line:

```
source("GateMeClass.R")

```

## EXECUTION

### 1: Classification using a manually defined marker table

In this first application of GateMeClass we will annotate the cells of a dataset of interest using a manual gated marker table.

Prior to this, it is necessary to create several variables that will be utilized throughout all sections of our tutorial.

We have to assign the dataset of interest to a variable `se` and to extract both the expression matrix `m` and the labels `lab` from it.

NOTE: both the dataset and the manually gated marker table files have to be present in the working directory!

```
se <- readRDS("Levine32.Rds")    
m <- se@assays@data$exprs    
lab <- se$labels
```
 
To use the annotation module of GateMeClass, you have to set correctly the fundamental input variables of the *GateMeClass_annotate* annotation function.

Assign to the `testing_set` variable the matrix of expression `m` and to the variable `gate` the manually created marker table file.

```
testing_set <- m
gate <- read_xls("Levine32_marker_table.xls")
```

After this, you can run the *GateMeClass_annotate* annotation function:


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
reject_option       : default = "T"
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Varing) or "E" (Equal), default = "V"
RSS                 : RSS (Ranked-Set-Sampling) observations, used to have a better resolution: default = "T"
k                   : k-NN (k-Nearest-Neighbour), used to refine the uncertain labels to the most similar already annotated
sampling            : default = "0.1"
verbose             : default = "T"
seed                : default = "1"
```

The output of this function will be saved in a variable `res`. 

This object is relatively complex, it contains different informations accessible with the `$` character:.


1) labels           (class: character)
2) marker_table     (class: tbl_df ; tbl ; data.frame)
3) cell_signatures  (class: data.frame)


### 2: Classification extracting the marker table from a reference annotated dataset

In this second application of GateMeClass, we will annotate the cells of a dataset of interest by utilizing the training module and a training set to construct a marker table.

We assign to the variable `training_set` the expression matrix `m`, and to the variable `training_set_lab` the labels fron the se object.

```
training_set <- m
training_set_lab <- lab
```

The next step involves utilizing the *GateMeClass_train* training function to construct the `gate` variable, which encompasses the gating strategy employed for the training set.

```
gate <- GateMeClass_train(training_set,
                          training_set_lab,
                          RSS = T,
                          GMM_parameterization = "V",
                          verbose = T, 
                          seed = 1)
```
Parameters description:
```
training_set        : Matrix of expression
training_set_lab    : Labels from the dataset
RSS                 : RSS (Ranked-Set-Sampling) observations, used to have a better resolution: default = "T"
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Variance) or "E" (Equal), default = "V"
verbose             : default = "T"
seed                : default = "1"
```

We can display the output of the function, which is the gating strategy of each cell type, using the command:

```
print(gate)

```
To keep the guided tutorial simple, we have chosen to apply the *GateMeClass_annotate* annotation function to the same matrix `m`:

```
testing_set <- m

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
marker_table        : gate (from the training module)
reject_option       : default = F
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Varing) or "E" (Equal), default = "V"
RSS                 : RSS (Ranked-Set-Sampling) observations, used to have a better resolution: default = "T"
k                   : k-NN (k-Nearest-Neighbour), used to refine the uncertain labels to the most similar already annotated
sampling            : default = "0.1"
verbose             : default = T
seed                : default = "1"
```


### 3: Training and classification in one step

```
res <- GateMeClass_annotate(m,
                            marker_table = NULL,
                            train_parameters = list(reference = m, labels = lab),
                            GMM_parameterization = "V",
                            sampling = 0.1,
                            verbose = T,
                            seed = 1)
```


Parameters description:
```
m                   : Matrix of expression
gate_table          : gate (from the training module)
train_parameters    : 
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Variance) or "E" (Equal), default = "V"
sampling            : default = "0.1"
verbose             : default = "T"
seed                : default = "1"
```


