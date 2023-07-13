# GateMeClass - Manual

## Description
<p align="justify">
GateMeClass, short for Gate-Mining-and-Classification is a tool implemented in the R programming language that automates the process of annotation of a cytometry dataset using marker tables. 
GateMeClass consists of two main modules: the training module and the annotation module. These modules can be used independently or in a sequence depending on the user's requirements. 
GateMeClass has the capability to learn the marker table from an externally annotated dataset (<i>GateMeClass training module</i>) and apply the corresponding cell labels to a different dataset of interest (<i>GateMeClass annotation module</i>). In addition, GateMeClass provides the option to manually define a marker table without relying on a already annotated dataset. This marker table can then be directly applied to the dataset of interest for cell classification. In the following sections we will show how to install and use GateMeClass. 

    
For the technical details of GateMeClass refer to the pubblication : <i>link to pubblication</i>
</p> 

## Preparation

### Step 1. Installation of the required R libraries.

This is a list of all the R packages required to run GateMeClass.

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


### Step 2. Load the libraries and GateMeClass code.

In order to run GateMeClass, you have to import the necessary R packages and the file called *GateMeClass.R* containing the complete architecture of GateMeClass:

```
library("stringi")
library("data.table")
library("moments")
library("mclust")
library("stringr")
library("scales")
library("recipes")
library("dplyr")
library("RcppAlgos")
library("caret")
library("BiocManager")
library("batchelor")

source("GateMeClass.R")
```

## Execution of GateMeClass

### Section 1. Annotation of cytometry data using a manually defined marker table
<p align="justify">
As an example we will use GateMeClass to annotate the cells of the Levine32 dataset (<i>Levine et al., 2015</i>) using the manually defined marker table used in ACDC (<i>Lee et al., 2017</i>). At first, we download the Levine32 dataset using the R package HDCytoData (<i>Weber and Soneson, 2019</i>):
</p>

```
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install("HDCytoData")

library(HDCytoData)

d_SE <- Levine_32dim_SE()
d_sub <- assay(d_SE[, colData(d_SE)$marker_class == "type"])
population <- rowData(d_SE)$population_id
cofactor <- 5
d_sub <- asinh(d_sub / cofactor)

exp_matrix <- t(d_sub[population != "unassigned", ])
population <- population[population != "unassigned"]
```

Then, we read the excel file with the marker table:

```
library(readxl) # install.packages("readxl")

gate <- as.data.frame(read_excel("Levine32.xlsx"))
colnames(gate)[which(colnames(gate) == "HLA-DR")] <- "HLA_DR"      # To avoid naming problems with marker names that contain keywords of GateMeClass
gate[is.na(gate)] <- "*"                                           # required for markers not set
```

Next, we can execute the *GateMeClass_annotate* function for annotating the dataset:
```
res <- GateMeClass_annotate(exp_matrix,
                            marker_table = gate,
                            reject_option = F,
                            GMM_parameterization = "V",
                            RSS = T,
                            k = 20,				
                            sampling = 0.1,
                            verbose = T,
                            seed = 1)
```
<p align="justify">
In this example we executed GateMeClass using GMM  with varying variance (GMM_parameterization = "V"), no reject_option (reject_option = F) in order to do not care of cells potentially not defined in the marker table, sampling 10% (sampling = 0.1) of cells, ranked set sampling (RSS = T) and k parameter set to 20 (k = 20) for k-NN and MNN algorithms for label refining. For a comphensive list of GateMeClass_annotate parameters refer to the following:
</p>

List of parameters
```
exp_matrix          : Expression matrix, class = matrix, (mandatory)
marker_table        : Manually curated marker table, class = data.frame, 
reject_option       : This parameter tries to detect cell types not defined in the marker table using MNN algorithm, class = logical, default = T
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), class = character, default = "E"
RSS                 : RSS (Ranked Set Sampling). It is particularly advised in combination with GMM_parameterization = "V" to have a better resolution of the marker distribution, class = logical, default = T
k                   : k parameter of k-NN (k-Nearest-Neighbour) used to refine uncertain labels to the most similar already annotated, class = numeric, default = 20
sampling            : Perform a sampling of the cells annotating the rest with k-NN, class = numeric, default = 0.2
verbose             : Show output information, class = logical, default = T
train_parameters    : A list of parameters to pass to the *GateMeClass_train* function, class = list, default = NULL
seed                : class = numeric, default = 1
```

The output of *GateMeClass_train* is a list with the following elements: 

1) labels            : The labels returned by GateMeClass, class = character
2) marker_table      : The marker table used by GateMeClass, class = data.frame
3) cell_signatures   : The marker signatures of each cell attributed by GateMeClass, class = data.frame

The output can be showed using the following:

```
table(res$labels)
print(res$marker_table)      
print(res$cell_signatures)
```

### Section 2. Annotation of cytometry data extracting the marker table from an annotated reference dataset

<p align="justify">
In this section, we will use GateMeClass to extract a marker table from an already annotated reference cytometry dataset (training set) to annotate our current dataset (control set). For a matter of semplicity, we will use the same dataset (Levine32) as training and control set. The <u>GateMeClass_train</u> function takes in input two main parameters, an expression matrix of the reference dataset and the corresponding labels and returns in output a marker table:
</p>
    
```
training_set <- exp_matrix        # Levine32 expression matrix
training_set_lab <- population    # Labels of Levine32
```

The next step involves the use the <u>GateMeClass_train</u> training function to obtain the `gate` variable, which encompasses the pseudo gating strategy employed for the training set:

```
new_gate <- GateMeClass_train(training_set,
                          training_set_lab,
                          RSS = T,
                          GMM_parameterization = "V",
                          verbose = T, 
                          seed = 1)
```
Some of the parameters are shared between <u>GateMeClass_train</u> and <u>GateMeClass_annotate</u> modules:

List of parameters

```
reference             : Expression matrix of the reference annotated dataset, class = matrix, (mandatory)
labels                : Labels of the reference dataset, class = character, (mandatory)
RSS                   : RSS (Ranked Set Sampling). It is particularly advised in combination with GMM_parameterization = "V" to have a better resolution of the marker distribution, class = logical, default = T
GMM_parameterization  : GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), class = character, default = "E"
verbose               : Show output information, class = logical, default = T
seed                  : class = numeric, default = 1
```

We can display the output of the function, which is the gating strategy of each cell type, using the command:

```
print(gate)

```
Next, <u>GateMeClass_annotate</u> can be executed to the same dataset with the following:

```
res <- GateMeClass_annotate(exp_matrix,
                            marker_table = new_gate,
                            reject_option = F,
                            GMM_parameterization = "V",
                            RSS = T,
                            k = 20,				
                            sampling = 0.1,
                            verbose = T,
                            seed = 1)
```

The output can be showed using the following:

```
table(res$labels)
print(res$marker_table)
print(res$cell_signatures)
```


### Section 3. Training and classification in one step

<p align="justify">
Training and classification can be performed in one step using the <u>GateMeClass_annotate</u> function using the parameter </p> `train_parameters`.

```
res <- GateMeClass_annotate(exp_matrix,
                            marker_table = NULL,
                            train_parameters = list(reference = exp_matrix, labels = training_set_lab),
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
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), default = "V"
sampling            : default = "0.1"
verbose             : default = "T"
seed                : default = "1"
```


