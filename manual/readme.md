# GateMeClass - Manual

 <img width="230" height="250" src="../logo2.jpg" alt = "Logo GateMeClass" align = "right">

## Description
<p align="justify">
GateMeClass (Gate Mining and Classification) is a tool implemented in the R programming language that automates the process of annotation of a cytometry dataset using marker tables. 
GateMeClass consists of two main modules: the training module and the annotation module. These modules can be used independently or in a sequence depending on the user's requirements. 
GateMeClass has the capability to learn the marker table from an externally annotated dataset (<i>GateMeClass training module</i>) and apply the corresponding cell labels to a different dataset of interest (<i>GateMeClass annotation module</i>). However, you can use GateMeClass with a manually defined a marker table without relying on a already annotated dataset. This marker table can then be directly applied to the dataset of interest for cell classification. In the following sections we will show how to install and use GateMeClass. For the technical details of GateMeClass refer to the pubblication : <i>link to pubblication</i>
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
source("GateMeClass.R")    # GateMeClass code
```

## Execution of GateMeClass

### Method 1. Annotation of cytometry data using a manually defined marker table
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
                            narrow_marker_table = F,
                            seed = 1)
```
<p align="justify">
In this example we executed GateMeClass using GMM  with varying variance (GMM_parameterization = "V"), no reject_option (reject_option = F) in order to do not care of cells potentially not defined in the marker table, sampling 10% (sampling = 0.1) of cells, ranked set sampling (RSS = T) and k parameter set to 20 (k = 20) for k-NN and mutual nearest neighbor (MNN) algorithms for label refining. For a comphensive list of <i>GateMeClass_annotate</i> parameters refer to the following:
</p>

List of parameters:
```
exp_matrix          : Expression matrix, class = matrix, (mandatory)
marker_table        : Manually curated marker table, class = data.frame, 
reject_option       : This parameter tries to detect cell types not defined in the marker table using MNN algorithm, class = logical, default = T
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Varying) or "E" (Equal), class = character, default = "E"
RSS                 : RSS (Ranked Set Sampling). It is particularly advised in combination with GMM_parameterization = "V" to have a better resolution of the marker distribution, class = logical, default = T
k                   : k parameter of k-NN (k-Nearest-Neighbour) used to refine uncertain labels to the most similar already annotated, class = numeric, default = 20
sampling            : Perform a sampling of the cells annotating the rest with k-NN, class = numeric, default = 0.2
narrow_marker_table : format of marker table. TRUE for using a marker table with two columns, "Cell" with the name of cell type and "Gate" with the gating strategy (e.g., CD3+CD4+). FALSE for using a column for each marker (as in the above example), class = logical, default = T. 
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
<p align="justify">
If we want to detect cell types not specified in the marker table, we can set the parameter reject_option = T. As an example, we define a "narrow" (narrow_marker_table = T) marker table and specify only the CD3 marker of T cells:
</p>

```
gate <- data.frame(Cell = c("T cells"), Gate = ("CD3+"))

res <- GateMeClass_annotate(exp_matrix,
                            marker_table = gate,
                            reject_option = T,          # we want to discriminate cell types not in marker table 
                            GMM_parameterization = "V",
                            RSS = T,
                            k = 20,				
                            sampling = 0.1,
                            verbose = T,
                            narrow_marker_table = T,    # we specified a narrow marker table
                            seed = 1)
table(res$labels)
```
<p align="justify">
GateMeClass supports the use of the or (|) and XOR (^) logical operators in the marker table. For example, if we want to identify cells that expresses CD11b, CD11c or both we can create the following marker table:
</p>

```
gate <- data.frame(Cell = c("Cells"), Gate = ("CD11b+|CD11c+|"))
```

To identify cells that expresses CD11b or CD11c but not both we can create the following marker table:

```
gate <- data.frame(Cell = c("Cells"), Gate = ("CD11b+^CD11c+^"))
```


### Method 2. Annotation of cytometry data extracting the marker table from an annotated reference dataset

<p align="justify">
In this section, we will use GateMeClass to extract a marker table from an already annotated reference cytometry dataset (training set) to annotate our current dataset (control set). For a matter of semplicity, we will use the same dataset (Levine32) as training and control set. The <i>GateMeClass_train</i> function takes in input two main parameters, an expression matrix of the reference dataset and the corresponding labels and returns in output a marker table:
</p>
    
```
training_set <- exp_matrix        # Levine32 expression matrix
training_set_lab <- population    # Labels of Levine32
```

The next step involves the use the *GateMeClass_train* training function to obtain the `gate` variable, which encompasses the pseudo gating strategy employed for the training set:

```
new_gate <- GateMeClass_train(training_set,
                          training_set_lab,
                          RSS = T,
                          GMM_parameterization = "V",
                          verbose = T, 
                          seed = 1)
```
Some of the parameters are shared between *GateMeClass_train* and *GateMeClass_annotate* modules:

List of parameters:

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
Next, *GateMeClass_annotate* can be executed to the same dataset with the following:

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
<p align="justify">
Suppose we want to extract a marker table from a dataset in which we know that monocytic subsets were gated considering the medium expression of CD11b. In that case, GateMeClass needs the presence of three subsets with the following labels: monocytes_CD11bhi, monocytes_CD11bmid and monocytes_CD11b-. The syntax is important and must be "cell_type_name"_"marker_name[hi|mid|-]". If we have these labels in the reference dataset we can take into account medium expression of that marker in <i>GateMeClass_train</i>.
</p>

### Method 3. Training and classification in one step


Training and classification can be performed in one step using the *GateMeClass_annotate* function using the parameter 'train_parameters'. In this case, 
It is sufficient to specify in the 'train_parameters' the expression matrix of the reference matrix and the corresponding labels:



```
res <- GateMeClass_annotate(exp_matrix,
                            marker_table = NULL,
                            train_parameters = list(reference = exp_matrix,
                                                    labels = training_set_lab),
                            GMM_parameterization = "V",
                            sampling = 0.1,
                            verbose = T,
                            seed = 1)
```

### Method 4. Annotation of clusters using GateMeClass


GateMeClass can be used also to annotate clusters obtained using other techniques. In this case *GateMeClass_train* can be executed setting the parameter 'labels' with the corresponding cluster identity obtained using an external clustering algorithm. In this way the marker table obtained will contain the marker signatures of each cluster that can be used for better interpretability and for subsequently annotation using *GateMeClass_annotate* function:

```
gate_clusters <- GateMeClass_train(exp_matrix,
                                   labels = cluster_labels,            # This labels are obtained executing clustering with an external tool (e.g., FlowSOM)
                                   RSS = T,
                                   GMM_parameterization = "V",
                                   verbose = T, 
                                   seed = 1)
```


