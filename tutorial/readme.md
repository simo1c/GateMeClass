# GateMeClass - Tutorial

## Description
 <img width="200" height="220" src="../logo2.jpg" alt = "Logo GateMeClass" align = "right">

<p align="justify"> 
 GateMeClass (Gate Mining and Classification) is a tool implemented in the R programming language that automates the process of annotation of a cytometry dataset using marker tables. 
GateMeClass consists of two main modules: the training module and the annotation module. These modules can be used independently or in a sequence depending on the user's requirements. 
GateMeClass has the capability to learn the marker table from an externally annotated dataset (<i>GateMeClass training module</i>) and apply the corresponding cell labels to a different dataset of interest (<i>GateMeClass annotation module</i>). However, you can use GateMeClass with a manually defined marker table without relying on a already annotated dataset. This marker table can then be directly applied to the dataset of interest for cell classification. In the following sections we will show you how to install and use GateMeClass. 
 
 <!--For the technical details of GateMeClass refer to the pubblication : <i>link to pubblication</i> -->
</p> 

## Preparation

### Installation and loading of GateMeClass package

To install and load GateMeClass please execute the following lines of code:
```
install.packages("devtools")
devtools::install_github("https://github.com/simo1c/GateMeClass")
library(GateMeClass)
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
library("tidyr")

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
res <- GateMeClass_annotate(exp_matrix = exp_matrix,
                            marker_table = gate,
                            GMM_parameterization = "V",
                            reject_option = F,
                            sampling = 0.1,
                            RSS = T,
                            k = 20,				
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
GMM_parameterization: GMM (Gaussian-Mixture-Model) parameter: "V" (Variable) or "E" (Equal), class = character, default = "E"
RSS                 : RSS (Ranked Set Sampling). It is particularly advised in combination with GMM_parameterization = "V" to have a better resolution of the marker distribution, class = logical, default = T
k                   : k parameter of k-NN (k-Nearest-Neighbour) used to refine uncertain labels to the most similar already annotated, class = numeric, default = 20
sampling            : Perform a sampling of the cells annotating the rest with k-NN, class = numeric, default = 0.2
narrow_marker_table : format of marker table. TRUE for using a marker table with two columns, "Cell" with the name of cell type and "Gate" with the gating strategy (e.g., CD3+CD4+). FALSE for using a column for each marker (as in the above example), class = logical, default = T. 
verbose             : Show output information, class = logical, default = T
train_parameters    : A list of parameters to pass to the GateMeClass_train function, class = list, default = NULL
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
To select the parameter GMM_parameterization (V or E), it can be useful to explore the distribution of each marker to investigate the variance of the expeceted Gaussian components: 
</p>

```
m <- as.data.frame(t(exp_matrix))

data_long <- m %>%
  pivot_longer(colnames(m)) %>% 
  as.data.frame()
head(data_long) 

ggp2 <- ggplot(data_long, aes(x = value)) +
  geom_density() + 
  facet_wrap(~ name, scales = "free") +
  theme_classic()
ggp2

```


<p align="justify">
If we want to detect cell types not specified in the marker table, we can set the parameter reject_option = T. As an example, we define a "narrow" (narrow_marker_table = T) marker table and specify only the CD3 marker of T cells:
</p>

```
gate <- data.frame(Cell = c("T cells"), Gate = ("CD3+"))

res <- GateMeClass_annotate(exp_matrix = exp_matrix,
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
GateMeClass supports the use of the OR (|) and XOR (^) logical operators in the marker table. For example, if we want to identify cells that expresses CD11b, CD11c or both we can create the following marker table:
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
new_gate <- GateMeClass_train(reference = training_set,
                          labels = training_set_lab,
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
GMM_parameterization  : GMM (Gaussian-Mixture-Model) parameter: "V" (Variable) or "E" (Equal), class = character, default = "E"
verbose               : Show output information, class = logical, default = T
seed                  : class = numeric, default = 1
```

We can display the output of the function, which is the gating strategy of each cell type, using the command:

```
print(gate)
```
Next, *GateMeClass_annotate* can be executed to the same dataset with the following:

```
res <- GateMeClass_annotate(exp_matrix = exp_matrix,
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
res <- GateMeClass_annotate(exp_matrix = exp_matrix,
                            marker_table = NULL,
                            train_parameters = list(reference = exp_matrix,
                                                    labels = training_set_lab),
                            GMM_parameterization = "V",
                            sampling = 0.1,
                            verbose = T,
                            seed = 1)
```

### Method 4. Annotation of clusters using GateMeClass

<p align="justify">
GateMeClass can be used also to annotate clusters obtained using other techniques. In this case <i>GateMeClass_train</i> can be executed setting the parameter 'labels' with the corresponding cluster identities obtained using an external clustering algorithm. In this way the marker table obtained will contain the marker signatures of each cluster that can be used for better interpretability and for subsequently annotation using <i>GateMeClass_annotate</i> function:
</p>

```
gate_clusters <- GateMeClass_train(exp_matrix = exp_matrix,
                                   labels = cluster_labels,            # This labels are obtained executing clustering with an external tool (e.g., FlowSOM)
                                   RSS = T,
                                   GMM_parameterization = "V",
                                   verbose = T, 
                                   seed = 1)
```
### Attached packages:

```
batchelor_1.14.1
SingleCellExperiment_1.20.0
BiocManager_1.30.20
caret_6.0-93
lattice_0.20-45
ggplot2_3.4.1
RcppAlgos_2.7.2
recipes_1.0.5
dplyr_1.1.0
scales_1.2.1               
stringr_1.5.0
mclust_6.0.0
moments_0.14.1
data.table_1.14.8
stringi_1.7.12             
readxl_1.4.2
HDCytoData_1.18.0
flowCore_2.10.0
SummarizedExperiment_1.28.0
Biobase_2.58.0
GenomicRanges_1.50.2
GenomeInfoDb_1.34.9
IRanges_2.32.0
S4Vectors_0.36.2
MatrixGenerics_1.10.0
matrixStats_0.63.0
ExperimentHub_2.6.0
AnnotationHub_3.6.0
BiocFileCache_2.6.1
dbplyr_2.3.1               
BiocGenerics_0.44.0 
```


