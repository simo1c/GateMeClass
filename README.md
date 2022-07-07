# GateMeClass: Gate Mining and Classification of flow and mass cytometry data

Example:

```
source("scGateMe.R")

gates <- data.frame(read_excel("gates.xlsx"))
colnames(gates) <- gsub("[.]", "-", colnames(gates))
gates

# Extracting the gate table
gate <- GateMeClass_train(m,                         # reference dataset
                       sce2$labels,                  # labels of the dataset
                       sampling = "none",            # Type of sampling ("all" or "class"), all = no sampling, class = use SMOTE sampling
                       imp_feature_thr = "GMM",      # Used only for dataset with two labels, criteria for marker selection by importance
                       gmm_parameterization = "V",   # Parameterization of GMM
                       sampling_feature_pre = 1000,  # Dataset size for Boruta
                       sampling_imp_vars = 1000,     # Dataset size for regression tree
                       thr_perc = -1,                # minimum percentage of + and - cells for a marker to be chosed 
                       seed = 1)
                       
gate

# Classification with the extracted gate table
res <- GateMeClass_annotate(m,                                # dataset to classify
                         gate_table = gate,                   # gate table
                         GMM_parameterization = "E",          # Parameterization of GMM
                         train_parameters = list(             # list of parameters for the training function
                           reference = NULL
                         ),
                         refine = T,                          # KNN algorithm to refine labels after classification
                         sampling = 0.25,                     # Sampling to classify a subset of the dataset
                         k = NULL,                            # k parameter for KNN
                         verbose = T,                         # Show or not output program
                         narrow_gate_table = T,               # format of gate table 
                         seed = 1)

# ---------------------------------------------------------------------------------------------------------------

# Training and classification in a step:
# Artificial example with the same dataset as reference and control set 
colnames(m) <- sce2$labels # column names of the reference are the labels
res <- GateMeClass_annotate(m,
                         gate_table = NULL,
                         train_parameters = list(
                           reference = m,
                           labels = colnames(m)
                         ),
                         gmm_parameterization = "V",
                         refine = T,
                         sampling = 0.25,
                         k = NULL,
                         verbose = T,
                         narrow_gate_table = T, 
                         seed = 1)


  
  
```
