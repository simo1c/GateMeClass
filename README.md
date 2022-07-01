# scGateMe
scGateMe is an R package for the classification of flow and mass cytometry data

Example:

```
source("scGateMe.R")

gates <- data.frame(read_excel("gates.xlsx"))
colnames(gates) <- gsub("[.]", "-", colnames(gates))
gates

# Extracting the gate table
gate <- scGateMe_train(m,                            # reference dataset
                       sce2$labels,                  # labels of the dataset
                       sampling = "all",             # Type of sampling ("all" or "class"), all = no sampling, class = use SMOTE sampling
                       Boruta = F,                   # Pre-processing to remove non important markers
                       imp_feature_thr = "GMM",      # Used only for dataset with two labels, criteria for marker selection by importance
                       gmm_parameterization = "V",   # Parameterization of GMM
                       sampling_feature_pre = 1000,  # Dataset size for Boruta
                       sampling_imp_vars = 1000,     # Dataset size for regression tree
                       thr_perc = -1,                # minimum percentage of + and - cells for a marker to be chosed 
                       seed = 1)
                       
gate

# Classification with the extracted gate table
res <- scGateMe(m,                                   # dataset to classify
                gates = gate,                        # gate table
                reference = NULL,                    # reference dataset to use for creating the gate table
                sampling_feature_pre = 1000,
                sampling_imp_vars = 1000,
                imp_feature_thr = "GMM",
                sampling = "all",
                thr_perc = -1,
                gmm_parameterization = "V",
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
res <- scGateMe(m,
                gates = gate,
                reference = m, 
                sampling_feature_pre = 1000,
                sampling_imp_vars = 1000,
                imp_feature_thr = "GMM",
                sampling_train = "all",
                thr_perc = -1,
                gmm_parameterization = "V",
                refine = T,
                sampling = 0.25,
                k = NULL,
                verbose = T,
                narrow_gate_table = T, 
                seed = 1)
  
  
```
