# scGateMe
scGateMe is an R package for the classification of flow and mass cytometry data

Example:

```
source("scGateMe.R")

gates <- data.frame(read_excel("gates.xlsx"))
colnames(gates) <- gsub("[.]", "-", colnames(gates))
gates

gate <- scGateMe_train(m,
                       sce2$labels,
                       sampling = "all",
                       Boruta = F,
                       imp_feature_thr = "GMM",
                       gmm_parameterization = "V",
                       sampling_feature_pre = 1000,
                       sampling_imp_vars = 1000,
                       thr_perc = -1,
                       seed = 1)
                       
gate

res <- scGateMe(m,
                gates = gate,
                # reference = m, 
                sampling_feature_pre = 1000,
                sampling_imp_vars = 1000,
                imp_feature_thr = "GMM",
                sampling = "all",
                thr_perc = -1,
                gmm_parameterization = "V",
                refine = T,
                sampling = 0.25,
                k = NULL,
                verbose = T,
                narrow_gate_table = T, 
                seed = 1)
  
```
