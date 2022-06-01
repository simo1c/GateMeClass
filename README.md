# scGateMe
scGateMe is an R package for the classification of Flow Cytometry data

Example:

```
source("scGateMe.R")

gates <- data.frame(read_excel("gates.xlsx"))
colnames(gates) <- gsub("[.]", "-", colnames(gates))
gates

gate <- scGateMe_train(m, sce2$labels, gmm_criteria = "BIC", sampling = 0.1)
gate

res <- scGateMe(m,
                train = F,
                gates = gate, 
                gmm_criteria = "BIC",
                refine = T,
                sampling = 0.1,
                k = NULL,
                verbose = T,
                narrow_gate_table = T, 
                seed = 1,
                marker_seq_eval = F,
                combine_seq_eval_res = F)
  
```
