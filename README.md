# GateMeClass: Gate Mining and Classification of cytometry data

Description:




Example:

```
source("GateMeClass.R")

   # Training
   gate <- GateMeClass_train(training_set,
                             training_set_lab,
                             RSS = T,
                             GMM_parameterization = "V",
                             verbose = T, 
                             seed = 1)


    # Annotation
    res <- GateMeClass_annotate(testing_set,
                                marker_table = gate,
                                reject_option = F,
                                GMM_parameterization = "V",
                                RSS = T,
                                k = 20,				
                                sampling = 0.1,
                                verbose = T,
                                seed = 1)


# Training and classification in one step:
# Artificial example with the same dataset as reference and control set 
colnames(m) <- sce2$labels # column names of the reference are the labels
res <- GateMeClass_annotate(m,
                         gate_table = NULL,
                         train_parameters = list(
                           reference = m,
                           labels = colnames(m)
                         ),
                         GMM_parameterization = "V",
                         sampling = 0.25,
                         verbose = T,
                         seed = 1)


  
  
```
