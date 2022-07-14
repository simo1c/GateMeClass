# options(warn=1)

## This function generates the possible value combinations of the markers of a cell type 
generate_set_values <- function(v, cell){
  
  sets_all <- list()
  sets_to_filter <- list()
  sets_to_filter2 <- list()
  
  special1 <- grep("[\\|]", v)
  special2 <- grep("[\\^]", v)
  
  if(length(special1) == 1 | length(special2) == 1){
    stop("Special characters (^, |) must be set in more than 2 markers!")
  }
  
  for(m in v){
    switch(m,
           `+` = {
             set <- list(c("+"))
             set_not <- list()
             set_not2 <- list()
           },
           `-` = {
             set <- list(c("-"))
             set_not <- list()
             set_not2 <- list()
           },
           `+|`={
             set <- list(c("+", "-"))
             set_not <- list(c("-"))
             set_not2 <- list()
           },
           `-|`={
             set <- list(c("+", "-"))
             set_not <- list(c("+"))
             set_not2 <- list()
           },
           `+^`={
             set <- list(c("+", "-"))
             set_not <- list()
             set_not2 <- list(c("+"), c("-"))
             # set_not2 <- list(c("-"))
             # set_not <- list()
           },
           `-^`={
             set <- list(c("+", "-"))
             set_not <- list()
             set_not2 <- list(c("-"), c("+"))
             # set_not2 <- list(c("+"))
           },
           `*`={
             set <- list(c("*"))
             # set <- list(c("*", "+", "-"))
             set_not <- list()
             set_not2 <- list()
           }
    )
    sets_all <- c(sets_all, set)
    sets_to_filter <- c(sets_to_filter, set_not)
    sets_to_filter2 <- c(sets_to_filter2, set_not2)
  }
  
  names(sets_all) <- names(v)
  
  sets_all <- expand.grid(sets_all)
  
  if(length(sets_to_filter) > 0){
    special_to_filter <- grep("[\\|]", v)
    sets_to_filter <- data.frame(matrix(as.character(sets_to_filter), ncol = length(special_to_filter)))
    colnames(sets_to_filter) <- names(v)[special_to_filter]
  }
  
  if(length(sets_to_filter2) > 0){
    special_to_filter2 <- grep("[\\^]", v)
    sets_to_filter2 <- data.frame(matrix(as.character(sets_to_filter2), ncol = length(special_to_filter2)))
    colnames(sets_to_filter2) <- names(v)[special_to_filter2]
  }
  
  if(length(sets_to_filter) > 0){
    sets_all_filtered <- anti_join(sets_all, sets_to_filter, by = names(sets_to_filter))
  }else{
    sets_all_filtered <- sets_all
  }
  
  if(length(sets_to_filter2) > 0){
    sets_all_filtered <- anti_join(sets_all_filtered, sets_to_filter2, by = names(sets_to_filter2))
  }else{
    sets_all_filtered <- sets_all_filtered
  }
  
  sets_all_filtered$Cell <- cell
  sets_all_filtered_temp <- sets_all_filtered[, -which(colnames(sets_all_filtered) == "Cell"), drop = F]
  markers <- colnames(sets_all_filtered_temp)
  gate <- apply(sets_all_filtered_temp, 1, stri_c, collapse = "")
  sets_all_filtered$Gate <- gate
  sets_all_filtered <- sets_all_filtered[, c("Cell", "Gate")]
  
  return(sets_all_filtered)
}

## Read the gate table and generate the possbile marker signature of ceach cell type
parse_gate_table <- function(gate_table, narrow_gate_table, extended_gate_table){
  
  # gate_table <- gate
  # narrow_gate_table = T
  # extended_gate_table = T
  
  if(any(duplicated(gate_table$Cell))){
    stop("The gate table must contains uniquely defined cell types!")
  }
  
  if(narrow_gate_table){
    temp <- sapply(gate_table$Gate, function(x){
      str <- strsplit(x, "[-]|-[\\|\\^]|[+]|\\+[\\|\\^]|[*]")
      str[[1]] <- stri_remove_empty(str[[1]])
      return(str)
    })
    
    names(temp) <- gate_table$Cell
    markers <- unique(unlist(temp))
    df_gates <- data.frame(matrix(nrow = nrow(gate_table), ncol = length(markers)))
    rownames(df_gates) <- gate_table$Cell
    colnames(df_gates) <- markers
    
    gate_table_exploded <- strsplit(perl=T, gate_table$Gate, '(?![\\+\\||\\-\\||\\+\\^|\\-\\^])')
    signs <- lapply(gate_table_exploded, function(x){ x[x %in% c("-^", "+^", "-|", "+|", "+", "-", "*")]})
    names(signs) <- names(temp)
    
    l_temp <- 1:length(temp)
    
    signs <- lapply(1:length(signs), function(i){
      el <- signs[[i]]
      marker <- temp[[i]]
      tt <- rep(NA, length(markers))
      names(tt) <- markers
      tt[temp[[i]]] <- el
      tt[is.na(tt)] <- "*"
      return(tt)
      
    })
    
    df_gates <- data.frame(matrix(unlist(signs), nrow=length(signs), ncol = length(markers), byrow=TRUE))
    colnames(df_gates) <- markers
    rownames(df_gates) <- names(temp)
    
    df_gates[is.na(df_gates)] <- "*"
    df_gates <- cbind(Cell = rownames(df_gates), df_gates)
    gate_table <- df_gates
  }
  
  celltypes <- gate_table$Cell
  
  if(length(grep("[\\*|\\||\\^]", celltypes)) > 0){
    stop("Cell names cannot contains special characters (e.g., *, ^)!")
  }
    
  if(extended_gate_table){
    df_list <- apply(gate_table, 1, function(v){
      
      # v <- as.character(gate_table[2,-1])
      # names(v) <- colnames(gate_table[2,-1])
      # cell <- gate_table[2,1]
      # 
      # v <- as.character(v)
      # cell <- v[1]
      # v <- v[-1] 
      # names(v) <- colnames(gate_table)[-1]
      # print(class(v))
      # print(cell)
      
      to_add <- generate_set_values(v[-1], v[1])
      return(to_add)
    })
    
    extended_gate_table <- as.data.frame(rbindlist(df_list))
    to_delete <- stri_c(colnames(gate_table[, -1]), "*", collapse = "")
    extended_gate_table <- extended_gate_table[extended_gate_table$Gate != to_delete, ]
    
    return(list(gate_table = gate_table, extended_gate_table = extended_gate_table))
  }

  return(list(gate_table = gate_table))
}

set_marker_expression_GMM <- function(X, GMM_parameterization){
  
  # X <- m["FS-A", ]
  # indexes = first$indexes
  # GMM_parameterization = "E"
  
  if(length(X) >= 100){

    sample <- X
    
    if(skewness(X) < 0){
      sk <- "min"
    }else{
      sk <- "max"
    }

    ############ Ranked Set Sampling (RSS) #########################
    sample <- X
    cycles <- floor(length(sample) / 4)
    n_unit <- 2
    n <- length(sample)
    samples <- cycles * n_unit
    indexes <- matrix(sample(1:n, samples * n_unit), nrow = n_unit)
    sel_samples <- matrix(sample[indexes], nrow = n_unit)
    ################################################################

    if(skewness(X) < 0){
      test <- apply(sel_samples, 2, min)
    }else{
      test <- apply(sel_samples, 2, max)
    }
  }else if(length(X) >= 2){
    test <- X
  }else{
    return(c("*"))
  }
  
  # test <- m["CD8", ]
  # GMM_parameterization <- "V"
  # X <- test
  
  # top <- mclustBIC(test, G = 2)
  # model_temp <- unlist(str_split(names(summary(top)[1]), ","))
  # type_model <- model_temp[1]

  # test <- m["CD10", ]
  
  # GMM_parameterization <- "V"
  
  cl <- Mclust(test, G = 2, verbose = F, modelNames = GMM_parameterization)
  
  max <- which.max(cl$parameters$mean)
  min <- which.min(cl$parameters$mean)
  
  min_cl <- min(test[cl$classification == min])
  max_cl <- max(test[cl$classification == max])
  
  pred <- predict.Mclust(cl, X)
  temp <- pred$classification
  
  temp <- ifelse(temp == max, "+", "-")
  
  # temp[X < min_cl] <- "-"
  # temp[X > max_cl] <- "+"
  
  # if(any(temp[X <= min_cl] == "+") | any(temp[X >= max_cl] == "-")){
  # 
  #   # print("************")
  #   # plot(cl, what = "classification")
  # 
  #   GMM_parameterization <- "E"
  #   cl <- Mclust(test, G = 2, verbose = F, modelNames = GMM_parameterization)
  #   max <- which.max(cl$parameters$mean)
  #   pred <- predict.Mclust(cl, X)
  #   temp <- pred$classification
  #   temp <- ifelse(temp == max, "+", "-")
  # }
  
  # temp[temp == names(means)[1]] <- "-"
  # temp[temp == names(means)[2]] <- "+"
  # 
  # test <- m["CD3", ]
  # X <- test
  # cl = normalmixEM(test, k = 2, verb = F, fast = T)
  # max <- which.max(cl$mu)
  # cl <- apply(cl$posterior, 1, which.max)
  # temp <- ifelse(cl == max, "+", "-")
  
  return(temp)
}

## This function set the marker signature of each cell
set_marker_expression <- function(exp_matrix, 
                                  markers,
                                  expr_markers, 
                                  verbose, 
                                  GMM_parameterization){


  # exp_matrix <- exp_matrix_2
  # markers <- colnames(new_gates$gate_table)[-1]
  # # # expr_markers
  # gates <- new_gates$gate_table
  # # verbose = F
  # marker_seq_eval = F
  # mm <- 2
  # rr <- 0.05
  # gmm_criteria <- "ICL"

  ## GMM probabilities initialization
  # prob <- rep(1, ncol(exp_matrix))

  queue <- list(list(indexes = 1:ncol(exp_matrix), markers = markers))

  while(length(queue) > 0){
    bimodal_markers <- c()
    not_bimodal_markers <- c()

    ## Pop operation
    first <- queue[[1]]
    queue <- queue[-1]

    for(m in first$markers){
      X <- exp_matrix[m, first$indexes]
      marker_expr <- set_marker_expression_GMM(X, GMM_parameterization)
      if(length(table(marker_expr)) > 1){
        bimodal_markers <- c(bimodal_markers, m)

        if(verbose){
          message(stri_c(" - ", stri_c(bimodal_markers, collapse = " ", sep = ""), collapse = " ", sep =))
        }

        expr_markers[m, first$indexes] <- marker_expr
        
      }else{
        not_bimodal_markers <- c(not_bimodal_markers, m)
      }
    }

    if(length(not_bimodal_markers) > 0 & length(bimodal_markers) == 0){
      expr_markers[not_bimodal_markers, first$indexes] <- "*"
      next
    }else if(length(not_bimodal_markers) == 0 & length(bimodal_markers) > 0){
      next
    }else{
      r_temp <- expr_markers[bimodal_markers, first$indexes]
      comb_markers <-  r_temp[!duplicated(as.list(r_temp))]
      comb_list <- as.list(comb_markers)
      names(comb_list) <- NULL
      comb_list2 <- sapply(comb_list, stri_c, collapse = "", sep = "")
      l <- expr_markers[bimodal_markers, first$indexes, drop = F]
      l <- sapply(l, stri_c, collapse = "", sep = "")
      
      to_add <- lapply(comb_list2, function(c){
        w <- which(sapply(l, function(c2){ 
          return(identical(c, c2))
        }))
      
        ## Push operation
        el <- list(indexes = first$indexes[w], markers = not_bimodal_markers)
        return(el)
      })
      
      queue <- c(queue, to_add)
    }
  }
  return(expr_markers)
}

## This function performs the cell classification
cell_classification <- function(marker_table, gates){
  
  # marker_table <- new_cells2
  # gates <- new_gates$extended_gate_table
  
  df_gate <- marker_table
  colnames(df_gate)[1] <- "Cell_ID"
  
  pos <- sapply(gates$Gate, function(g){
    gs <- unlist(str_split(g, pattern = ""))
    w <- which(gs != "*")
    w <- stri_c(w, collapse = "_", sep = "")
    return(w)
  })
  
  gates$Pos <- pos
  
  clean <- sapply(gates$Gate, function(g){
    gs <- unlist(str_split(g, pattern = ""))
    w <- which(gs != "*")
    clean_gs <- stri_c(gs[w], collapse = "", sep = "")
    return(clean_gs)
  })
  
  gates$Gate <- clean
  
  gates$N <- sapply(gates$Gate, function (g){
    return(sum(unlist(str_split(g, pattern = "")) == "+"))
  })
  
  gates$N2 <- sapply(gates$Gate, function (g){
    return(sum(unlist(str_split(g, pattern = "")) == "-"))
  })
  
  gates <- gates[order(gates$N2, decreasing = F), ]
  gates <- gates[order(gates$N, decreasing = F), ]
  
  labels <- rep("Unclassified", nrow(df_gate))
  multiple_cl <- rep(0, nrow(df_gate))
  
  ###  Optimized case case
  if(length(table(gates$Pos)) == 1){
    
    index <- as.numeric(unlist(str_split(gates$Pos[1], pattern = "_")))
    
    gates_split <- sapply(1:length(gates$Gate), function(j){
      gate_split <- unlist(str_split(gates$Gate[j], pattern = ""))
      x <- stri_c(gate_split[index], collapse = "", sep = "")
      return(x)
    })
    
    cell_split <- sapply(1:length(df_gate$Gate), function(j){
      gate_split <- unlist(str_split(df_gate$Gate[j], pattern = ""))
      x <- stri_c(gate_split[index], collapse = "", sep = "")
      return(x)
    })
    
    for(i in 1:length(cell_split)){
      cell_gate <- cell_split[i]
      
      test <- sapply(1:length(gates_split), function(j){
        return(gates_split[j] == cell_gate)
      })
      
      if(sum(test) > 1){
        labels[i] <- "Unclassified"
      }else if(sum(test) == 1){
        labels[i] <- gates$Cell[test]
      }
    }
  }else{
    gates_split <- lapply(gates$Pos, function(g){
      pos <- as.numeric(unlist(str_split(g, pattern = "_")))
      return(pos)
    })
    
    cell_split <- lapply(df_gate$Gate, function(g){
      gate <- unlist(str_split(g, pattern = ""))
      return(gate)
    })
    
    for(i in 1:length(df_gate$Gate)){
      
      # i <- 1
      
      c_split <- cell_split[[i]]
      
      test <- sapply(1:length(gates$Gate), function(j){
        index_split <- gates_split[[j]]
        x <- stri_c(c_split[index_split], collapse = "", sep = "")
        y <- gates$Gate[j] == x
        names(y) <- paste0(gates[j, c("N", "N2")], collapse = "_")
        return(y)
      })
      
      if(sum(test) == 1){
        labels[i] <- gates$Cell[test]
      }else if(sum(test) > 1){
        
        w <- which(test == T)
        
        tmp1 <- as.numeric(sapply(strsplit(names(w),"_"), `[`, 1))
        tmp2 <- as.numeric(sapply(strsplit(names(w),"_"), `[`, 2))
        
        w_max_tmp1 <- which(tmp1 == max(tmp1))
        
        if(length(w_max_tmp1) == 1){
          
          ind <- w[w_max_tmp1]
          labels[i] <- gates$Cell[ind]
          
        }else{
          
          w_max_tmp2 <- which(tmp2[w_max_tmp1] == max(tmp2[w_max_tmp1]))
          
          if(length(w_max_tmp2) == 1){
            ind <- w[w_max_tmp2]
            labels[i] <- gates$Cell[ind]
          }else{
            labels[i] <- "Unclassified"
          }
        }
      }
      
      # if(sum(test) > 1){
      #   labels[i] <- "Unclassified"
      # }else if(sum(test) == 1){
      #   labels[i] <- gates$Cell[test]
      # }
      
    }
  }
  
  rownames(gates) <- NULL
  marker_table$Celltype <- labels
  
  cl_res <- list(labels = labels, marker_table = gates, cell_signatures = marker_table)
  return(cl_res)
}

check_marker_names <- function(marker_names){
  to_modify <- grep("^[0-9]", marker_names)
  
  if(length(to_modify) > 0){
    marker_names[to_modify] <- stri_c("P_", marker_names[to_modify], sep = "")
  }
  
  to_modify2 <- grep("[\\+]|[\\-]|[\\|]|[\\^]|[\\*]", marker_names)
  
  if(length(to_modify2) > 0){
    marker_names[to_modify2] <- gsub("[\\+]|[\\-]|[\\|]|[\\^]|[\\*]", "_", marker_names[to_modify2])
  }
  
  modified <- c(to_modify, to_modify2)
  return(list(marker_names = marker_names, modified = modified))
}

# This is the train module of GateMeClass
GateMeClass_train <- function(reference = NULL,
                              labels = NULL, 
                              GMM_parameterization = NULL,
                              sampling = "class",
                              sampling_perc = 0.01,
                              perc.over = 100,
                              perc.under = 0,
                              sampling_k = min(5, table(labels)-1),
                              sampling_imp_vars = 1000,
                              seed = 1,
                              verbose = T){

# reference <- m
# labels <- lab
# GMM_parameterization = "V"
# sampling = "class"
# seed = 1
# rr = 0.1
# sampling_feature_pre = 1000
# sampling_imp_vars = 1000
# thr_perc = -1
# verbose = T
# sampling_feature_method ="all"
# imp_feature_thr = "all"
# sampling_perc = 0.01
# sampling_k = 5
# perc.over = 100
# perc.under = 0
  
  set.seed(seed)
  
  if(!is.null(reference)){
    old_names <- rownames(reference)
    check <- check_marker_names(rownames(reference))
    rownames(reference) <- check$marker_names
    
    if(length(check$modified) > 0){
      if(verbose){
        message("GateMeClass train - Renaming markers of reference dataset: ", 
                stri_c(old_names[check$modified], collapse = " ", sep = ""), " --> ", 
                stri_c(rownames(reference)[check$modified], collapse = " ", sep = ""))
      }
    }
  }else{
    stop("The expression matrix is mandatory!")
  }
  
  if(is.null(labels)){
    stop("The labels of the reference are mandatory!")
  }
  
  if(is.null(GMM_parameterization)){
    stop("The parameter 'GMM_parameterization' is mandatory!")
  }
  
  if(sampling == "class"){
    if(verbose){
      message("GateMeClass train - Executing oversampling to balance the training set...")
    }
    
    new_reference <- reference
    new_reference <- t(new_reference)
    new_reference <- data.frame(new_reference)
    new_reference$labels <- labels 
      
    new_labels <- labels
    t_lab <- table(labels)
    max <- max(t_lab)
    w <- names(which.max(t_lab))
    t_lab_perc <- t_lab / max
    w_sample <- names(which(t_lab_perc < sampling_perc))
    
    for(w2 in w_sample){
      temp <- as.data.frame(t(reference[, labels %in% c(w, w2)]))
      temp$labels <- labels[labels %in% c(w, w2)]
      temp$labels <- factor(temp$labels)

      newData <- SMOTE(labels ~ ., temp, perc.over = perc.over, k = sampling_k, perc.under = perc.under)
      new_reference <- rbind(new_reference, newData)
    }
    
    labels <- new_reference$labels 
    new_reference <- new_reference[, -which(colnames(new_reference) == "labels")]
    new_reference <- t(new_reference)
    reference <- new_reference
  }
  
  markers <- rownames(reference)
  celltypes <- factor(unique(labels))
  reference_2 <- reference[markers, , drop = F]
  
  if(verbose){
    message("GateMeClass train - Determining the marker signature of each cell...")
  }
  
  expr_markers <- data.frame(matrix(ncol = ncol(reference_2), nrow = nrow(reference_2)))
  rownames(expr_markers) <- markers
  
  reference_2 <- reference[markers, , drop = F]
  expr_markers <- data.frame(matrix(ncol = ncol(reference_2), nrow = nrow(reference_2)))
  rownames(expr_markers) <- markers
  
  ## Obtain the signature for the cells of the dataset to be annotated
  expr_markers <- set_marker_expression(reference_2,
                                        markers, 
                                        expr_markers, 
                                        # new_gates$gate_table, 
                                        verbose = verbose, 
                                        GMM_parameterization = GMM_parameterization)
  
  expr_markers <- expr_markers[markers, ]
  expr_markers <- as.data.frame(t(expr_markers))
  
  cells <- cbind(rownames(expr_markers), expr_markers)
  colnames(cells)[1] <- "Cell"
  
  expr_markers <- apply(cells[,- 1], 1, stri_c, collapse = "")
  new_cells2 <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
  
  res <- list(cell_signatures = new_cells2)
  signatures <- res$cell_signatures

  gate_ext <- sapply(signatures$Gate, function(v){
    split <- unlist(str_split(v, ""))
    sig <- stri_c(markers, split, collapse = "", sep = "")
    return(sig)
  })
  
  res$cell_signatures[, "Gate"] <- gate_ext
  res <- res$cell_signatures
  
  if(verbose){
    message("GateMeClass train - Pairwise comparison between cell types...")
  }
  
  cell_df <- res
  cell_df$labels <- labels
  cell_df <- cell_df[, c("Gate", "labels")]
  cell_df <- na.omit(cell_df)
  
  gate_table <- cell_df
  gate_table$Cell <- stri_c(gate_table$labels, as.character(1:nrow(gate_table)), sep = "__")
  
  gate_parsed <- parse_gate_table(gate_table, T, F)
  gate_parsed <- gate_parsed$gate_table
  gate_parsed$Cell <- sapply(str_split(gate_parsed$Cell,"__"), `[`, 1)
  
  new_gate_table <- data.frame(Cell = unique(gate_parsed$Cell), Gate = "")
  rownames(new_gate_table) <- new_gate_table$Cell
  
  markers <- colnames(gate_parsed[-1])
  celltypes <- unique(gate_parsed$Cell)
  
  marker_df <- data.frame(Marker = markers, Pos = NA, Neg = NA)
  rownames(marker_df) <- marker_df$Marker
  
  data <- gate_parsed
  data <- data.frame(lapply(data, factor), check.names = F)
  
  cell_markers <- vector("list", length(celltypes))
  names(cell_markers) <- celltypes
  pairs <- as.matrix(comboGrid(as.character(celltypes), as.character(celltypes), repetition = F))
  
  for(i in 1:nrow(pairs)){
    
    sig <- c()
    signs <- c()
    
    pair <- pairs[i, ]
    
    # pair <- pairs[1, ]
    
    
    c <- pair[1]
    c2 <- pair[2]
    
    # c <- "Plasmacytoid_DC_cells"
    # c2 <- "NK_cells"
    
    
    if(verbose){
      message(stri_c(" - (", c,", ", c2, ")", sep = ""))
    }
    
    to_exclude <- c()
    
    gates3 <- gate_parsed
    gates3 <- gates3[gates3$Cell %in% c(c,c2), ]
    
    data <- gates3
    data <- data[, c("Cell", markers)]
    
    s1 <- w1 <- which(data$Cell == c)
    s2 <- w2 <- which(data$Cell == c2)
    
    if(length(w1) > sampling_imp_vars){
      s1 <- sample(w1, sampling_imp_vars)
    }
    
    if(length(w2) > sampling_imp_vars){
      s2 <- sample(w2, sampling_imp_vars)
    }
    
    data <- data[c(s1,s2), ]
    data <- data.frame(lapply(data, factor), check.names = F)
    
    mas <- c()
    
    nzv <- nearZeroVar(data[, -1], saveMetrics = TRUE)
    to_exclude <- rownames(nzv)[nzv$zeroVar == T]
    
    if(length(to_exclude) > 0){
      data <- data[, -which(colnames(data) %in% to_exclude)]
    }
    
    flag <- -1
    
    tryCatch({
      control <- trainControl(method = "repeatedcv", number = 5)
      model <- train(Cell ~ ., data = data, method = "rpart", trControl = control)
      importance <- varImp(model, useModel = F)
      plot(importance)
      imp <- importance$importance
      mm <- rownames(imp)
      imp <- imp[, 1]
      names(imp) <- mm
      imp <- imp[order(imp, decreasing = T)]
      mas <- names(imp)
      flag <- 1
    },
    error = function(e){
      stop(e)
    })
    
    if(flag < 0){
      mas <- markers
    }
    
    t <- prop.table(table(gate_parsed[gate_parsed$Cell == c, mas[1]]))
    w1 <- names(t)[which.max(t)]
    
    t2 <- prop.table(table(gate_parsed[gate_parsed$Cell == c2, mas[1]]))
    w2 <- names(t2)[which.max(t2)]
    
    if(w1 != w2){
      
      signs <- w1
      
      # w2 <- which.max(t2)
      # signs2 <- names(t2)[w2]
      # 
      # if(signs == signs2){
      #   if(t[signs] < t2[signs2]){
      #     signs <- ifelse(signs == "+", "-", "+")
      #   }
      # }
    
      sig <- mas[1]
      
      if(!is.null(sig)){
        top_marker <- stri_c(sig, signs, sep = "")
        int_pos <- top_marker %in% cell_markers[[c]]
        sg <- ifelse(signs == "+", "-", "+")
        
        # Top discriminant marker is not present in gate table
        if(!int_pos){
          cell_markers[[c]] <- c(cell_markers[[c]], top_marker)
        }
        
        ## Check if the top marker is present in c2 
        top_marker_opp <- gsub(stri_c("\\", signs, sep = ""), sg, top_marker)
        int_c2 <- top_marker_opp %in% cell_markers[[c2]]
        
        if(!int_c2){
          cell_markers[[c2]] <- c(cell_markers[[c2]], top_marker_opp)
        }
      }
    }
  }

  # g_temp <- cell_markers
  # g_temp <- lapply(g_temp, unique)
  # g_temp <- lapply(g_temp, unique)
  # g_temp <- lapply(g_temp, gsub, pattern = "\\+|-", replacement = "")
  # 
  # to_remove <- sapply(1:length(g_temp), function(i){
  #   ns <- names(g_temp[i])
  #   el <- g_temp[[i]]
  #   w <- which(duplicated(el))
  #   to_remove <- el[w]
  #   l <- list(which(el == to_remove))
  #   names(l) <- ns
  #   return(l)
  # })
  # 
  # for(i in 1:length(cell_markers)){
  #   ns <- names(cell_markers[i])
  #   w <- -to_remove[[ns]]
  #   if(length(to_remove[[ns]])){
  #     cell_markers[[ns]] <- cell_markers[[ns]][w]
  #   }
  # }
  
  cell_markers <- lapply(cell_markers, sort, decreasing = F)
  g <- sapply(cell_markers, stri_c, collapse = "", sep = "")
  new_gate_table <- data.frame(Cell = names(g), Gate = g)
  rownames(new_gate_table) <- NULL
  new_gate_table <- new_gate_table[new_gate_table$Gate != "", ]
  
  return(new_gate_table)
}

## This is the annotation module of GateMeClass
GateMeClass_annotate <- function(exp_matrix = NULL,
                                 gate_table = NULL, 
                                 GMM_parameterization = NULL,
                                 train_parameters = list(
                                   reference = NULL
                                 ),
                                 refine = F,
                                 k = NULL,
                                 sampling = 1,
                                 narrow_gate_table = T,
                                 verbose = T,
                                 seed = 1){
  
  gate_table <- NULL
  refine = T
  seed = 1
  exp_matrix <- testing_set
  verbose = T
  narrow_gate_table = T
  sampling <- 1
  k = NULL
  # reference <- NULL
  # labels <- colnames(m)
  GMM_parameterization = "E"
  train_parameters = list(
    reference = training_set,
    labels = training_set_lab,
    GMM_parameterization = "E"
  )

  set.seed(seed)
  
  if(!is.null(exp_matrix)){
    old_names <- rownames(exp_matrix)
    check <- check_marker_names(rownames(exp_matrix))
    rownames(exp_matrix) <- check$marker_names
    
    if(length(check$modified) > 0){
      if(verbose){
        message("GateMeClass annotate - Renaming markers of the dataset: ", 
                stri_c(old_names[check$modified], collapse = " ", sep = ""), " --> ", 
                stri_c(rownames(exp_matrix)[check$modified], collapse = " "), sep = "")
      }
    }
  }else{
    stop("The expression matrix is mandatory!")
  }
  
  if(!is.null(gate_table)){
    temp <- sapply(gate_table$Gate, function(x){
      str <- strsplit(x, "[-]|-[\\|\\^]|[+]|\\+[\\|\\^]|[*]")
      str[[1]] <- stri_remove_empty(str[[1]])
      return(str)
    })
    
    if(length(setdiff(unique(unlist(temp)), rownames(exp_matrix))) > 0){
      for(i in 1:length(old_names[check$modified])){
        x <- old_names[check$modified][i]
        gate_table$Gate <- gsub(x, rownames(exp_matrix)[check$modified][i], gate_table$Gate)
      }
    }
  }
  
  if(sampling < 1){
    n <- floor(ncol(exp_matrix) * sampling)
    s <- sample(1:ncol(exp_matrix), n)
    exp_matrix_pre_sampling <- exp_matrix
    exp_matrix <- exp_matrix[, s]
  }else{
    exp_matrix_pre_sampling <- exp_matrix
    s <- 1:ncol(exp_matrix_pre_sampling)
  }
  
  if(!is.null(train_parameters$reference) & is.null(gate_table)){
    reference <- train_parameters$reference
    old_names <- rownames(reference)
    check <- check_marker_names(rownames(reference))
    rownames(reference) <- check$marker_names
    markers <- check$marker_names
    
    if(length(check$modified) > 0){
      if(verbose){
        message("GateMeClass annotate - Renaming markers of reference dataset: ", 
                stri_c(old_names[check$modified], collapse = " ", sep = ""), " --> ", 
                stri_c(rownames(reference)[check$modified], collapse = " "), sep = "")
      }
    }
    
    common <- intersect(rownames(exp_matrix), rownames(reference))
    
    if(length(common) > 0){
      reference <- reference[rownames(reference) %in% common, ]
      train_parameters$reference <- reference
      train_parameters$seed <- seed
      train_parameters$verbose <- verbose
      gate_table <- do.call("GateMeClass_train", train_parameters)
      new_gates <- parse_gate_table(gate_table, T, T)
      markers <- colnames(new_gates$gate_table)[-1]
    }else{
      stop("No common markers between actual and reference datasets!")
    }
  }else{
    if(is.null(gate_table)){
      stop("Please, specify a gate table or a reference dataset!")
    }
    
    if(verbose){
      message("GateMeClass annotate - Parsing gate table...")
    }
    new_gates <- parse_gate_table(gate_table, narrow_gate_table, T)
    markers <- colnames(new_gates$gate_table)[-1]
  }
  
  exp_matrix_2 <- exp_matrix[markers, , drop = F]
  
  if(verbose){
    message("GateMeClass annotate - Determining the marker signature of each cell...")
  }
  
  expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  rownames(expr_markers) <- markers
  
  exp_matrix_2 <- exp_matrix[markers, , drop = F]
  expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  rownames(expr_markers) <- markers
  
  ## Obtain the signature for the cells of the dataset to be annotated
  expr_markers <- set_marker_expression(exp_matrix_2,
                                        markers, 
                                        expr_markers, 
                                        # new_gates$gate_table, 
                                        verbose = verbose, 
                                        GMM_parameterization = GMM_parameterization)
  
  expr_markers <- expr_markers[markers, ]
  expr_markers <- as.data.frame(t(expr_markers))
  
  cells <- cbind(rownames(expr_markers), expr_markers)
  colnames(cells)[1] <- "Cell"
  
  expr_markers <- apply(cells[,- 1, drop = F], 1, stri_c, collapse = "")
  new_cells2 <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
  
  if(verbose){
    message("GateMeClass annotate - Cell annotation...")
  }
  
  res <- cell_classification(new_cells2, new_gates$extended_gate_table)
  
  if(sampling < 1){
    res_temp <- rep("Unclassified", ncol(exp_matrix_pre_sampling))
    res_temp[s] <- res$labels
    cell_signatures <- data.frame(Cell = rep(NA, ncol(exp_matrix_pre_sampling)), 
                                  Gate = rep(NA, ncol(exp_matrix_pre_sampling)), 
                                  Celltype = rep(NA, ncol(exp_matrix_pre_sampling)))
    
    
    cell_signatures[s, ] <- res$cell_signatures
    cell_signatures[setdiff(1:ncol(exp_matrix_pre_sampling), s), "Celltype"] <- "Unclassified"
    
    res_temp <- list(labels = res_temp, marker_table = res$marker_table, cell_signatures = cell_signatures)
    res <- res_temp
  }
  
  real_uncl <- which(res$labels == "Unclassified" & (1:ncol(exp_matrix_pre_sampling)) %in% s)
  not_real_uncl <- which(res$labels == "Unclassified" & !(1:ncol(exp_matrix_pre_sampling)) %in% s)
  tot_uncl <- c(real_uncl, not_real_uncl)
  not_uncl <- which(res$labels != "Unclassified")
  real_not_uncl <- c(real_uncl, not_uncl)
  
  uncl_prec <- length(res$labels[tot_uncl])
  
  ## Refinement of the unclassified cells using K-NN classification
  if(refine & uncl_prec > 0 & uncl_prec < ncol(exp_matrix_pre_sampling)){
    if(verbose){
      message("GateMeClass annotate - Refinement of the labels...")
    }
    
    if(is.null(k)){
      t <- table(res$labels)
      min <- t[which.min(t)]
      if(min > 1){
        k <- min - 1
      }else{
        k <- 1
      }
      # k <- floor(sqrt(tt))
    }
    
    # w <- which(res$labels == "T cells")
    # s <- sample(w, 2000)
    # res$labels[s] <- "Unclassified"
    
    # train_labels <- factor(res$labels[not_uncl])
    # training_set <- data.frame(labels = train_labels, t(exp_matrix_pre_sampling[, not_uncl]))
    
    training_set <- t(exp_matrix_pre_sampling)[real_not_uncl, ]
    control <- t(exp_matrix_pre_sampling)[tot_uncl, ]
  
    nn <- FNN::knnx.index(training_set, control, k = k+1)
    nn[1:length(real_uncl), 1] <- 0
    
    if(length(real_uncl) < nrow(nn)){
      nn[(length(real_uncl)+1):nrow(nn), k+1] <- 0
    }
    
    knn_res <- apply(nn, 1, function(r){
      t <- table(res$labels[real_not_uncl][r])
      max <- max(t)
      if(length(max) > 1){
        return("Unclassified")
      }else{
        return(names(t)[which.max(t)])
      }
    })
  
    new_labels <- rep(NA, length(res$labels))
    new_labels[tot_uncl] <- knn_res
    knn_res <- new_labels[res$labels == "Unclassified"]
    
    # training_set <- data.frame(labels = res$labels[not_uncl], t(exp_matrix_pre_sampling[, not_uncl]))
    # control <- t(exp_matrix_pre_sampling[, res$labels == "Unclassified"])
    # control <- t(exp_matrix_pre_sampling[, res$labels == "Unclassified"])
    
    # ctrl <- trainControl(method="none")
    # knnFit <- train(labels ~ ., data = training_set, method = "knn", trControl = ctrl, tuneGrid = expand.grid(k = k))
    # knn_res <- predict(knnFit, newdata = control, type = "prob")
    # knn_res <- predict(knnFit, newdata = control)
    
    ####### TO CHECK ############################################
    # perc_to_remove <- 1 / k
    # # knn_res$Unclassified <- knn_res$Unclassified - perc_to_remove
    # knn_res <- apply(knn_res, 1, function(r){
    #   return(colnames(knn_res)[which.max(r)])
    # })
    #############################################################
    
    # Class R package:
    # system.time(knn_res <- knn(training_set[, -1],
    #                control,
    #                cl = factor(res$labels[!res$labels %in% c("Unclassified")]),
    #                k = k,
    #                prob = T))
    
    res$labels[res$labels == "Unclassified"] <- as.character(knn_res)
    res$cell_signatures[res$cell_signatures$Celltype == "Unclassified", "Celltype"] <- as.character(knn_res)
  }
  
  signatures <- res$cell_signatures

  gate_ext <- sapply(signatures$Gate, function(v){
    split <- unlist(str_split(v, ""))
    sig <- stri_c(markers, split, collapse = "", sep = "")
    return(sig)
  })
  
  res$cell_signatures[, "Gate"] <- gate_ext
  res <- list(labels = res$labels,
              gate_table = gate_table,
              cell_signatures = res$cell_signatures)
  
  return(res)
}

















