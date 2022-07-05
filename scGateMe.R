# options(warn=1)

## This function generates the possible value combinations of the markers of a cell type 
generate_set_values <- function(v, cell){
  
  sets_all <- list()
  sets_to_filter <- list()
  sets_to_filter2 <- list()
  
  special1 <- grep("[\\|]", v)
  special2 <- grep("[\\^]", v)
  
  if(length(special1) == 1 | length(special2) == 1){
    stop("Error! Special characters (^, |) must be set in more than 2 markers!")
  }
  
  for(m in v){
    switch(m,
           `+` = {
             set <- list(c("+"))
             set_not <- list()
             set_not2 <- list()
           },
           # `++` = {
           #   set <- list(c("*", "+", "++"))
           #   set_not <- list()
           # },
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
             set <- list(c("+", "-"))   #### to check ######
             set_not <- list()
             set_not2 <- list(c("+"), c("-"))
             # set_not2 <- list(c("-"))
             # set_not <- list()
           },
           `-^`={
             set <- list(c("+", "-"))   #### to check ######
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
  gate <- apply(sets_all_filtered_temp, 1, stringi::stri_c, collapse = "")
  sets_all_filtered$Gate <- gate
  sets_all_filtered <- sets_all_filtered[, c("Cell", "Gate")]
  
  return(sets_all_filtered)
}

## Read the gate table and generate the possbile marker signature of ceach cell type
parse_gate_table <- function(gate_table, narrow_gate_table, extended_gate_table){
  
  # gate_table <- gates
  # narrow_gate_table = T
  # extended_gate_table = T
  
  if(any(duplicated(gate_table$Cell))){
    stop("Error! The gate table must contains uniquely defined cell types!")
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
    stop("Error! Cell names cannot contains special characters (e.g., *, ^)!")
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
      
      # to_add <- generate_set_values(v, cell)
      
      
      return(to_add)
    })
    
    extended_gate_table <- as.data.frame(rbindlist(df_list))
    to_delete <- stringi::stri_c(colnames(gate_table[, -1]), "*", collapse = "")
    extended_gate_table <- extended_gate_table[extended_gate_table$Gate != to_delete, ]
    
    return(list(gate_table = gate_table, extended_gate_table = extended_gate_table))
  }

  return(list(gate_table = gate_table))
}

set_marker_expression_GMM <- function(X, indexes, mm, GMM_parameterization){
  
  # X <- m["CD19", ]
  # gmm_criteria = "BIC"
  # mm = 2
  # rr = 0.05
  # indexes = first$indexes
  
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
  
  cl <- Mclust(test, G = 2, verbose = F, modelNames = GMM_parameterization)
  pred <- predict.Mclust(cl, X)
  temp <- pred$classification
  means <- cl$parameters$mean
  means <- sort(means)
  temp[temp == names(means)[1]] <- "-"
  temp[temp == names(means)[2]] <- "+"
  return(temp)
}

## This function set the marker signature of each cell
set_marker_expression <- function(exp_matrix, markers,
                                  expr_markers, 
                                  gates, 
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

  queue <- list(list(indexes = 1:ncol(exp_matrix), markers = markers))

  while(length(queue) > 0){
    bimodal_markers <- c()
    not_bimodal_markers <- c()

    ## Pop operation
    first <- queue[[1]]
    queue <- queue[-1]

    for(m in first$markers){
      X <- exp_matrix[m, first$indexes]
      marker_expr <- set_marker_expression_GMM(X, indexes = first$indexes, m, GMM_parameterization)
      
      if(length(table(marker_expr)) > 1){
        bimodal_markers <- c(bimodal_markers, m)

        if(verbose){
          message(stringi::stri_c(" - ", stringi::stri_c(bimodal_markers, collapse = " ", sep = ""), collapse = " ", sep =))
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
      comb_list2 <- sapply(comb_list, stringi::stri_c, collapse = "", sep = "")
      l <- expr_markers[bimodal_markers, first$indexes, drop = F]
      l <- sapply(l, stringi::stri_c, collapse = "", sep = "")
      
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
  # gates <- new_gates2$extended_gate_table
  
  df_gate <- marker_table
  colnames(df_gate)[1] <- "Cell_ID"
  
  pos <- sapply(gates$Gate, function(g){
    gs <- unlist(str_split(g, pattern = ""))
    w <- which(gs != "*")
    w <- stringi::stri_c(w, collapse = "_", sep = "")
    return(w)
  })
  
  gates$Pos <- pos
  
  clean <- sapply(gates$Gate, function(g){
    gs <- unlist(str_split(g, pattern = ""))
    w <- which(gs != "*")
    clean_gs <- stringi::stri_c(gs[w], collapse = "", sep = "")
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
      x <- stringi::stri_c(gate_split[index], collapse = "", sep = "")
      return(x)
    })
    
    cell_split <- sapply(1:length(df_gate$Gate), function(j){
      gate_split <- unlist(str_split(df_gate$Gate[j], pattern = ""))
      x <- stringi::stri_c(gate_split[index], collapse = "", sep = "")
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
      
      c_split <- cell_split[[i]]
      
      test <- sapply(1:length(gates$Gate), function(j){
        index_split <- gates_split[[j]]
        x <- stringi::stri_c(c_split[index_split], collapse = "", sep = "")
        return(gates$Gate[j] == x)
      })
      
      if(sum(test) > 1){
        labels[i] <- "Unclassified"
      }else if(sum(test) == 1){
        labels[i] <- gates$Cell[test]
      }
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
    marker_names[to_modify] <- stringi::stri_c("P_", marker_names[to_modify], sep = "")
  }
  
  to_modify2 <- grep("[\\+]|[\\-]|[\\|]|[\\^]|[\\*]", marker_names)
  
  if(length(to_modify2) > 0){
    marker_names[to_modify2] <- gsub("[\\+]|[\\-]|[\\|]|[\\^]|[\\*]", "_", marker_names[to_modify2])
  }
  
  modified <- c(to_modify, to_modify2)
  return(list(marker_names = marker_names, modified = modified))
}

scGateMe_train <- function(reference,
                           labels, 
                           imp_feature_thr = "all",
                           GMM_parameterization = "V",
                           sampling = "none",
                           sampling_perc = 0.01,
                           perc.over = 100,
                           perc.under = 0,
                           sampling_k = min(5, table(labels)-1),
                           sampling_feature_method = "all",
                           sampling_imp_vars = 1000,
                           thr_perc = -1, 
                           seed = 1,
                           verbose = T){

# reference <- m
# labels <- colnames(m)
# GMM_parameterization = "E"
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
# sampling_k = 2
  
  # options(warn=1)
  set.seed(seed)
  
  if(!is.null(reference)){
    old_names <- rownames(reference)
    check <- check_marker_names(rownames(reference))
    rownames(reference) <- check$marker_names
    
    if(length(check$modified) > 0){
      if(verbose){
        message("scGateMe train - Renaming markers of reference dataset: ", 
                stringi::stri_c(old_names[check$modified], collapse = " ", sep = ""), 
                " --> ", 
                stringi::stri_c(rownames(reference)[check$modified], collapse = " ", sep = ""))
      }
    }
  }else{
    message("scGateMe train - Error! The expression matrix is mandatory!")
    stop()
  }
  
  if(sampling == "class"){
    if(verbose){
      message("scGateMe train - Executing SMOTE algorithm to balance the training set...")
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

      newData <- SMOTE(labels ~ ., temp, perc.over = 100, k = sampling_k, perc.under = 0)
      new_reference <- rbind(new_reference, newData)
    }
    
    labels <- new_reference$labels 
    new_reference <- new_reference[, -which(colnames(new_reference) == "labels")]
    new_reference <- t(new_reference)
    reference <- new_reference
  }
  
  markers <- rownames(reference)
  celltypes <- factor(unique(labels))
  gate_table <- data.frame(Cell = c("Unknown"), Gate = stringi::stri_c(markers, "+", collapse = "", sep = ""))
  
  new_gates <- parse_gate_table(gate_table, T, T)
  reference_2 <- reference[colnames(new_gates$gate_table)[-1], , drop = F]
  
  if(verbose){
    message("scGateMe train - Determining the marker signature of each cell...")
  }
  
  expr_markers <- data.frame(matrix(ncol = ncol(reference_2), nrow = nrow(reference_2)))
  rownames(expr_markers) <- colnames(new_gates$gate_table)[-1]
  
  reference_2 <- reference[colnames(new_gates$gate_table)[-1], , drop = F]
  expr_markers <- data.frame(matrix(ncol = ncol(reference_2), nrow = nrow(reference_2)))
  rownames(expr_markers) <- colnames(new_gates$gate_table)[-1]
  
  ## Obtain the signature for the cells of the dataset to be annotated
  expr_markers <- set_marker_expression(reference_2,
                                        colnames(new_gates$gate_table)[-1], 
                                        expr_markers, 
                                        new_gates$gate_table, 
                                        verbose = verbose, 
                                        GMM_parameterization = GMM_parameterization)
  
  expr_markers <- expr_markers[colnames(new_gates$gate_table)[-1], ]
  expr_markers <- as.data.frame(t(expr_markers))
  
  cells <- cbind(rownames(expr_markers), expr_markers)
  colnames(cells)[1] <- "Cell"
  
  expr_markers <- apply(cells[,- 1], 1, stringi::stri_c, collapse = "")
  new_cells2 <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
  
  res <- list(cell_signatures = new_cells2)
  
  signatures <- res$cell_signatures
  markers <- colnames(new_gates$gate_table)[-1]
  
  gate_ext <- sapply(signatures$Gate, function(v){
    split <- unlist(str_split(v, ""))
    sig <- stringi::stri_c(markers, split, collapse = "", sep = "")
    return(sig)
  })
  
  res$cell_signatures[, "Gate"] <- gate_ext
  res <- res$cell_signatures
  
  if(verbose){
    message("scGateMe train - Pairwise comparison between cell types...")
  }
  
  cell_df <- res
  cell_df$labels <- labels
  cell_df <- cell_df[, c("Gate", "labels")]
  cell_df <- na.omit(cell_df)
  
  gate_table <- cell_df
  gate_table$Cell <- stringi::stri_c(gate_table$labels, as.character(1:nrow(gate_table)), sep = "__")
  
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
    c <- pair[1]
    c2 <- pair[2]
    
    if(verbose){
      message(stringi::stri_c(" - (", c,", ", c2, ")", sep = ""))
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
    to_exclude <- rownames(nzv)[nzv$nzv == T]
    
    if(length(to_exclude) > 0){
      data <- data[, -which(colnames(data) %in% to_exclude)]
    }
    
    flag <- -1
    
    tryCatch({
      control <- trainControl(method = "repeatedcv", number = 5)
      model <- train(Cell ~ ., data = data, method = "rpart", trControl = control)
      importance <- varImp(model, useModel = F)
      # plot(importance)
      imp <- importance$importance
      mm <- rownames(imp)
      mm <- gsub("`", "", mm)
      mm <- gsub("\\+", "", mm)
      mm <- gsub("\\-", "", mm)
      imp <- imp[, 1]
      names(imp) <- mm
      imp <- imp[order(imp, decreasing = T)]
      
      if(imp_feature_thr == "all"){
        mas <- names(imp[imp > 0])
      }else if(imp_feature_thr == "median"){
        mas <- names(imp[imp >= median(imp)])
      }else if(imp_feature_thr == "GMM"){
        cl <- Mclust(imp, G = 2, modelNames = "E", verbose = F)
        mas <- names(imp[cl$classification == "2"])
      }else{
        messagae("Error! The specified value of parameter 'imp_feature_thr' does not exist!")
        stop()
      }
      flag <- 1
    },
    error = function(e){
      message(e)
      stop()
    })
    
    if(flag < 0){
      mas <- markers
    }
    
    for(k in mas){
      
      # k <- "CD14"
      
      t <- prop.table(table(gate_parsed[gate_parsed$Cell == c, k]))
      max <- max(t)
      
      if(flag > 0){
        if(max >= thr_perc){
          sign <- names(t)[which.max(t)]
          signs <- c(signs, sign)
          sig <- c(sig, k)
        }
      }else{
        if(max > 0.9){
          sign <- names(t)[which.max(t)]
          signs <- c(signs, sign)
          sig <- c(sig, k)
        }
      }
    }
    
    if(!is.null(sig)){
      top_marker <- stringi::stri_c(sig[1], signs[1], sep = "")
      int_pos <- top_marker %in% cell_markers[[c]]
      sg <- ifelse(signs[1] == "+", "-", "+")
      
      # Top discriminant marker is not present in gate table
      if(!int_pos){
        cell_markers[[c]] <- c(cell_markers[[c]], top_marker)
      }
      
      ## Check if the top marker is present in c2 
      top_marker_opp <- gsub(stringi::stri_c("\\", signs[1], sep = ""), sg, top_marker)
      int_c2 <- top_marker_opp %in% cell_markers[[c2]]
      
      if(!int_c2){
        cell_markers[[c2]] <- c(cell_markers[[c2]], top_marker_opp)
      }
    }
  }

  g <- sapply(cell_markers, stringi::stri_c, collapse = "", sep = "")
  new_gate_table <- data.frame(Cell = names(g), Gate = g)
  rownames(new_gate_table) <- NULL
  new_gate_table <- new_gate_table[new_gate_table$Gate != "", ]
  
  return(new_gate_table)
}

## This is the core function for the classification of the single cells
scGateMe_classify <- function(exp_matrix,
                              gate_table = NULL, 
                              GMM_parameterization = "V",
                              train_parameters = list(
                                reference = NULL
                              ),
                              refine = F,
                              k = NULL,
                              sampling = 1,
                              narrow_gate_table = T,
                              verbose = T,
                              seed = 1){
  
  # gate_table <- gate
  # refine = T
  # seed = 1
  # exp_matrix <- m
  # verbose = T
  # narrow_gate_table = T
  # sampling <- 0.25
  # k = NULL
  # reference <- NULL
  # labels <- colnames(m)
  # GMM_parameterization = "V"
  # train_parameters = list(
  #   reference = NULL
  # )
  
  set.seed(seed)
  
  if(!is.null(exp_matrix)){
    old_names <- rownames(exp_matrix)
    check <- check_marker_names(rownames(exp_matrix))
    rownames(exp_matrix) <- check$marker_names
    
    if(length(check$modified) > 0){
      if(verbose){
        message("scGateMe classify - Renaming markers of the dataset: ", 
                stringi::stri_c(old_names[check$modified], collapse = " ", sep = ""), 
                " --> ", 
                stringi::stri_c(rownames(exp_matrix)[check$modified], collapse = " "), sep = "")
      }
    }
  }else{
    message("scGateMe classify - Error! The expression matrix is mandatory!")
    stop()
  }
  
  if(!is.null(gate_table)){
    temp <- sapply(gate_table$Gate, function(x){
      str <- strsplit(x, "[-]|-[\\|\\^]|[+]|\\+[\\|\\^]|[*]")
      str[[1]] <- stri_remove_empty(str[[1]])
      return(str)
    })
    
    if(length(setdiff(unique(unlist(temp)), rownames(exp_matrix))) > 0){
      for(i in 1:length(old_names[modified])){
        x <- old_names[modified][i]
        gate_table$Gate <- gsub(x, rownames(exp_matrix)[modified][i], gate_table$Gate)
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
  }
  
  if(!is.null(train_parameters$reference)){
    reference <- train_parameters$reference
    old_names <- rownames(reference)
    check <- check_marker_names(rownames(reference))
    rownames(reference) <- check$marker_names
    
    if(length(check$modified) > 0){
      if(verbose){
        message("scGateMe classify - Renaming markers of reference dataset: ", 
                stringi::stri_c(old_names[check$modified], collapse = " ", sep = ""), 
                " --> ", 
                stringi::stri_c(rownames(reference)[check$modified], collapse = " "), sep = "")
      }
    }
    
    common <- intersect(rownames(exp_matrix), rownames(reference))
    
    if(length(common) > 0){
      reference <- reference[rownames(reference) %in% common, ]
      train_parameters$reference <- reference
      train_parameters$seed <- seed
      train_parameters$verbose <- verbose
      gate_table <- do.call("scGateMe_train", train_parameters)
      new_gates <- parse_gate_table(gate_table, T, T)
    }else{
      message("Error! No common markers between actual and reference datasets!")
      stop()
    }
  }else{
    if(is.null(gate_table)){
      message("scGateMe classify - Error! Please, specify a gate table or a reference dataset!")
      stop()
    }
    
    if(verbose){
      message("scGateMe classify - Parsing gate table...")
    }
    new_gates <- parse_gate_table(gate_table, narrow_gate_table, T)
  }
  
  exp_matrix_2 <- exp_matrix[colnames(new_gates$gate_table)[-1], , drop = F]
  
  if(verbose){
    message("scGateMe classify - Determining the marker signature of each cell...")
  }
  
  expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  rownames(expr_markers) <- colnames(new_gates$gate_table)[-1]
  
  exp_matrix_2 <- exp_matrix[colnames(new_gates$gate_table)[-1], , drop = F]
  expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  rownames(expr_markers) <- colnames(new_gates$gate_table)[-1]
  
  ## Obtain the signature for the cells of the dataset to be annotated
  expr_markers <- set_marker_expression(exp_matrix_2,
                                        colnames(new_gates$gate_table)[-1], 
                                        expr_markers, 
                                        new_gates$gate_table, 
                                        verbose = verbose, 
                                        GMM_parameterization = GMM_parameterization)
  
  expr_markers <- expr_markers[colnames(new_gates$gate_table)[-1], ]
  expr_markers <- as.data.frame(t(expr_markers))
  
  cells <- cbind(rownames(expr_markers), expr_markers)
  colnames(cells)[1] <- "Cell"
  
  expr_markers <- apply(cells[,- 1, drop = F], 1, stringi::stri_c, collapse = "")
  new_cells2 <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
  
  if(verbose){
    message("scGateMe classify - Cell annotation...")
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
  
  uncl_prec <- length(res$labels[res$labels == "Unclassified"])
  
  ## Refinement of the unclassified cells using K-NN classification
  if(refine & uncl_prec > 0 & uncl_prec < ncol(exp_matrix_pre_sampling)){
    if(verbose){
      message("scGateMe classify - Refinement of the labels using K-NN classification...")
    }
    
    training_set <- t(exp_matrix_pre_sampling[, !res$labels %in% c("Unclassified")])
    control <- t(exp_matrix_pre_sampling[, res$labels %in% c("Unclassified")])
    
    if(is.null(k)){
      t <- table(res$labels)
      tt <- t[which.min(t)]
      k <- floor(sqrt(tt))
    }
    
    knn_res <- knn(training_set,
                   control,
                   cl = factor(res$labels[!res$labels %in% c("Unclassified")]),
                   k = k,
                   l = 0,
                   prob = T)
    
    res$labels[res$labels == "Unclassified"] <- as.character(knn_res)
    res$cell_signatures[res$cell_signatures$Celltype == "Unclassified", "Celltype"] <- as.character(knn_res)
  }
  
  signatures <- res$cell_signatures
  markers <- colnames(new_gates$gate_table)[-1]
  
  gate_ext <- sapply(signatures$Gate, function(v){
    split <- unlist(str_split(v, ""))
    sig <- stringi::stri_c(markers, split, collapse = "", sep = "")
    return(sig)
  })
  
  res$cell_signatures[, "Gate"] <- gate_ext
  res <- list(labels = res$labels,
              gate_table = gate_table,
              cell_signatures = res$cell_signatures)
  
  return(res)
}

















