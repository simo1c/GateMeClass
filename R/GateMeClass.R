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
           `m` = {
             set <- list(c("m"))
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
  gate <- apply(sets_all_filtered_temp, 1, paste0, collapse = "")
  sets_all_filtered$Gate <- gate
  sets_all_filtered <- sets_all_filtered[, c("Cell", "Gate")]

  return(sets_all_filtered)
}

## Read the gate table and generate the possbile marker signature of ceach cell type
parse_marker_table <- function(marker_table, narrow_marker_table, extended_marker_table){

  if(any(duplicated(marker_table$Cell))){
    stop("The gate table must contains uniquely defined cell types!")
  }

  if(narrow_marker_table){
    temp <- sapply(marker_table$Gate, function(x){
      str <- strsplit(x, "[m]|-[\\|\\^]|[-]|-[\\|\\^]|[+]|\\+[\\|\\^]|[*]")
      str[[1]] <- stri_remove_empty(str[[1]])
      return(str)
    })

    names(temp) <- marker_table$Cell
    markers <- unique(unlist(temp))
    df_gates <- data.frame(matrix(nrow = nrow(marker_table), ncol = length(markers)))
    rownames(df_gates) <- marker_table$Cell
    colnames(df_gates) <- markers

    marker_table_exploded <- strsplit(perl=T, marker_table$Gate, '(?![m\\||\\+\\||\\-\\||m\\^|\\+\\^|\\-\\^])')
    signs <- lapply(marker_table_exploded, function(x){ x[x %in% c("m^", "-^", "+^", "m|", "-|", "+|", "m", "+", "-", "*")]})
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
    marker_table <- df_gates
  }

  celltypes <- marker_table$Cell

  if(length(grep("[\\*|\\||\\^]", celltypes)) > 0){
    stop("Cell names cannot contains special characters (e.g., *, ^)!")
  }

  if(extended_marker_table){
    df_list <- apply(marker_table, 1, function(v){
      to_add <- generate_set_values(v[-1], v[1])
      return(to_add)
    })

    extended_marker_table <- as.data.frame(rbindlist(df_list))
    to_delete <- paste0(colnames(marker_table[, -1]), "*", collapse = "")
    extended_marker_table <- extended_marker_table[extended_marker_table$Gate != to_delete, ]

    w <- apply(marker_table, 2, function(c){
      return(any(c == "m"))
    })

    w <- w[-1]

    return(list(marker_table = marker_table, extended_marker_table = extended_marker_table, bimodal = !w))
  }


  w <- apply(marker_table, 2, function(c){
    return(any(c == "m"))
  })

  w <- w[-1]

  return(list(marker_table = marker_table, bimodal = !w))
}

set_marker_expression_GMM <- function(X, GMM_parameterization, type, RSS){

  if(length(X) >= 100){

    sample <- X

    if(RSS){
      ############ Ranked Set Sampling (RSS) #########################
      sample <- X
      cycles <- floor(length(sample) / 4)
      n_unit <- 2
      n <- length(sample)
      samples <- cycles * n_unit
      indexes <- matrix(sample(1:n, samples * n_unit), nrow = n_unit)
      sel_samples <- matrix(sample[indexes], nrow = n_unit)
      ################################################################

      if(moments::skewness(X) < 0){
        test <- apply(sel_samples, 2, base::min)

        index_sel <- sapply(1:ncol(sel_samples), function(i){
          w <- base::which.min(sel_samples[, i])
          return(indexes[w, i])
        })

      }else{
        test <- apply(sel_samples, 2, base::max)

        index_sel <- sapply(1:ncol(sel_samples), function(i){
          w <- base::which.max(sel_samples[, i])
          return(indexes[w, i])
        })
      }
    }else{
      test <- X
    }
  }else if(length(X) >= 2){
    test <- X
  }else{
    return(c("*"))
  }

  if(type){
    cl <- Mclust(test, G = 2, verbose = F, modelNames = GMM_parameterization)
    max <- which.max(cl$parameters$mean)
    min <- which.min(cl$parameters$mean)
    temp <- rep(NA, length(X))
    if(RSS){
      pred <- predict.Mclust(cl, X[-index_sel])
      temp[index_sel] <- ifelse(cl$classification == max, "+", "-")
      temp[-index_sel] <- ifelse(pred$classification == max, "+", "-")
      
      ########## Correction for V ##########################
      temp[X < cl$parameters$mean[min] & temp == "+"] <- "-"
      temp[X > cl$parameters$mean[max] & temp == "-"] <- "+"
      #######################################################
    }else{
      temp <- ifelse(cl$classification == max, "+", "-")
    }
  }else{

    cl <- Mclust(test, G = 3, verbose = F, modelNames = GMM_parameterization)
    n_cl <- length(table(cl$classification))

    if(n_cl == 2){
      cl <- Mclust(test, G = 2, verbose = F, modelNames = GMM_parameterization)
      max <- which.max(cl$parameters$mean)
      min <- which.min(cl$parameters$mean)

      if(RSS){
        pred <- predict.Mclust(cl, X[-index_sel])
        temp[index_sel] <- ifelse(cl$classification == max, "+", "-")
        temp[-index_sel] <- ifelse(pred$classification == max, "+", "-")
        
        ########## Correction for V ##########################
        temp[X < cl$parameters$mean[min] & temp == "+"] <- "-"
        temp[X > cl$parameters$mean[max] & temp == "-"] <- "+"
        #######################################################
      }else{
        temp <- ifelse(cl$classification == max, "+", "-")
      }
    }else{

      cl2 <- sort(cl$parameters$mean)
      min <- names(cl2)[1]
      mid <- names(cl2)[2]
      max <- names(cl2)[3]

      temp <- rep(NA, length(X))

      if(RSS){
        pred <- predict.Mclust(cl, X[-index_sel])
        temp[index_sel] <- cl$classification
        temp[-index_sel] <- pred$classification
        temp <- plyr::mapvalues(temp, c(min, mid, max), c("-", "m", "+"))
        
        ########## Correction for V ##########################
        temp[X < cl$parameters$mean[min] & temp == "+"] <- "-"
        temp[X < cl$parameters$mean[mid] & temp == "+"] <- "m"
        
        temp[X > cl$parameters$mean[max] & temp == "-"] <- "+"
        temp[X > cl$parameters$mean[mid] & temp == "-"] <- "m"
        
        temp[X < cl$parameters$mean[min] & temp == "m"] <- "-"
        temp[X > cl$parameters$mean[max] & temp == "m"] <- "+"
        ######################################################
        
      }else{
        temp <- plyr::mapvalues(cl$classification, c(min, mid, max), c("-", "m", "+"))
      }
    }
  }
  return(temp)
}

## This function set the marker signature of each cell
set_marker_expression <- function(exp_matrix,
                                  markers,
                                  type,
                                  expr_markers,
                                  verbose,
                                  GMM_parameterization,
                                  RSS){

  queue <- list(list(indexes = 1:ncol(exp_matrix), markers = markers, type = type))

  while(length(queue) > 0){
    bimodal_markers <- c()
    not_bimodal_markers <- c()

    ## Pop operation
    first <- queue[[1]]
    queue <- queue[-1]

    for(m in first$markers){
      X <- exp_matrix[m, first$indexes]
      
      if(length(GMM_parameterization) > 1){
        marker_expr <- set_marker_expression_GMM(X, GMM_parameterization[m], first$type[m], RSS[m])        
      }else{
        marker_expr <- set_marker_expression_GMM(X, GMM_parameterization, first$type[m], RSS)
      }
      
      if(length(table(marker_expr)) > 1){
        bimodal_markers <- c(bimodal_markers, m)

        if(verbose){
          message(paste0(" - ", paste0(bimodal_markers, collapse = " ", sep = ""), collapse = " ", sep = ""))
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
      comb_list2 <- sapply(comb_list, paste0, collapse = "", sep = "")
      l <- expr_markers[bimodal_markers, first$indexes, drop = F]
      l <- sapply(l, paste0, collapse = "", sep = "")

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

  df_gate <- marker_table
  colnames(df_gate)[1] <- "Cell_ID"

  pos <- sapply(gates$Gate, function(g){
    gs <- unlist(strsplit(perl=T, g, '(?![m])'))
    w <- which(gs != "*")
    w <- paste0(w, collapse = "_", sep = "")
    return(w)
  })

  gates$Pos <- pos

  clean <- sapply(gates$Gate, function(g){
    gs <- unlist(strsplit(perl=T, g, '(?![m])'))
    w <- which(gs != "*")
    clean_gs <- paste0(gs[w], collapse = "", sep = "")
    return(clean_gs)
  })

  gates$Gate <- clean

  gates$N <- sapply(gates$Gate, function (g){
    return(sum(unlist(str_split(g, pattern = "")) == "+"))
  })

  gates$N2 <- sapply(gates$Gate, function (g){
    return(sum(unlist(str_split(g, pattern = "")) == "-"))
  })

  m_check <- length(grep("m", gates$Gate))

  if(m_check > 0){
    gates$N3 <- sapply(gates$Gate, function (g){
      return(sum(unlist(str_split(g, pattern = "")) == "m"))
    })
  }

  gates <- gates[order(gates$N2, decreasing = F), ]

  if(m_check > 0){
    gates <- gates[order(gates$N3, decreasing = F), ]
  }

  gates <- gates[order(gates$N, decreasing = F), ]

  labels <- rep("Unclassified", nrow(df_gate))
  multiple_cl <- rep(0, nrow(df_gate))

  ###  Optimized base case
  if(length(table(gates$Pos)) == 1){

    index <- as.numeric(unlist(str_split(gates$Pos[1], pattern = "_")))

    gates_split <- sapply(1:length(gates$Gate), function(j){
      gate_split <- unlist(strsplit(perl=T, gates$Gate[j], '(?![m])'))
      x <- paste0(gate_split[index], collapse = "", sep = "")
      return(x)
    })

    cell_split <- sapply(1:length(df_gate$Gate), function(j){
      gate_split <- unlist(strsplit(perl=T, df_gate$Gate[j], '(?![m])'))
      x <- paste0(gate_split[index], collapse = "", sep = "")
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
      # sgate <- unlist(str_split(g, pattern = ""))

      gate <- unlist(strsplit(perl=T, g, '(?![m])'))


      return(gate)
    })

    for(i in 1:length(df_gate$Gate)){

      c_split <- cell_split[[i]]

      test <- sapply(1:length(gates$Gate), function(j){
        index_split <- gates_split[[j]]
        x <- paste0(c_split[index_split], collapse = "", sep = "")
        y <- gates$Gate[j] == x
        names(y) <- paste0(gates[j, c("N", "N2")], collapse = "_")
        return(y)
      })

      if(sum(test) == 1){
        labels[i] <- gates$Cell[test]
      }else if(sum(test) > 1){
        t <- table(names(test)[test])
        if(any(t > 1)){
          labels[i] <- "Unclassified"
        }else{
          ind <- which(test)
          ind <- ind[length(ind)]
          labels[i] <- gates$Cell[ind]
        }
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
    marker_names[to_modify] <- paste0("P_", marker_names[to_modify], sep = "")
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
                              GMM_parameterization = "E",
                              sampling_imp_vars = 0.05,
                              seed = 1,
                              verbose = T){

  set.seed(seed)

  if(!is.null(reference)){
    old_names <- rownames(reference)
    check <- check_marker_names(rownames(reference))
    rownames(reference) <- check$marker_names

    if(length(check$modified) > 0){
      if(verbose){
        message("GateMeClass train - Renaming markers of reference dataset: ",
                paste0(old_names[check$modified], collapse = " ", sep = ""), " --> ",
                paste0(rownames(reference)[check$modified], collapse = " ", sep = ""))
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

  markers <- rownames(reference)
  celltypes <- factor(unique(labels))
  reference_2 <- reference[markers, , drop = F]

  if(verbose){
    message(paste0("GateMeClass annotate - Determining the marker signature of each cell..."))
  }

  expr_markers <- data.frame(matrix(ncol = ncol(reference_2), nrow = nrow(reference_2)))
  rownames(expr_markers) <- markers

  reference_2 <- reference[markers, , drop = F]
  expr_markers <- data.frame(matrix(ncol = ncol(reference_2), nrow = nrow(reference_2)))
  rownames(expr_markers) <- markers

  ####### Mid expression?? ########
  m_neg <- paste0(markers, "-")
  m_mid <- paste0(markers, "mid")
  m_hi <- paste0(markers, "hi")


  type_neg <- sapply(m_neg, function(x){
    g <- grep(x, unique(labels), ignore.case = T)
    return(length(g) == 0)
  })

  type_mid <- sapply(m_mid, function(x){
    g <- grep(x, unique(labels), ignore.case = T)
    return(length(g) == 0)
  })

  type_hi <- sapply(m_hi, function(x){
    g <- grep(x, unique(labels), ignore.case = T)
    return(length(g) == 0)
  })

  type <- type_neg | type_mid | type_hi
  names(type) <- gsub("mid|[-]|hi", "", names(type))
  ################################


  ## Obtain the signature for the cells of the dataset to be annotated
  expr_markers <- set_marker_expression(reference_2,
                                        markers,
                                        type,
                                        expr_markers,
                                        # new_gates$marker_table,
                                        verbose = verbose,
                                        GMM_parameterization = GMM_parameterization,
                                        RSS = ifelse(GMM_parameterization == "E", F, T))

  expr_markers <- expr_markers[markers, ]
  expr_markers <- as.data.frame(t(expr_markers))

  cells <- cbind(rownames(expr_markers), expr_markers)
  colnames(cells)[1] <- "Cell"

  expr_markers <- apply(cells[,- 1], 1, paste0, collapse = "")
  new_cells2 <- data.frame(Cell = names(expr_markers), Gate = expr_markers)

  res <- list(cell_signatures = new_cells2)
  signatures <- res$cell_signatures

  gate_ext <- sapply(signatures$Gate, function(v){
    split <- unlist(str_split(v, ""))
    sig <- paste0(markers, split, collapse = "", sep = "")
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

  marker_table <- cell_df
  marker_table$Cell <- paste0(marker_table$labels, as.character(1:nrow(marker_table)), sep = "__")

  gate_parsed <- parse_marker_table(marker_table, T, F)
  gate_parsed <- gate_parsed$marker_table
  gate_parsed$Cell <- sapply(str_split(gate_parsed$Cell,"__"), `[`, 1)
  gate_parsed$labels <- marker_table$labels

  new_marker_table <- data.frame(Cell = unique(gate_parsed$Cell), Gate = "")
  rownames(new_marker_table) <- new_marker_table$Cell

  markers <- colnames(gate_parsed[-which(colnames(gate_parsed) %in% c("Cell", "labels"))])

  celltypes <- unique(gate_parsed$labels)

  marker_df <- data.frame(Marker = markers, Pos = NA, Neg = NA)
  rownames(marker_df) <- marker_df$Marker

  data <- gate_parsed
  data <- data.frame(lapply(data, factor), check.names = F)

  cell_markers <- vector("list", length(celltypes))
  names(cell_markers) <- celltypes

  pairs <- as.matrix(comboGrid(as.character(celltypes), as.character(celltypes), repetition = F))
  sampling_imp_vars <- floor(sampling_imp_vars * ncol(reference))

  for(i in 1:nrow(pairs)){
    sig <- c()
    signs <- c()

    pair <- pairs[i, ]
    c <- pair[1]
    c2 <- pair[2]

    if(verbose){
      message(paste0(" - (", c,", ", c2, ")", sep = ""))
    }

    to_exclude <- c()

    gates3 <- gate_parsed
    gates3 <- gates3[gates3$labels %in% c(c,c2), ]

    data <- gates3
    data <- data[, c("labels", markers)]

    s1 <- w1 <- which(data$labels == c)
    s2 <- w2 <- which(data$labels == c2)

    if(length(w1) > sampling_imp_vars){
      s1 <- sample(w1, sampling_imp_vars)
    }

    if(length(w2) > sampling_imp_vars){
      s2 <- sample(w2, sampling_imp_vars)
    }

    data <- data[c(s1,s2), ]
    data <- data.frame(lapply(data, factor), check.names = F)
    data[, -1] <- data.frame(lapply(data[, -1], factor, levels = c("-", "m", "+"), ordered = T), check.names = F)

    mas <- c()

    nzv <- nearZeroVar(data[, -1], saveMetrics = TRUE)
    to_exclude <- rownames(nzv)[nzv$zeroVar == T]

    if(length(to_exclude) > 0){
      data <- data[, -which(colnames(data) %in% to_exclude)]
    }

    flag <- -1

    tryCatch({
      control <- trainControl(method = "none")
      model <- caret::train(labels ~ ., data = data, method = "rpart", trControl = control)
      importance <- varImp(model, useModel = F)
      imp <- importance$importance
      mm <- rownames(imp)
      imp <- imp[, 1]
      names(imp) <- mm
      imp <- imp[order(imp, decreasing = T)]
      mas <- names(imp)
      flag <- 1
    },
    error = function(e){
      #stop(e)
    })

    if(flag < 0){
      mas <- markers
    }

    t <- prop.table(table(gate_parsed[gate_parsed$labels == c, mas[1]]))
    w1 <- names(t)[which.max(t)]

    t2 <- prop.table(table(gate_parsed[gate_parsed$labels == c2, mas[1]]))
    w2 <- names(t2)[which.max(t2)]


    if(w1 != w2){

      signs <- w1
      sig <- mas[1]

      if(!is.null(sig)){
        top_marker <- paste0(sig, signs, sep = "")
        int_pos <- top_marker %in% cell_markers[[c]]

        sg <- w2

        # Top discriminant marker is not present in gate table
        if(!int_pos){
          cell_markers[[c]] <- c(cell_markers[[c]], top_marker)
        }

        ## Check if the top marker is present in c2
        top_marker_opp <- gsub(paste0("\\", signs, sep = ""), sg, top_marker)
        int_c2 <- top_marker_opp %in% cell_markers[[c2]]

        if(!int_c2){
          cell_markers[[c2]] <- c(cell_markers[[c2]], top_marker_opp)
        }
      }
    }
  }

  cell_markers <- lapply(cell_markers, sort, decreasing = F)
  g <- sapply(cell_markers, paste0, collapse = "", sep = "")
  new_marker_table <- data.frame(Cell = names(g), Gate = g)
  rownames(new_marker_table) <- NULL
  new_marker_table <- new_marker_table[new_marker_table$Gate != "", ]

  return(new_marker_table)
}

## This is the annotation module of GateMeClass
GateMeClass_annotate <- function(exp_matrix = NULL,
                                 marker_table = NULL,
                                 GMM_parameterization = "E",
                                 train_parameters = list(
                                   reference = NULL
                                 ),
                                 reject_option = F,
                                 k = 20,
                                 sampling = 0.2,
                                 narrow_marker_table = T,
                                 verbose = T,
                                 seed = 1){

  set.seed(seed)

  if(!is.null(exp_matrix)){
    old_names <- rownames(exp_matrix)
    check <- check_marker_names(rownames(exp_matrix))
    rownames(exp_matrix) <- check$marker_names

    if(length(check$modified) > 0){
      if(verbose){
        message("GateMeClass annotate - Renaming markers of the dataset: ",
                paste0(old_names[check$modified], collapse = " ", sep = ""), " --> ",
                paste0(rownames(exp_matrix)[check$modified], collapse = " "), sep = "")
      }
    }
  }else{
    stop("The expression matrix is mandatory!")
  }

  if(!is.null(marker_table)){
    temp <- sapply(marker_table$Gate, function(x){
      str <- strsplit(x, "[m]|m[\\|\\^]|[-]|-[\\|\\^]|[+]|\\+[\\|\\^]|[*]")
      str[[1]] <- stri_remove_empty(str[[1]])
      return(str)
    })

    if(length(setdiff(unique(unlist(temp)), rownames(exp_matrix))) > 0){
      for(i in 1:length(old_names[check$modified])){
        x <- old_names[check$modified][i]
        marker_table$Gate <- gsub(x, rownames(exp_matrix)[check$modified][i], marker_table$Gate)
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

  if(!is.null(train_parameters$reference) & is.null(marker_table)){
    reference <- train_parameters$reference
    old_names <- rownames(reference)
    check <- check_marker_names(rownames(reference))
    rownames(reference) <- check$marker_names
    markers <- check$marker_names

    if(length(check$modified) > 0){
      if(verbose){
        message("GateMeClass annotate - Renaming markers of reference dataset: ",
                paste0(old_names[check$modified], collapse = " ", sep = ""), " --> ",
                paste0(rownames(reference)[check$modified], collapse = " "), sep = "")
      }
    }

    common <- intersect(rownames(exp_matrix), rownames(reference))

    if(length(common) > 0){
      reference <- reference[rownames(reference) %in% common, ]
      train_parameters$reference <- reference
      train_parameters$seed <- seed
      train_parameters$verbose <- verbose
      marker_table <- do.call("GateMeClass_train", train_parameters)
      new_gates <- parse_marker_table(marker_table, T, T)
      markers <- colnames(new_gates$marker_table)[-1]
    }else{
      stop("No common markers between actual and reference datasets!")
    }
  }else{
    if(is.null(marker_table)){
      stop("Please, specify a gate table or a reference dataset!")
    }

    if(verbose){
      message("GateMeClass annotate - Parsing marker table...")
    }
    new_gates <- parse_marker_table(marker_table, narrow_marker_table, T)
    markers <- colnames(new_gates$marker_table)[-1]
  }

  exp_matrix_2 <- exp_matrix[markers, , drop = F]

  if(verbose){
    message(paste0("GateMeClass annotate - Determining the marker signature of each cell..."))
  }

  expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  rownames(expr_markers) <- markers

  exp_matrix_2 <- exp_matrix[markers, , drop = F]
  expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  rownames(expr_markers) <- markers

  ## Obtain the signature for the cells of the dataset to be annotated
  expr_markers <- set_marker_expression(exp_matrix_2,
                                        markers,
                                        type = new_gates$bimodal,
                                        expr_markers,
                                        verbose = verbose,
                                        GMM_parameterization = GMM_parameterization,
                                        RSS = ifelse(GMM_parameterization == "E", F, T))

  expr_markers <- expr_markers[markers, ]
  expr_markers <- as.data.frame(t(expr_markers))

  cells <- cbind(rownames(expr_markers), expr_markers)
  colnames(cells)[1] <- "Cell"

  expr_markers <- apply(cells[,- 1, drop = F], 1, paste0, collapse = "")
  new_cells2 <- data.frame(Cell = as.character(s), Gate = expr_markers)

  if(verbose){
    message("GateMeClass annotate - Cell annotation...")
  }


  ######################## Slow part of code ##########################
  res <- cell_classification(new_cells2, new_gates$extended_marker_table)
  #####################################################################


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
  if(uncl_prec > 0 & uncl_prec < ncol(exp_matrix_pre_sampling)){
    if(verbose){
      message(paste0("GateMeClass annotate - Refinement of the labels ","[k = ", k, ", Reject option = ", as.character(reject_option),"]..."))
    }

    if(!reject_option){
      training_set <- data.frame(labels = c(res$labels[not_uncl]), t(exp_matrix_pre_sampling[, not_uncl]))
      control <- t(exp_matrix_pre_sampling[, tot_uncl])
      ctrl <- trainControl(method="none")
      knnFit <- train(labels ~ ., data = training_set, method = "knn", trControl = ctrl, tuneGrid = expand.grid(k = k))
      knn_res <- predict(knnFit, newdata = control)
      res$labels[tot_uncl] <- as.character(knn_res)
      res$cell_signatures$Celltype <- res$labels
      res$cell_signatures$Cell[tot_uncl] <- tot_uncl
      res$cell_signatures$Gate[tot_uncl] <- NA
    }else{

      if(length(real_uncl) > 0){
        ## High confident "Unclassified" are cells with not MNN in the set of annotated cells
        d1 <- t(exp_matrix_pre_sampling[, real_uncl])
        d2 <- t(exp_matrix_pre_sampling[, not_uncl])
        out <- findMutualNN(d1, d2, k1=k)
        not_mnn <- setdiff(1:nrow(d1), unique(out$first))

        ## Annotation errors are, probably, annotated cells with many MNN in the set of "Unclassified" cells
        d1 <- t(exp_matrix_pre_sampling[, not_uncl])
        tot <- c(not_uncl, real_uncl)
        d2 <- t(exp_matrix_pre_sampling[, tot])
        out <- findMutualNN(d1, d2, k1=k)

        label_mapped <- res$labels[tot[out$second]]
        names(label_mapped) <- out$first

        t <- data.table(table(names(label_mapped), label_mapped))
        colnames(t) <- c("Index", "Cell", "N")
        t$Index <- as.factor(t$Index)
        t$N <- as.numeric(t$N)

        sel_max <- t %>% group_by(Index) %>% slice_max(N)
        sel_max$Index <- as.numeric(as.character(sel_max$Index))
        dup <- sel_max[which(duplicated(sel_max$Index)), "Index"]$Index
        sel_max <- sel_max[!sel_max$Index %in% dup, ]
        w_uncl <- c(sel_max[sel_max$Cell == "Unclassified", "Index"]$Index)
        to_reannotate <- not_uncl[w_uncl]
        res$labels[to_reannotate] <- "Unclassified"

        tot_uncl <- c(setdiff(tot_uncl, real_uncl[not_mnn]), to_reannotate)
        not_uncl <- setdiff(not_uncl, to_reannotate)

        training_set <- data.frame(labels = c(res$labels[not_uncl], rep("Unclassified", length(not_mnn))), t(exp_matrix_pre_sampling[, c(not_uncl, real_uncl[not_mnn])]))
        control <- t(exp_matrix_pre_sampling[, tot_uncl])
      }else{
        training_set <- data.frame(labels = res$labels[not_uncl], t(exp_matrix_pre_sampling[, not_uncl]))
        control <- t(exp_matrix_pre_sampling[, not_real_uncl])
      }

      ctrl <- trainControl(method="none")
      knnFit <- train(labels ~ ., data = training_set, method = "knn", trControl = ctrl, tuneGrid = expand.grid(k = k))
      knn_res <- predict(knnFit, newdata = control)
      res$labels[tot_uncl] <- as.character(knn_res)

      res$cell_signatures$Celltype <- res$labels
      res$cell_signatures$Cell[not_real_uncl] <- not_real_uncl
      res$cell_signatures$Gate[tot_uncl] <- NA
    }
  }

  signatures <- res$cell_signatures

  gate_ext <- sapply(signatures$Gate, function(v){
    if(!is.na(v)){
      split <- unlist(str_split(v, ""))
      sig <- paste0(markers, split, collapse = "", sep = "")
      return(sig)
    }else{
      return(NA)
    }
  })

  res$cell_signatures[, "Gate"] <- gate_ext
  res <- list(labels = res$labels,
              marker_table = marker_table,
              cell_signatures = res$cell_signatures)

  return(res)
}
