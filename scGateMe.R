# options(warn=1)

## This function generates the possible value combinations of the markers of a cell type 
generate_set_values <- function(v, cell){
  
  sets_all <- list()
  sets_to_filter <- list()
 
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
           },
           `-` = {
             set <- list(c("-"))
             set_not <- list()
           },
           # `+|`={
           #   set <- list(c("+", "-"))
           #   set_not <- list(c("-"))
           # },
           # `-|`={
           #   set <- list(c("*","+", "-"))
           #   set_not <- list(c("+"))
           # },
           # `+^`={
           #   set <- list(c("*","+", "-"))   #### to check ######
           #   # set_not <- list(c("int", "-"))
           #   set_not <- list()
           # },
           # `-^`={
           #   set <- list(c("*","+", "-"))   #### to check ######
           #   set_not <- list(c("-"))
           # },
           `*`={
             # set <- list(c("*"))
             set <- list(c("*", "+", "-"))
             set_not <- list()
           }
    )
    sets_all <- c(sets_all, set)
    sets_to_filter <- c(sets_to_filter, set_not)
  }
  
  names(sets_all) <- names(v)
  sets_all <- do.call(expand.grid, sets_all)
  
  if(length(sets_to_filter) > 0){
    special_to_filter <- grep("[\\|]", v)
    names(sets_to_filter) <- names(v)[special_to_filter]
    sets_to_filter <- do.call(expand.grid, sets_to_filter)
  }
  
  if(length(sets_to_filter) > 0){
    sets_all_filtered <- anti_join(sets_all, sets_to_filter, by = names(sets_to_filter))
  }else{
    sets_all_filtered <- sets_all
  }
  
  sets_all_filtered$Cell <- cell
  
  gate <- apply(sets_all_filtered, 1, function(r){ paste0(colnames(sets_all_filtered)[-which(colnames(sets_all_filtered) == "Cell")], 
                                                          r[-which(colnames(sets_all_filtered) == "Cell")], collapse = "") })
  
  sets_all_filtered$Gate <- gate
  sets_all_filtered <- sets_all_filtered[, c("Cell", "Gate")]
  return(sets_all_filtered)
}

## Read the gate table and generate the possbile marker signature of ceach cell type
parse_gate_table <- function(gate_table, narrow_gate_table){
  
  #gate_table <- gates
  
  if(any(duplicated(gate_table$Cell))){
    stop("Error! The gate table must contains uniquely defined cell types!")
  }
  
  if(narrow_gate_table){
    gates2 <- gate_table$Gate
    temp <- sapply(gates2, function(x){ str <- strsplit(x, "[+]|[-]|[*]"); return(str)})
    names(temp) <- gates$Cell
    markers <- unique(unlist(temp))
    df_gates <- data.frame(matrix(nrow = nrow(gates), ncol = length(markers)))
    rownames(df_gates) <- gates$Cell
    colnames(df_gates) <- markers
    gates2_exploded <- strsplit(gates2, split = "")
    signs <- lapply(gates2_exploded, function(x){ x[x %in% c("+", "-", "*")]})
    names(signs) <- names(temp)
    
    for(j in 1:length(temp)){
      el <- temp[j]
      for(i in 1:length(unlist(el))){
        df_gates[names(el), unlist(el)[i]] <- unlist(signs[names(el)])[i]
      }
    }
    
    df_gates[is.na(df_gates)] <- "*"
    df_gates <- cbind(Cell = rownames(df_gates), df_gates)
    gate_table <- df_gates
  }

  celltypes <- gate_table$Cell
  
  if(length(grep("[\\*|\\||\\^]", celltypes)) > 0){
    stop("Error! Cell names cannot contains special characters (e.g., *, ^)!")
  }
  
  extended_gate_table <- data.frame(Cell = NA, Gate = NA)
      
  for(i in 1:nrow(gate_table)){
        
    v <- gate_table[i, -1, drop = F]
    cell <- gate_table[i, "Cell"]
    to_add <- generate_set_values(v, cell)
    
    if(is.null(to_add)){
      return(NULL)
    }
    extended_gate_table <- rbind(extended_gate_table, to_add)
  }
  
  extended_gate_table <- extended_gate_table[-1, ]
  to_delete <- paste0(colnames(gate_table[, -1]), "*", collapse = "")
  extended_gate_table <- extended_gate_table[extended_gate_table$Gate != to_delete, ]
  
  return(list(gate_table = gate_table, extended_gate_table = extended_gate_table))
}

## This function parse the marker signatures for each cell types and prioritize the choice of 
## cell type in the case of different cell types have the same marker signature 
prioritize_gate_table <- function(new_gates){
  
  # gates
  # new_gates
  # gates = new_gates$gate_table
  # new_gates
  
  gates <- new_gates$gate_table
  extended_gate_table <- new_gates$extended_gate_table
  star <- apply(gates, 1, function(r){return(sum(r == "*"))})
  star_check <- table(star)
  
  star_df <- data.frame(Cell = gates$Cell, Star = star)
  gate_dup <- extended_gate_table[duplicated(extended_gate_table$Gate), "Gate"]
  
  if(length(gate_dup) > 0){
    
    gate_not_dup <- extended_gate_table[!extended_gate_table$Gate %in% gate_dup, ]
    extended_gate_table_dup <- extended_gate_table[extended_gate_table$Gate %in% gate_dup, ]
    extended_gate_table_dup <- merge(extended_gate_table_dup, star_df, by = "Cell")
    extended_gate_table_dup$id <- 1:nrow(extended_gate_table_dup)
        
    ############ TO OPTIMIZE!!!!!! ##########
    to_keep <- c()
    for(g in unique(extended_gate_table_dup$Gate)){
            
      #g <- extended_gate_table_dup$Gate[1]
      
      temp <- extended_gate_table_dup[extended_gate_table_dup$Gate == g, ]
      temp <- temp[order(temp$Star, decreasing = F), ]
      mins <- which(temp$Star == min(temp$Star))
      
      if(length(mins) > 1){
        to_keep <- c(to_keep, temp[mins[1], "id"])
        extended_gate_table_dup[extended_gate_table_dup$Gate == temp$Gate[1] & extended_gate_table_dup$Cell %in% temp$Cell, "Cell"] <- paste0(temp$Cell, collapse = " | ")
      }else{
        to_keep <- c(to_keep, temp[mins, "id"])
      }
    }
    
    extended_gate_table_dup <- extended_gate_table_dup[to_keep, ]
    extended_gate_prior <- rbind(extended_gate_table_dup[, -which(colnames(extended_gate_table_dup) %in% c("Star", "id"))], gate_not_dup)
  }else{
    return(list(gate_table = new_gates$gate_table, extended_gate_table = new_gates$extended_gate_table))
  }
  
  return(list(gate_table = new_gates$gate_table, extended_gate_table = extended_gate_prior))
}


## This function set tha marker signature of each cell
set_marker_expression <- function(exp_matrix, markers, expr_markers, gates, verbose){
  
  # exp_matrix <- exp_matrix[not_bimodal_markers, w, drop=FALSE]
  # markers <- not_bimodal_markers
  # expr_markers <- expr_markers[not_bimodal_markers, w, drop=FALSE]
  # exp_matrix <- exp_matrix_2
  # markers <- colnames(gates)[-1]
  
  r <- data.frame(matrix(nrow = length(markers), ncol = ncol(exp_matrix)))
  rownames(r) <- markers
  
  bimodal_markers <- c()
  not_bimodal_markers <- c()
  
  if(length(markers) > 0){
    for(m in markers){
      
      #m <- "HLA-DR"
      
      X <- exp_matrix[m, ]
      
      ## Test for multimodality
      p <- suppressWarnings(dip.test(X))

      if(p$p.value < 0.05){
        
        bimodal_markers <- c(bimodal_markers, m)
        
        if(verbose){
          print(paste0("Evaluate bimodal markers ", paste0(bimodal_markers, collapse = " "), collapse = " "))
        }
        
        x <- tryCatch({
          sink(tempfile())
          out <-  normalmixEM(X, k = 2, fast = T, maxit = 5000)
          closeAllConnections()
          cluster <- apply(out$posterior, 1, function(row){ which.max(row) })
          
          mu <- order(out$mu, decreasing = T)
          temp <- cluster
          temp[temp == mu[2]] <- "-"
          temp[temp == mu[1]] <- "+"
          r[m,] <- temp
          
          p <- suppressWarnings(dip.test(X[temp == "-"]))
          
          if(p$p.value < 0.05){
            sink(tempfile())
            out <- normalmixEM(X, k = 3, fast = F, maxit = 5000)
            closeAllConnections()
            #plot(out, which = 2)
            cluster <- apply(out$posterior, 1, function(row){ which.max(row) })
            mu <- order(out$mu, decreasing = T)
            temp <- cluster
            temp[temp == mu[3]] <- "-"
            temp[temp == mu[2]] <- "+"
            temp[temp == mu[1]] <- "+"
            r[m,] <- temp
          }
        },
        error = function(e){
          expr_markers[markers, ] <- "*"
        },
        warning = function(e){
          print(e)})
      }else{
        #print(m)
        not_bimodal_markers <- c(not_bimodal_markers, m)
      }
    }
    
    if(length(bimodal_markers) > 0 & length(not_bimodal_markers) > 0){
      
      #print(bimodal_markers)
      
      expr_markers[bimodal_markers, ] <- r[bimodal_markers, ]
      
      r_temp <- r[bimodal_markers, ]
      comb_markers <-  r_temp[!duplicated(as.list(r_temp))]
      
      expr_markers2 <- data.frame(as.matrix(rep(NA, length(not_bimodal_markers))))
      colnames(expr_markers2) <- "Test"        
      rownames(expr_markers2) <- not_bimodal_markers
      
      
      
      for(cc in 1:ncol(comb_markers)){
        w <- which(apply(r[bimodal_markers, , drop = F], 2, function(c){ names(c)<- NULL;identical(c, comb_markers[, cc])}))

        if(verbose){
          print(paste0("Evaluate unimodal markers ", paste0(not_bimodal_markers, collapse = " "), " in: ", paste0(rownames(comb_markers), comb_markers[, cc], collapse = "")))
        }
        
        temp <- set_marker_expression(exp_matrix[not_bimodal_markers, w, drop=FALSE],
                                      not_bimodal_markers,
                                      expr_markers[not_bimodal_markers, w, drop=FALSE],
                                      verbose = verbose)
        
        if(!identical(rownames(expr_markers2), rownames(temp))){
          print("Error in the combination of data frames!")
          print(rownames(expr_markers2))
          print(rownames(temp))
          stop()
        }
        
        expr_markers2 <- cbind(expr_markers2, temp)
      }
      
      temp2 <- rbind(expr_markers[bimodal_markers, ], expr_markers2[, -1])
      temp2 <- temp2[markers, ]
      return(temp2)
      
      
    }else if(length(bimodal_markers) > 0 & length(not_bimodal_markers) == 0){
      expr_markers[bimodal_markers, ] <- r[bimodal_markers, ]
      return(expr_markers)
    }else{
      expr_markers[markers, ] <- "*"
      return(expr_markers)
    }
  }else{
    return(expr_markers)
  }
}

## This function performs the cell classification
cell_classification <- function(marker_table, gates){
  
  # marker_table <- new_cells
  # gates <- new_gates$extended_gate_table
  # new_cells, new_gates$extended_gate_table, exact
  
  gates <- gates[!duplicated(gates$Gate), ]
  rownames(gates) <- gates$Gate
  labels <- gates[marker_table$Gate, "Cell"]
  labels[is.na(labels)] <- "Unclassified"
  rownames(gates) <- NULL
  cl_res <- list(labels = labels, marker_table = gates)
  return(cl_res)
}

extract_gate_table <- function(exp_matrix, markers, clusters){
  
  # clusters <- res$labels_post_clustering
  # markers <- colnames(gates)[-1]
  # exp_matrix <- m
  
  exp_matrix <- m[, clusters != "Unclassified"]
  clusters <- clusters[clusters != "Unclassified"]
  
  message("Extracting gate table from data...")
  useful_clusters <- unique(clusters)
  
  thresholds <- select_thresholds_from_clusters(exp_matrix,
                                                markers, 
                                                clusters, 
                                                prop_cluster_min = 0)
  
  gate_table <- data.frame(matrix(nrow = length(useful_clusters), ncol = length(markers) + 1))
  colnames(gate_table) <- c("Cell", markers)
  gate_table$Cell <- useful_clusters
  rownames(gate_table) <- gate_table$Cell
  
  for(mm in markers){
    
    #mm <- "CD45"
    
    mds <- c()
    q1s <- c()
    q3s <- c()
    
    for(c in useful_clusters){
      d <- exp_matrix[mm, as.character(clusters) %in% c]
      md <- median(d)
      q1 <- summary(d)["1st Qu."] 
      q3 <- summary(d)["3rd Qu."]
      names(md) <- c
      names(q1) <- c
      names(q3) <- c  
      mds <- c(mds, md)
      q1s <- c(q1s, q1)
      q3s <- c(q3s, q3)
    }
    
    pos <- which(thresholds[mm] < q1s)
    neg <- which(thresholds[mm] > q3s)
    star <- which(thresholds[mm] > q1s & thresholds[mm] < q3s)
    
    exp <- rep(NA, length(useful_clusters))
    names(exp) <- useful_clusters
    exp[names(star)] <- "*"
    exp[names(pos)] <- "pos"
    exp[names(neg)] <- "neg"
    
    gate_table[names(exp), mm] <- exp
    
  }
  return(gate_table)
}

## This is the core function for the classification of the single cells
scGateMe <- function(exp_matrix,
                     gates, 
                     ignore_markers = NULL,
                     refine = F,
                     narrow_gate_table = T,
                     verbose = T,
                     seed = 1){
  
  # gates
  # refine = F
  # n_clusters = 20
  # cluster_level_labels = T
  # prop_cluster_min = 0.01
  # plots_thr = F
  # seed = 1
  # prop_cluster_min = 0.01
  # clusters = NULL
  # # clusters <- res$labels_post_clustering
  # exp_matrix <- m
  # ignore_markers = NULL
  # verbose = T
  # narrow_gate_table = T
  
  set.seed(seed)
  message("Loading gate table...")
  new_gates <- parse_gate_table(gates, narrow_gate_table)
  message("Prioritizing gate table...")
  new_gates2 <- prioritize_gate_table(new_gates)
  new_gates2$extended_gate_table <- rbind(new_gates2$extended_gate_table, 
                                          data.frame(Cell = "Unclassified", 
                                                     Gate = paste0(colnames(new_gates2$gate_table)[-1], "*", collapse = "")))
  
  exp_matrix_2 <- exp_matrix[colnames(new_gates$gate_table)[-1], , drop = F]
  
  if(!is.null(ignore_markers)){
    n <- length(ignore_markers)
    for(i in 1:n){
      el <- ignore_markers[i]
      cell <- names(el)
      markers <- unlist(el)
      gates[gates$Cell == cell, markers] <- "*"
    }
  }
    
  message("Determining the marker signature for each cell...")
  expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  rownames(expr_markers) <- colnames(new_gates2$gate_table)[-1]
  expr_markers <- set_marker_expression(exp_matrix_2, colnames(new_gates2$gate_table)[-1], expr_markers, new_gates2$gate_table, verbose = verbose)
  expr_markers <- expr_markers[colnames(new_gates2$gate_table)[-1], ]
  expr_markers <- as.data.frame(t(expr_markers))
  cells <- cbind(rownames(expr_markers), expr_markers)
  colnames(cells)[1] <- "Cell"
  expr_markers <- apply(cells, 1, function(r){ gate <- paste(colnames(cells)[-1], r[-1], sep = "", collapse = ""); return(gate)})
  new_cells <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
  
  message("Classification of the cells...")
  res <- cell_classification(new_cells, new_gates2$extended_gate_table)
  
  res_2 <- res
  uncl <- 0
  uncl_prec <- length(res$labels[res$labels == "Unclassified"])
  
  ## Refinement of the unclassified cells using K-NN classification
  if(refine & uncl_prec > 0 & uncl_prec < ncol(exp_matrix)){
    message("Refinement of the labels using k-nn classification...")
    train <- t(exp_matrix[, !res$labels %in% c("Unclassified")])
    control <- t(exp_matrix[, res$labels %in% c("Unclassified")])
    
    t <- table(res$labels)
    tt <- t[which.min(t)]
    k <- floor(sqrt(tt))
    
    knn_res <- knn(train,
                   control,
                   cl = factor(res$labels[!res$labels %in% c("Unclassified")]),
                   k = k,
                   l = 1,
                   prob = T)
    
    res$labels[res$labels == "Unclassified"] <- as.character(knn_res)
  }
  
  res <- list(labels = res$labels,
              marker_table = res$marker_table,
              thresholds = antimodes)
  
  return(res)
}
