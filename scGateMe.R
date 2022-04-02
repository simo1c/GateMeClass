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

parse_gate_table <- function(gate_table, narrow_gate_table){
  
  #gate_table <- gates
  
  if(any(duplicated(gate_table$Cell))){
    stop("Error! The gate table must contains uniquely defined cell types!")
    #return(NULL)
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
    message("Error! Cell names cannot contains special characters (e.g., *, ^)!")
    return(NULL)
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
        # message("Warning! Two or more cell types have the same gating strategy, random choice!")
        # min <- sample(1:length(mins), 1)
        # min <- mins[1]
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

# options(warn=1)

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
  
  # # Bimodal and not bimodal markers
  # for(m in markers){
  #   X <- data.matrix(exp_matrix[m, ])
  #   if(is.multimodal(X, min.size=0.1)){
  #     bimodal_markers <- c(bimodal_markers, m)
  #   }else{
  #     not_bimodal_markers <- c(not_bimodal_markers, m)
  #   }
  # }
  # 
  # print(bimodal_markers)
  
  if(length(markers) > 0){
    for(m in markers){
      
      #m <- "HLA-DR"
      
      X <- exp_matrix[m, ]
      
      ## Test for multimodality
      # x <- tryCatch({
      p <- suppressWarnings(dip.test(X))
      # }, warning = function(e){
      # })
      
      if(p$p.value < 0.05){
        
        bimodal_markers <- c(bimodal_markers, m)
        
        if(verbose){
          print(paste0("Evaluate bimodal markers ", paste0(bimodal_markers, collapse = " "), collapse = " "))
        }
        
        x <- tryCatch({
          
          #if(method == "mixtools"){
          sink(tempfile())
          out <-  normalmixEM(X, k = 2, fast = T, maxit = 5000)
          closeAllConnections()
          cluster <- apply(out$posterior, 1, function(row){ which.max(row) })
          
          mu <- order(out$mu, decreasing = T)
          temp <- cluster
          temp[temp == mu[2]] <- "-"
          temp[temp == mu[1]] <- "+"
          r[m,] <- temp
          
          #x <- tryCatch({
          p <- suppressWarnings(dip.test(X[temp == "-"]))
          # }, warning = function(e){
          # })
          
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
        
        #cc <- 2
        
        w <- which(apply(r[bimodal_markers, , drop = F], 2, function(c){ names(c)<- NULL;identical(c, comb_markers[, cc])}))
        
        
        
        
        
        #comb_markers[, cc]
        if(verbose){
          print(paste0("Evaluate unimodal markers ", paste0(not_bimodal_markers, collapse = " "), " in: ", paste0(rownames(comb_markers), comb_markers[, cc], collapse = "")))
        }
        # expr_markers2 <- cbind(expr_markers2, set_marker_expression(exp_matrix[not_bimodal_markers, w, drop=FALSE],
        #                                                             not_bimodal_markers,
        #                                                             expr_markers[not_bimodal_markers, w, drop=FALSE], gates))
        
        
        
        
        # expr_markers2 <- cbind(expr_markers2, set_marker_expression(exp_matrix[not_bimodal_markers, w, drop=FALSE],
        #                                                             not_bimodal_markers,
        #                                                             expr_markers[not_bimodal_markers, w, drop=FALSE]))
        
        temp <- set_marker_expression(exp_matrix[not_bimodal_markers, w, drop=FALSE],
                                      not_bimodal_markers,
                                      expr_markers[not_bimodal_markers, w, drop=FALSE],
                                      verbose = verbose)
        
        
        # print(dim(expr_markers2))
        # print(dim(temp))
        # 
        # print(rownames(expr_markers2))
        # print(rownames(temp))
        
        if(!identical(rownames(expr_markers2), rownames(temp))){
          print("*******************ERROR**********************")
          print(rownames(expr_markers2))
          print(rownames(temp))
          stop()
        }
        
        expr_markers2 <- cbind(expr_markers2, temp)
        
        
        
        # if(ncol(expr_markers2) > 1 & nrow(expr_markers2) > 1){
        #   print(expr_markers2[1:2, 1:2])
        # }
        # 
        # print(set_marker_expression(exp_matrix[not_bimodal_markers, w, drop=FALSE],
        #                             not_bimodal_markers,
        #                             expr_markers[not_bimodal_markers, w, drop=FALSE], gates)[1:2, 1:2])
        
        # 
        # expr_markers2 <- merge(expr_markers2, set_marker_expression(exp_matrix[not_bimodal_markers, w, drop=FALSE],
        #                                                             not_bimodal_markers,
        #                                                             expr_markers[not_bimodal_markers, w, drop=FALSE], gates), by='row.names', all=TRUE)
        
        # rownames(expr_markers2) <- expr_markers2$Row.names
        #expr_markers2 <- expr_markers2[, -1]
        
        
      }
      
      
      # print(dim(expr_markers[bimodal_markers, , drop = F]))
      # print(dim(expr_markers2))
      
      temp2 <- rbind(expr_markers[bimodal_markers, ], expr_markers2[, -1])
      temp2 <- temp2[markers, ]
      return(temp2)
      
      
    }else if(length(bimodal_markers) > 0 & length(not_bimodal_markers) == 0){
      expr_markers[bimodal_markers, ] <- r[bimodal_markers, ]
      
      #print(dim(expr_markers))
      
      return(expr_markers)
    }else{
      expr_markers[markers, ] <- "*"
      
      
      #print(dim(expr_markers))
      
      
      return(expr_markers)
    }
  }else{
    
    #print(expr_markers[1:2, 1:2])
    
    
    #expr_markers[markers, ] <- "*"
    return(expr_markers)
  }
}

# expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix), nrow = length(markers)))
# rownames(expr_markers) <- markers
# test <- set_marker_expression(exp_matrix_2, colnames(gates)[-1], expr_markers)
# test[, 1:10]

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

locate_peaks <- function(d, dec = F){
  modes <- Modes(d)$modes
  if(length(modes) > 1){
    modes <- modes[order(modes, decreasing = dec)]
  }
  return(modes)
}

cutoff_between_peaks <- function(d, peaks){
  # 
  # d <- density(exp_matrix[mm, ])
  # peaks <- locate_peaks(exp_matrix[mm, ])[1:2]
  
  cut <- optimize(approxfun(d$x,d$y), interval=c(peaks[1], peaks[2]))$minimum
  
  return(cut)
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

check_overlap_intervals <- function(a, b){
  
  # a <- a[a > summary(a)["Min."] & a < summary(a)["Max."]]
  # b <- b[b > summary(b)["Min."] & b < summary(b)["Max."]]
  
  a_min <- min(a)
  a_max <- max(a)
  b_min <- min(b)
  b_max <- max(b)
  
  if (a_min > b_max || a_max < b_min) {
    return(F)
  }
  else {
    return(T)
  }
} 

find_dist_intersection <- function(a, b){
  
  
  #a <- c(0, 1,2,2,3,4,5,6)
  #b <- c(0, 25,26,27,28)
  
  
  c1 <- "a"
  c2 <- "b"
  
  if(check_overlap_intervals(a, b)){
    df <- data.frame(value = c(a,b), cluster = c(rep(c1, length(a)), rep(c2, length(b))))
    c1.density <- density(subset(df, cluster == c1)$value, from = min(df$value), to = max(df$value), n = 2^10, kernel = "gaussian")
    c2.density <- density(subset(df, cluster == c2)$value, from = min(df$value), to = max(df$value), n = 2^10, kernel = "gaussian")
    intersection.point <- c1.density$x[which(diff((c1.density$y - c2.density$y) > 0) != 0) + 1]
    intersection.point <- sort(intersection.point)
    
    if(length(intersection.point) > 1){
      med1 <- median(a)
      med2 <- median(b)
      min_mean_distance <- which.min(abs(intersection.point - med1) + abs(intersection.point - med2) / 2)
      x <- intersection.point[min_mean_distance]
    }else{
      x <- intersection.point
    }
  }else{
    return(NULL)
  }
  
  return(x)
  
}

select_thresholds_from_clusters <- function(exp_matrix,
                                            markers, 
                                            clusters, 
                                            min.size=0.1,
                                            prop_cluster_min = 0.01,
                                            iter = 0){
  
  
  # exp_matrix_3, colnames(gates)[-1], clusters2, min.size=0.1, prop_cluster_min = prop_cluster_min, count  # 
  
  # prop_cluster_min <- 0.01
  # clusters <- clusters2
  # markers <- colnames(gates)[-1]
  # exp_matrix <- exp_matrix_3
  
  # useful_clusters <- unique(res$labels_post_clustering)
  
  
  useful_clusters <- table(as.character(clusters))
  perc <- useful_clusters / sum(useful_clusters)
  useful_clusters <- names(useful_clusters)[useful_clusters > 1 & perc > prop_cluster_min]
  
  # if(return_gate_table){
  #   gate_table <- data.frame(matrix(nrow = length(unique(useful_clusters)), ncol = length(markers) + 1))
  #   colnames(gate_table) <- c("Cell", markers)
  #   gate_table$Cell <- unique(useful_clusters)
  #   rownames(gate_table) <- gate_table$Cell
  # }
  
  if(length(useful_clusters) == 0){
    message("Error! Cannot select thresholds! Try to relax the parameter 'prop_cluster_min'!")
    # thr <- rep(-1, length(markers))
    # names(thr) <- markers
    return(NULL)
  }else{
    thr <- c()
    for(mm in markers){
      
      #mm <- "CD19"
      
      ## If the distribution of the marker is multimodal
      if(length(locate_peaks(exp_matrix[mm, ])) > 1){
        th <- cutoff_between_peaks(density(exp_matrix[mm, ], kernel = "gaussian"), locate_peaks(exp_matrix[mm, ])[1:2])
        names(th) <- mm
        thr <- c(thr, th)
        ## If the distribution of the marker is unimodal
      }else{
        
        peaks_min <- c()
        peaks_max <- c()
        
        ## For each cluster
        for(c in useful_clusters){
          
          d <- exp_matrix[mm, as.character(clusters) %in% c]
          p_min <- locate_peaks(d, dec = F)
          
          if(length(p_min) == 0){
            p_min <- max(d)
          }
          
          p_min <- list(p_min)
          names(p_min) <- c
          peaks_min <- c(peaks_min, p_min)
          p_max <- locate_peaks(d, dec = T)
          
          if(length(p_max) == 0){
            p_max <- max(d)
          }
          
          p_max <- list(p_max)
          names(p_max) <- c
          peaks_max <- c(peaks_max, p_max)
        }
        
        peaks_max <- sort(sapply(peaks_max, `[[`, 1))
        peaks_min <- sort(sapply(peaks_min, `[[`, 1))
        
        #### By default take into consideration clusters > 1%
        # perc <- table(clusters) / sum(table(clusters))
        # to_remove <- names(perc)[which(perc < prop_cluster_min)]
        # 
        # peaks_min <- peaks_min[!names(peaks_min) %in% to_remove]
        # peaks_max <- peaks_max[!names(peaks_max) %in% to_remove]
        
        inter <- c()
        diff_peaks <- c()
        
        if(length(peaks_max) > 1){
          
          
          
          peaks <- locate_peaks(exp_matrix[mm, clusters == names(peaks_max[1]) | clusters == names(peaks_max[length(peaks_max)])])
          
          if(length(peaks) > 1){
            th <- cutoff_between_peaks(density(exp_matrix[mm, clusters == names(peaks_max[1]) |
                                                            clusters == names(peaks_max[length(peaks_max)])], 
                                               kernel = "gaussian"), peaks)[1:2]
            names(th) <- mm
          }else{
            
            th <- peaks
            names(th) <- mm
            
          }
          
          
          
          
        }else{
          th <- peaks_max
          names(th) <- mm
          inter <- c(inter, th)
        }
        thr <- c(thr, th)
      }
    }
  }
  
  
  perc <- table(clusters) / sum(table(clusters))
  useful_clusters <- names(perc)[perc > 0.01]
  
  # for(mm in markers){
  #   df <- data.frame(cluster = clusters[clusters %in% useful_clusters], value = exp_matrix[mm, clusters %in% useful_clusters])
  #   pdf(paste0(mm, "_thr_",iter, ".pdf"))
  #   gg <- ggplot(df, aes(value, fill = cluster), aes_string(y = "..ndensity..")) + geom_density(alpha = 0.6, kernel = "gaussian") +
  #     geom_vline(xintercept = thr[mm], color = "red") +
  #     theme_classic()
  # 
  #   print(gg)
  # 
  #   dev.off()
  # }
  
  return(thr)
}

#### DA IMPLEMENTARE:
## 1) Annotazione gerarchica. Alcuni marker devono essere controllati dopo aver fatto il 
## subsetting di specifiche sotto-popolazioni. Per fare questo, per ciascun marker di questo 
## tipo si deve specificare in quale popolazione andare a vederlo. Si può aggiungere una colonna nella
## tabella dei marcatori dove, per ogni, popolazione si specificano i marker da guardare considerando 
## solo quel subset.
## 2) Possibilità di specificare che l'espressione di un marker deve essere la più alta di tutti ("Hi")
## 3) Provare a vedere se possibile mettere una thr dinamica per ogni cluster, es. 

scGateMe <- function(exp_matrix,
                     gates, 
                     ignore_markers = NULL,
                     refine = T,
                     narrow_gate_table = T,
                     verbose = F,
                     seed = 1){
  
  gates
  refine = F
  n_clusters = 20
  cluster_level_labels = T
  prop_cluster_min = 0.01
  plots_thr = F
  seed = 1
  prop_cluster_min = 0.01
  clusters = NULL
  # clusters <- res$labels_post_clustering
  exp_matrix <- m
  ignore_markers = NULL
  verbose = T
  narrow_gate_table = T
  
  set.seed(seed)
  
  message("Loading gate table...")
  new_gates <- parse_gate_table(gates, narrow_gate_table)
  message("Prioritizing gate table...")
  new_gates2 <- prioritize_gate_table(new_gates)
  new_gates2$extended_gate_table <- rbind(new_gates2$extended_gate_table, data.frame(Cell = "Unclassified", Gate = paste0(colnames(new_gates2$gate_table)[-1], "*", collapse = "")))
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
  
  # if(is.null(clusters)){
  #   message("Clustering...")
  #   sc <- SingleCellExperiment(assays=list(exp_matrix = exp_matrix, exprs = exp_matrix))
  #   sc <- CATALYST::cluster(sc, features = NULL, verbose = F, maxK = n_clusters)
  #   meta_cl <- paste0("meta", n_clusters)
  #   clusters <- cluster_ids(sc, meta_cl)
  # }else{
  #   clusters <- factor(clusters)
  # }
  
  
  # ### Chose number of clusters ##
  # purities <- rep(NA, 18)
  # names(purities) <- as.character(3:20)
  # message("Calculating the number of clusters...")
  # for(n in 3:3){
  #   message(paste(n-2, " "), appendLF = F)
  #   meta_cl <- paste0("meta", n)
  #   clusters <- cluster_ids(sc, meta_cl)
  #   antimodes <- select_thresholds_from_clusters(exp_matrix_2, colnames(gates)[-1], clusters = clusters, min.size=0.1, prop_cluster_min = prop_cluster_min)
  #   new_cells <- set_marker_expression(exp_matrix_2, colnames(gates)[-1], antimodes)
  #   res <- cell_classification(new_cells, new_gates$extended_gate_table)
  #   p <- sum(apply(table(res$labels, clusters), 2, max)) / length(clusters)
  #   purities[as.character(n)] <- p
  # }
  # message("")
  # 
  # meta_cl <- paste0("meta", names(purities)[which(purities == max(purities))])
  # clusters <- cluster_ids(sc, meta_cl)
  
  
  #message("Choosing thresholds for each marker...")
  #antimodes <- select_thresholds_from_clusters(exp_matrix_2, colnames(gates)[-1], clusters = clusters, min.size=0.1, prop_cluster_min = prop_cluster_min)
  #antimodes <- -1
  
  #if(length(antimodes) > 0){
    
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
    
    # count <- 1
    # while(uncl != uncl_prec & any(res$labels == "Unclassified")){
    #   message(paste(count, " "), appendLF = F)
    #   uncl_prec <- sum(res$labels == "Unclassified")
    #   exp_matrix_3 <- matrix(exp_matrix_2[, res$labels == "Unclassified"], nrow(exp_matrix_2), uncl_prec)
    #   rownames(exp_matrix_3) <- rownames(exp_matrix_2)
    #   colnames(exp_matrix_3) <- colnames(exp_matrix_2)[res$labels == "Unclassified"]
    #   
    #   clusters2 <- clusters[res$labels == "Unclassified"]
    #   
    #   # sc2 <- CATALYST::cluster(sc[, res$labels == "Unclassified"], features = NULL, verbose = F, maxK = n_clusters)
    #   # meta_cl <- paste0("meta", n_clusters)
    #   # clusters2 <- cluster_ids(sc2, meta_cl)
    #   
    #   useful_clusters <- table(as.character(clusters2))
    #   perc <- useful_clusters / sum(useful_clusters)
    #   useful_clusters <- names(useful_clusters)[useful_clusters > 1 & perc > prop_cluster_min]
    #   
    #   # if(length(useful_clusters) >= 1){
    #   #   antimodes <- select_thresholds_from_clusters(exp_matrix_3, colnames(gates)[-1], clusters2, min.size=0.1, prop_cluster_min = prop_cluster_min, count)
    #   # }
    #   
    #   antimodes <- -1
    #   
    #   expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_3), nrow = nrow(exp_matrix_3)))
    #   rownames(expr_markers) <- colnames(gates)[-1]
    #   
    #   expr_markers_2 <- set_marker_expression(exp_matrix_3, colnames(gates)[-1], expr_markers, gates, method = method)
    #   expr_markers_2 <- expr_markers_2[colnames(gates)[-1], ]
    #   expr_markers_2 <- as.data.frame(t(expr_markers_2))
    #   cells <- cbind(rownames(expr_markers_2), expr_markers_2)
    #   colnames(cells)[1] <- "Cell"
    #   expr_markers_2 <- apply(cells, 1, function(r){ gate <- paste(colnames(cells)[-1], r[-1], sep = "", collapse = ""); return(gate)})
    #   new_cells <- data.frame(Cell = names(expr_markers_2), Gate = expr_markers_2)
    #   res_2 <- cell_classification(new_cells, new_gates$extended_gate_table)
    #   
    #   res$labels[res$labels == "Unclassified"] <- res_2$labels
    #   
    #   uncl <- length(res$labels[res$labels == "Unclassified"])
    #   count <- count + 1
    # }
    # if(count > 1){
    #   message("")
    # }
  }
  # if(cluster_level_labels){
  #   
  #   message("Refinement of the labels using clustering...")
  #   # tb <- prop.table(table(res$labels, clusters), 2)
  #   # top_lab_index <- unlist(apply(tb, 2, function(c){
  #   #   return(which.max(c))
  #   # }))
  #   # 
  #   # top_lab_cluster <- rep(NA, length(levels(clusters)))
  #   # top_lab_cluster <- rownames(tb)[top_lab_index]
  #   # labels_cl <- plyr::mapvalues(clusters, levels(clusters), top_lab_cluster)
  #   
  #   train <- t(exp_matrix[, !res$labels %in% c("Unclassified")])
  #   control <- t(exp_matrix[, res$labels %in% c("Unclassified")])
  #   
  #   t <- table(res$labels)
  #   tt <- t[which.min(t)]
  #   k <- floor(sqrt(tt))
  #   
  #   knn_res <- knn(train,
  #            control,
  #            cl = factor(res$labels[!res$labels %in% c("Unclassified")]),
  #            k = k,
  #            prob = T)
  #   
  #   labels_cl <- res$labels
  #   labels_cl[labels_cl == "Unclassified"] <- as.character(knn_res)
  #   
  #   res <- list(labels = res$labels,
  #               labels_post_clustering = labels_cl,
  #               marker_table = res$marker_table,
  #               thresholds = antimodes)
  
  #}else{
  res <- list(labels = res$labels,
              marker_table = res$marker_table,
              thresholds = antimodes)
  #}
#}
  
  return(res)
}
