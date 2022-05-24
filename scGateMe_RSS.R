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
             set <- list(c("*", "+"))
             set_not <- list()
           },
           # `++` = {
           #   set <- list(c("*", "+", "++"))
           #   set_not <- list()
           # },
           `-` = {
             set <- list(c("*", "-"))
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
  
  sets_all_filtered_temp <- sets_all_filtered[, -which(colnames(sets_all_filtered) == "Cell"), drop = F]
  markers <- colnames(sets_all_filtered_temp)
  
  gate <- apply(sets_all_filtered_temp, 1, function(r){
    paste0(markers, r, collapse = "")
  })
  
  sets_all_filtered$Gate <- gate
  sets_all_filtered <- sets_all_filtered[, c("Cell", "Gate")]
  return(sets_all_filtered)
}

## Read the gate table and generate the possbile marker signature of ceach cell type
parse_gate_table <- function(gate_table, narrow_gate_table){
  
  # gate_table <- gates
  
  if(any(duplicated(gate_table$Cell))){
    stop("Error! The gate table must contains uniquely defined cell types!")
  }
  
  if(narrow_gate_table){
    
    gates2 <- gate_table
    temp <- sapply(gates2$Gate, function(x){ str <- strsplit(x, "[+]+|[-]|[*]"); return(str)})
    
    # print(temp)
    
    names(temp) <- gates2$Cell
    markers <- unique(unlist(temp))
    df_gates <- data.frame(matrix(nrow = nrow(gates2), ncol = length(markers)))
    rownames(df_gates) <- gates2$Cell
    colnames(df_gates) <- markers
    
    
    
    # gates2_exploded <- strsplit(gates2$Gate, split = "")
    gates2_exploded <- strsplit(perl=T, gates2$Gate, '(?!\\++)');
    
    #print(gates2_exploded)
    
    
    signs <- lapply(gates2_exploded, function(x){ x[x %in% c("+", "-", "*")]})
    
    # print(signs)
    
    names(signs) <- names(temp)
    
    l_temp <- 1:length(temp)
    
    for(j in l_temp){
      el <- temp[j]
      l_temp2 <- 1:length(unlist(el))
      for(i in l_temp2){
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
    
  df_list <- apply(gate_table, 1, function(v){
    to_add <- generate_set_values(v[-1], v[1])
    return(to_add)
  })

  extended_gate_table <- as.data.frame(rbindlist(df_list))
  to_delete <- paste0(colnames(gate_table[, -1]), "*", collapse = "")
  extended_gate_table <- extended_gate_table[extended_gate_table$Gate != to_delete, ]
  
  return(list(gate_table = gate_table, extended_gate_table = extended_gate_table))
}

## This function parse the marker signatures for each cell types and prioritize the choice of 
## cell type in the case of different cell types have the same marker signature 
# prioritize_gate_table <- function(new_gates){
#   
#   # gates
#   # new_gates
#   # gates = new_gates$gate_table
#   # new_gates
#   
#   gates <- new_gates$gate_table
#   extended_gate_table <- new_gates$extended_gate_table
#   star <- apply(gates, 1, function(r){return(sum(r == "*"))})
#   star_check <- table(star)
#   
#   star_df <- data.frame(Cell = gates$Cell, Star = star)
#   
#   # gate_dup <- extended_gate_table[duplicated(extended_gate_table$Gate), "Gate"]
#   
#   extended_gate_table <- unique(extended_gate_table)
#   gate_dup <- any(duplicated(extended_gate_table$Gate))
#   
#   if(gate_dup){
#     
#     gate_not_dup <- extended_gate_table[!duplicated(extended_gate_table$Gate), ]
#     extended_gate_table_dup <- extended_gate_table[duplicated(extended_gate_table$Gate), ]
#     
#     extended_gate_table_dup <- merge(extended_gate_table_dup, star_df, by = "Cell")
#     extended_gate_table_dup$id <- 1:nrow(extended_gate_table_dup)
#     
#     ############ TO OPTIMIZE!!!!!! ##########
#     to_keep <- c()
#     unique_gates <- unique(extended_gate_table_dup$Gate)
#     
#     # for(g in unique_gates){
#     #   
#     #   #g <- extended_gate_table_dup$Gate[1]
#     #   
#     #   temp <- extended_gate_table_dup[extended_gate_table_dup$Gate == g, ]
#     #   temp <- temp[order(temp$Star, decreasing = F), ]
#     #   mins <- which(temp$Star == min(temp$Star))
#     #   
#     #   if(length(mins) > 1){
#     #     to_keep <- c(to_keep, temp[mins[1], "id"])
#     #     extended_gate_table_dup[extended_gate_table_dup$Gate == temp$Gate[1] & extended_gate_table_dup$Cell %in% temp$Cell, "Cell"] <- paste0(temp$Cell, collapse = " | ")
#     #   }else{
#     #     to_keep <- c(to_keep, temp[mins, "id"])
#     #   }
#     # }
#     
#     to_keep <- sapply(unique_gates, function(g){{
#       
#       temp <- extended_gate_table_dup[extended_gate_table_dup$Gate == g, ]
#       temp <- temp[order(temp$Star, decreasing = F), ]
#       mins <- which(temp$Star == min(temp$Star))
#       
#       # if(length(mins) > 1){
#       #   to_keep <- c(to_keep, temp[mins[1], "id"])
#       #   # p <- paste0(temp$Cell, collapse = " | ")
#       #   extended_gate_table_dup[extended_gate_table_dup$Gate == temp$Gate[1] & extended_gate_table_dup$Cell %in% temp$Cell, "Cell"] <- "Unclassified"
#       # }else{
#       #   to_keep <- c(to_keep, temp[mins, "id"])
#       # }
#       
#       if(length(mins) == 1){
#         to_keep <- c(to_keep, temp[mins, "id"])
#       }
#       
#       return(to_keep)
#     }})
#     
#     to_keep <- unlist(to_keep)
#     extended_gate_table_dup <- extended_gate_table_dup[to_keep, ]
#     extended_gate_prior <- rbind(extended_gate_table_dup[, -which(colnames(extended_gate_table_dup) %in% c("Star", "id"))], gate_not_dup)
#   }else{
#     return(list(gate_table = new_gates$gate_table, extended_gate_table = new_gates$extended_gate_table))
#   }
#   
#   return(list(gate_table = new_gates$gate_table, extended_gate_table = extended_gate_prior))
# }

# set_marker_expression_GMM <- function(X, trimodality){
#   
#   
#   print(dim(X))
#   
#   # X <- m["CD34",]
#   # X <- X[sample(1:100)]
#   # plot(density(X))
#   
#   # icl <- mclustICL(X, G = 1:2, verbose = F)
#   # cl <- Mclust(X, G = 2)
#   # plot(cl, what = "classification")
#   
#   # rss(X, m = 2, sets = F, r = 500)
#   
#   s <- obsno.Mrss(X, m = 2, r = 500, type = "r")
#   index <- as.numeric(gsub("Obs.  ", "", s[, 2]))
#   
#   test <- as.numeric(X[index])
#   # test <- s$sample[, 2]
#   
#   # plot(density(test))
#   
#   # test <- as.numeric(s$sample[, 2])
#   
#   icl <- mclustICL(test, G = 1:2, verbose = F)
#   # cl <- Mclust(test, G = 2)
#   # plot(cl, what = "classification")
#   
#   # table(pred$classification, sce2$labels)
#   
#   
#   
#   model_temp <- unlist(strsplit(names(summary(icl)[1]), ","))
#   type_model <- model_temp[1]
#   model <- as.numeric(model_temp[2])
#   
#   if(!is.na(model) & model > 1){
#     
#     cl <- Mclust(test, G = model, verbose = F, modelNames = type_model)
#     
#     temp <- predict.Mclust(cl, X)
#     # temp <- cl$classification
#     
#     means <- cl$parameters$mean
#     means <- sort(means)
#     
#     if(model == 2){
#       temp[temp == names(means)[1]] <- "-"
#       temp[temp == names(means)[2]] <- "+"
#     }else{
#       
#       temp[temp == names(means)[1]] <- "-"
#       
#       if(any(names(trimodality) == "++")){
#         temp[temp == names(means)[2]] <- "+"
#         temp[temp == names(means)[3]] <- "++"
#       }else{
#         temp[temp == names(means)[2]] <- "+"
#         temp[temp == names(means)[3]] <- "+"
#       }
#  
#     }
#     return(temp)
#   }else{
#     temp <- rep("*", length(X)) 
#     return(temp)
#   }
# }

set_marker_expression_GMM <- function(X, trimodality, m, r){
  
  
  # X <- m["CD45", ]
  # plot(density(X))

  s <- obsno.Mrss(X, m = 2, r = 1000, type = "r")
  
  if(skewness(X) < 0){
    index <- as.numeric(gsub("Obs.  ", "", s[, 1]))
  }else{
    index <- as.numeric(gsub("Obs.  ", "", s[, 2]))
  }
  
  test <- as.numeric(X[index])
  
  # plot(density(test))
  # cl <- Mclust(test, G = 2, modelNames = "E")
  # plot(cl, what = "classification")
  
  icl <- mclustICL(test, G = 1:2, verbose = F)
  model_temp <- unlist(strsplit(names(summary(icl)[1]), ","))
  type_model <- model_temp[1] 
  model <- as.numeric(model_temp[2])
  
  if(!is.na(model) & model > 1){
    
    cl <- Mclust(test, G = model, verbose = F, modelNames = type_model)
    #temp <- cl$classification
    temp <- predict.Mclust(cl, X)$classification
    means <- cl$parameters$mean
    means <- sort(means)
    
    if(model == 2){
      temp[temp == names(means)[1]] <- "-"
      temp[temp == names(means)[2]] <- "+"
    }else{
      temp[temp == names(means)[1]] <- "-"
      
      if(trimodality == "strictly_pos"){
        temp[temp == names(means)[2]] <- "-"
      }else{
        temp[temp == names(means)[2]] <- "+"
      }
      
      temp[temp == names(means)[3]] <- "+"
    }
    return(temp)
  }else{
    temp <- rep("*", length(X)) 
    return(temp)
  }
}

## This function set the marker signature of each cell
set_marker_expression <- function(exp_matrix, markers, 
                                  expr_markers, 
                                  gates, 
                                  verbose, 
                                  marker_seq_eval, 
                                  trimodality,
                                  m,
                                  r){


 #  exp_matrix <- exp_matrix_2
 # markers <- colnames(new_gates$gate_table)[-1]
 #  # # expr_markers
 #  gates <- new_gates$gate_table
 #  # verbose = F
 #  marker_seq_eval = F

  ##################################
  ## Sequential marker evaluation ##
  ##################################
  # if(marker_seq_eval){
  #   
  #   queue <- list(list(indexes = 1:ncol(exp_matrix), markers = markers))
  #   
  #   while(length(queue) > 0){
  #     
  #     bimodal_markers <- c()
  #     not_bimodal_markers <- c()
  #     
  #     ## Pop operation
  #     first <- queue[[1]]
  #     queue <- queue[-1]
  #     
  #     for(m in first$markers){
  #       
  #       X <- exp_matrix[m, first$indexes]
  #       
  #       # if(length(X) < 1){
  #       #   next
  #       # }
  #       t <- table(gates[, m])
  #       
  #       marker_expr <- set_marker_expression_GMM(X, t)
  #       
  #       if(length(table(marker_expr)) > 1){
  #         bimodal_markers <- c(bimodal_markers, m)
  #         other_markers <- first$markers[first$markers != bimodal_markers]
  #         break
  #       }
  #     }
  #     
  #     if(length(bimodal_markers) > 0){
  #       expr_markers[bimodal_markers, first$indexes] <- marker_expr
  #       
  #       if(length(other_markers) > 0){
  #         
  #         w1 <- which(marker_expr == "+")
  #         w2 <- which(marker_expr == "-")
  #         
  #         el1 <- list(indexes = first$indexes[w1], markers = other_markers)
  #         el2 <- list(indexes = first$indexes[w2], markers = other_markers)
  #         
  #         queue <- c(queue, list(el1), list(el2))
  #       }
  #     }else{
  #       expr_markers[first$markers, first$indexes] <- "*"
  #     }
  #   }
  #   return(expr_markers)
  # 
  #   ######################################
  #   ## Not sequential marker evaluation ##
  #   ######################################
  # }else{

      queue <- list(list(indexes = 1:ncol(exp_matrix), markers = markers))

      while(length(queue) > 0){

        bimodal_markers <- c()
        not_bimodal_markers <- c()

        ## Pop operation
        first <- queue[[1]]
        queue <- queue[-1]

        for(m in first$markers){

          X <- exp_matrix[m, first$indexes]
          
          
          # ### C'è un bug
          # if(length(X) < 1){
          #   next
          # }
          
          t <- table(gates[, m])
          
          marker_expr <- set_marker_expression_GMM(X, t, m, r)

          if(length(table(marker_expr)) > 1){
            bimodal_markers <- c(bimodal_markers, m)

            if(verbose){
              message(paste0(" - ", paste0(bimodal_markers, collapse = " "), collapse = " "))
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
          
          ### Ci sono NA values in comb_markers
          comb_list <- as.list(comb_markers)
          names(comb_list) <- NULL
          
          to_add <- lapply(comb_list, function(c){
            w <- which(sapply(expr_markers[bimodal_markers, first$indexes, drop = F], function(c2){ return(all(c == c2))}))
            ## Push operation
            el <- list(indexes = first$indexes[w], markers = not_bimodal_markers)
            return(el)
          })
          queue <- c(queue, to_add)
        }
      }
    return(expr_markers)
  #}
}

## This function performs the cell classification
cell_classification <- function(marker_table, gates){
  
  # marker_table <- new_cells
  # gates <- new_gates2$extended_gate_table
  
  # ################ TO CHECK ###############
  # gates <- gates[!duplicated(gates$Gate), ]
  # #########################################
  
  rownames(gates) <- gates$Gate
  labels <- gates[marker_table$Gate, "Cell"]
  labels[is.na(labels)] <- "Unclassified"
  rownames(gates) <- NULL
  marker_table$Cell_type <- labels
  
  cl_res <- list(labels = labels, marker_table = gates, cell_signatures = marker_table)
  return(cl_res)
}

# extract_gate_table <- function(exp_matrix, markers, clusters){
#   
#   # clusters <- res$labels_post_clustering
#   # markers <- colnames(gates)[-1]
#   # exp_matrix <- m
#   
#   exp_matrix <- m[, clusters != "Unclassified"]
#   clusters <- clusters[clusters != "Unclassified"]
#   
#   message("Extracting gate table from data...")
#   useful_clusters <- unique(clusters)
#   
#   thresholds <- select_thresholds_from_clusters(exp_matrix,
#                                                 markers, 
#                                                 clusters, 
#                                                 prop_cluster_min = 0)
#   
#   gate_table <- data.frame(matrix(nrow = length(useful_clusters), ncol = length(markers) + 1))
#   colnames(gate_table) <- c("Cell", markers)
#   gate_table$Cell <- useful_clusters
#   rownames(gate_table) <- gate_table$Cell
#   
#   for(mm in markers){
#     
#     #mm <- "CD45"
#     
#     mds <- c()
#     q1s <- c()
#     q3s <- c()
#     
#     for(c in useful_clusters){
#       d <- exp_matrix[mm, as.character(clusters) %in% c]
#       md <- median(d)
#       q1 <- summary(d)["1st Qu."] 
#       q3 <- summary(d)["3rd Qu."]
#       names(md) <- c
#       names(q1) <- c
#       names(q3) <- c  
#       mds <- c(mds, md)
#       q1s <- c(q1s, q1)
#       q3s <- c(q3s, q3)
#     }
#     
#     pos <- which(thresholds[mm] < q1s)
#     neg <- which(thresholds[mm] > q3s)
#     star <- which(thresholds[mm] > q1s & thresholds[mm] < q3s)
#     
#     exp <- rep(NA, length(useful_clusters))
#     names(exp) <- useful_clusters
#     exp[names(star)] <- "*"
#     exp[names(pos)] <- "pos"
#     exp[names(neg)] <- "neg"
#     
#     gate_table[names(exp), mm] <- exp
#     
#   }
#   return(gate_table)
# }

## DA IMPLEMENTARE:
## 1) Possibilità di settare il parametro "hi", in questo caso "+" diventa il valore intermedio nel GMM con 3 
## componenti. Possibilità di settare il valore intermedio come "+" oppure come "-".
## 2) Restituire nei risultati la signature delle cellule "Unclassified"
## 3) Ottimizzazione del codice

## This is the core function for the classification of the single cells
scGateMe <- function(exp_matrix,
                     gates, 
                     refine = F,
                     k = NULL,
                     m = 2,
                     r = floor(0.01*ncol(exp_matrix)),
                     sampling = 1,
                     narrow_gate_table = T,
                     marker_seq_eval = F,
                     combine_seq_eval_res = T,
                     verbose = T,
                     seed = 1){
  
  # gates
  # refine = F
  # n_clusters = 20
  # cluster_level_labels = T
  # plots_thr = F
  # seed = 1
  # clusters = NULL
  # # clusters <- res$labels_post_clustering
  # exp_matrix <- m
  # ignore_markers = NULL
  # verbose = T
  # narrow_gate_table = T
  # marker_seq_eval = F
  # #trimodality = "pos"
  # sampling <- 1
  # combine_seq_eval_res = F
  
  set.seed(seed)
  
  if(sampling < 1){
    n <- floor(ncol(exp_matrix) * sampling)
    s <- sample(1:ncol(exp_matrix), n)
    exp_matrix_pre_sampling <- exp_matrix
    exp_matrix <- exp_matrix[, s]
  }else{
    exp_matrix_pre_sampling <- exp_matrix
  }
  
  # if(marker_seq_eval){
  #   
  #   res_all <- list()
  #   celltypes <- gates$Cell
  #   
  #   for(celltype in celltypes){
  #     
  #     if(verbose){
  #       message(paste0("Loading gate table for ", celltype, "..."))
  #     }
  #     
  #     new_gates <- parse_gate_table(gates[gates$Cell == celltype, , drop = F], narrow_gate_table)
  #     new_gates$extended_gate_table <- rbind(new_gates$extended_gate_table,
  #                                            data.frame(Cell = "Unclassified", 
  #                                                       Gate = paste0(colnames(new_gates$gate_table)[-1], "*", collapse = "")))
  #     
  #     exp_matrix_2 <- exp_matrix[colnames(new_gates$gate_table)[-1], , drop = F]
  #     
  #     if(verbose){
  #       message("Determining the marker signature for each cell...")
  #     }
  #     
  #     expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
  #     rownames(expr_markers) <- colnames(new_gates$gate_table)[-1]
  #     
  #     expr_markers <- set_marker_expression(exp_matrix_2, 
  #                                           colnames(new_gates$gate_table)[-1], 
  #                                           expr_markers, new_gates$gate_table, 
  #                                           verbose = verbose, 
  #                                           marker_seq_eval = marker_seq_eval,
  #                                           trimodality = trimodality)
  #     
  #     expr_markers <- expr_markers[colnames(new_gates$gate_table)[-1], ]
  #     expr_markers <- as.data.frame(t(expr_markers))
  #     cells <- cbind(rownames(expr_markers), expr_markers)
  #     colnames(cells)[1] <- "Cell"
  #     expr_markers <- apply(cells, 1, function(r){ gate <- paste(colnames(cells)[-1], r[-1], sep = "", collapse = ""); return(gate)})
  #     new_cells <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
  #     
  #     if(verbose){
  #       message("Classification of the cells...")
  #     }
  #     
  #     res <- cell_classification(new_cells, new_gates$extended_gate_table)
  #     res_all <- c(res_all, res)
  #     
  #     if(length(gates$Cell) > 1){
  #       
  #       if(verbose){
  #         message("------------------------------------------")
  #       }
  #     }
  #   }
  #   
  #   res_all <- res_all[names(res_all) == "labels"]
  #   names(res_all) <- gates$Cell
  #   
  #   if(length(res_all) > 1 & combine_seq_eval_res == T){
  #     res_all_df <- matrix(as.character(unlist(res_all)), ncol=length(res_all), byrow=F)
  #     colnames(res_all_df) <- gates$Cell
  #     
  #     if(verbose){
  #       message("Combine all labels...")
  #     }
  #     
  #     # res_all <- unlist(apply(res_all_df, 1, function(r){
  #     #   r <- as.character(r)
  #     #   t <- table(r)
  #     #   l <- names(t)
  #     #   ll <- l[l != "Unclassified"]
  #     #   if(length(ll) >= 1){
  #     #     return(paste(ll, collapse = " | "))
  #     #   }else{
  #     #     return("Unclassified")
  #     #   }
  #     # }))
  #     
  #     res_all <- rep(NA, nrow(res_all_df))
  #     
  #     w <- apply(res_all_df, 1, function(x){
  #       return(any(x != "Unclassified"))
  #     })
  #     
  #     if(sum(w) != nrow(res_all_df)){
  #       res_all_df_1 <- res_all_df[w, ]
  #       p <- apply(res_all_df_1, 1, function(x){
  #         return(paste(x[x != "Unclassified"], collapse = " | "))
  #       })
  #       w2 <- which(w)
  #       res_all[w2] <- p
  #       res_all[-w2] <- "Unclassified"
  #     }else{
  #       res_all <- "Unclassified"
  #     }
  #     
  #     res <- list(labels = res_all, marker_table = gates, cell_signatures = res$cell_signatures)
  #     
  #   }else{
  #     res <- list(labels = res_all, marker_table = gates, cell_signatures = res$cell_signatures)
  #   }
  # }else{
    
    if(verbose){
      message("Loading gate table...")
    }
    
    new_gates <- parse_gate_table(gates, narrow_gate_table)
    
    ### Cells with the same gate are "Unclassified" (Check the classification function!!!!!!!)
    not_dup <- !duplicated(new_gates$extended_gate_table$Gate)
    new_gates$extended_gate_table <- new_gates$extended_gate_table[not_dup, ]
    
    # if(verbose){
    #   message("Prioritizing gate table...")
    # }
    # 
    # new_gates2 <- prioritize_gate_table(new_gates)
    
    new_gates2 <- new_gates
    new_gates2$extended_gate_table <- rbind(new_gates2$extended_gate_table,
                                            data.frame(Cell = "Unclassified", 
                                                       Gate = paste0(colnames(new_gates2$gate_table)[-1], "*", collapse = "")))
    
    exp_matrix_2 <- exp_matrix[colnames(new_gates2$gate_table)[-1], , drop = F]
    
    
    if(verbose){
      message("Determining the marker signature for each cell...")
    }
    
    expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
    rownames(expr_markers) <- colnames(new_gates2$gate_table)[-1]
    
    expr_markers <- set_marker_expression(exp_matrix_2,
                                          colnames(new_gates2$gate_table)[-1], 
                                          expr_markers, 
                                          new_gates2$gate_table, 
                                          verbose = verbose, 
                                          marker_seq_eval = marker_seq_eval,
                                          trimodality = trimodality,
                                          m,
                                          r)
    
    expr_markers <- expr_markers[colnames(new_gates2$gate_table)[-1], ]
    expr_markers <- as.data.frame(t(expr_markers))
    cells <- cbind(rownames(expr_markers), expr_markers)
    colnames(cells)[1] <- "Cell"
    
    markers <- colnames(cells)[-1]
    
    expr_markers <- apply(cells, 1, function(r){
      return(paste0(markers, r[-1], collapse = ""))
      # return(gate)
    })
    
    new_cells <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
    
    if(verbose){
      message("Classification of the cells...")
    }
    
    res <- cell_classification(new_cells, new_gates2$extended_gate_table)
  #}
  
  
  ###### DA RIVEDERE PER L'ANNOTAZIONE SEQUENZIALE ##########################
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
      message("Refinement of the labels using k-nn classification...")
    }
    
    train <- t(exp_matrix_pre_sampling[, !res$labels %in% c("Unclassified")])
    control <- t(exp_matrix_pre_sampling[, res$labels %in% c("Unclassified")])
    
    if(is.null(k)){
      t <- table(res$labels)
      tt <- t[which.min(t)]
      k <- floor(sqrt(tt))
    }
    
    knn_res <- knn(train,
                   control,
                   cl = factor(res$labels[!res$labels %in% c("Unclassified")]),
                   k = k,
                   l = 0,
                   prob = T)
    
    res$labels[res$labels == "Unclassified"] <- as.character(knn_res)
    res$cell_signatures[res$cell_signatures$Celltype == "Unclassified", "Celltype"] <- as.character(knn_res)
  }
  
  res <- list(labels = res$labels,
              marker_table = res$marker_table, 
              cell_signatures = res$cell_signatures)
  
  return(res)
}
