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
    sets_to_filter <- data.frame(matrix(as.character(sets_to_filter), ncol = 2))
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
  
  # print(class(sets_all_filtered))
  
  sets_all_filtered_temp <- sets_all_filtered[, -which(colnames(sets_all_filtered) == "Cell"), drop = F]
  markers <- colnames(sets_all_filtered_temp)
  
  # system.time(gate <- apply(sets_all_filtered_temp, 1, function(r){
  #   # paste0(markers, r, collapse = "")
  #   # names(r) <- NULL
  #   # rle <- rle(r)
  #   # rle <- paste(rle$values, rle$lengths, collapse = "", sep = "")
  #   # return(rle)
  #   paste0(r, collapse = "")
  #   # stringi::stri_c(r, collapse = "")
  # }))
  
  gate <- apply(sets_all_filtered_temp, 1, stringi::stri_c, collapse = "")
  
  
  sets_all_filtered$Gate <- gate
  sets_all_filtered <- sets_all_filtered[, c("Cell", "Gate")]
  return(sets_all_filtered)
}

## Read the gate table and generate the possbile marker signature of ceach cell type
parse_gate_table <- function(gate_table, narrow_gate_table){
  
  # gate_table <- gates
  # narrow_gate_table = T
  
  if(any(duplicated(gate_table$Cell))){
    stop("Error! The gate table must contains uniquely defined cell types!")
  }
  
  if(narrow_gate_table){
    
    gates2 <- gate_table
    temp <- sapply(gates2$Gate, function(x){ str <- strsplit(x, "[-]|-[\\|\\^]|[+]|\\+[\\|\\^]|[*]"); return(str)})
    
    # print(temp)
    
    names(temp) <- gates2$Cell
    markers <- unique(unlist(temp))
    df_gates <- data.frame(matrix(nrow = nrow(gates2), ncol = length(markers)))
    rownames(df_gates) <- gates2$Cell
    colnames(df_gates) <- markers
    
    # gates2_exploded <- strsplit(gates2$Gate, split = "")
    gates2_exploded <- strsplit(perl=T, gates2$Gate, '(?![\\+\\||\\-\\||\\+\\^|\\-\\^])');
    signs <- lapply(gates2_exploded, function(x){ x[x %in% c("-^", "+^", "-|", "+|", "+", "-", "*")]})
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
    
  system.time(df_list <- apply(gate_table, 1, function(v){
    
    # v <- as.character(gate_table[2,-1])
    # names(v) <- colnames(gate_table[2,-1])
    # cell <- gate_table[2,1]
    
    to_add <- generate_set_values(v[-1], v[1])
    return(to_add)
  }))

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

set_marker_expression_GMM <- function(X, indexes, gmm_criteria, mm, rr){
  
  
  # X
  # gmm_criteria = "BIC"
  # mm = 2
  # rr = 0.05
  # indexes = first$indexes
  
  
  # 
  # X <- m["CD123", sce2$labels %in% c("CD34+CD38+CD123-_HSPCs", "CD34+CD38+CD123+_HSPCs")]
  # skewness(X)
  # plot(density(X))
  # 
  # bic <- mclustBIC(X, G = 1:2, verbose = F)
  # icl <- mclustICL(X, G = 1:2, verbose = F)
  # 
  # 
  # 
  # is.multimodal(test)
  # 
  # s <- obsno.Mrss(X, m = 2, r = 1000, type = "r")
  # index <- as.numeric(gsub("Obs.  ", "", s[, 2]))
  # test <- as.numeric(X[index])
  # plot(density(test))
  # 
  # icl <- mclustBIC(test, G = 1:2, verbose = F)
  # model_temp <- unlist(strsplit(names(summary(icl)[1]), ","))
  # type_model <- model_temp[1]
  # model <- as.numeric(model_temp[2])
  # cl <- Mclust(test, G = model, verbose = F, modelNames = "E")
  # plot(cl, what = "classification")
  # 
  # cl <- Mclust(X, G = 2, verbose = F)
  # plot(cl, what = "classification")
  # 
  
  # X <- m["CD3", ]
  
  if(length(X) >= 100){
    
    rr <- ceiling(rr*length(X))
    
    if(rr < 2){
      rr <- 100
    }
    
    # s <- obsno.Mrss(X, m = mm, r = rr, type = "r")
    
    ############ Ranket Set Sampling (RSS) ##############
    sample = X
    cycles = length(sample) / 4
    n_unit = 2
    n <- length(sample)
    samples <- cycles * n_unit
    indexes <- matrix(sample(1:n, samples * n_unit), nrow = n_unit)
    sel_samples <- matrix(sample[indexes], nrow = n_unit)
    #####################################################
    
    
    # if(skewness(X) < 0){
    #   index <- as.numeric(gsub("Obs.  ", "", s[, 1]))
    # }else{
    #   index <- as.numeric(gsub("Obs.  ", "", s[, 2]))
    # }
    
    if(skewness(X) < 0){
      test <- apply(sel_samples, 2, min)
    }else{
      test <- apply(sel_samples, 2, max)
    }
    
    # test <- as.numeric(X[index])
    
    # plot(density(test))
    
  }else if(length(X) >= 2){
    test <- X
  }else{
    return(c("*"))
  }
  
  # plot(density(test))
  # cl <- Mclust(test, G = 2, modelNames = "E")
  # plot(cl, what = "classification")
  
  switch(gmm_criteria,
         # best = {
         #   bic <- mclustBIC(test, G = 1:2, verbose = F)
         #   icl <- mclustICL(test, G = 1:2, verbose = F)
         #
         #   if(summary(bic)[1] >= summary(icl)[1]){
         #     crit <- mclustBIC(test, G = 1:2, verbose = F)
         #   }else{
         #     crit <- mclustICL(test, G = 1:2, verbose = F)
         #   }
         # },
         BIC = {
           crit <- mclustBIC(test, G = 1:2, verbose = F)
         },
         ICL = {
           crit <- mclustICL(test, G = 1:2, verbose = F)
         },
  )
  
  model_temp <- unlist(strsplit(names(summary(crit)[1]), ","))
  type_model <- model_temp[1]
  model <- as.numeric(model_temp[2])
  
  if(!is.na(model) & model > 1){
    
    # cl <- Mclust(test, G = model, verbose = F, modelNames = type_model)
    cl <- Mclust(test, G = 2, verbose = F, modelNames = "E")
    
    #plot(cl, what = "classification")
    
    
    #temp <- cl$classification
    pred <- predict.Mclust(cl, X)
    temp <- pred$classification
    means <- cl$parameters$mean
    means <- sort(means)
    
    # if(model == 2){
    
    temp[temp == names(means)[1]] <- "-"
    temp[temp == names(means)[2]] <- "+"
      
    # pred_df <- data.frame(Pred = temp, pred$z)
    # pred_df$indexes <- indexes
    # colnames(pred_df)[c(2,3)] <- c("1", "2")
    # 
    # pred_df_minus <- pred_df[, c("Pred", names(means)[1], "indexes")]
    # index_minus <- pred_df_minus[order(pred_df_minus[, names(means)[1]], decreasing = T), "indexes"][1:10]
    # 
    # pred_df_plus <- pred_df[, c("Pred", names(means)[2], "indexes")]
    # index_plus <- pred_df_plus[order(pred_df_plus[, names(means)[2]], decreasing = T), "indexes"][1:10]
    
    # View(pred_df)
      
    #}else{
      #temp[temp == names(means)[1]] <- "-"
      
      # if(trimodality == "strictly_pos"){
      #   temp[temp == names(means)[2]] <- "-"
      # }else{
      #   temp[temp == names(means)[2]] <- "+"
      # }
      
      #temp[temp == names(means)[3]] <- "+"
    #}
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
                                  gmm_criteria,
                                  mm,
                                  rr){


  # exp_matrix <- exp_matrix_2
  # markers <- colnames(new_gates$gate_table)[-1]
  # # # expr_markers
  # gates <- new_gates$gate_table
  # # verbose = F
  # marker_seq_eval = F
  # mm <- 2
  # rr <- 0.05
  # gmm_criteria <- "ICL"

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
      # 
      # markers[1] <- "CD274"
      # markers[8] <- "FS_A"

      queue <- list(list(indexes = 1:ncol(exp_matrix), markers = markers))

      while(length(queue) > 0){


        bimodal_markers <- c()
        not_bimodal_markers <- c()

        ## Pop operation
        first <- queue[[1]]
        queue <- queue[-1]

        
        # first$markers[1] <- "CD274"
        # first$markers[8] <- "FS_A"
      
        
        for(m in first$markers){

          # 
          # if(m == "CD274"){
          #   stop()
          # }
  
          
          #m <- "CD274"
          
          X <- exp_matrix[m, first$indexes]
          
          
          # ### C'è un bug
          # if(length(X) < 1){
          #   next
          # }
          
          # t <- table(gates[, m])
          
          marker_expr <- set_marker_expression_GMM(X, indexes = first$indexes, gmm_criteria, mm, rr)
          


          if(length(table(marker_expr)) > 1){
            
            bimodal_markers <- c(bimodal_markers, m)

            # if(bimodal_markers == "CD44"){
            #   stop()
            # }
            
            
            if(verbose){
              message(paste0(" - ", paste0(bimodal_markers, collapse = " "), collapse = " "))
            }

            expr_markers[m, first$indexes] <- marker_expr
            
            
          }else{
            not_bimodal_markers <- c(not_bimodal_markers, m)
          }
          
          
          # print(table(expr_markers["CD274", ], nr_t0$labels))
          
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
          
          # if(m == "CXCR4"){
          #   stop()
          # }
      
          
          comb_list2 <- sapply(comb_list, paste0, collapse = "")
          l <- expr_markers[bimodal_markers, first$indexes, drop = F]
          l <- sapply(l, paste0, collapse = "")
          
          to_add <- lapply(comb_list2, function(c){
            
            # c <- comb_list[[1]]
            
            
            # w <- which(sapply(expr_markers[bimodal_markers, first$indexes, drop = F], function(c2){ 
            #   return(identical(c, c2))
            # }))
            
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
  #}
}

## This function performs the cell classification
cell_classification <- function(marker_table, gates){
  
  # marker_table <- new_cells2
  # gates <- new_gates2$extended_gate_table
  
  # 
  # 
  # intersect(marker_table$Gate, gates$Gate)
  
  # ################ TO CHECK ###############
  # gates <- gates[!duplicated(gates$Gate), ]
  # #########################################
  
  # rownames(gates) <- gates$Gate
  
  
  df_gate <- marker_table
  colnames(df_gate)[1] <- "Cell_ID"
  
  
  pos <- sapply(gates$Gate, function(g){
    gs <- unlist(strsplit(g, split = ""))
    w <- which(gs != "*")
    w <- paste(w, collapse = "_", sep = "")
    return(w)
  })
  
  gates$Pos <- pos
  
  clean <- sapply(gates$Gate, function(g){
    gs <- unlist(strsplit(g, split = ""))
    w <- which(gs != "*")
    clean_gs <- paste(gs[w], collapse = "", sep = "")
    return(clean_gs)
  })
  
  gates$Gate <- clean
  
  gates$N <- sapply(gates$Gate, function (g){
    return(sum(unlist(strsplit(g, split = "")) == "+"))
  })
  
  gates <- gates[order(gates$N, decreasing = F), ]
  
  labels <- rep("Unclassified", nrow(df_gate))
  
  for(i in 1:nrow(gates)){
    
    pos <- as.numeric(unlist(strsplit(gates$Pos[i], split = "_")))
    
    test <- sapply(df_gate$Gate, function(g){
      

      
      x <- paste0(unlist(strsplit(g, split = ""))[pos], collapse = "")
      
      
      return(x == gates$Gate[i])
      
      
    })
    
    labels[test] <- gates$Cell[i]
    
  }
  
  

  
  
  
  # 
  # gates$Gate <- clean
  # 
  # clean <- sapply(df_gate$Gate, function(g){
  #   gs <- unlist(strsplit(g, split = ""))
  #   # w <- which(gs != "*")
  #   clean_gs <- paste(w, gs, collapse = "", sep = "")
  # })
  # 
  # df_gate$Gate <- clean
  # 
  # l <- lapply(gates$Gate, function(g){
  #   
  #   # g <- gates$Gate[1]
  #   
  #   w <- str_which(df_gate$Gate, coll(g))
  #   
  #   
  #   #w <- which(d)
  #   return(w)
  # })
  
  # labels <- rep("Unclassified", nrow(df_gate))
  # 
  # for(i in 1:length(l)){
  # 
  #   el <- l[i]
  #   index <- unlist(el)
  #   
  #   if(length(index) > 0){
  #     labels[index] <- names(el)
  #     
  #   }
  #   
  # }

  # labels <- plyr::join(df_gate, gates, by='Gate')
  # # labels <- merge(df_gate, gates, by = "Gate", all.x = T)
  # labels <- as.character(labels$Cell)
  # labels[is.na(labels)] <- "Unclassified"
  
  rownames(gates) <- NULL
  marker_table$Celltype <- labels
  
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
                     train = F,
                     gates, 
                     gmm_criteria = "BIC",
                     refine = F,
                     k = NULL,
                     mm = 2,
                     rr = 0.05,
                     sampling = 1,
                     narrow_gate_table = T,
                     marker_seq_eval = F,
                     combine_seq_eval_res = T,
                     verbose = T,
                     seed = 1){
  
  # gates <- gate_table
  # refine = F
  # n_clusters = 20
  # cluster_level_labels = T
  # plots_thr = F
  # seed = 2
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
  # mm = 2
  # rr = 0.05
  # gmm_criteria = "BIC"
  # train = T
  # k = NULL

  # print(dim(exp_matrix))
  
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
    
    #if(is.null(reference)){
      if(verbose){
        message("Loading gate table...")
      }
      
      new_gates <- parse_gate_table(gates, narrow_gate_table)
      
      ### Cells with the same gate are "Unclassified" (Check the classification function!!!!!!!)
      dup <- new_gates$extended_gate_table$Gate[duplicated(new_gates$extended_gate_table$Gate)]
      new_gates$extended_gate_table <- new_gates$extended_gate_table[!new_gates$extended_gate_table$Gate %in% dup, ]
      
      # if(verbose){
      #   message("Prioritizing gate table...")
      # }
      # 
      # new_gates2 <- prioritize_gate_table(new_gates)
      
      new_gates2 <- new_gates
      # new_gates2$extended_gate_table <- rbind(new_gates2$extended_gate_table,
      #                                         data.frame(Cell = "Unclassified", 
      #                                                    Gate = paste0(rep("*", ncol(new_gates2$gate_table)), collapse = "")))
      
      exp_matrix_2 <- exp_matrix[colnames(new_gates2$gate_table)[-1], , drop = F]
      
      if(verbose){
        message("Determining the marker signature for each cell...")
      }
      
      expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
      rownames(expr_markers) <- colnames(new_gates2$gate_table)[-1]
      
    #}
      # else{
      # markers <- colnames(reference)[-which(colnames(reference) == "Cell")]
      # 
      # if(length(intersect(rownames(m), markers)) == 0){
      #   stop("No common markers between current expression matrix and reference dataset!")
      # }
      # 
      # gate_table <- data.frame(Cell = c("Unknown"), Gate = paste0(markers, "+", collapse = ""))
      # new_gates2 <- parse_gate_table(gate_table, narrow_gate_table)
      # 
      # ref <- reference
      # labels <- ref$Cell
      # ref <- ref[, -1]
      # ref <- t(ref)
      # s <- sample(1:ncol(ref), 10000)
      # ref <- ref[, s]
      # labels <- labels[s]
      # 
      # ref <- ref[colnames(new_gates2$gate_table)[-1], , drop = F]
      # expr_markers <- data.frame(matrix(ncol = ncol(ref), nrow = nrow(ref)))
      # rownames(expr_markers) <- colnames(new_gates2$gate_table)[-1]
      # 
      # ## Obtain the signature for the cells of the reference dataset 
      # expr_markers <- set_marker_expression(ref,
      #                                       colnames(new_gates2$gate_table)[-1], 
      #                                       expr_markers, 
      #                                       new_gates2$gate_table, 
      #                                       verbose = verbose, 
      #                                       marker_seq_eval = marker_seq_eval,
      #                                       gmm_criteria = gmm_criteria,
      #                                       mm,
      #                                       rr)
      # 
      # expr_markers <- expr_markers[colnames(new_gates2$gate_table)[-1], ]
      # expr_markers <- as.data.frame(t(expr_markers))
      # cells <- cbind(rownames(expr_markers), expr_markers)
      # colnames(cells)[1] <- "Cell"
      # 
      # markers <- colnames(cells)[-1]
      # 
      # expr_markers <- apply(cells, 1, function(r){
      #   return(paste0(markers, r[-1], collapse = ""))
      #   # return(gate)
      # })
      # 
      # new_cells <- data.frame(Cell = names(expr_markers), Gate = expr_markers)
      # new_cells$Cell <- labels
      # new_cells <- na.omit(new_cells)
      
    #}
  
    exp_matrix_2 <- exp_matrix[colnames(new_gates2$gate_table)[-1], , drop = F]
    expr_markers <- data.frame(matrix(ncol = ncol(exp_matrix_2), nrow = nrow(exp_matrix_2)))
    rownames(expr_markers) <- colnames(new_gates2$gate_table)[-1]
  
    ## Obtain the signature for the cells of the dataset to be annotated
    expr_markers <- set_marker_expression(exp_matrix_2,
                                          colnames(new_gates2$gate_table)[-1], 
                                          expr_markers, 
                                          new_gates2$gate_table, 
                                          verbose = verbose, 
                                          marker_seq_eval = marker_seq_eval,
                                          gmm_criteria = gmm_criteria,
                                          mm,
                                          rr)
    
    expr_markers <- expr_markers[colnames(new_gates2$gate_table)[-1], ]
    expr_markers <- as.data.frame(t(expr_markers))
    
    cells <- cbind(rownames(expr_markers), expr_markers)
    colnames(cells)[1] <- "Cell"
    
    # markers <- colnames(cells)[-1]
    
    # expr_markers <- apply(cells, 1, function(r){
    #   return(paste0(markers, r[-1], collapse = ""))
    #   # return(gate)
    # })
    
    expr_markers <- apply(cells[,- 1], 1, stringi::stri_c, collapse = "")
    new_cells2 <- data.frame(Cell = names(expr_markers), Gate = expr_markers)

    if(!train){
      message("Classification of the cells...")
      res <- cell_classification(new_cells2, new_gates2$extended_gate_table)
    }else{
      res <- list(cell_signatures = new_cells2)
    }
    
    # else{
    #   res <- cell_classification(new_cells2, new_cells)
    # }
    
  #}
  
  
  ###### DA RIVEDERE PER L'ANNOTAZIONE SEQUENZIALE ##########################
  
  if(!train){
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
  }
  
  # not_na <- which(!is.na(res$cell_signatures))
  # signatures <- res$cell_signatures[not_na, ]
  
  signatures <- res$cell_signatures
  
  markers <- colnames(new_gates$gate_table)[-1]
  
  gate_ext <- sapply(signatures$Gate, function(v){
    split <- unlist(strsplit(v, ""))
    sig <- stringi::stri_c(markers, split, collapse = "", sep = "")
    return(sig)
  })
  
  res$cell_signatures[, "Gate"] <- gate_ext
  
  if(train){
    res <- res$cell_signatures
  }else{
    res <- list(labels = res$labels,
                marker_table = res$marker_table, 
                cell_signatures = res$cell_signatures)
  }

  return(res)
}

scGateMe_train <- function(reference, labels, gmm_criteria = "BIC", sampling, thr_pos = 0.75, thr_neg = 0.02, rr = 0.1, seed = 1){
  
  # reference <- m
  # labels <- sce2$labels
  # gmm_criteria = "BIC"
  # thr = 0.95
  # sampling = 1
  # seed = 300
  # rr = 0.1
  
  set.seed(seed)
  
  s <- c()
  
  if(sampling < 1){
    for(c in unique(labels)){
      w <- which(labels == c)
      n <- floor(length(w) * sampling)
      s <- c(s, sample(w, n))
    }
    
    reference <- reference[, s]
    labels <- labels[s] 
    
  }
  
  # ms <- sample(1:32, 20)
  # reference <- reference[ms, s]

  
  markers <- rownames(reference)
  celltypes <- factor(unique(labels))
  gate_table <- data.frame(Cell = c("Unknown"), Gate = paste0(markers, "+", collapse = ""))
  
  res <- scGateMe(reference,
                  train = T,
                  gate_table, 
                  gmm_criteria = gmm_criteria,
                  refine = F,
                  sampling = 1,
                  k = NULL,
                  verbose = T,
                  narrow_gate_table = T, 
                  marker_seq_eval = F,
                  combine_seq_eval_res = F,
                  rr = rr,
                  seed = seed)
  
  message("-----------------------")
  cell_df <- res
  cell_df$labels <- labels
  cell_df <- cell_df[, c("Gate", "labels")]
  cell_df <- na.omit(cell_df)
  
  gate_table <- table(cell_df$labels, cell_df$Gate)
  gate_table <- as.data.frame(data.table(gate_table))
  colnames(gate_table) <- c("Cell", "Gate", "N")
  
  # gate_table <- na.omit(gate_table)
  gate_table <- gate_table[gate_table$N > 1, ]
  # labels <- gate_table$Cell
  
  labels2 <- gate_table$Cell
  
  # gate_table <- gate_table[gate_table$Cell == "CD4_T_cells", ]
  gate_table$Cell <- paste(gate_table$Cell, as.character(1:nrow(gate_table)), sep = "_")
  gate_table$Cell_2 <- labels2
  
  gates2 <- parse_gate_table(gate_table, T)
  
  #### Check order !!!!
  gates2$gate_table <- gates2$gate_table[rep(seq_len(nrow(gates2$gate_table)), times = gate_table$N), ]
  
  gates2 <- gates2$gate_table
  gates2$Cell <- rep(gate_table$Cell_2, times = gate_table$N)
  # gate3 <- reshape2::melt(gates2, id.vars = "Cell")
  
  new_gate_table <- data.frame(Cell = unique(gates2$Cell), Gate = "")
  rownames(new_gate_table) <- new_gate_table$Cell
  
  markers <- colnames(gates2[-1])
  celltypes <- unique(gates2$Cell)
  
  # n_comparisons <- length(markers) * length(celltypes)
  
  for(c in celltypes){
  #   
    message(paste0(" - ", c))
  #   
  # c <- "Basophils"
  #   
  #   sig <- ""
    
    #for(c2 in celltypes[celltypes != c]){
      
      # ma <- "CD274"
      
      # c2 <- "Monocytes"
      
      gates3 <- gates2
      
      gates3$Cell[gates3$Cell != c] <- "Other_type"
      # gates3 <- gates3[gates3$Cell %in% c(c,c2), ]
      
      gates3$Cell <- factor(as.character(gates3$Cell), levels = c(c, "Other_type"))
      
      # gates3$Cell <- factor(as.character(gates3$Cell), levels = c(c, c2))
      
      
      data <- gates3
      # data$Cell[data$Cell != c] <- "Other_type"
      # data <- data[data$Cell == c, ]
      
      s1 <- w1 <- which(data$Cell == c)
      s2 <- w2 <- which(data$Cell != c)
      
      s <- 2000
      
      if(length(w1) > s){
        s1 <- sample(w1, s)
      }
      
      if(length(w2) > s){
        s2 <- sample(w2, s)
      } 

      # s1 <- sample(w1, 2000)
      # s2 <- sample(w2, 2000)
      
      #if(nrow(data) > 100){
      data <- data[c(s1,s2), ]
      #}
      
      data <- as.data.frame(lapply(data, factor))
      
      # data$Cell <- factor(data$Cell)
      # data$FS_A <- factor(data$FS_A)
      # data$SS_A <- factor(data$SS_A) 
      # data$CD10 <- factor(data$CD10)
      # data$HLA_DR <- factor(data$HLA_DR)
      # data$CD14 <- factor(data$CD14)
      # data$CD38 <- factor(data$CD38)
      # data$CD279 <- factor(data$CD279)
      # data$CD274 <- factor(data$CD274)
      # data$CD45 <- factor(data$CD45)
      # data$CD3 <- factor(data$CD3)
      # data$CD16 <- factor(data$CD16)
      # data$CD4 <- factor(data$CD4)
      
      b <- Boruta(Cell ~ ., data = data, getImp=getImpFerns)
      mas <- names(b$finalDecision)[b$finalDecision  == "Confirmed"]
      # table(data$Cell, data$CD16)
      
      
      # rf_mod <- randomForest(Cell ~ ., data=data, ntree=100)
      # 
      # imp <- data.frame(rf_mod$importance)
      # imp <- imp[order(imp$MeanDecreaseGini, decreasing = T), , drop = F]
      # 
      # mas <- rownames(imp)
      # 
      # for(i in length(mas):1){
      # 
      #   if(sum(imp$MeanDecreaseGini[1:i-1]) / sum(imp$MeanDecreaseGini) >= 0.85){
      #     mas <- mas[mas != rownames(imp)[i]]
      #   }
      # }
      
      
  

      
      
      #to.remove<-c(which(data.frame(iris.rf$importance)$MeanDecreaseAccuracy==min(data.frame(iris.rf$importance)$MeanDecreaseAccuracy)))
      
      
      # varImpPlot(rf_mod)
      
      
      #control <- rfeControl(functions = rfFuncs, method = "repeatedcv")
      

      # cl <- makeCluster(detectCores())
      # registerDoParallel(cl)

      

      # result_rfe1 <- caret::rfe(x = data[, -1],
      #                    y = factor(data$Cell),
      #                    sizes = c(1:12),
      #                    rfeControl = rfeControl(functions = rfFuncs,
      #                                            method = "repeatedcv"))
      
      #stopCluster(cl)

      # mas <- result_rfe1$optVariables
      
      sig <- ""
      
      for(k in mas){
        
        # k <- "CD16"
        
        t <- prop.table(table(data[data$Cell == c, k]))
        
        # print(k)
        # print(t)
        
        max <- max(t)
        
        if(max >= 0.9){
          sign <- names(t)[which.max(t)]
          sig <- paste0(sig, k, sign)
        }
      }
      
      
      
    #   df <- as.data.frame.matrix(table(gates3$Cell, gates3[, ma]))
    #   df$Cell <- rownames(df)
    #   
    #   
    #   # df <- df + 1
    #   star <- which(colnames(df) == "*")
    #   
    #   if(length(star) > 0){
    #     df <- df[, -which(colnames(df) == "*")]
    #   }
    #   
    #   df[, -3] <- df[, -3] + 1
    #   df
    #   
    #   # f <- fisher.test(df)
    #   # 
    #   # if(f$p.value < 0.05 & f$estimate > 1){
    #   #   sig <- paste0(sig, ma, "-")
    #   # }else if(f$p.value < 0.05 & f$estimate < 1){
    #   #   sig <- paste0(sig, ma, "+")
    #   # }
    #   
    #   df$Cell <- rownames(df)
    #   
    #   #df$Cell[df$Cell != "Monocytes"] <- "Other"
    #   
    #   colnames(df) <- c("neg", "pos", "Cell")
    #   
    #   
    #   # df <- reshape2::melt(df, id.vars = "Cell")
    #   # df$Other <- sum(df$value) - df$value
    #   df$Cell <- factor(df$Cell, levels = c("Other", c))
    #   # df$value[df$value == 0] <- 1
    #   
    #   model0 <- glm(
    #     formula = cbind(pos, neg) ~ Cell,
    #     family = binomial(link = 'logit'),
    #     data = df
    #   )
    #   
    #   emm1 <- emmeans(model0, specs = revpairwise ~ Cell)
    #   emm1$contrasts %>%
    #     summary(infer = TRUE, type = 'response') %>%
    #     rbind() %>%
    #     as.data.frame() -> c_results
    #   
    #   c_results
    #   
    #   if(c_results$p.value < (0.05 / n_comparisons)){
    #     
    #     p <- prop.test(df[c, "pos"], df[c, "pos"] + df[c, "neg"])
    #     
    #     if(c_results$odds.ratio > 1 & p$p.value < 0.05 & p$estimate > thr_pos){
    #       sig <- paste0(sig, ma, "+")
    #     }else if(c_results$odds.ratio < 1 & p$p.value < 0.05 & p$estimate < thr_neg){
    #       sig <- paste0(sig, ma, "-")
    #     }
    #   }
    # #}
    
    new_gate_table[c, "Gate"] <- sig
    
  

    
    # w_plus <- which(c_results$p.value < 0.05 & c_results$odds.ratio > 1 & c_results$odds.ratio > 2)
    # w_minus <- which(c_results$p.value < 0.05 & c_results$odds.ratio < 1  & c_results$odds.ratio < 0.1)
    # 
    # if(length(w_minus) != length(celltypes) & length(w_plus) != length(celltypes)){
    #   sig_plus <- rep(paste0(ma, "+"), length(w_plus))
    #   sig_minus <- rep(paste0(ma, "-"), length(w_minus))
    # 
    #   new_gate_table[as.character(c_results$Cell[w_plus]), "Gate"] <- paste0(new_gate_table[as.character(c_results$Cell[w_plus]), "Gate"], as.character(sig_plus))
    #   new_gate_table[as.character(c_results$Cell[w_minus]), "Gate"] <- paste0(new_gate_table[as.character(c_results$Cell[w_minus]), "Gate"], as.character(sig_minus))
    # }
  }
  
  # for(cell in unique(gates2$Cell)){
  #   
  #   # cell <- "Mature_B_cells"
  #   
  #   test <- gates2[gates2$Cell == cell, ]
  #   sig <- c()
  #   
  #   acc <- c()
  #   
  #   for(s in c("+", "-")){
  #     
  #     # s <- "+"
  #   
  #     for(j in 1:100){
  #       
  #       # j <- 1
  #       
  #       message(paste0(cell, " iteration ", j))
  #       
  #       n <- ceiling(0.1 * nrow(test))
  #       
  #       #test2 <- test[sample(1:nrow(test), n), ]
  #       
  #       
  #       if(n > 100){
  #         test2 <- test[sample(1:nrow(test), 10), ]
  #       }else{
  #         test2 <- test
  #       }
  #       
  #       # test2 <- test2[, -1]
  #       
  #       
  #       for(i in 2:nrow(test2)){
  #         
  #         #i <- 2
  #         
  #         if(i == 2){
  #           
  #           w1 <- which(test2[i, -1] == s)
  #           w2 <- which(test2[1, -1] == s)
  #           
  #           a <- paste(colnames(test2)[-1][w1], test2[i, -1][w1], sep = "")
  #           b <-  paste(colnames(test2)[-1][w2], test2[1, -1][w2], sep = "")
  #           prec <- LCS(a,b)
  #           
  #           
  #         }else{
  #           w1 <- which(test2[i, -1] == s)
  #           a <- paste(colnames(test2)[-1][w1], test2[i, -1][w1], sep = "")
  #           b <- prec$LCS
  #           prec <- LCS(a, b)
  #         }
  #       }
  #       
  #       signature <- paste0(sort(prec$LCS), collapse = "")
  #       # names(signature) <- n
  #       
  #       if(signature != ""){
  #         acc <- unique(c(acc, prec$LCS))
  #       }
  #     }
  #     
  #     sig <- paste0(sort(acc), collapse = "")
  #     
  #     
  #     # acc <- table(acc)
  #     # ord <- order(acc, decreasing = T)
  #     # acc <- acc[ord]
  #     # acc
  #     
  #     # acc <- acc[order(acc, decreasing = T)]
  #     # sig <- paste0(sig, names(acc[1]), sep = "")
  #   }
  #   
  #   new_gate_table[cell, "Gate"] <- sig
  #   
  # }
  
  rownames(new_gate_table) <- NULL
  new_gate_table <- new_gate_table[new_gate_table$Gate != "", ]
  return(new_gate_table)
  
}
















