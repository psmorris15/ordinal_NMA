#generate code and JAGS data object; format of JAGS code follows from that given 
#   for the data example presented in Schmid et al. (2014)

gen_inputs_adj_cat <- function(Study, Treat, Outcomes, N, hyperpar_values){
  
  # determine the outcome categorizations for each arm
  unique_cats <- unique(!is.na(Outcomes))*1
  outcome_cat <- numeric(length = nrow(Outcomes))
  
  for(i in 1:nrow(Outcomes)){
    outcome_cat[i] <- which(apply(unique_cats, 1, function(x) return(all(x == !is.na(Outcomes[i,])))))
  }
  
  #order Study, Treat, Outcomes, and N by outcome_cat, Study, and Treat; JAGS code requires that
  #   data objects be ordered by outcome categorization
  idx <- order(outcome_cat, Study, Treat)
  
  Study <- Study[idx]
  Treat <- Treat[idx]
  Outcomes <- Outcomes[idx,]
  N <- N[idx]
  outcome_cat <- outcome_cat[idx]
  
  # calculate number of arms per outcome categorization
  n_arms_cat <- as.numeric(table(outcome_cat))
  
  # initialize lists to hold outcome data for each separate outcome categorization
  data_list <- list()
  rep_data_list <- list()
  
  #initialize list to hold character vectors that indicate columns from Outcomes
  #   that are combined in each outcome categorization
  groups <- list()
  
  # for each outcome categorization:
  for(i in 1:nrow(unique_cats)){
    # keep rows with given outcome categorization
    Outcomes_subset <- Outcomes[which(outcome_cat == i),]
    # remove columns with NAs
    Outcomes_subset <- Outcomes_subset[,colSums(is.na(Outcomes_subset))<nrow(Outcomes_subset)]
    Outcomes_subset <- unname(Outcomes_subset)
    
    # populate lists with appropriate outcome data
    data_list[[i]] <- Outcomes_subset
    rep_data_list[[i]] <- Outcomes_subset
    
    # calculate groups: columns of Outcomes that match up with individual or combined levels
    group_starts <- as.numeric(which(unique_cats[i,] == 1))
    
    group_char <- character(length = length(group_starts))
    
    for(j in 1:length(group_starts)){
      
      if(j < length(group_starts)){
        group_char[j] <- ifelse(group_starts[j+1] - group_starts[j] == 1, 
                                group_starts[j], paste0(group_starts[j], ":", group_starts[j+1] - 1))
      }
      if(j == length(group_starts)){
        group_char[j] <- ifelse(group_starts[j] == ncol(unique_cats), group_starts[j],
                                paste0(group_starts[j], ":", ncol(unique_cats)))
      }
    }
    
    groups[[i]] <- group_char
    
  }
  
  # name data lists
  names(data_list) <- paste0("D", 1:length(data_list))
  names(rep_data_list) <- paste0("DD", 1:length(data_list))
  
  # calculate baseline intervention for each arm
  comb <- data.frame(Study = Study, Treat = Treat)
  comb <- comb %>% group_by(Study) %>% summarize(baseline = min(Treat)) %>% as.data.frame()
  baseline.Treat <- numeric(length = length(Treat))
  for(i in 1:length(baseline.Treat)){
    baseline.Treat[i] <- comb[which(comb$Study == Study[i]), "baseline"]
  }
  
  ### generate jags code ###
  
  code <- paste0("model",
                 "\n{")
  
  # multinomial likelihood by outcome categorization
  
  for(i in 1:nrow(unique_cats)){
    
    code <- paste0(code, 
                   "\n\tfor (i in ", 1, ":", n_arms_cat[i], ") {",
                   "\n\t\tD", i, "[i,1:", length(groups[[i]]),"] ~ dmulti(P", i, "[i,],N", i, "[i])",
                   "\n\t\tN", i, "[i] <- sum(DD", i, "[i,1:", length(groups[[i]]), "])")
    
    
    group <- groups[[i]]
    
    for(j in 1:length(group)){
      
        if(grepl(pattern = ":", x = group[j]) == FALSE){
          code <- paste0(code, 
                         "\n\t\tP", i, "[i,", j, "] <- P[i", 
                         ifelse(i == 1, "", 
                          paste0("+", sum(n_arms_cat[1:(i-1)]))), ",", group[j], "]")
        }
      
      if(grepl(pattern = ":", x = group[j]) == TRUE){
        code <- paste0(code, 
                       "\n\t\tP", i, "[i,", j, "] <- sum(P[i", 
                       ifelse(i == 1, "", 
                              paste0("+", sum(n_arms_cat[1:(i-1)]))), ",", group[j], "])")
      }
        
    }
    
    
    
    code <- paste0(code, "\n\t}")
  }
  
  # response probabilities model
  
  code <- paste0(code,
                 "\n\tfor (i in 1:", sum(n_arms_cat), "){",
                 "\n\t\tfor (k in 1:", ncol(Outcomes), ") {",
                 "\n\t\t\tP[i,k] <- prods[i,k]/sum(prods[i,])",
                 "\n\t\t\tprods[i,k] <- prod(theta[i, 1:k])",
                 "\n\t\t\tlog(theta[i,k]) <- Eta[Study[i],k] + Delta[i,k]*(1-equals(Treat[i],baseline.Treat[i]))",
                 "\n\t}",
                 "\n\tDelta[i,1] <- 0",
                 "\n\tDelta[i,2:", ncol(Outcomes), "] <- d[Treat[i],] - d[baseline.Treat[i],]",
                 "\n\t}")
  
  for(i in 1:length(unique(Study))){
    
    code <- paste0(code, 
                   "\n\tEta[", i, ", 1] <- 0",
                   "\n\tEta[", i, ", 2:", ncol(Outcomes), "] ~ dmnorm(mean.Eta", baseline.Treat[i], 
                   "[], Prec.Eta", baseline.Treat[i], "[,])")
    
  }
                 
  code <- paste0(code, "\n\tfor (k in 1:", ncol(Outcomes) - 1, ") {",
                       "\n\t\t d[1,k] <- 0",
                       "\n\t}",
                       "\n\tfor (j in 2:", length(unique(Treat)), ") {",
                       "\n\t\td[j, 1:", ncol(Outcomes) - 1, "] ~ dmnorm(mean.d[], Prec.d[,])",
                       "\n\t}",
                       "\n}")
  
  
  # put necessary objects into jags.data list
  
  jags.data <- list()
  
  for(i in 1:length(data_list)){
    
    jags.data[[paste0("D", i)]] <- data_list[[i]]
  
  }
  
  for(i in 1:length(rep_data_list)){
    
    jags.data[[paste0("DD", i)]] <- rep_data_list[[i]]
    
  }
  
  jags.data[["Study"]] <- Study
  jags.data[["Treat"]] <- Treat
  jags.data[["baseline.Treat"]] <- baseline.Treat
  jags.data[["mean.d"]] <- rep(0, ncol(Outcomes) - 1) 
  jags.data[["Prec.d"]] <- diag(1e-06, ncol(Outcomes) - 1) 
  
  for(i in 1:length(unique(baseline.Treat))){
    
    jags.data[[paste0("mean.Eta", unique(baseline.Treat)[i])]] <- hyperpar_values[[i]]$mean.Eta 
    jags.data[[paste0("Prec.Eta", unique(baseline.Treat)[i])]] <- hyperpar_values[[i]]$Prec.Eta 
    
  }
  
  # JAGS code
  jags.code <- code
  
  # output jags.data and jags.code
  output <- list(jags.data = jags.data, jags.code = jags.code)  
  
  return(output)
}
