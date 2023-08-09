adj_cat_baseline_priors <- function(Outcomes, Treat, Study, nadapt, nburnin, niter, n.chains, seeds){
  
  require(dplyr)
  require(rjags)
  
  # determine the baseline intervention for each trial and put into integer vector
  baseline_Treat <- as.data.frame(cbind(Study, Treat)) %>% 
    group_by(Study) %>% summarise(baseline_Treat = min(Treat)) %>% 
    select(baseline_Treat)
  baseline_Treat <- baseline_Treat$baseline_Treat
  
  # get the unique values from baseline_Treat
  baseline_Treats <- unique(baseline_Treat)
  
  # add 0.5 to each cell that contains a 0 event count to ensure our empirical
  #   log-odds fall on the real line
  for(i in 1:nrow(Outcomes)){
    
    if(0 %in% Outcomes[i,]){
      
      Outcomes[i,] <- Outcomes[i,] + 0.5
      
    }

  }
  
  # JAGS model we will use to obtain hyperparameter values for the mean and precision
  #   of the Eta parameters; code modified from baseline effects model code in
  #   Hu et al. (2020)
  model <- 
    "
  ## get prior value for baseline treatment in logit scale         
  model{                          # *** PROGRAM STARTS

    for (i in 1:ns){                # LOOP THROUGH STUDIES
      eta[i] ~ dnorm(m,tau.m)       
    }


    m ~ dnorm(0, hyperVarInNormMean)              # vague prior for mean
    sd.m ~ dunif(0, hyperSDInUnif)               # vague prior for between-trial SD
    tau.m <- pow(sd.m,-2)           # between-trial precision = (1/between-trial variance)

  }
  "
  
  # initialize list to store results
  hyperpar_values <- list()

  # for each unique trial-specific baseline treatment in the network:
  for(i in 1:length(baseline_Treats)){
    
    # keep rows of outcome matrix corresponding to the given baseline intervention
    Outcomes_baseline <- Outcomes[which(Treat == baseline_Treats[i]),]
    
    # initialize m vector
    post_means_m <- numeric(length = ncol(Outcomes) - 1)
    names(post_means_m) <- paste0("level_pair", 1:(ncol(Outcomes) - 1))
    
    # initialize sd vector; will later be transformed to precision matrix
    post_means_sd <- numeric(length = ncol(Outcomes) - 1)
    names(post_means_sd) <- paste0("level_pair", 1:(ncol(Outcomes) - 1))
    
    # initialize matrix for the Gelman-Rubin diagnostic values
    gelman_diags <- matrix(nrow = ncol(Outcomes) - 1, ncol = 2, data = NA)
    row.names(gelman_diags) <- paste0("level_pair", 1:(ncol(Outcomes) - 1))
    colnames(gelman_diags) <- c("m", "sd.m")
    
    # for each level comparison:
    for(j in 2:ncol(Outcomes)){
      
      # keep only rows that observed complete data for that level comparison
      if(j < ncol(Outcomes)){
        
        Outcomes_baseline_full <- Outcomes_baseline[!is.na(Outcomes_baseline[, j-1]) & !is.na(Outcomes_baseline[, j]) & !is.na(Outcomes_baseline[, j + 1]),]  
        
      }
      
      if(j == ncol(Outcomes)){
        
        Outcomes_baseline_full <- Outcomes_baseline[!is.na(Outcomes_baseline[, j-1]) & !is.na(Outcomes_baseline[, j]),]   
        
      }
      
      # calculate the vector of empirical log-odds
      empirical_LO <- log(Outcomes_baseline_full[,j]/Outcomes_baseline_full[,j-1])
      
      # if seed is specified:
      if(!is.null(seeds)){
        
        # add seeds to the inits list
        inits <- list()
          
        for(k in 1:n.chains){
          inits[[k]] <- sapply(c(".RNG.seed", ".RNG.name"),function(x) NULL)
          inits[[k]]$.RNG.seed <- seeds[[i]][k, j-1]
          inits[[k]]$.RNG.name <- "base::Wichmann-Hill" 
        
        }
        
        # initialize the model
        mod = rjags::jags.model(textConnection(model),
                                data = list(eta = empirical_LO,
                                            ns = length(empirical_LO),
                                            hyperSDInUnif = 5,
                                            hyperVarInNormMean = 1e-4),
                                inits = inits,
                                n.chains=n.chains,
                                n.adapt=nadapt)
        
      }
      
      # if seed is not specified:
      if(is.null(seeds)){
        
        # initialize the model
        mod = rjags::jags.model(textConnection(model),
                                data = list(eta = empirical_LO,
                                            ns = length(empirical_LO),
                                            hyperSDInUnif = 5,
                                            hyperVarInNormMean = 1e-4),
                                n.chains=n.chains,
                                n.adapt=nadapt)        
        
      }
      
      # run the burn-in phase
      update(mod, n.iter = nburnin)
      
      # run the sampling phase
      res = coda.samples(mod, variable.names = c("m", "sd.m"), niter)
      
      # input appropriate results
      gelman_diags[j-1, 1] <- gelman.diag(res, multivariate = FALSE)[[1]][1,1]
      gelman_diags[j-1, 2] <- gelman.diag(res, multivariate = FALSE)[[1]][2,1]
      post_means_m[j-1] <- summary(res)[[1]][,1][1]
      post_means_sd[j-1] <- summary(res)[[1]][,1][2]
      
    }
    
    # convert standard deviations to precision
    prec.m <- solve(diag(post_means_sd^2))
    colnames(prec.m) <- row.names(prec.m) <- names(post_means_sd)
    
    # put results into list element
    hyperpar_values[[i]] <- list(mean.Eta = post_means_m,
                                 Prec.Eta = prec.m,
                                 gelman_diags = gelman_diags) 
    
  }  

  return(hyperpar_values)
  
}
