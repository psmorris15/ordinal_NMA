### Following code is from the network.inits function of the BNMA R package ###
### Modification: empirical logits based on adjacent-categories framework ###

adj_cat_inits <- function(Study, Treat, Outcomes, N, n.chains){
    
    Dimputed = Outcomes + 0.5
    
    logits <- matrix(data = NA, nrow = nrow(Dimputed),
                     ncol = ncol(Dimputed) - 1)
    
    ### modification: empirical logits based on adjacent-categories framework
    for(i in 2:ncol(Dimputed)){
      logits[,i-1] <- as.matrix(log(Dimputed[, i]) - log(Dimputed[, i-1]))
    }
    
    
    se.logits <- matrix(data = NA, nrow = nrow(Dimputed),
                     ncol = ncol(Dimputed) - 1)
    
    for(i in 2:ncol(Dimputed)){
      se.logits[,i-1] <- as.matrix(sqrt(1/Dimputed[, i] + 1/Dimputed[, i-1]))
    }
    
    
    nstudy <- length(unique(Study)) 
    ncat <- ncol(Dimputed)
    
    df <- data.frame(Study, Treat) %>% group_by(Study) %>% 
      summarize(baseline_treat = min(Treat)) %>% as.matrix()
    
    treat_mins <- merge(cbind(Study, Treat), df, by = "Study", 
                  sort = FALSE)
    
    b.id <- treat_mins$Treat == treat_mins$baseline_treat
    
    Eta <- se.Eta <- matrix(NA, nstudy, ncat)
    Eta[,2:ncat] <- logits[b.id,]
    se.Eta[,2:ncat] <- se.logits[b.id,]
    
    df <- data.frame(Study) %>% group_by(Study) %>% 
      summarise(n = n()) %>% as.matrix()
    
    na <- merge(data.frame(Study = unique(Study)), df, by = "Study", sort = FALSE)$n
    
    delta <- logits - apply(as.matrix(Eta[, -1]), 2, rep, times = na)
    rows.of.basetreat <- seq(dim(as.matrix(delta))[1])*as.numeric(b.id)
    delta <- delta[-rows.of.basetreat,,drop=F]   # Eliminate base treatment arms
    
    ###################### Using delta, Eta, and se.Eta make initial values
    
    ntreat <- length(unique(Treat))
    
    y <- delta            # dependent variable for regression (part of Delta)
    d <- se.d <- matrix(NA, length(unique(Treat)), ncat - 1)
    resid.var <- rep(NA, ncat -1)
    base.tx <- Treat[b.id]    # base treatment for N studies
    end.Study <- c(0, cumsum(na))  # end row number of each trial
    rows <- end.Study - seq(0, nstudy)   # end number of each trial not including base treatment arms
    design.mat <- matrix(0, sum(na) - nstudy, ntreat) # no. non-base arms x #txs
    
    for (i in seq(nstudy)){
      studytx <- Treat[(end.Study[i]+1):end.Study[i+1]]  #treatments in ith Study
      nonbase.tx <- studytx[studytx!=base.tx[i]]    #non-baseline treatments for ith Study
      design.mat[(1+rows[i]):rows[i+1],base.tx[i]] <- -1
      for (j in seq(length(nonbase.tx)))
        design.mat[j+rows[i],nonbase.tx[j]] <- 1
    }
    design.mat <- design.mat[,-1,drop=F]
    
    for(k in 1:(ncat - 1)){
      fit <- summary(lm(y[,k] ~ design.mat - 1))
      
      if(length( coef(fit)[1:(ntreat-1), 1]) == (ntreat-1)){
        d[-1,k] <- coef(fit)[1:(ntreat-1), 1]
        se.d[-1,k] <- coef(fit)[1:(ntreat-1), 2]
      } else{
        d[-1,k] <- rep(0, ntreat - 1)
        se.d[-1,k] <- rep(0.1, ntreat -1)
      }
      resid.var[k] <- fit$sigma^2
    }
    
   
    
    ################################################
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    
    for(i in 1:n.chains){
      random.Eta <- matrix(rnorm(dim(Eta)[1]*dim(Eta)[2]),dim(Eta)[1],dim(Eta)[2])
      initial.values[[i]][["Eta"]] <- Eta + se.Eta * random.Eta
    }
    
    if(!is.nan(fit$fstat[1])){
      for(i in 1:n.chains){
        random.d = matrix(rnorm(dim(d)[1]*dim(d)[2]),dim(d)[1],dim(d)[2])
        initial.values[[i]][["d"]] = d + se.d * random.d
        
      }
    }

    return(initial.values)
}