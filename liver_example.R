library(tidyverse)
library(rjags)

rm(list = ls())

set.seed(641394)

source("adj_cat_inits.R")
source("adj_cat_baseline_priors.R")
source("gen_inputs_adj_cat.R")

# read in data
liver_data <- read.csv("liver_data.csv")

###### create objects to be used in functions; Each row/element corresponds to a single trial arm ######

## Data format follows that used when running functions in the BNMA R package (see readme file) ##

# Study: integer value indicating which trial each arm belongs to
Study <- liver_data$Study

# Treat: integer value indicating which intervention was applied to each arm
Treat <- liver_data$Treat

# Outcomes: matrix of event counts for the categories/levels measured for each arm
#   If two or more adjoining levels are reported together, indicate the combined count for those levels in the 
#   left-most column belonging to the combined levels and insert NAs in the cells for the other levels.
#   FUNCTIONS WILL NOT WORK AS INTENDED IF THIS FORMATTING IS NOT FOLLOWED; 
#   COMBINED LEVELS SHOULD BE ADJOINING UNDER ORDINAL OUTCOME FRAMEWORK.

#   Following from the liver dataset with categories H, A-, A, and A+: 
#   Study 1 reported each of H, A-, A, and A+ separately
#   Study 4 reported H and (A-, A, A+) 
#   Study 6 reported H and (A-, A) and A+
Outcomes <- as.matrix(liver_data[, c("H", "Aminus", "A", "Aplus")])

# N: Total sample size for each arm
N <- rowSums(Outcomes, na.rm = TRUE)


###### determine outcome categorization for each arm ######
# This will be used when filling in the outcome data which is necessary to 
# generate starting values

# matrix showing the unique categorizations present in the Outcomes matrix 
unique_cats <- unique(!is.na(Outcomes))*1

# vector indicating the outcome categorization for each arm; integer value that
#   corresponds to the row number of unique_cats
outcome_cat <- numeric(length = nrow(Outcomes))

for(i in 1:nrow(Outcomes)){
  outcome_cat[i] <- which(apply(unique_cats, 1, function(x) return(all(x == !is.na(Outcomes[i,])))))
}


###### determine starting values ######

### Must fill in outcome data to generate starting values ###
# The best way to do this can depend on the dataset and preferences of the analyst
# Therefore we do not provide a function to do this, but provide an example of how
# it was done for our dataset.

# use only complete studies for simplicity 
complete_data <- liver_data[which(outcome_cat == 1), ]

# calculate observed proportions of counts that fall into A- vs A 
complete_data$Aminus_2_prop <- with(complete_data, Aminus/(Aminus + A))
complete_data$A_2_prop <- with(complete_data, A/(Aminus + A))

# calculate observed proportions of counts that fall into A- vs A vs A+
complete_data$Aminus_3_prop <- with(complete_data, Aminus/(Aminus + A + Aplus))
complete_data$A_3_prop <- with(complete_data, A/(Aminus + A + Aplus))
complete_data$Aplus_3_prop <- with(complete_data, Aplus/(Aminus + A + Aplus))

# by intervention, calculate average proportions for A- vs A and A- vs A vs A+
imputatation_means <- complete_data %>% group_by(Treat) %>% 
  summarise(Aminus_2_prop = mean(Aminus_2_prop),
            A_2_prop = mean(A_2_prop),
            Aminus_3_prop = mean(Aminus_3_prop),
            A_3_prop = mean(A_3_prop),
            Aplus_3_prop = mean(Aplus_3_prop)) %>% as.data.frame()

# fill in outcome data using proportions calculated above 
liver_data_imputed <- liver_data

for(i in 1:nrow(liver_data_imputed)){
  if(outcome_cat[i] == 3){
    total <- liver_data_imputed[i, "Aminus"]
    Aminus <- round(imputatation_means[which(imputatation_means$Treat == liver_data_imputed[i, "Treat"]),
                                       "Aminus_2_prop"]*total)
    A <- total - Aminus
    liver_data_imputed[i, "Aminus"] <- Aminus
    liver_data_imputed[i, "A"] <- A
  }
  if(outcome_cat[i] == 2){
    total <- liver_data_imputed[i, "Aminus"]
    
    Aminus <- round(imputatation_means[which(imputatation_means$Treat == liver_data_imputed[i, "Treat"]),
                                       "Aminus_3_prop"]*total)
    A <- round(imputatation_means[which(imputatation_means$Treat == liver_data_imputed[i, "Treat"]),
                                  "A_3_prop"]*total)
    
    Aplus <- total - Aminus - A
    
    liver_data_imputed[i, c("Aminus")] <- Aminus
    liver_data_imputed[i, c("A")] <- A
    liver_data_imputed[i, c("Aplus")] <- Aplus
    
  }
}

# create the filled in Outcomes matrix
Outcomes_imputed <- as.matrix(liver_data_imputed[, c("H", "Aminus", "A", "Aplus")])

# use the adj_cat_inits function to generate starting values for trial-specific 
# baseline parameters (Eta) and the comparative effect parameters (d)
# Returns a list  consisting of matrices of starting values; number of elements
# in list corresponds to number of chains
inits <- adj_cat_inits(Study = Study, 
                       Treat = Treat, 
                       Outcomes = Outcomes_imputed, 
                       N = N, 
                       n.chains = 4)


###### determine hyperparameter values for trial-specific baseline parameters ######

### determine number of trial-specific baseline interventions ###

# will be used when determining values for hyperparameters, as we need separate 
#   values for each trial-specific baseline; trial-specific baseline is the
#   lowest-valued intervention in the trial as specific by the integer given in Treat 
 
num_baseline_Treats <- as.data.frame(cbind(Study, Treat)) %>% 
  group_by(Study) %>% summarise(baseline_Treat = min(Treat)) %>%
  as.data.frame() %>% 
  select(baseline_Treat) %>% unique() %>% length()


### optional : set seeds for JAGS samplers used in procedure to determine 
#   hyperparameter values for trial-specific baseline parameters ###

# number of chains
n.chains <- 4

hyperpar_seeds <- list()

# each element of list corresponds to a trial-specific baseline intervention
for(i in 1:length(num_baseline_Treats)){
  
  # matrix with nrow = number of chains and ncol = number of level pairings (3 for our dataset)  
  hyperpar_seeds[[i]] <- matrix(nrow = n.chains, ncol = ncol(Outcomes) - 1,
         data = sample(1:1000000, size = n.chains*(ncol(Outcomes) - 1),
                       replace = FALSE))
  
}  

# adj_cat_baseline_priors returns a list where the number of elements corresponds to
# number of trial-specific baseline parameters
# Each element of list consists of: 
#   vector mea.Eta of length ncol(Outcomes) - 1 consisting of calculated values for the means 
#     of the priors of the Eta parameters; first element is for Eta corresponding to 
#     category 2 vs category 1, second element for category 3 vs category 2, etc.
#   matrix Prec.Eta corresponding to the precision matrix of the prior for Eta of 
#     dimension ncol(Outcomes) - 1
#   Gelman-Rubin convergence diagnostics for each JAGS sampler ran to calculate values
#     for the m and prec.m
hyperpar_values <- adj_cat_baseline_priors(Outcomes, Treat, Study, nadapt = 10000, 
                                           nburnin = 10000, niter = 10000, n.chains = 4, 
                                           seeds = hyperpar_seeds)

###### fit NMA model using JAGS ######

# given the data objects and the hyperpar_values list, gen_inputs_adj_cat returns:
#   a list containing the data needed to run the JAGS sampler
#   a JAGS code file that will run the ordinal NMA model
# JAGS code is formatted based on that provided for the example dataset in 
#   Schmid et al. (2014)
jags_inputs <- gen_inputs_adj_cat(Study = Study, Treat = Treat, 
                                     Outcomes = Outcomes, N = N, 
                                     hyperpar_values = hyperpar_values)


### optional : set seeds for JAGS sampler used in estimating model parameters ###
#   Number of elements corresponding to number of chains
seeds <- sample(1:1000000, size = n.chains, replace = FALSE)

for(k in 1:n.chains){
  
  inits[[k]]$.RNG.seed <- seeds[k]
  inits[[k]]$.RNG.name <- "base::Wichmann-Hill"
  
}

# initialize JAGS model
load.module("glm")
mod_adj_cat = rjags::jags.model(textConnection(jags_inputs$jags.code),
                                data = jags_inputs$jags.data,
                                inits = inits,
                                n.chains=4,
                                n.adapt=0)


# adaptation phase; this chunk of code from the network.run function of the BNMA R package
adapted <- FALSE
count <- 0
setsize <- 10000
while(!adapted){
  adapted <- rjags::adapt(mod_adj_cat, setsize, end.adaptation = FALSE)
  count <- count + 1
  print(count)
  if(count == 100){
    stop("algorithm has not adapted")
  }
}

# end adaptation phase
rjags::adapt(mod_adj_cat, n.iter = 0, end.adaptation = TRUE)

# burn-in phase
update(mod_adj_cat, n.iter = 50000)

# sampling phase
res_adj_cat = coda.samples(mod_adj_cat, c("d", "Eta"), 50000)

#Gelman-Rubin convergence diagnostic
gelman.diag(res_adj_cat, multivariate = FALSE)

#examine summary statistics from MCMC run
summary(res_adj_cat)
