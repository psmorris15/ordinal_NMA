This directory contains R functions used to run the analysis presented in "Network meta-analysis for an ordinal outcome when outcome categorization varies across trials", as well as a workflow example using the data for the network presented in the paper. The functions supplied here follow heavily from the code given in Schmid et al. (2014) for the unordered outcome case and the BNMA R package (Seo and Schmid, 2022), as detailed below. 

liver_example.R: Runs through the analysis for the liver network from beginning to end. Shows how data should be formatted for use with the R functions provided here and works through the selection of starting values for use in the sampler, specification of values for the hyperparameters for the priors of the trial-specific baseline parameters, and the generation and running of the JAGS code used to estimate the model parameters.    

adj_cat_inits.R: Function to generate starting values to be used in JAGS sampler. Modified from the multinomial section of the network.inits function provided in the source code of the BNMA R package. Modification consists of substituting empirical logits based on the baseline-category logit for those based on the adjacent-categories logit.

adj_cat_baseline_priors.R: Function to generate values for the hyperparameters of the priors for the study-specific baseline parameters. Follows the procedure outlined in the manuscript. JAGS code used in this function modified from that used for the baseline effects model in Hu et al. (2020).  

gen_inputs_adj_cat.R: Function to generate the data list used as an input to run the JAGS sampler and the JAGS code used to specify the ordinal NMA model. Format of the output JAGS code follows that of the code given in Schmid et al. (2014). 


*** How to format outcome data for a particular network ***
If two or more adjoining levels are reported together, indicate the combined count for those levels in the left-most column belonging to the combined levels and insert NAs in the cells for the other levels. FUNCTIONS WILL NOT WORK AS INTENDED IF THIS FORMATTING IS NOT FOLLOWED; COMBINED LEVELS SHOULD BE ADJOINING UNDER ORDINAL OUTCOME FRAMEWORK.

Following from the liver dataset with categories H, A-, A, and A+: 
	Study 1 reported each of H, A-, A, and A+ separately
	Study 4 reported H and (A-, A, A+) 
	Study 6 reported H and (A-, A) and A+


*** References ***
Schmid, C.H., Trikalinos, T.A., Olkin, I.: Bayesian network meta-analysis for unordered categorical outcomes with incomplete data. Research synthesis methods 5(2), 162–185 (2014) 

Seo, M., Schmid, C.: Bnma: Bayesian Network Meta-Analysis Using ’JAGS’.
(2022). R package version 1.5.0. https://CRAN.R-project.org/package=bnma

Hu, D., O’Connor, A.M., Wang, C., Sargeant, J.M., Winder, C.B.: How to conduct a bayesian network meta-analysis. Frontiers in veterinary science 7, 271 (2020)