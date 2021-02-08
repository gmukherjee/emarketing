rm(list=ls())

# Please set the path to the input data and sub-routines here.
# IT takes the following input files: starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx)
# and output from the learning code for the focal model signified as MODELNAME_learning_output.txt 
# It creates the coefficients for open, click, and purchase stages (MODELNAME_parameter_stagename.xlsx) as well as the correlation matrix (MODELNAME_parameter_correlation.xlsx)

##############################################################################

param_fun <- function(predname_o, beta_o, betas_o){
  
  ave_beta_o = apply(beta_o, 2, FUN=mean);
  
  ave_betas_o = apply(betas_o, 2, FUN=mean);
  
  ave_beta_o_sd = apply(beta_o, 2, FUN=sd);
  
  ave_betas_o_sd = apply(betas_o, 2, FUN=sd);
  
  ave_beta_o_low = apply(beta_o, 2, FUN=function(x){quantile(x, 0.025)});
  
  ave_betas_o_low = apply(betas_o, 2, FUN=function(x){quantile(x, 0.025)});
  
  ave_beta_o_up = apply(beta_o, 2, FUN=function(x){quantile(x, 0.975)});
  
  ave_betas_o_up = apply(betas_o, 2, FUN=function(x){quantile(x, 0.975)});
  
  all_model4_open = data.frame(Coeffs=predname_o, value=append(ave_beta_o, ave_betas_o), 
                               sd=append(ave_beta_o_sd, ave_betas_o_sd), low=append(ave_beta_o_low, ave_betas_o_low), 
                               up=append(ave_beta_o_up, ave_betas_o_up));
  
  return(all_model4_open);
  
  
}

library(readxl)


###################################################################
params = read.table("NCL_learning_output.txt");

basics <- read_excel("initial_coeff_NCL_open.xlsx",sheet=1);
predname_o = as.character(basics$Coefficients);
hh = basics$Estimate;
beta_o = hh[c(1:13)];
betas_o = hh[c(14:dim(basics)[1])]; 

beta_o = params[,1:length(beta_o)];
betas_o = params[,(length(beta_o)+1):(length(beta_o)+length(betas_o))];

gg = param_fun(predname_o, beta_o, betas_o);

########### for click stage ############
basics <- read_excel("initial_coeff_NCL_click.xlsx",sheet=1);
predname_c = as.character(basics$Coefficients);
hh = basics$Estimate;
beta_c = hh[c(1:13)];
betas_c = hh[c(14:dim(basics)[1])];

beta_c = params[,(length(beta_o)+length(betas_o)+1):(length(beta_o)+length(betas_o)+length(beta_c))];
betas_c = params[,(length(beta_o)+length(betas_o)+length(beta_c)+1):(length(beta_o)+length(betas_o)+
                                                                       length(beta_c)+length(betas_c))];

gg1 = param_fun(predname_c, beta_c, betas_c);

####### for purchase stage ###########
basics <- read_excel("initial_coeff_NCL_purchase.xlsx",sheet=1);
predname_p = as.character(basics$Coefficients);
hh = basics$Estimate;
beta_p = hh[c(1:13)];
betas_p = hh[c(14:dim(basics)[1])]; 

beta_p = params[,(length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+1):(length(beta_o)+length(betas_o)+
                                                                                      length(beta_c)+length(betas_c)+length(beta_p))];
betas_p = params[,(length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+length(beta_p)+1):(length(beta_o)+length(betas_o)+
                                                                                                      length(beta_c)+length(betas_c)+length(beta_p)
                                                                                                    +length(betas_p))];


gg2 = param_fun(predname_p, beta_p, betas_p);

sig_sim = params[,(length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+length(beta_p)
                   +length(betas_p)+1):dim(params)[2]];
sig_est = matrix(apply(sig_sim, 2, mean), 3, 3);

write.csv(gg, "NCL_parameter_open.csv", row.names=F, quote=F);

write.csv(gg1, "NCL_parameter_click.csv", row.names=F, quote=F);

write.csv(gg2, "NCL_parameter_purchase.csv", row.names=F, quote=F);

write.csv(sig_est, "NCL_parameter_correlation.csv", row.names=F, quote=F);

######################### end of the code for generating components from parallel runs ####################################33

