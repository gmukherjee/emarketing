rm(list=ls())
setwd("C:\\Users\\wkar\\Downloads\\EmailPaper_SWG_Final\\trialFinal\\SCRON")
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


################## load the rejected regression sites ##############
reject = read.table("reject_sites_open_parallel.txt");
reject = as.vector(reject[,1]);

reject2 = read.table("reject_sites_click_parallel.txt");
reject2 = as.vector(reject2[,1]);

reject3 = read.table("reject_sites_purchase_parallel.txt");
reject3 = as.vector(reject3[,1]);

reject = reject + 12+1; reject2 = reject2 + 12+1; reject3 = reject3 + 12+1;

###################################################################
params = read.table("SCRON_learning_output.txt");

basics <- read_excel("initial_coeff_SCRON_open.xlsx",sheet=1);
basics = basics[c(1:37,54:dim(basics)[1],38:53),];

hh = basics$Estimate;

beta_o1 = hh[c(1:87)];
beta_o = beta_o1[-reject];
rm(beta_o1);

betas_o = hh[c(88:dim(basics)[1])]; 

basics = basics[-(reject),];
predname_o = as.character(basics$Coefficients);

beta_o = params[,1:length(beta_o)];
betas_o = params[,(length(beta_o)+1):(length(beta_o)+length(betas_o))];

gg = param_fun(predname_o, beta_o, betas_o);

########### for click stage ############
basics <- read_excel("initial_coeff_SCRON_click.xlsx",sheet=1);
basics = basics[c(1:37,54:dim(basics)[1],38:53),];

hh = basics$Estimate;

beta_c1 = hh[c(1:87)];
beta_c = beta_c1[-reject2];
rm(beta_c1);

betas_c = hh[c(88:dim(basics)[1])]; 

basics = basics[-(reject2),];
predname_c = as.character(basics$Coefficients);

beta_c = params[,(length(beta_o)+length(betas_o)+1):(length(beta_o)+length(betas_o)+length(beta_c))];
betas_c = params[,(length(beta_o)+length(betas_o)+length(beta_c)+1):(length(beta_o)+length(betas_o)+
                                                                       length(beta_c)+length(betas_c))];

gg1 = param_fun(predname_c, beta_c, betas_c);

####### for purchase stage ###########
basics <- read_excel("initial_coeff_SCRON_purchase.xlsx",sheet=1);
basics = basics[c(1:37,46:dim(basics)[1],38:45),];

hh = basics$Estimate;

beta_p1 = hh[c(1:87)];
beta_p = beta_p1[-reject3];
rm(beta_p1);

betas_p = hh[c(88:dim(basics)[1])]; 

basics = basics[-(reject3),];
predname_p = as.character(basics$Coefficients);

beta_p = params[,(length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+1):(length(beta_o)+length(betas_o)+
                                                                                      length(beta_c)+length(betas_c)+length(beta_p))];
betas_p = params[,(length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+length(beta_p)+1):(length(beta_o)+length(betas_o)+
                                                                                                      length(beta_c)+length(betas_c)+length(beta_p)
                                                                                                    +length(betas_p))];


gg2 = param_fun(predname_p, beta_p, betas_p);

sig_sim = params[,(length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+length(beta_p)
                   +length(betas_p)+1):dim(params)[2]];
sig_est = matrix(apply(sig_sim, 2, mean), 3, 3);

write.csv(gg, "SCRON_parameter_open.csv", row.names=F, quote=F);

write.csv(gg1, "SCRON_parameter_click.csv", row.names=F, quote=F);

write.csv(gg2, "SCRON_parameter_purchase.csv", row.names=F, quote=F);

write.csv(sig_est, "SCRON_parameter_correlation.csv", row.names=F, quote=F);

######################### end of the code for generating components from parallel runs ####################################33

