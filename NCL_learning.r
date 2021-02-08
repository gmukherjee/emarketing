rm(list=ls())

# Please set the path to the input data and sub-routines here.
# IT takes the following input files: training_data.txt, starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx) of the focal model
# and locations of sub-sampling for the purchase stage locations_subsampling_purchase.csv
# it creates the output MODELNAME_learning_output.txt, which becomes an input in the coeff generation code, test and training prediction codes.

####### initial runs were set and output was stored in results_save2/. TMCMC formula is used to calculate scales from the
#initial output
library(matrixsampling)
library(emdbook)
library(SimDesign)

library(readxl)
library(mvtnorm)
#library(LaplacesDemon)

library(foreach)
library(doParallel)
library(doSNOW)

part_split = 1;

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
registerDoSNOW(cl)

pb <- txtProgressBar(max = part_split, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

###########################################################################
start_time <- Sys.time()



##############################
options(digits=15);
source("joint_subroutine_random_mod.R");

eps=0.00000000016;
degree = 2; #### degree of spline
########## fix initial probabilities ###################
aprob = rep(0, 2); aprob_c = rep(0, 2); aprob_p = rep(0, 2);


aprob_p[1] = 0.5; aprob_p[2] = 1.0;


############################ loading data #######################################

############ subsampling locations for fitting the model #################
locs_purchase = read.csv("locations_subsampling_purchase.csv", T);

# Training Set
total_data = read.table("training_data.txt", T);
total_data = data.frame(total_data);

pos_c = which(total_data$opened==1);
total_data_c = total_data[pos_c,];

total_data_p = total_data_c[locs_purchase$location,];

total_data = total_data[order(total_data$recipient_id),];

################## match click with purchase locations ##########################
serial_c = data.frame(recipient_id=total_data_c$recipient_id, location=1:length(total_data_c[,1]));

################ declaring output details ########################################
total_customers = unique(total_data$recipient_id);
total_customers = data.frame(recipient_id=total_customers, num=1:length(total_customers));

size = round(length(total_customers$recipient_id)/part_split, digits=0);

cutoff = seq(1, length(total_customers$recipient_id), by=size);
if(cutoff[length(cutoff)] < length(total_customers$recipient_id)){
  cutoff = append(cutoff, length(total_customers$recipient_id));
}

part_split = length(cutoff) - 1;

cutoff[2:length(cutoff)] = cutoff[2:length(cutoff)] + 1;

###########################################################################

finalMatrix <- foreach(rr=1:part_split, .combine='+', .options.snow = opts) %dopar% {
  
  set.seed(11);
  
  library(matrixsampling)
  library(emdbook)
  library(SimDesign)
  
  library(readxl)
  library(mvtnorm)
  
  options(digits=15);
  source("joint_subroutine_random_mod.R");
  
  ########## initialising parameters ###############
  mcmc_run = 8000; ##### number of MCMC run
  burn = 3000; ###### number of burn-in 
  
  ID_details = data.frame(recipient_id=total_customers$recipient_id[(cutoff[rr]:(cutoff[rr+1]-1))]);
  
  data = merge(ID_details, total_data, by=c("recipient_id"));
  data = data[order(data$recipient_id),];
  
  pos_c = which(data$opened==1);
  data_c = data[pos_c,];
  
  data = data[order(data$recipient_id),];
  
  ################ random components ##################################
  customers = unique(data$recipient_id);
  customers = data.frame(recipient_id=customers, num=1:length(customers));
  
  data = merge(data, customers, by=c("recipient_id"));
  data = data[order(data$recipient_id),];
  
  locs = data.frame(recipient_id=data$recipient_id, num=data$num);
  data$num = NULL;
  
  ################## match with purchase ############################
  pick_serial_c = merge(ID_details, serial_c, by=c("recipient_id"));
  pick_serial_c = pick_serial_c[order(pick_serial_c$recipient_id),];
  pick_serial_c$loc_num = 1:length(pick_serial_c[,1]);
  
  pick_serial_p = merge(pick_serial_c, locs_purchase, by=c("location"));
  pick_serial_p = pick_serial_p[order(pick_serial_p$recipient_id),];
  
  ################ random components ##################################
  
  customers_c = unique(data_c$recipient_id);
  customers_c = data.frame(recipient_id=customers_c, num=1:length(customers_c));
  
  data_c = merge(data_c, customers_c, by=c("recipient_id"));
  data_c = data_c[order(data_c$recipient_id),];
  
  locs_c = data.frame(recipient_id=data_c$recipient_id, num=data_c$num);
  data_c$num = NULL;
  
  ################ random components ##################################
  data_p = merge(ID_details, total_data_p, by=c("recipient_id")); 
  data_p = data_p[order(data_p$recipient_id),];
  
  customers_p = unique(data_p$recipient_id);
  customers_p = data.frame(recipient_id=customers_p, num=1:length(customers_p));
  
  
  data_p = merge(data_p, customers_p, by=c("recipient_id"));
  data_p = data_p[order(data_p$recipient_id),];
  
  locs_p = data.frame(recipient_id=data_p$recipient_id, num=data_p$num);
  data_p$num = NULL;
  
  #######################################################################
  
  N =  dim(data)[1];   #### number of data points
  N_click = dim(data_c)[1];
  N_purchase = dim(data_p)[1];
  
  ########## defining knots #########################
  knots_open = c(1, 3, 10, 15, 30, 50);  knots_click = c(1, 3, 10, 20, 30, 50); 
  knots_purchase = c(1, 3, 10, 15, 30, 50);
  
  ############# open stage #########################################
  
  income = factor(data$income);
  hh=model.matrix(~income);
  
  data = as.data.frame(cbind(data, hh[,2:dim(hh)[2]]));
  
  ############### covariate structure non-spline #####################
  
  covar_all = as.matrix(cbind(rep(1, dim(data)[1]), data[,c(7:10,21:dim(data)[2])]));
  
  ######### for click stage #############################
  
  income = factor(data_c$income);
  hh=model.matrix(~income);
  
  data_c = as.data.frame(cbind(data_c, hh[,2:dim(hh)[2]]));
  
  covar_all_c = as.matrix(cbind(rep(1, dim(data_c)[1]), data_c[,c(7:10,21:dim(data_c)[2])]));
  
  ########### for purchase stage ####################
  income = factor(data_p$income);
  hh=model.matrix(~income);
  
  data_p = as.data.frame(cbind(data_p, hh[,2:dim(hh)[2]]));
  
  covar_all_p = as.matrix(cbind(rep(1, dim(data_p)[1]), data_p[,c(7:10,21:dim(data_p)[2])]));
  
  ################ prepare the covariance structure spline ########
  ##### for open stage #############
  open_covars_o = covar_bspline_open45(trans="log", degree, knots_open, knots_purchase, data$daysSincePurchased, data$daysSinceOpened);
  
  ######## for click stage ###########
  click_covars_o = covar_bspline_open45(trans="log", degree, knots_click, knots_purchase, data_c$daysSincePurchased, data_c$daysSinceClicked);
  
  ########## for purchase stage ###########
  purchase_covars_o = covar_b.spline_logpurchase(trans="log", degree, knots_purchase, data_p$daysSincePurchased);
  
  ###################### structures for producing training probs ##############
  train_c = covar_bspline_open45(trans="log", degree, knots_open, knots_purchase, data$daysSincePurchased, data$daysSinceClicked);
  
  train_p = covar_b.spline_logpurchase(trans="log", degree, knots_purchase, data$daysSincePurchased);
  
  ####### logit parameters ##########################
  len_o = dim(covar_all)[2]; ####### number of covariates excluding the constant term in x #####
  beta_o = rep(1, len_o);
  beta_c = rep(1, len_o);
  beta_p = rep(1, len_o);
  
  ########## for open ##############
  spline_len_open = dim(open_covars_o)[2]; 
  betas_o = rep(1, spline_len_open);
  
  ########### for click ############
  spline_len_click = dim(click_covars_o)[2];
  betas_c = rep(1, spline_len_click);
  
  ######## for purchase ###########
  spline_len_purchase = dim(purchase_covars_o)[2];
  betas_p = rep(1, spline_len_purchase);
  
  ########## defining the probability terms #########
  prob_o = rep(0, N); 
  prob_c = rep(0, N_click);
  prob_p = rep(0, N_purchase);
  
  ################### initialising beta ##############
  beta_rand = matrix(0, nrow=length(customers$recipient_id), ncol=3);
  
  ########### for open stage ############
  basics <- read_excel("initial_coeff_NCL_open.xlsx",sheet=1);
  predname_o = as.character(basics$Coefficients);
  hh = basics$Estimate;
  beta_o = hh[c(1:13)];
  betas_o = hh[c(14:dim(basics)[1])]; 
  
  scales_o = basics$Std..Error;
  sig_betao = (basics$Std..Error)^2;
  mean_betao = basics$Estimate;
  
  ########### for click stage ############
  basics <- read_excel("initial_coeff_NCL_click.xlsx",sheet=1);
  predname_c = as.character(basics$Coefficients);
  hh = basics$Estimate;
  beta_c = hh[c(1:13)];
  betas_c = hh[c(14:dim(basics)[1])]; 
  
  scales_c = basics$Std..Error;
  sig_betac = (basics$Std..Error)^2;
  mean_betac = basics$Estimate;
  
  ####### for purchase stage ###########
  basics <- read_excel("initial_coeff_NCL_purchase.xlsx",sheet=1);
  predname_p = as.character(basics$Coefficients);
  hh = basics$Estimate;
  beta_p = hh[c(1:13)];
  betas_p = hh[c(14:dim(basics)[1])]; 
  
  scales_p = basics$Std..Error;
  sig_betap = (basics$Std..Error)^2;
  mean_betap = basics$Estimate;
  
  #######################################################################################
  
  mean_rand = rep(0, 3);
  sig_rand = matrix(c(5,0,0,0,5,0,0,0,5), nrow=3, ncol=3);
  scales_rand = c(0.005, 0.005, 0.005);
  
  beta_rand_exp = beta_rand[as.integer(locs$num),];
  
  beta_rand_exp_c = beta_rand_exp[pos_c,2];
  
  beta_rand_exp_pick = beta_rand_exp[pos_c,];
  beta_rand_exp_p = beta_rand_exp_pick[pick_serial_p$loc_num,3];
  
  w_param1 = 6; w_param2 = matrix(c(20,0,0,0,20,0,0,0,20), nrow=3, ncol=3);
  
  ##################################################################
  eps = 0.000000000016;
  
  scales_o = rep(0.001, length(beta_o)+length(betas_o));
  scales_c = rep(0.001, length(beta_c)+length(betas_c));
  scales_p = rep(0.001, length(beta_p)+length(betas_p));
  
  
  ############## loading the indicators ######################
  OO = data$opened; CC = data_c$clicked; PP = data_p$ordered;
  
  ######### initialise probability ########################
  ########## open stage ####################
  prob_o = logit_calc_spline_rand(N, beta_o, betas_o, beta_rand_exp[,1], covar_all, open_covars_o, eps);
  
  ########## click stage ####################
  prob_c = logit_calc_spline_rand(N_click, beta_c, betas_c, beta_rand_exp_c, covar_all_c, click_covars_o, eps);
  
  ########## click stage ####################
  prob_p = logit_calc_spline_rand(N_purchase, beta_p, betas_p, beta_rand_exp_p, covar_all_p, purchase_covars_o, eps);
  
  ############################################
  ############## loading the indicators ######################
  OO = data$opened; CC = data_c$clicked; PP = data_p$ordered;
  
  OO1 = OO; CC1 = CC; PP1 = PP;
  
  ############## TMCMC parameters determination ##################
  c_tmcmc = 0.95; l_opt = 1.715;
  
  ############# declaring output file ###########################
  total_output = rep(0, length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+length(beta_p)+length(betas_p)+length(sig_rand));
  
  
  ############## MCMC loop #######################
  for(ll in 1:mcmc_run){
    cat("sim no: ", ll, "\n");
    ############ updating beta ###################################
    
    ############### beta for X #############
    
    ######## open stage ################
    
    #set.seed(11)
    result = update_beta_x_tmcmc_model2_parallel(method="tmcmc", part_split, N, len_o, OO, covar_all, open_covars_o, beta_o, betas_o, beta_rand_exp[,1], eps, mean_betao, sig_betao, scales_o, prob_o);
    
    beta_o = result[1:length(beta_o),1];
    betas_o = result[1:length(betas_o),2];
    prob_o = result[,3];
    
    cat("end open MCMC\n");
    
    ################ click stage ##################
    
    #set.seed(11)
    result = update_beta_x_tmcmc_model2_parallel(method="tmcmc", part_split, N_click, len_o, CC, covar_all_c, click_covars_o, beta_c, betas_c, beta_rand_exp_c, eps, mean_betac, sig_betac, scales_c, prob_c);
    
    beta_c = result[1:length(beta_c),1];
    betas_c = result[1:length(betas_c),2];
    prob_c = result[,3];
    
    
    cat("end click MCMC\n");
    
    ################# purchase stage ##################
    
    #set.seed(11)
    result = update_beta_x_tmcmc_model2_parallel(method="tmcmc", part_split, N_purchase, len_o, PP, covar_all_p, purchase_covars_o, beta_p, betas_p, beta_rand_exp_p, eps, mean_betap, sig_betap, scales_p, prob_p);
    
    beta_p = result[1:length(beta_p),1];
    betas_p = result[1:length(betas_p),2];
    prob_p = result[,3];
    
    cat("end purchase MCMC\n");
    
    ##################### random components update ###################################
    
    #set.seed(11)
    result = update_rand_new2(method="tmcmc", N, N_click, N_purchase, locs, locs_c, locs_p, pos_c, pick_serial_p$loc_num, OO, CC, PP, covar_all, open_covars_o, covar_all_c, click_covars_o,
                              covar_all_p, purchase_covars_o, beta_o, betas_o, beta_c, betas_c, beta_p, betas_p, beta_rand, eps, mean_rand,
                              sig_rand, scales_rand, prob_o, prob_c, prob_p);
    
    
    beta_rand = result[1:(dim(customers)[1]*3),1];
    beta_rand = matrix(beta_rand, nrow=dim(customers)[1], ncol=3, byrow=TRUE);
    prob_o = result[1:length(prob_o),2];
    prob_c = result[1:length(prob_c),3];
    prob_p = result[1:length(prob_p),4];
    
    ###############################################################3
    beta_rand_exp = beta_rand[as.integer(locs$num),];
    
    beta_rand_exp_c = beta_rand_exp[pos_c,2];
    
    beta_rand_exp_pick = beta_rand_exp[pos_c,];
    beta_rand_exp_p = beta_rand_exp_pick[pick_serial_p$loc_num,3];
    
    ################### updating Wishart matrix #############################
    param1 = ((w_param1+3+1)/(part_split))+((dim(beta_rand)[1]-3-1));  
    
    param_prepare = apply(beta_rand, 1, function(x){dd=matrix(x, nrow=length(x), ncol=1); hh=dd%*%t(dd); return(hh)});
    add_param = apply(param_prepare, 1, FUN=sum);
    add_param = matrix(add_param, nrow=3, ncol=3);
    
    param2 = (w_param2/part_split) + add_param;     
    sig_rand_sim = rinvwishart(n=1, param1, param2);
    sig_rand = sig_rand_sim[,,1];
    
    
    ############ saving the MCMC simulations ############
    if(ll>burn){
      
      total_beta = append(append(beta_o, betas_o), append(beta_c, betas_c));
      total_beta = append(total_beta, append(beta_p, betas_p));
      total_beta = append(total_beta, sig_rand);
      
      total_output = rbind(total_output, total_beta);
      
    }
    
  } ### mcmc loop
  
  total_output = total_output[-1,];
  
  total_output
  
} #### do parallel loop


close(pb)
#stop cluster
stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

finalMatrix = finalMatrix/part_split;

write.table(finalMatrix, "NCL_learning_output.txt", col.names=F, row.names=F, quote=F);
save(finalMatrix, file = "NCL_learning_output.rda");


