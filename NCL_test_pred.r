rm(list=ls())

# Please set the path to the input data and sub-routines here.
# IT takes the following input files: training_data.txt, starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx) of the focal model,
# locations of sub-sampling for the purchase stage locations_subsampling_purchase.csv, expanded form of the test data predict_log.csv and 
# output from the learning code MODELNAME_learning_output.txt
# it creates the test prediction for each of the three stages and stores them in MODELNAME_test_pred.txt.

cojoin <- function(X, Y){
  
  ZZ = list(rbind(X[[1]], Y[[1]]), rbind(X[[2]], Y[[2]]), rbind(X[[3]], Y[[3]]));
  
  return(ZZ);
}


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

eps = 0.000000000016;

############# assigning first stage MCMC values ####################

first_mcmc = read.table("NCL_learning_output.txt");

sep_beta_o = first_mcmc[,1:(length(beta_o)+length(betas_o))];

sep_beta_c = first_mcmc[,(length(beta_o)+length(betas_o)+1):(length(beta_o)+length(betas_o)+
                                                               length(beta_c)+length(betas_c))];

sep_beta_p = first_mcmc[,(length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+1):
                          (length(beta_o)+length(betas_o)+length(beta_c)+length(betas_c)+length(beta_p)+length(betas_p))];

sep_sigma = first_mcmc[,(length(beta_o)+length(betas_o)+length(beta_c)+
                           length(betas_c)+length(beta_p)+length(betas_p)+1):dim(first_mcmc)[2]];

###################################################
beta_o_sim = sep_beta_o[,1:length(beta_o)];
betas_o_sim = sep_beta_o[,(length(beta_o)+1):dim(sep_beta_o)[2]];

beta_c_sim = sep_beta_c[,1:length(beta_c)];
betas_c_sim = sep_beta_c[,(length(beta_c)+1):dim(sep_beta_c)[2]];

beta_p_sim = sep_beta_p[,1:length(beta_p)];
betas_p_sim = sep_beta_p[,(length(beta_p)+1):dim(sep_beta_p)[2]];

rm(sep_beta_c); rm(sep_beta_o); rm(sep_beta_p);

############## TMCMC parameters determination ##################
c_tmcmc = 0.95; l_opt = 1.715;


############ loading the data for prediction #######################
data_pred_all = read.csv("predict_log.csv", T);

##########################################################
pred_id = unique(data_pred_all$recipient_id);
pred_id = data.frame(recipient_id=pred_id, num=1:length(pred_id));

pred_param = data.frame(recipient_id=data_pred_all$recipient_id);

locs_pred = merge(pred_param, pred_id, by=c("recipient_id"));
locs_pred = locs_pred[order(locs_pred$recipient_id),];

################ declaring output details ########################################
total_customers_pred = unique(locs_pred$recipient_id);
total_customers_pred = data.frame(recipient_id=total_customers_pred, num=1:length(total_customers_pred));

size_pred = round(length(total_customers_pred$recipient_id)/part_split, digits=0);

cutoff_pred = seq(1, length(total_customers_pred$recipient_id), by=size_pred);
if(cutoff_pred[length(cutoff_pred)] < length(total_customers_pred$recipient_id)){
  cutoff_pred = append(cutoff_pred, length(total_customers_pred$recipient_id));
}

cutoff_pred[2:length(cutoff_pred)] = cutoff_pred[2:length(cutoff_pred)] + 1;

###########################################################################
start_time <- Sys.time()

finalMatrix <- foreach(rr=1:part_split, .combine='rbind', .options.snow = opts) %dopar% {
  
  set.seed(11);
  
  ########## initialising parameters ###############
  mcmc_run = 8000; ##### number of MCMC run
  burn = 3000; ###### number of burn-in 
  
  stage2_run = 1;
  
  library(matrixsampling)
  library(emdbook)
  library(SimDesign)
  
  library(readxl)
  library(mvtnorm)
  
  options(digits=15);
  source("joint_subroutine_random_mod.R");
  
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
  
  ############## loading the indicators ######################
  OO = data$opened; CC = data_c$clicked; PP = data_p$ordered;
  
  OO1 = OO; CC1 = CC; PP1 = PP;
  
  
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
  
  ############### making the prediction data set ###########################
  ID_pred = data.frame(recipient_id=total_customers_pred$recipient_id[(cutoff_pred[rr]:(cutoff_pred[rr+1]-1))]);
  
  data_pred = merge(ID_pred, data_pred_all, by=c("recipient_id"));
  data_pred = data_pred[order(data_pred$recipient_id),];
  
  data_pred$facID = as.factor(data_pred$ID);
  
  ff = model.matrix(~facID:type_high + facID:type_mid, data=data_pred);
  
  data_pred = as.data.frame(cbind(data_pred, as.data.frame(ff)));
  
  
  data_pred_o = data_pred[,predname_o];
  data_pred_o = as.matrix(data_pred_o);
  
  data_pred_c = data_pred[,predname_c];
  data_pred_c = as.matrix(data_pred_c);
  
  data_pred_p = data_pred[,predname_p];
  data_pred_p = as.matrix(data_pred_p);
  
  ##########################################################
  pred_id = unique(data_pred$recipient_id);
  pred_id = data.frame(recipient_id=pred_id, num=1:length(pred_id));
  
  pred_param = data.frame(recipient_id=data_pred$recipient_id);
  
  locs_pred = merge(pred_param, pred_id, by=c("recipient_id"));
  locs_pred = locs_pred[order(locs_pred$recipient_id),];
  
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
  
  
  ###########################################
  OO1 = OO; CC1 = CC; PP1 = PP;
  
  ############## TMCMC parameters determination ##################
  c_tmcmc = 0.95; l_opt = 1.715;
  
  ############## MCMC loop #######################
  
  beta_join_sim = cbind(beta_o_sim, betas_o_sim);
  prob_o_pred = (data_pred_o%*%t(beta_join_sim));
  
  beta_join_sim = cbind(beta_c_sim, betas_c_sim);
  prob_c_pred = (data_pred_c%*%t(beta_join_sim));
  
  beta_join_sim = cbind(beta_p_sim, betas_p_sim);
  prob_p_pred = (data_pred_p%*%t(beta_join_sim));
  
  rm(beta_join_sim);
  
  ############# declaring output file ###########################
  
  mat_pred_o = 0;
  mat_pred_c = 0;
  mat_pred_p = 0;
  
  for(ll in 1:(mcmc_run-burn)){
    cat("sim no: ", ll, "\n");
    ############ loading simulated beta ###################################
    
    sig_rand = matrix(as.numeric(sep_sigma[ll,]), nrow=3, ncol=3);
    
    ############## for predictions ###########################################
    beta_rand_pred = rmvnorm(n=dim(pred_id)[1], mean=mean_rand, sigma=sig_rand);
    
    beta_rand_pred_exp = beta_rand_pred[locs_pred$num,];
    
    ############################################################
    
    join_mat = cbind(prob_o_pred[,ll], prob_c_pred[,ll], prob_p_pred[,ll]);
    
    ##########################################################
    
    beta_rand_pred_exp[,3] = beta_rand_pred_exp[,3] - log(aprob_p[2]/aprob_p[1]);
    
    join_mat1 = beta_rand_pred_exp + join_mat;
    
    join_mat1 = exp(join_mat1)/(1+exp(join_mat1));
    
    mat_pred_o = mat_pred_o + join_mat1[,1];
    
    mat_pred_c = mat_pred_c + join_mat1[,2];
    
    mat_pred_p = mat_pred_p + join_mat1[,3];
    
  } ### mcmc loop
  
  mat_pred_o = mat_pred_o/(mcmc_run-burn); mat_pred_c = mat_pred_c/(mcmc_run-burn);
  mat_pred_p = mat_pred_p/(mcmc_run-burn);
  
  
  whole_result = cbind(data_pred$opened, data_pred$clicked, data_pred$ordered,
                       mat_pred_o, mat_pred_c, mat_pred_p);
  
  
  whole_result
  
  
} #### do parallel loop


close(pb)
#stop cluster
stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

finalMatrix = data.frame(opened=finalMatrix[,1], clicked=finalMatrix[,2], ordered=finalMatrix[,3],
                         pred_est_o=finalMatrix[,4], pred_est_c=finalMatrix[,5], pred_est_p=finalMatrix[,6]);


write.table(finalMatrix, "NCL_test_pred.txt", col.names=T, row.names=F, quote=F);


