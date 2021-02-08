
options(digits=15)
library(mvtnorm)
#pos_p = locs_purchase$location; xx=OO; yy=CC;zz=PP; covar_x=covar_all; covars_x = open_covars_o;
#covar_y=covar_all_c; covars_y=click_covars_o; covar_z=covar_all_p; covars_z=purchase_covars_o;
#beta_x=beta_o; betas_x=betas_o; beta_y=beta_c; betas_y=betas_c; beta_z=beta_p; betas_z = betas_p;
#prob_x=prob_o; prob_y=prob_c; prob_z=prob_p;

update_rand_new2_check <- function(method, N, N_click, N_purchase, locs, locs_c, 
                             locs_p, pos_c, pos_p, xx, yy, zz, covar_x, covars_x, covar_y, covars_y,
                             covar_z, covars_z, beta_x, betas_x, beta_y, betas_y, beta_z, betas_z, beta_rand, eps, mean_rand, sig_rand,
                             scales_rand, prob_x, prob_y, prob_z){
  
  #cat(summary(locs$num)); cat("\n");
  #cat(summary(locs_c$num)); cat("\n");
  
  #cat(sig_rand, "\n");
  #cat(dim(covar_x), "\n");
  scales_exp = rep(scales_rand, dim(beta_rand)[1]);
  scales_exp = matrix(scales_exp, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
  
  ########### create the summation frame ################################################################
  frame_open = data.frame(num=locs$num, opart1=(xx*log(prob_x)), opart2=((1-xx)*log(1-prob_x)));
  frame_open$cpart1 = 0; frame_open$cpart2=0;  frame_open$ppart1=0;  frame_open$ppart2 = 0;
  
  frame_open$cpart1[pos_c] = (yy*log(prob_y)); frame_open$cpart2[pos_c] = ((1-yy)*log(1-prob_y));
  frame_open$ppart1[pos_c[pos_p]] = (zz*log(prob_z)); 
  frame_open$ppart2[pos_c[pos_p]] = ((1-zz)*log(1-prob_z));
  
  #frame_click = data.frame(num=locs_c$num, cpart1=(yy*log(prob_y)), cpart2=((1-yy)*log(1-prob_y)));
  #frame_purchase = data.frame(num=locs_p$num, ppart1=(zz*log(prob_z)), ppart2=((1-zz)*log(1-prob_z)));
  
  agg1 = aggregate(frame_open, by=list(num1=frame_open$num), FUN=sum, na.rm=T);
  agg1$num = NULL;
  agg1 = agg1[order(agg1$num1),];
  
  prior_calc = apply(beta_rand, 1, 
                     FUN=function(x)((matrix((x-mean_rand), nrow=1, ncol=length(x))%*%solve(sig_rand))%*%matrix((x-mean_rand), nrow=length(x), ncol=1)));
  prior_calc = -(0.5*prior_calc)-log((2*pi)^(3/2)*sqrt(det(sig_rand)));
  
  total_frame = data.frame(num1=1:length(agg1$num1), 
                           prior=prior_calc);
  
  total_frame = merge(total_frame, agg1, by=c("num1"));
  total_frame = total_frame[order(total_frame$num1),];
  
  
  sum1 = total_frame$opart1+total_frame$opart2+total_frame$cpart1+total_frame$cpart2+total_frame$ppart1+
    total_frame$ppart2+total_frame$prior;
  
  if(method=="ordinary"){
    ########## new component ############
    mrkv = rmvnorm(n=dim(beta_rand)[1], mean=c(0,0,0), sigma=matrix(c(0.05,0,0,0,0.05,0,0,0,0.05), nrow=3, ncol=3));
    
    beta_rand_new = beta_rand + mrkv;
    ####################################################
  }
  
  if(method=="tmcmc"){
    ########## new component ############
    beta_rand_new = matrix(0, dim(beta_rand)[1], ncol=3);
    
    #change_ep = rep(0, length(beta_rand));
    change_ep = rnorm(n=dim(beta_rand)[1], mean=0, sd=1);
    change_ep = abs(change_ep);
    
    #z_change = rep(0, length(change_ep));
    z_change = rbinom(n=(dim(beta_rand)[1]*3), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    z_change = matrix(z_change, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
    
    beta_rand_new = beta_rand + (scales_exp*z_change*change_ep);
    ####################################################
  }
  
  beta_rand_exp_new = beta_rand_new[as.integer(locs$num),];
  
  beta_rand_exp_new_c = beta_rand_exp_new[as.integer(pos_c),2];
  
  beta_rand_exp_new_pick = beta_rand_exp_new[as.integer(pos_c),];
  beta_rand_exp_new_p = beta_rand_exp_new_pick[as.integer(pos_p),3];
  
  ########## open stage ####################
  prob_x_new = logit_calc_spline_rand(N, beta_x, betas_x, beta_rand_exp_new[,1], covar_x, covars_x, eps);
  pos1 = which(prob_x_new<eps);
  prob_x_new[pos1] = prob_x_new[pos1] + eps;
  
  ########## click stage ####################
  prob_y_new = logit_calc_spline_rand(N_click, beta_y, betas_y, beta_rand_exp_new_c, covar_y, covars_y, eps);
  pos2 = which(prob_y_new<eps);
  prob_y_new[pos2] = prob_y_new[pos2] + eps;
  
  ########## purchase stage ####################
  prob_z_new = logit_calc_spline_rand(N_purchase, beta_z, betas_z, beta_rand_exp_new_p, covar_z, covars_z, eps);
  pos3 = which(prob_z_new<eps);
  prob_z_new[pos3] = prob_z_new[pos3] + eps;
  
  pos_check = which(prob_x_new==1); pos_check2 = which(prob_y_new==1); pos_check3 = which(prob_z_new==1);
  
  if((length(pos_check)==0) & (length(pos_check2)==0) & (length(pos_check3)==0)){
    ########## new likelihood ##########
    frame_open$opart1 = (xx*log(prob_x_new));  frame_open$opart2 = ((1-xx)*log(1-prob_x_new));
    frame_open$cpart1[pos_c] = (yy*log(prob_y_new)); frame_open$cpart2[pos_c] = ((1-yy)*log(1-prob_y_new));
    frame_open$ppart1[pos_c[pos_p]] = (zz*log(prob_z_new)); 
    frame_open$ppart2[pos_c[pos_p]] = ((1-zz)*log(1-prob_z_new));
    
    agg1 = aggregate(frame_open, by=list(num1=frame_open$num), FUN=sum);
    agg1$num = NULL;
    agg1 = agg1[order(agg1$num1),];
    
    ##################################################################
    prior_calc_new = apply(beta_rand_new, 1, 
                           FUN=function(x)((matrix((x-mean_rand), nrow=1, ncol=length(x))%*%solve(sig_rand))%*%matrix((x-mean_rand), nrow=length(x), ncol=1)));
    prior_calc_new = -(0.5*prior_calc_new)-log((2*pi)^(3/2)*sqrt(det(sig_rand)));
    
    total_frame = data.frame(num1=1:length(agg1$num1), 
                             prior=prior_calc_new);
    
    total_frame = merge(total_frame, agg1, by=c("num1"));
    total_frame = total_frame[order(total_frame$num1),];
    
    sum_new = total_frame$opart1+total_frame$opart2+total_frame$cpart1+total_frame$cpart2+total_frame$ppart1+
      total_frame$ppart2+total_frame$prior;
    
    diff = sum_new - sum1;
    uu = runif(n=length(diff), 0, 1);  
    #cat("sum=", sum1, "sum_new=", sum_new, "\n");
    #cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    pos_diff = which(log(uu)<diff);
    
    if(length(pos_diff)!=0){
      beta_rand[pos_diff,] = beta_rand_new[pos_diff,];
      prob_x[pos_diff] = prob_x_new[pos_diff];
      prob_y[pos_diff] = prob_y_new[pos_diff];
      prob_z[pos_diff] = prob_z_new[pos_diff];
      
      prior_calc = prior_calc_new;
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=5);
  result[1:(dim(beta_rand)[1]*3),1] = as.vector(t(beta_rand));
  result[,2] = prob_x;
  result[1:length(prob_y),3] = prob_y;
  result[1:length(prob_z),4] = prob_z;
  result[1:dim(beta_rand)[1],5] = prior_calc;
  
  return(result);
  
}


update_rand_new2_joint <- function(method, N, N_click, N_purchase, locs, locs_c, 
                                   locs_p, pos_c, pos_p, xx, yy, zz, covar_x, covars_x, covar_y, covars_y,
                                   covar_z, covars_z, beta_x, betas_x, beta_y, betas_y, beta_z, betas_z, beta_rand, eps, mean_rand, sig_rand,
                                   scales_rand, prob_x, prob_y, prob_z){
  
  #cat(summary(locs$num)); cat("\n");
  #cat(summary(locs_c$num)); cat("\n");
  
  #cat(sig_rand, "\n");
  #cat(dim(covar_x), "\n");
  scales_exp = rep(scales_rand, dim(beta_rand)[1]);
  scales_exp = matrix(scales_exp, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
  
  ########### create the summation frame ################################################################
  old = sum(xx*log(prob_x)) + sum((1-xx)*log(1-prob_x)) + sum(yy*log(prob_y)) + sum((1-yy)*log(1-prob_y));
  old = old + sum(zz*log(prob_z)) + sum((1-zz)*log(1-prob_z));
    
  old = old + sum(dmvnorm(beta_rand, mean_rand, sigma=matrix(sig_rand, nrow=3, ncol=3), log=FALSE));
  
  
  if(method=="ordinary"){
    ########## new component ############
    mrkv = rmvnorm(n=dim(beta_rand)[1], mean=c(0,0,0), sigma=matrix(c(0.05,0,0,0,0.05,0,0,0,0.05), nrow=3, ncol=3));
    
    beta_rand_new = beta_rand + mrkv;
    ####################################################
  }
  
  if(method=="tmcmc"){
    ########## new component ############
    beta_rand_new = matrix(0, dim(beta_rand)[1], ncol=3);
    
    #change_ep = rep(0, length(beta_rand));
    change_ep = rnorm(n=dim(beta_rand)[1], mean=0, sd=1);
    change_ep = abs(change_ep);
    
    #z_change = rep(0, length(change_ep));
    z_change = rbinom(n=(dim(beta_rand)[1]*3), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    z_change = matrix(z_change, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
    
    beta_rand_new = beta_rand + (scales_exp*z_change*change_ep);
    ####################################################
  }
  
  beta_rand_exp_new = beta_rand_new[as.integer(locs$num),];
  
  beta_rand_exp_new_c = beta_rand_exp_new[as.integer(pos_c),2];
  
  beta_rand_exp_new_pick = beta_rand_exp_new[as.integer(pos_c),];
  beta_rand_exp_new_p = beta_rand_exp_new_pick[as.integer(pos_p),3];
  
  ########## open stage ####################
  prob_x_new = logit_calc_spline_rand(N, beta_x, betas_x, beta_rand_exp_new[,1], covar_x, covars_x, eps);
  pos1 = which(prob_x_new<eps);
  prob_x_new[pos1] = prob_x_new[pos1] + eps;
  
  ########## click stage ####################
  prob_y_new = logit_calc_spline_rand(N_click, beta_y, betas_y, beta_rand_exp_new_c, covar_y, covars_y, eps);
  pos2 = which(prob_y_new<eps);
  prob_y_new[pos2] = prob_y_new[pos2] + eps;
  
  ########## purchase stage ####################
  prob_z_new = logit_calc_spline_rand(N_purchase, beta_z, betas_z, beta_rand_exp_new_p, covar_z, covars_z, eps);
  pos3 = which(prob_z_new<eps);
  prob_z_new[pos3] = prob_z_new[pos3] + eps;
  
  pos_check = which(prob_x_new==1); pos_check2 = which(prob_y_new==1); pos_check3 = which(prob_z_new==1);
  
  if((length(pos_check)==0) & (length(pos_check2)==0) & (length(pos_check3)==0)){
    ########## new likelihood ##########
    new_sum = sum(xx*log(prob_x_new)) + sum((1-xx)*log(1-prob_x_new)) + sum(yy*log(prob_y_new)) +
      sum((1-yy)*log(1-prob_y_new)) + sum(zz*log(prob_z_new)) + sum((1-zz)*log(1-prob_z_new));
    
    new_sum = new_sum + sum(dmvnorm(beta_rand_new, mean_rand, sigma=matrix(sig_rand, nrow=3, ncol=3), log=FALSE));
    
    diff = new_sum - old;
    uu = runif(n=1, 0, 1);  
    #cat("sum=", sum1, "sum_new=", sum_new, "\n");
    #cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    
    
    if(log(uu)<diff){
      beta_rand = beta_rand_new;
      prob_x = prob_x_new;
      prob_y = prob_y_new;
      prob_z = prob_z_new;
      
      
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=5);
  result[1:(dim(beta_rand)[1]*3),1] = as.vector(t(beta_rand));
  result[,2] = prob_x;
  result[1:length(prob_y),3] = prob_y;
  result[1:length(prob_z),4] = prob_z;
  
  return(result);
  
}

