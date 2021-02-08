##############################
options(digits=22);

covar_bspline_open45 <- function(trans, degree, knots, knots_p, purchase, covar){
  
  no_var = length(covar);
  no_knots = length(knots);
  no_knots_p = length(knots_p);
  
  prepare = rep(1, length(covar));
  if(trans=="log"){
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (log(covar+1))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots){
      ff = (log(covar+1)-log(knots[ll]+1));
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
    
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (log(purchase+1))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = (log(purchase+1)-log(knots_p[ll]+1));
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="ordinary"){
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (covar)^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots){
      ff = (covar-knots[ll]);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
    
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (purchase)^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = (purchase-knots_p[ll]);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="std10"){
    std_div1 = 10;
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((covar/std_div1))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots){
      ff = ((covar-knots[ll])/std_div1);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
    
    std_div2 = 10;
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = ((purchase-knots_p[ll])/std_div2);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="stdmax"){
    std_div1 = max(covar);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((covar/std_div1))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots){
      ff = ((covar-knots[ll])/std_div1);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
    
    std_div2 = max(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = ((purchase-knots_p[ll])/std_div2);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="stdsd"){
    std_div1 = sd(covar);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((covar/std_div1))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots){
      ff = ((covar-knots[ll])/std_div1);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
    
    std_div2 = sd(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = ((purchase-knots_p[ll])/std_div2);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  
  prepare = prepare[,-1];
  
  return(prepare);
}

covar_b.spline_model3 <- function(trans, degree, purchase, covar){
  
  no_var = length(covar);
  
  prepare = rep(1, length(covar));
  if(trans=="log"){
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (log(covar+1))^kk);
      kk = kk+1;
    }
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (log(purchase+1))^kk);
      kk = kk+1;
    }
  }
  if(trans=="ordinary"){
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (covar)^kk);
      kk = kk+1;
    }
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (purchase)^kk);
      kk = kk+1;
    }
  }
  if(trans=="stdmax"){
    std_div1 = max(covar);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((covar/std_div1))^kk);
      kk = kk+1;
    }
    std_div2 = max(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
  }
  if(trans=="stdsd"){
    std_div1 = sd(covar);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((covar/std_div1))^kk);
      kk = kk+1;
    }
    std_div2 = sd(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
  }
  
  prepare = prepare[,-1];
  
  return(prepare);
}

covar_b.spline_model3_purchase <- function(trans, degree, purchase){
  
  no_var = length(purchase);
  
  prepare = rep(1, length(purchase));
  if(trans=="log"){
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (log(purchase+1))^kk);
      kk = kk+1;
    }
  }
  if(trans=="ordinary"){
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (purchase)^kk);
      kk = kk+1;
    }
  }
  if(trans=="stdmax"){
    std_div2 = max(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
  }
  if(trans=="stdsd"){
    std_div2 = sd(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
  }
  
  prepare = prepare[,-1];
  
  return(prepare);
}

covar_b.spline_logpurchase <- function(trans, degree, knots_p, purchase){
  
  no_knots_p = length(knots_p);
  
  prepare = rep(1, length(purchase));
  if(trans=="log"){
    
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (log(purchase+1))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = (log(purchase+1)-log(knots_p[ll]+1));
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="ordinary"){
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, (purchase)^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = (purchase-knots_p[ll]);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="std10"){
    std_div2 = 10;
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = ((purchase-knots_p[ll])/std_div2);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="stdmax"){
    std_div2 = max(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = ((purchase-knots_p[ll])/std_div2);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  if(trans=="stdsd"){
    std_div2 = sd(purchase);
    kk=1;
    while(kk<=degree){
      prepare = cbind(prepare, ((purchase/std_div2))^kk);
      kk = kk+1;
    }
    for(ll in 1:no_knots_p){
      ff = ((purchase-knots_p[ll])/std_div2);
      pos = which(ff<0);
      ff[pos] = 0;
      
      prepare = cbind(prepare, (ff)^degree);
    }
  }
  
  prepare = prepare[,-1];
  
  return(prepare);
}

logit_calc_spline <- function(N, beta_v, betas_v, cov, covs, eps){
  
  sum = regression_spline(beta_v, betas_v, cov, covs); 
  
  prob = exp(sum)/(1+exp(sum));         
  cat("beta_v=", max(beta_v), "betas_v=", max(betas_v), "\n");
  cat("check for logit function=", max(sum), "logit=", max(prob), "\n");
  
  pos = which(prob<eps);
  prob[pos] = prob[pos] + eps;
  
  return(prob);
}

logit_calc_spline_rand <- function(N, beta_v, betas_v, rand_v, cov, covs, eps){
  
  sum = regression_spline(beta_v, betas_v, cov, covs);
  sum = sum + rand_v;
  
  prob = exp(sum)/(1+exp(sum));         
  cat("beta_v=", max(beta_v), "betas_v=", max(betas_v), "\n");
  cat("check for logit function=", max(sum), "logit=", max(prob), "\n");
  
  pos = which(prob<eps);
  prob[pos] = prob[pos] + eps;
  
  return(prob);
}

regression_spline <-function(beta_x, betas_x, covar, covars, eps){
  
  ##beta_all = append(beta_o, beta_s);
  sum = (covar%*%beta_x)+(covars%*%betas_x);
  
  #cat("all part=", max(covar%*%beta_x), "spline part=", max(covars%*%betas_x), "\n");
  
  return(sum);
}

update_beta_x_tmcmc <- function(N, len_x, XX, covar_x, covars_x, beta_x, betas_x, eps, mean_betax, sig_betax, scales_r, prob_x){
  
  betax = append(beta_x, betas_x);
  ######## old component #############
  s1 = XX*log(prob_x); 
  s2 = (1-XX)*log(1-prob_x);
  
  sum = sum(s1) + sum(s2);
  sum = sum - sum((0.5*(betax-mean_betax)^2/sig_betax));
  
  ########## new component ############
  betax_new = rep(0, length(betax));
  beta_x_new = rep(0, length(beta_x));
  betas_x_new = rep(0, length(betas_x));
  
  change_ep = rep(0, length(beta_x));
  change_ep = rnorm(n=1, mean=0, sd=1);
  change_ep = abs(change_ep);
  
  z_change = rep(0, length(change_ep));
  z_change = rbinom(n=length(beta_x), size=1, prob=0.5);
  z_change[z_change==0] = -1;
  beta_x_new = beta_x + (scales_r[1:length(beta_x)]*z_change*change_ep);
  ####################################################
  
  change_ep = rep(0, length(betas_x));
  change_ep = rnorm(n=1, mean=0, sd=1);
  change_ep = abs(change_ep);
  
  z_change = rep(0, length(change_ep));
  z_change = rbinom(n=length(betas_x), size=1, prob=0.5);
  z_change[z_change==0] = -1;
  betas_x_new = betas_x + (scales_r[(length(beta_x)+1):length(betax)]*z_change*change_ep);
  
  betax_new = append(beta_x_new, betas_x_new);
  
  prob_x_new = logit_calc_spline(N, beta_x_new, betas_x_new, covar_x, covars_x, eps);
  pos1 = which(prob_x_new<eps);
  prob_x_new[pos1] = prob_x_new[pos1] + eps;
  
  pos_check = which(prob_x_new==1);
  if(length(pos_check)==0){
    ########## new likelihood ##########
    s1_new = XX*log(prob_x_new);
    s2_new = (1-XX)*log(1-prob_x_new);
    
    sum_new = sum(s1_new) + sum(s2_new);
    sum_new = sum_new - sum((0.5*(betax_new-mean_betax)^2/sig_betax));
    
    diff = sum_new - sum;
    uu = runif(n=1, 0, 1);  
    cat("sum=", sum, "sum_new=", sum_new, "\n");
    cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    if(log(uu)<diff){
      beta_x = beta_x_new;
      betas_x = betas_x_new;
      prob_x = prob_x_new;
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=3);
  result[1:length(beta_x),1] = beta_x;
  result[1:length(betas_x),2] = betas_x;
  result[,3] = prob_x;
  
  return(result)
}


update_beta_x_tmcmc_model2_parallel <- function(method, S_part, N, len_x, XX, covar_x, covars_x, beta_x, betas_x, rand, eps, mean_betax, sig_betax, scales_r, prob_x){
  
  betax = append(beta_x, betas_x);
  ######## old component #############
  s1 = XX*log(prob_x); 
  s2 = (1-XX)*log(1-prob_x);
  
  sum = sum(s1) + sum(s2);
  sum = sum - (sum((0.5*(betax-mean_betax)^2/sig_betax))/S_part);
  
  if(method=="ordinary"){
    betax_new = rep(0, length(betax));
    beta_x_new = rep(0, length(beta_x));
    betas_x_new = rep(0, length(betas_x));
    
    #beta_x_new = rnorm(n=length(beta_x_new), mean=mean_betax[1:length(beta_x_new)], sd=sqrt(sig_betax[1:length(beta_x_new)]));
    #betas_x_new = rnorm(n=length(betas_x_new), mean=mean_betax[(length(beta_x_new)+1):length(betax)], sd=sqrt(sig_betax[(length(beta_x_new)+1):length(betax)]));
    
    beta_x_new = rnorm(n=length(beta_x_new), mean=mean_betax[1:length(beta_x_new)], sd=scales_r[1:length(beta_x_new)]);
    betas_x_new = rnorm(n=length(betas_x_new), mean=mean_betax[(length(beta_x_new)+1):length(betax)], sd=scales_r[(length(beta_x_new)+1):length(betax)]);
    
  }
  
  if(method=="tmcmc"){
    ########## new component ############
    betax_new = rep(0, length(betax));
    beta_x_new = rep(0, length(beta_x));
    betas_x_new = rep(0, length(betas_x));
    
    change_ep = rep(0, length(beta_x));
    change_ep = rnorm(n=1, mean=0, sd=1);
    change_ep = abs(change_ep);
    
    z_change = rep(0, length(change_ep));
    z_change = rbinom(n=length(beta_x), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    beta_x_new = beta_x + (scales_r[1:length(beta_x)]*z_change*change_ep);
    ####################################################
    
    change_ep = rep(0, length(betas_x));
    change_ep = rnorm(n=1, mean=0, sd=1);
    change_ep = abs(change_ep);
    
    z_change = rep(0, length(change_ep));
    z_change = rbinom(n=length(betas_x), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    betas_x_new = betas_x + (scales_r[(length(beta_x)+1):length(betax)]*z_change*change_ep);
    
  }
  
  betax_new = append(beta_x_new, betas_x_new);
  
  prob_x_new = logit_calc_spline_rand(N, beta_x_new, betas_x_new, rand, covar_x, covars_x, eps);
  pos1 = which(prob_x_new<eps);
  prob_x_new[pos1] = prob_x_new[pos1] + eps;
  
  pos_check = which(prob_x_new==1);
  if(length(pos_check)==0){
    ########## new likelihood ##########
    s1_new = XX*log(prob_x_new);
    s2_new = (1-XX)*log(1-prob_x_new);
    
    sum_new = sum(s1_new) + sum(s2_new);
    sum_new = sum_new - (sum((0.5*(betax_new-mean_betax)^2/sig_betax))/S_part);
    
    diff = sum_new - sum;
    uu = runif(n=1, 0, 1);  
    cat("sum=", sum, "sum_new=", sum_new, "\n");
    cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    if(log(uu)<diff){
      beta_x = beta_x_new;
      betas_x = betas_x_new;
      prob_x = prob_x_new;
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=3);
  result[1:length(beta_x),1] = beta_x;
  result[1:length(betas_x),2] = betas_x;
  result[,3] = prob_x;
  
  return(result)
}

update_beta_x_tmcmc_model2 <- function(method, N, len_x, XX, covar_x, covars_x, beta_x, betas_x, rand, eps, mean_betax, sig_betax, scales_r, prob_x){
  
  betax = append(beta_x, betas_x);
  ######## old component #############
  s1 = XX*log(prob_x); 
  s2 = (1-XX)*log(1-prob_x);
  
  sum = sum(s1) + sum(s2);
  sum = sum - sum((0.5*(betax-mean_betax)^2/sig_betax));
  
  if(method=="ordinary"){
    betax_new = rep(0, length(betax));
    beta_x_new = rep(0, length(beta_x));
    betas_x_new = rep(0, length(betas_x));
    
    #beta_x_new = rnorm(n=length(beta_x_new), mean=mean_betax[1:length(beta_x_new)], sd=sqrt(sig_betax[1:length(beta_x_new)]));
    #betas_x_new = rnorm(n=length(betas_x_new), mean=mean_betax[(length(beta_x_new)+1):length(betax)], sd=sqrt(sig_betax[(length(beta_x_new)+1):length(betax)]));
    
    beta_x_new = rnorm(n=length(beta_x_new), mean=mean_betax[1:length(beta_x_new)], sd=scales_r[1:length(beta_x_new)]);
    betas_x_new = rnorm(n=length(betas_x_new), mean=mean_betax[(length(beta_x_new)+1):length(betax)], sd=scales_r[(length(beta_x_new)+1):length(betax)]);
    
  }
  
  if(method=="tmcmc"){
    ########## new component ############
    betax_new = rep(0, length(betax));
    beta_x_new = rep(0, length(beta_x));
    betas_x_new = rep(0, length(betas_x));
    
    change_ep = rep(0, length(beta_x));
    change_ep = rnorm(n=1, mean=0, sd=1);
    change_ep = abs(change_ep);
    
    z_change = rep(0, length(change_ep));
    z_change = rbinom(n=length(beta_x), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    beta_x_new = beta_x + (scales_r[1:length(beta_x)]*z_change*change_ep);
    ####################################################
    
    change_ep = rep(0, length(betas_x));
    change_ep = rnorm(n=1, mean=0, sd=1);
    change_ep = abs(change_ep);
    
    z_change = rep(0, length(change_ep));
    z_change = rbinom(n=length(betas_x), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    betas_x_new = betas_x + (scales_r[(length(beta_x)+1):length(betax)]*z_change*change_ep);
    
  }
  
  betax_new = append(beta_x_new, betas_x_new);
  
  prob_x_new = logit_calc_spline_rand(N, beta_x_new, betas_x_new, rand, covar_x, covars_x, eps);
  pos1 = which(prob_x_new<eps);
  prob_x_new[pos1] = prob_x_new[pos1] + eps;
  
  pos_check = which(prob_x_new==1);
  if(length(pos_check)==0){
    ########## new likelihood ##########
    s1_new = XX*log(prob_x_new);
    s2_new = (1-XX)*log(1-prob_x_new);
    
    sum_new = sum(s1_new) + sum(s2_new);
    sum_new = sum_new - sum((0.5*(betax_new-mean_betax)^2/sig_betax));
    
    diff = sum_new - sum;
    uu = runif(n=1, 0, 1);  
    cat("sum=", sum, "sum_new=", sum_new, "\n");
    cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    if(log(uu)<diff){
      beta_x = beta_x_new;
      betas_x = betas_x_new;
      prob_x = prob_x_new;
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=3);
  result[1:length(beta_x),1] = beta_x;
  result[1:length(betas_x),2] = betas_x;
  result[,3] = prob_x;
  
  return(result)
}

pick_subsample <- function(aprob, dataframe){
  zprob = aprob[as.integer(dataframe$opened+1)];
  
  #cat(table(dataframe$opened), "\n")
  #cat(aprob[1], aprob[2], "\n")
  
  zsim = rbinom(n=length(zprob), size=1, prob=zprob);
  
  pos = which(zsim==1);
  pos = data.frame(location=pos);
  
  return(pos);
}

update_rand <- function(method, N, N_click, N_purchase, xx, yy, zz, covar_x, covars_x, covar_y, covars_y,
    covar_z, covars_z, beta_x, betas_x, beta_y, betas_y, beta_z, betas_z, beta_rand, eps, mean_rand, sig_rand,
    scales_rand, prob_x, prob_y, prob_z){
  
  #cat(sig_rand, "\n");
  
  s2 = sum(xx*log(prob_x)); 
  s2 = s2 + sum((1-xx)*log(1-prob_x));
  s2 = s2 + sum(yy*log(prob_y));
  s2 = s2 + sum((1-yy)*log(1-prob_y));
  s2 = s2 + sum(zz*log(prob_z));
  sum1 = s2 + sum((1-zz)*log(1-prob_z));
  sum1 = sum1 + log(dmvnorm(beta_rand, mean=mean_rand, sigma=sig_rand));
  
  if(method=="tmcmc"){
    ########## new component ############
    beta_rand_new = rep(0, length(beta_rand));
    
    change_ep = rep(0, length(beta_rand));
    change_ep = rnorm(n=1, mean=0, sd=1);
    change_ep = abs(change_ep);
    
    z_change = rep(0, length(change_ep));
    z_change = rbinom(n=length(beta_rand), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    beta_rand_new = beta_rand + (scales_rand*z_change*change_ep);
    ####################################################
  }
  
  ########## open stage ####################
  prob_x_new = logit_calc_spline_rand(N, beta_x, betas_x, beta_rand_new[1], covar_x, covars_x, eps);
  pos1 = which(prob_x_new<eps);
  prob_x_new[pos1] = prob_x_new[pos1] + eps;
  
  ########## click stage ####################
  prob_y_new = logit_calc_spline_rand(N_click, beta_y, betas_y, beta_rand_new[2], covar_y, covars_y, eps);
  pos2 = which(prob_y_new<eps);
  prob_y_new[pos2] = prob_y_new[pos2] + eps;
  
  ########## purchase stage ####################
  prob_z_new = logit_calc_spline_rand(N_purchase, beta_z, betas_z, beta_rand_new[3], covar_z, covars_z, eps);
  pos3 = which(prob_z_new<eps);
  prob_z_new[pos3] = prob_z_new[pos3] + eps;
  
  pos_check = which(prob_x_new==1); pos_check2 = which(prob_y_new==1); pos_check3 = which(prob_z_new==1);
  
  if((length(pos_check)==0) & (length(pos_check2)==0) & (length(pos_check3)==0)){
    ########## new likelihood ##########
    s2_new = sum(xx*log(prob_x_new));
    s2_new = s2_new + sum((1-xx)*log(1-prob_x_new));
    s2_new = s2_new + sum(yy*log(prob_y_new));
    s2_new = s2_new + sum((1-yy)*log(1-prob_y_new));
    s2_new = s2_new + sum(zz*log(prob_z_new));
    sum_new = s2_new + sum((1-zz)*log(1-prob_z_new));
    sum_new = sum_new + log(dmvnorm(beta_rand_new, mean=mean_rand, sigma=sig_rand));
    
    diff = sum_new - sum1;
    uu = runif(n=1, 0, 1);  
    cat("sum=", sum1, "sum_new=", sum_new, "\n");
    cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    if(log(uu)<diff){
      beta_rand = beta_rand_new;
      prob_x = prob_x_new;
      prob_y = prob_y_new;
      prob_z = prob_z_new;
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=4);
  result[1:length(beta_rand),1] = beta_rand;
  result[,2] = prob_x;
  result[1:length(prob_y),3] = prob_y;
  result[1:length(prob_z),4] = prob_z;
  
  return(result);
  
}

update_rand_new <- function(method, N, N_click, N_purchase, locs_num, pos_c, pos_p, xx, yy, zz, covar_x, covars_x, covar_y, covars_y,
                        covar_z, covars_z, beta_x, betas_x, beta_y, betas_y, beta_z, betas_z, beta_rand, eps, mean_rand, sig_rand,
                        scales_rand, prob_x, prob_y, prob_z){
  
  #cat(sig_rand, "\n");
  #cat(dim(covar_x), "\n");
  scales_exp = rep(scales_rand, dim(beta_rand)[1]);
  scales_exp = matrix(scales_exp, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
  
  s2 = sum(xx*log(prob_x)); 
  s2 = s2 + sum((1-xx)*log(1-prob_x));
  s2 = s2 + sum(yy*log(prob_y));
  s2 = s2 + sum((1-yy)*log(1-prob_y));
  s2 = s2 + sum(zz*log(prob_z));
  sum1 = s2 + sum((1-zz)*log(1-prob_z));
  sum1 = sum1 + sum(log(dmvnorm(beta_rand, mean=mean_rand, sigma=sig_rand)));
  
  if(method=="tmcmc"){
    ########## new component ############
    beta_rand_new = matrix(0, dim(beta_rand)[1], ncol=3);
    
    #change_ep = rep(0, length(beta_rand));
    change_ep = rnorm(n=1, mean=0, sd=1);
    change_ep = abs(change_ep);
    
    #z_change = rep(0, length(change_ep));
    z_change = rbinom(n=(dim(beta_rand)[1]*3), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    z_change = matrix(z_change, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
    
    beta_rand_new = beta_rand + (scales_exp*z_change*change_ep);
    ####################################################
  }
  
  beta_rand_exp_new = beta_rand_new[as.integer(locs_num),];
  
  beta_rand_exp_new_c = beta_rand_exp_new[as.integer(pos_c),2];
  beta_rand_exp_new_p = beta_rand_exp_new_c[as.integer(pos_p)];
  
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
    s2_new = sum(xx*log(prob_x_new));
    s2_new = s2_new + sum((1-xx)*log(1-prob_x_new));
    s2_new = s2_new + sum(yy*log(prob_y_new));
    s2_new = s2_new + sum((1-yy)*log(1-prob_y_new));
    s2_new = s2_new + sum(zz*log(prob_z_new));
    sum_new = s2_new + sum((1-zz)*log(1-prob_z_new));
    sum_new = sum_new + sum(log(dmvnorm(beta_rand_new, mean=mean_rand, sigma=sig_rand)));
    
    diff = sum_new - sum1;
    uu = runif(n=1, 0, 1);  
    cat("sum=", sum1, "sum_new=", sum_new, "\n");
    cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    if(log(uu)<diff){
      beta_rand = beta_rand_new;
      prob_x = prob_x_new;
      prob_y = prob_y_new;
      prob_z = prob_z_new;
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=4);
  result[1:(dim(beta_rand)[1]*3),1] = as.vector(t(beta_rand));
  result[,2] = prob_x;
  result[1:length(prob_y),3] = prob_y;
  result[1:length(prob_z),4] = prob_z;
  
  return(result);
  
}


update_rand_new_block <- function(block, method, N, N_click, N_purchase, locs_num, pos_c, pos_p, xx, yy, zz, covar_x, covars_x, covar_y, covars_y,
                            covar_z, covars_z, beta_x, betas_x, beta_y, betas_y, beta_z, betas_z, beta_rand, eps, mean_rand, sig_rand,
                            scales_rand, prob_x, prob_y, prob_z){
  
  #cat(sig_rand, "\n");
  #cat(dim(covar_x), "\n");
  scales_exp = rep(scales_rand, (length(block)*3));
  scales_exp = matrix(scales_exp, nrow=length(block), ncol=3, byrow=TRUE);
  
  s2 = sum(xx*log(prob_x)); 
  s2 = s2 + sum((1-xx)*log(1-prob_x));
  s2 = s2 + sum(yy*log(prob_y));
  s2 = s2 + sum((1-yy)*log(1-prob_y));
  s2 = s2 + sum(zz*log(prob_z));
  sum1 = s2 + sum((1-zz)*log(1-prob_z));
  sum1 = sum1 + sum(log(dmvnorm(beta_rand, mean=mean_rand, sigma=sig_rand)));
  
  
  if(method=="tmcmc"){
    ########## new component ############
    beta_rand_new = matrix(0, dim(beta_rand)[1], ncol=3);
    beta_rand_new = beta_rand;
    
    #change_ep = rep(0, length(beta_rand));
    change_ep = rnorm(n=1, mean=0, sd=1);
    change_ep = abs(change_ep);
    
    #z_change = rep(0, length(change_ep));
    z_change = rbinom(n=(length(block)*3), size=1, prob=0.5);
    z_change[z_change==0] = -1;
    z_change = matrix(z_change, nrow=length(block), ncol=3, byrow=TRUE);
    
    beta_rand_new[block,] = beta_rand[block,] + (scales_exp*z_change*change_ep);
    ####################################################
  }
  
  beta_rand_exp_new = beta_rand_new[as.integer(locs_num),];
  
  beta_rand_exp_new_c = beta_rand_exp_new[as.integer(pos_c),2];
  beta_rand_exp_new_p = beta_rand_exp_new_c[as.integer(pos_p)];
  
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
    s2_new = sum(xx*log(prob_x_new));
    s2_new = s2_new + sum((1-xx)*log(1-prob_x_new));
    s2_new = s2_new + sum(yy*log(prob_y_new));
    s2_new = s2_new + sum((1-yy)*log(1-prob_y_new));
    s2_new = s2_new + sum(zz*log(prob_z_new));
    sum_new = s2_new + sum((1-zz)*log(1-prob_z_new));
    sum_new = sum_new + sum(log(dmvnorm(beta_rand_new, mean=mean_rand, sigma=sig_rand)));
    
    diff = sum_new - sum1;
    uu = runif(n=1, 0, 1);  
    cat("sum=", sum1, "sum_new=", sum_new, "\n");
    cat("diff=", diff, "log(uu)=", log(uu), "\n");
    
    if(log(uu)<diff){
      beta_rand = beta_rand_new;
      prob_x = prob_x_new;
      prob_y = prob_y_new;
      prob_z = prob_z_new;
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=4);
  result[1:(dim(beta_rand)[1]*3),1] = as.vector(t(beta_rand));
  result[,2] = prob_x;
  result[1:length(prob_y),3] = prob_y;
  result[1:length(prob_z),4] = prob_z;
  
  return(result);
  
}

update_rand_new1 <- function(method, N, N_click, N_purchase, locs, locs_c, 
                             locs_p, pos_c, pos_p, xx, yy, zz, covar_x, covars_x, covar_y, covars_y,
                            covar_z, covars_z, beta_x, betas_x, beta_y, betas_y, beta_z, betas_z, beta_rand, eps, mean_rand, sig_rand,
                            scales_rand, prob_x, prob_y, prob_z){
  
  #cat(summary(locs$num)); cat("\n");
  #cat(summary(locs_c$num)); cat("\n");
  
  #cat(sig_rand, "\n");
  #cat(dim(covar_x), "\n");
  scales_exp = rep(scales_rand, dim(beta_rand)[1]);
  scales_exp = matrix(scales_exp, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
  
  ########### create the summation frame ########################
  frame_open = data.frame(num=locs$num, opart1=(xx*log(prob_x)), opart2=((1-xx)*log(1-prob_x)));
  frame_click = data.frame(num=locs_c$num, cpart1=(yy*log(prob_y)), cpart2=((1-yy)*log(1-prob_y)));
  frame_purchase = data.frame(num=locs_p$num, ppart1=(zz*log(prob_z)), ppart2=((1-zz)*log(1-prob_z)));
  
  agg1 = aggregate(frame_open, by=list(num1=frame_open$num), FUN=sum);
  agg1$num = NULL;
  agg1 = agg1[order(agg1$num1),];
  
  agg2 = aggregate(frame_click, by=list(num1=frame_click$num), FUN=sum);
  agg2$num = NULL;
  agg2 = agg2[order(agg2$num1),];
  
  agg3 = aggregate(frame_purchase, by=list(num1=frame_purchase$num), FUN=sum);
  agg3$num = NULL;
  agg3 = agg3[order(agg3$num1),];
  
  total_frame = data.frame(num1=1:length(agg1$num1), 
            prior=log(dmvnorm(beta_rand, mu=mean_rand, Sigma=sig_rand)));
  
  total_frame = merge(total_frame, agg1, by=c("num1"), all.x=T, all.y=T);
  total_frame = merge(total_frame, agg2, by=c("num1"), all.x=T, all.y=T);
  total_frame = merge(total_frame, agg3, by=c("num1"), all.x=T, all.y=T);
  total_frame = total_frame[order(total_frame$num1),];
  
  pos1_c = which(is.na(total_frame$cpart1)==T);  pos2_p = which(is.na(total_frame$ppart1)==T);
  if(length(pos1_c)!=0){
    total_frame$cpart1[pos1_c] = 0; total_frame$cpart2[pos1_c] = 0;
  }
  if(length(pos2_p)!=0){
    total_frame$ppart1[pos2_p] = 0; total_frame$ppart2[pos2_p] = 0;
  }
  
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
  beta_rand_exp_new_p = beta_rand_exp_new_c[as.integer(pos_p)];
  
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
    frame_click$cpart1 = (yy*log(prob_y_new));  frame_click$cpart2 = ((1-yy)*log(1-prob_y_new));
    frame_purchase$ppart1 = (zz*log(prob_z_new));  frame_purchase$ppart2 = ((1-zz)*log(1-prob_z_new));
    
    agg1 = aggregate(frame_open, by=list(num1=frame_open$num), FUN=sum);
    agg1$num = NULL;
    agg1 = agg1[order(agg1$num1),];
    
    agg2 = aggregate(frame_click, by=list(num1=frame_click$num), FUN=sum);
    agg2$num = NULL;
    agg2 = agg2[order(agg2$num1),];
    
    agg3 = aggregate(frame_purchase, by=list(num1=frame_purchase$num), FUN=sum);
    agg3$num = NULL;
    agg3 = agg3[order(agg3$num1),];
    
    total_frame = data.frame(num1=1:length(agg1$num1), 
                             prior=log(dmvnorm(beta_rand_new, mu=mean_rand, Sigma=sig_rand)));
    
    total_frame = merge(total_frame, agg1, by=c("num1"), all.x=T, all.y=T);
    total_frame = merge(total_frame, agg2, by=c("num1"), all.x=T, all.y=T);
    total_frame = merge(total_frame, agg3, by=c("num1"), all.x=T, all.y=T);
    total_frame = total_frame[order(total_frame$num1),];
    
    pos1_c = which(is.na(total_frame$cpart1)==T);  pos2_p = which(is.na(total_frame$ppart1)==T);
    if(length(pos1_c)!=0){
      total_frame$cpart1[pos1_c] = 0; total_frame$cpart2[pos1_c] = 0;
    }
    if(length(pos2_p)!=0){
      total_frame$ppart1[pos2_p] = 0; total_frame$ppart2[pos2_p] = 0;
    }
    
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
    }
  }
  
  result = matrix(0, nrow=length(prob_x), ncol=4);
  result[1:(dim(beta_rand)[1]*3),1] = as.vector(t(beta_rand));
  result[,2] = prob_x;
  result[1:length(prob_y),3] = prob_y;
  result[1:length(prob_z),4] = prob_z;
  
  return(result);
  
}


update_rand_new2 <- function(method, N, N_click, N_purchase, locs, locs_c, 
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

update_rand_new2_parallel <- function(method, N, N_click, N_purchase, locs, locs_c, 
                             locs_p, pos_c, pos_p, xx, yy, zz, pick_prob_o, pick_prob_c, pick_prob_p,
                             beta_rand, eps, mean_rand, sig_rand,
                             scales_rand, prob_x, prob_y, prob_z){
  
  
  scales_exp = rep(scales_rand, dim(beta_rand)[1]);
  scales_exp = matrix(scales_exp, nrow=dim(beta_rand)[1], ncol=3, byrow=TRUE);
  
  ########### create the summation frame ################################################################
  frame_open = data.frame(num=locs$num, opart1=(xx*log(prob_x)), opart2=((1-xx)*log(1-prob_x)));
  frame_open$cpart1 = 0; frame_open$cpart2=0;  frame_open$ppart1=0;  frame_open$ppart2 = 0;
  
  frame_open$cpart1[pos_c] = (yy*log(prob_y)); frame_open$cpart2[pos_c] = ((1-yy)*log(1-prob_y));
  frame_open$ppart1[pos_c[pos_p]] = (zz*log(prob_z)); 
  frame_open$ppart2[pos_c[pos_p]] = ((1-zz)*log(1-prob_z));
  
  
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
  
  prob_x_new = beta_rand_exp_new[,1] + pick_prob_o;
  prob_x_new = exp(prob_x_new)/(1+exp(prob_x_new));
  
  pos1 = which(prob_x_new<eps);
  prob_x_new[pos1] = prob_x_new[pos1] + eps;
  
  ########## click stage ####################
  prob_y_new = beta_rand_exp_new_c + pick_prob_c;
  prob_y_new = exp(prob_y_new)/(1+exp(prob_y_new));
  
  pos2 = which(prob_y_new<eps);
  prob_y_new[pos2] = prob_y_new[pos2] + eps;
  
  ########## purchase stage ####################
  prob_z_new = beta_rand_exp_new_p + pick_prob_p;
  prob_z_new = exp(prob_z_new)/(1+exp(prob_z_new));
  
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


