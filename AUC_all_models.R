rm(list=ls())

# Please set the path to the input data and sub-routines here.
# IT takes the following input files: training_data.txt, test_data.txt, output files from MODELNAME_test_pred.r, MODELNAME_training_pred.r
# Depending on stage and datatype (test/training), it creates overall AUC for the focal stage e.g. AUC_purchase_training.txt
# as well hit rate for different ranges of false rates e.g. NCG_model_training_roc_purchase.txt

####### specify if the exercise is for test or training data set ################
option_data = "training"; ## otherwise "training"

########## specify the option for producing ROC graph #######################
option = "purchase"; ### other options "open" or "click"

#################################################################################
if(option_data=="test"){
  model1 = read.table("NCG_test_pred.txt", T);
  
  model2 = read.table("NCL_test_pred.txt", T);
  
  model3 = read.table("CRON_test_pred.txt", T);
  
  model4 = read.table("SCRON_test_pred.txt", T);
  
}

if(option_data=="training"){
  model1 = read.table("NCG_training_pred.txt", T);
  
  model2 = read.table("NCL_training_pred.txt", T);
  
  model3 = read.table("CRON_training_pred.txt", T);
  
  model4 = read.table("SCRON_training_pred.txt", T);
}

####### initial runs were set and output was stored in results_save1/. TMCMC formula is used to calculate scales from the
#initial output

library(readxl)
##############################
options(digits=15);

eps=0.00000000016;
degree = 2; #### degree of spline
no_knots = 3 #### number of knots
########## initialising parameters ###############
mcmc_run = 8000; ##### number of MCMC run
burn = 3000; ###### number of burn-in 

############# load locations #####################

# Test and training Set
if(option_data=="test"){
  data = read.table("test_data.txt", T);
}
if(option_data=="training"){
  data = read.table("training_data.txt", T);
}
data = data.frame(data);


############## loading the indicators ######################
OO = data$opened; CC = data$clicked; PP = data$order;

OO1 = rep(0, length(OO)); CC1 = rep(0, length(CC)); PP1 = rep(0, length(PP));

pos_o1 = which(OO==0);
OO1[pos_o1] = rep(1, length(pos_o1));
pos_c1 = which(CC==0);
CC1[pos_c1] = rep(1, length(pos_c1)); 
pos_p1 = which(PP==0);
PP1[pos_p1] = rep(1, length(pos_p1));

######################################################
library(verification)


if(option=="open"){
  s1 =  roc.plot(model1$opened, model1$pred_est_o, plot.thres=NULL);
  
  s2 =  roc.plot(model2$opened, model2$pred_est_o, plot.thres=NULL);
  
  s3 =  roc.plot(model3$opened, model3$pred_est_o, plot.thres=NULL);
  
  s5 =  roc.plot(model4$opened, model4$pred_est_o, plot.thres=NULL);
}

if(option=="click"){
  s1 =  roc.plot(model1$clicked, model1$pred_est_c, plot.thres=NULL);
  
  s2 =  roc.plot(model2$clicked, model2$pred_est_c, plot.thres=NULL);
  
  s3 =  roc.plot(model3$clicked, model3$pred_est_c, plot.thres=NULL);
  
  s5 =  roc.plot(model4$clicked, model4$pred_est_c, plot.thres=NULL);
}

if(option=="purchase"){
  s1 =  roc.plot(model1$ordered, model1$pred_est_p, plot.thres=NULL);
  
  s2 =  roc.plot(model2$ordered, model2$pred_est_p, plot.thres=NULL);
  
  s3 =  roc.plot(model3$ordered, model3$pred_est_p, plot.thres=NULL);
  
  s5 =  roc.plot(model4$ordered, model4$pred_est_p, plot.thres=NULL);
}


ff1 = s1$plot.data;
gg1 = ff1[,,1];

ff2 = s2$plot.data;
gg2 = ff2[,,1];

ff3 = s3$plot.data;
gg3 = ff3[,,1];

ff5 = s5$plot.data;
gg5 = ff5[,,1];

########################################################################

first_area = data.frame(model=1, area=s1$roc.vol[2], pvalue=s1$roc.vol[3]);
second_area = data.frame(model=2, area=s2$roc.vol[2], pvalue=s2$roc.vol[3]);
third_area = data.frame(model=3, area=s3$roc.vol[2], pvalue=s3$roc.vol[3]);
fifth_area = data.frame(model=5, area=s5$roc.vol[2], pvalue=s5$roc.vol[3]);

add = rbind(first_area, second_area);
add = rbind(add, third_area);
add = rbind(add, fifth_area);
add = as.data.frame(add);

#######################################################################
if(option_data=="test"){
  if(option=="open"){
    first_model = data.frame(false_rate=gg1[,3], hit_rate=gg1[,2]);
    write.table(first_model, "NCG_test_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    second_model = data.frame(false_rate=gg2[,3], hit_rate=gg2[,2]);
    write.table(second_model, "NCL_model_test_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    third_model = data.frame(false_rate=gg3[,3], hit_rate=gg3[,2]);
    write.table(third_model, "CRON_model_test_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    fourth_model = data.frame(false_rate=gg5[,3], hit_rate=gg5[,2]);
    write.table(fourth_model, "SCRON_model_test_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    write.table(add, "AUC_open_test.txt", col.names=T, row.names=F, quote=F);
    
  }
  
  if(option=="click"){
    first_model = data.frame(false_rate=gg1[,3], hit_rate=gg1[,2]);
    write.table(first_model, "NCG_model_test_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    second_model = data.frame(false_rate=gg2[,3], hit_rate=gg2[,2]);
    write.table(second_model, "NCL_model_test_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    third_model = data.frame(false_rate=gg3[,3], hit_rate=gg3[,2]);
    write.table(third_model, "CRON_model_test_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    fourth_model = data.frame(false_rate=gg5[,3], hit_rate=gg5[,2]);
    write.table(fourth_model, "SCRON_model_test_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    write.table(add, "AUC_click_test.txt", col.names=T, row.names=F, quote=F);
    
  }
  
  if(option=="purchase"){
    first_model = data.frame(false_rate=gg1[,3], hit_rate=gg1[,2]);
    write.table(first_model, "NCG_model_test_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    second_model = data.frame(false_rate=gg2[,3], hit_rate=gg2[,2]);
    write.table(second_model, "NCL_model_test_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    third_model = data.frame(false_rate=gg3[,3], hit_rate=gg3[,2]);
    write.table(third_model, "CRON_test_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    fourth_model = data.frame(false_rate=gg5[,3], hit_rate=gg5[,2]);
    write.table(fourth_model, "SCRON_model_test_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    write.table(add, "AUC_purchase_test.txt", col.names=T, row.names=F, quote=F);
    
  }
  
}


if(option_data=="training"){
  if(option=="open"){
    first_model = data.frame(false_rate=gg1[,3], hit_rate=gg1[,2]);
    write.table(first_model, "NCG_model_training_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    second_model = data.frame(false_rate=gg2[,3], hit_rate=gg2[,2]);
    write.table(second_model, "NCL_model_training_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    third_model = data.frame(false_rate=gg3[,3], hit_rate=gg3[,2]);
    write.table(third_model, "CRON_model_training_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    fourth_model = data.frame(false_rate=gg5[,3], hit_rate=gg5[,2]);
    write.table(fourth_model, "SCRON_model_training_roc_open.txt", col.names=T, row.names=F, quote=F);
    
    write.table(add, "AUC_open_training.txt", col.names=T, row.names=F, quote=F);
    
  }
  
  if(option=="click"){
    first_model = data.frame(false_rate=gg1[,3], hit_rate=gg1[,2]);
    write.table(first_model, "NCG_model_training_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    second_model = data.frame(false_rate=gg2[,3], hit_rate=gg2[,2]);
    write.table(second_model, "NCL_model_training_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    third_model = data.frame(false_rate=gg3[,3], hit_rate=gg3[,2]);
    write.table(third_model, "CRON_model_training_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    fourth_model = data.frame(false_rate=gg5[,3], hit_rate=gg5[,2]);
    write.table(fourth_model, "SCRON_model_training_roc_click.txt", col.names=T, row.names=F, quote=F);
    
    write.table(add, "AUC_click_training.txt", col.names=T, row.names=F, quote=F);
    
  }
  
  if(option=="purchase"){
    first_model = data.frame(false_rate=gg1[,3], hit_rate=gg1[,2]);
    write.table(first_model, "NCG_model_training_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    second_model = data.frame(false_rate=gg2[,3], hit_rate=gg2[,2]);
    write.table(second_model, "NCL_model_training_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    third_model = data.frame(false_rate=gg3[,3], hit_rate=gg3[,2]);
    write.table(third_model, "CRON_model_training_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    fourth_model = data.frame(false_rate=gg5[,3], hit_rate=gg5[,2]);
    write.table(fourth_model, "SCRON_model_training_roc_purchase.txt", col.names=T, row.names=F, quote=F);
    
    write.table(add, "AUC_purchase_training.txt", col.names=T, row.names=F, quote=F);
    
  }
  
}






