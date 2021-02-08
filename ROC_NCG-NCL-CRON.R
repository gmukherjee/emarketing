rm(list=ls())

# Please set the path to the input data and sub-routines here.
# IT takes the following input files: training_data.txt, test_data.txt, NCG_test_pred.txt, NCL_test_pred.txt, CRON_test_pred.txt,
# NCG_training_pred.txt, NCL_training_pred.txt, and CRON_training_pred.txt
# Depeding on stage and datatype (test/training), it creates ROC Curves comparing NCG, NCL and CRON e.g. ROC_compare_purchase_training_NCG-NCL-CRON.pdf

####### specify if the exercise is for test or training data set ################
option_data = "training"; ## otherwise "training"

########## specify the option for producing ROC graph #######################
option = "purchase"; ### other options "open" or "click"

###################################################################################
if(option_data=="test"){
  model1 = read.table("NCG_test_pred.txt", T);
  
  model2 = read.table("NCL_test_pred.txt", T);
  
  model3 = read.table("CRON_test_pred.txt", T);
  
}

if(option_data=="training"){
  model1 = read.table("NCG_training_pred.txt", T);
  
  model2 = read.table("NCL_training_pred.txt", T);
  
  model3 = read.table("CRON_training_pred.txt", T);
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

############################ loading data #######################################

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
  
}

if(option=="click"){
  s1 =  roc.plot(model1$clicked, model1$pred_est_c, plot.thres=NULL);
  
  s2 =  roc.plot(model2$clicked, model2$pred_est_c, plot.thres=NULL);
  
  s3 =  roc.plot(model3$clicked, model3$pred_est_c, plot.thres=NULL);
  
}

if(option=="purchase"){
  s1 =  roc.plot(model1$ordered, model1$pred_est_p, plot.thres=NULL);
  
  s2 =  roc.plot(model2$ordered, model2$pred_est_p, plot.thres=NULL);
  
  s3 =  roc.plot(model3$ordered, model3$pred_est_p, plot.thres=NULL);
  
}


ff1 = s1$plot.data;
gg1 = ff1[,,1];

ff2 = s2$plot.data;
gg2 = ff2[,,1];

ff3 = s3$plot.data;
gg3 = ff3[,,1];

if(option_data=="test"){
  if(option=="open"){ pdf("ROC_compare_open_test_NCG-NCL-CRON.pdf"); }
  
  if(option=="click"){ pdf("ROC_compare_click_test_NCG-NCL-CRON.pdf"); }
  
  if(option=="purchase"){ pdf("ROC_compare_purchase_test_NCG-NCL-CRON.pdf"); }
  
  plot(gg1[,3], gg1[,2], xlim=c(0, 1), ylim=c(0, 1), 'l', lty=1, col="green", xlab="False rate", ylab="Hit rate", main="");
  lines(gg2[,3], gg2[,2], col="blue")
  lines(gg3[,3], gg3[,2], col="red")
  legend(0.7, 0.3, legend=c("NCG","NCL","CRON"),col=c("green", "blue", "red"), lwd=c(1,1,1), cex=0.7)
  
  
  dev.off()
}

if(option_data=="training"){
  if(option=="open"){ pdf("ROC_compare_open_training_NCG-NCL-CRON.pdf"); }
  
  if(option=="click"){ pdf("ROC_compare_click_training_NCG-NCL-CRON.pdf"); }
  
  if(option=="purchase"){ pdf("ROC_compare_purchase_training_NCG-NCL-CRON.pdf"); }
  
  plot(gg1[,3], gg1[,2], xlim=c(0, 1), ylim=c(0, 1), 'l', lty=1, col="green", xlab="False rate", ylab="Hit rate", main="");
  lines(gg2[,3], gg2[,2], col="blue")
  lines(gg3[,3], gg3[,2], col="red")
  legend(0.7, 0.3, legend=c("NCG","NCL","CRON"),col=c("green", "blue", "red"), lwd=c(1,1,1), cex=0.7)
  
  dev.off()
}

