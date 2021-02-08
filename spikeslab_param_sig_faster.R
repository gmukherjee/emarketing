
################ This part of the code is for the nonparallel run ###################

aa = read.table("results_save2/beta_o_model7_spikeslab_part1.txt");

param = aa[,c(14:dim(aa)[2])];

mean_param = apply(param, 2, FUN=mean, na.rm=T);

sd_param = apply(param, 2, FUN=sd, na.rm=T);

ratio = abs(mean_param)/sd_param;

pos = which(ratio<2.56);

##############################################################
aa = read.table("results_save2/beta_c_model7_spikeslab_part1.txt");

param = aa[,c(14:dim(aa)[2])];

mean_param = apply(param, 2, FUN=mean, na.rm=T);

sd_param = apply(param, 2, FUN=sd, na.rm=T);

ratio = abs(mean_param)/sd_param;

pos = which(ratio<2.56);

##############################################################
aa = read.table("results_save2/beta_p_model7_spikeslab_part1.txt");

param = aa[,c(14:dim(aa)[2])];

mean_param = apply(param, 2, FUN=mean, na.rm=T);

sd_param = apply(param, 2, FUN=sd, na.rm=T);

ratio = abs(mean_param)/sd_param;

pos = which(ratio<2.56);


###########################################################

bb=read.table("results_save2/delta_spikeslab_model7_spikeslab_part1.txt");

prop = apply(bb, 2, mean);

rspike = length(which(prop<0.5))/length(which(is.na(prop)==F));

pick = which(prop<0.5);
######################################3

prop1 = prop[1:74];

pick1 = which(prop1<0.5);

prop2 = prop[75:148];

pick2 = which(prop2<0.5);

prop3 = prop[149:222];

pick3 = which(prop3<0.5);

write.table(pick1, "spikeslab/reject_sites_open.txt", col.names=F, row.names=F, quote=F);

write.table(pick2, "spikeslab/reject_sites_click.txt", col.names=F, row.names=F, quote=F);

write.table(pick3, "spikeslab/reject_sites_purchase.txt", col.names=F, row.names=F, quote=F);


########### This part of the code is for parallel run ####################################


len1 = 102;
len2 = 102;
len3 = 94;

mcmc_run = 8000;
burn = 3000;

bb = read.table("results_parallel/spikeslab_part1_beta_sigma.txt");

prop = bb[(len1+len2+len3+1):dim(bb)[1],1]/(mcmc_run-burn);

prop1 = prop[1:74];

pick1 = which(prop1<0.5);

prop2 = prop[75:148];

pick2 = which(prop2<0.5);

prop3 = prop[149:222];

pick3 = which(prop3<0.5);

write.table(pick1, "spikeslab/reject_sites_open_parallel.txt", col.names=F, row.names=F, quote=F);

write.table(pick2, "spikeslab/reject_sites_click_parallel.txt", col.names=F, row.names=F, quote=F);

write.table(pick3, "spikeslab/reject_sites_purchase_parallel.txt", col.names=F, row.names=F, quote=F);


