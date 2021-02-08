rm(list=ls())

# This code generates the rejected promo codes for each stage reject_sites_stagename.txt.
# SCRON_learning_outputA.txt which is the Output from SCRON_learning_partA.r becomes the input here

len1 = 102;
len2 = 102;
len3 = 94;

mcmc_run = 8000;
burn = 3000;

bb = read.table("SCRON_learning_outputA.txt");

prop = bb[(len1+len2+len3+1):dim(bb)[1],1]/(mcmc_run-burn);

prop1 = prop[1:74];

pick1 = which(prop1<0.5);

prop2 = prop[75:148];

pick2 = which(prop2<0.5);

prop3 = prop[149:222];

pick3 = which(prop3<0.5);

write.table(pick1, "reject_sites_open.txt", col.names=F, row.names=F, quote=F);

write.table(pick2, "reject_sites_click.txt", col.names=F, row.names=F, quote=F);

write.table(pick3, "reject_sites_purchase.txt", col.names=F, row.names=F, quote=F);


