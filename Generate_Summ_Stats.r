rm(list=ls())
# Please set the path of the input files here

# This code is responsible for creating the summary statistics mentioned in the paper
# namely Table S1, Table S2, Table S3, Figure S1.
# Input files required are:
# - EmailCampaign_MainData.csv
# - EmailSumm_Proj.csv

library(sqldf)

dataC = data.frame(read.csv("EmailCampaign_MainData.csv",header=TRUE))

temp1 = sqldf('select recipient_id, count(ID) promo_cnt from dataC group by recipient_id')
tableS1 = sqldf('select promo_cnt No_of_Promotions, count(recipient_id) No_of_Recipients from temp1 group by promo_cnt')
#write.csv(tableS1,"tableS1.csv", row.names=FALSE)

dataP = read.csv("EmailSumm_Proj.csv", header=TRUE)
dataCPT = data.frame(merge(dataC,dataP,by = "ID"))
tableS2 = sqldf("select NewPromoId, PromoType, count(recipient_id) Count, sum(opened) OpenRate, sum(clicked) ClickRate, sum(ordered) OrderRate, avg(daysSinceOpened) MeanLastOpen, stdev(daysSinceOpened) SDLastOpen, avg(daysSinceClicked) MeanLastClick, stdev(daysSinceClicked) SDLastClick, avg(daysSincePurchased) MeanLastPur, stdev(daysSincePurchased) SDLastPur, avg(Order_cnt) MeanOderCnt, stdev(Order_cnt) SDOrderCnt, avg(aov_retail) MeanAovRetail  , stdev(aov_retail) SDAOVRetail ,avg(aov_web) MeanAovWeb  , stdev(aov_web) SDAovWeb from dataCPT group by NewPromoId")
tableS2 = data.frame(tableS2)
tableS2$OpenRate = tableS2$OpenRate/tableS2$Count
tableS2$ClickRate = tableS2$ClickRate/tableS2$Count
tableS2$OrderRate = tableS2$OrderRate/tableS2$Count
#write.csv(tableS2,"tableS2.csv", row.names=FALSE)

tableS3  = sqldf("select Type, PromoType, count(recipient_id) Count, sum(opened) OpenRate, sum(clicked) ClickRate, sum(ordered) OrderRate, avg(daysSinceOpened) MeanLastOpen, stdev(daysSinceOpened) SDLastOpen, avg(daysSinceClicked) MeanLastClick, stdev(daysSinceClicked) SDLastClick, avg(daysSincePurchased) MeanLastPur, stdev(daysSincePurchased) SDLastPur, avg(Order_cnt) MeanOderCnt, stdev(Order_cnt) SDOrderCnt, avg(aov_retail) MeanAovRetail  , stdev(aov_retail) SDAOVRetail ,avg(aov_web) MeanAovWeb  , stdev(aov_web) SDAovWeb from dataCPT group by Type, PromoType")
tableS3 = data.frame(tableS3)
tableS3$OpenRate = tableS3$OpenRate/tableS3$Count
tableS3$ClickRate = tableS3$ClickRate/tableS3$Count
tableS3$OrderRate = tableS3$OrderRate/tableS3$Count
#write.csv(tableS3,"tableS3.csv", row.names=FALSE)

N = length(unique(dataC$recipient_id))
temp2 = sqldf('select NewPromoID, count(recipient_id) no_Cust from dataCPT group by NewPromoID')
figureS1 = data.frame(temp2)
figureS1$no_Cust = (figureS1$no_Cust/N)*100
figureS1 = merge(figureS1,dataP, merge = "NewPromoId" )
#write.csv(figureS1,"figureS1.csv", row.names=FALSE)