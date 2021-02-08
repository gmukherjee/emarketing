#####################################
Description of Raw Data:
#####################################

EmailCampaign_MainData.csv: 
This is the main data which has 1.62M observations.
These observations were created by cataloging responses (open, click, purchase) to 25 promotion emails sent to one or more of 77986 customers. 
Given below is the data dictionary for the above data:
	* recipient_id - customer indetifier number
	* sendDate - data at which the email was send
	* opened - whether the email was opened or not
	* clicked - whether the email was clicked or not
	* ordered - whether a purchase was made or not
	* total_orders_val - if a purchase was made, what was the amount
	* aov_retail - average retail spend
	* aov_web - average web spend
	* Order_cnt - previous number of purchases
	* age - age binned into five levels (description in Data section of main paper)
	* gender - male/female (not well populated, therefore was not used)
	* income - income binned into nine levels (description in Data section of main paper)
	* daysSinceAccOpen - days since account opned
	* daysSinceOpened - days since the last time the customer opened an email
	* daysSinceClicked - days since the last time the customer clicked on an email
	* daysSincePurchased - days since the last time the customer purchased from an email
	* ID - Email Promo ID (25 options). The true ID is masked. The serial ID could be learned by merging with EmailSumm_Proj.csv.
	* Type - Customer Type (3 levels)


EmailSumm_Proj.csv:
This dataset has the classification information on whether an email is a price promotion (P), non-price promotion (N) or a combination of both (C). 

training_data.txt:
This dataset is created by random sampling 67986 customers from EmailCampaign_MainData.csv.
This dataset is used to train models specified in codes denoted by MODELNAME_learning.r

test_data.txt:
This dataset is created by taking customers from EmailCampaign_MainData.csv who are not included in the training_data.txt.
This dataset is used to create test predictions specified in codes denoted by MODELNAME_test_pred.r

locations_subsampling_purchase.csv:
This data stores the sub-sampling locations to be used in the purchase stage.
This sub-sampling is done to correct for the imbalance in the data created by rare instances of purchase.

predict_log.csv:
This is also another version of the test data where the factors and interactions are created as separate columns.

initial_coeff_MODELNAME_stage.csv:
This input file is a generic name for initial estimates that are necessary for the learning models. Each stage (open, click, purchase) and each model (NCG, NCL, CRON, SCRON) have their own file.