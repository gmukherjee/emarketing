#####################################
Generic description of Codes for NCG, NCL, and CRON models:
#####################################

MODELNAME_learning.r:
This code produces MCMC simulations for fixed and random effects coefficients and for their variances.
Such output are used for generating predictions and parameter estimates for the three stages of conversion (namely, open, click, and purchase).
* It takes the following input files: 
	- training_data.txt
	- starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx) of the focal model
	- locations of sub-sampling for the purchase stage locations_subsampling_purchase.csv
* it creates the output MODELNAME_learning_output.txt, which becomes an input in the coeff generation code, test and training prediction codes.

MODELNAME_coeff_gen.r
This code is used to generate parameter estimates for three stages of conversion (namely, open, click, and purchase)
* It takes the following input files: 
	- starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx)
	- output from the learning code for the focal model signified as MODELNAME_learning_output.txt 
* It creates output files which have the parameter estimates for open, click, and purchase stages (MODELNAME_parameter_stagename.xlsx) as well as the correlation matrix (MODELNAME_parameter_correlation.xlsx)

MODELNAME_test_pred.r
This code creates the probability predictions from the test data for each stage of conversion (namely, open, click, and purchase).
* It takes the following input files: 
	- training_data.txt
	- starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx) of the focal model
	- locations of sub-sampling for the purchase stage locations_subsampling_purchase.csv
	- expanded form of the test data predict_log.csv
	- output from the learning code MODELNAME_learning_output.txt
# it creates the test prediction for each of the three stages and stores them in MODELNAME_test_pred.txt.

MODELNAME_training_pred.r
This code creates the probability predictions from the training data for each stage of conversion (namely, open, click, and purchase).
* It takes the following input files: 
	- training_data.txt
	- starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx) of the focal model
	- locations of sub-sampling for the purchase stage locations_subsampling_purchase.csv
	- output from the learning code MODELNAME_learning_output.txt
# it creates the probability predictions on the training data for each of the three stages and stores them in MODELNAME_training_pred.txt.

#####################################
Codes specific to SCRON models:
#####################################

SCRON_learning_partA.r
This code implements the spike and slab.
* It takes the following input files: 
	- training_data.txt
	- starting values of coefficients of each stage (initial_coeff_MODELNAME_stagename.xlsx) of the focal model
	- locations of sub-sampling for the purchase stage locations_subsampling_purchase.csv
* It creates the output SCRON_learning_outputA.txt, which becomes an input in the code that generates the rejection sides after variable selection (SCRON_learning_rejectSite.r).

SCRON_learning_rejectSite.r
This code generates the rejected promo codes for each stage.
* It takes the following input files: 
	- SCRON_learning_outputA.txt which is the Output from SCRON_learning_partA.r becomes the input here
* It creates three files (reject_sites_stagename.txt), one for each stage which stores the promo code which will not be included in the analysis for MCMC simulations for fixed and random effects coefficients and for their variances.

Following the above two codes, SCRON model uses the following codes, the descriptions of which are given above:
	* SCRON_learning_partB.r performs the same function as MODELNAME_learning.r described above
	* SCRON_coeff_gen.r
	* SCRON_test_pred.r
	* SCRON_training_pred.r

#####################################
AUC_all_models.r
#####################################
This code is respnsible for creating the AUC values.
*It takes the following input files: 
	- training_data.txt
	- test_data.txt, output files from MODELNAME_test_pred.r
	- MODELNAME_training_pred.r
* It also takes two arguments: 
	- dataype: test or training
	- stage type: open or click or purchase
* Depending on stage type and datatype (test/training), it creates overall AUC for the focal stage e.g. AUC_purchase_training.txt as well hit rate for different ranges of false rates e.g. NCG_model_training_roc_purchase.txt

#####################################
ROC_NCG-NCL-CRON.r
#####################################
This code is responsible for creating the ROC curves for each stage and data type for the models NCG, NCL, and CRON.
* It takes the following input files: 
	- training_data.txt
	- test_data.txt
	- NCG_test_pred.txt
	- NCL_test_pred.txt
	- CRON_test_pred.txt
	- NCG_training_pred.txt
	- NCL_training_pred.txt
	- CRON_training_pred.txt
* It also takes two arguments: 
	- dataype: test or training
	- stage type: open or click or purchase
* Depeding on stage type and datatype (test/training), it creates ROC Curves comparing NCG, NCL and CRON e.g. ROC_compare_purchase_training_NCG-NCL-CRON.pdf

#####################################
ROC_CRON-SCRON.r
#####################################
This code is responsible for creating the ROC curves for each stage and data type for the models CRON and SCRON.
* It takes the following input files:
	- training_data.txt
	- test_data.txt
	- CRON_test_pred.txt
	- SCRON_test_pred.txt,
	- CRON_training_pred.txt
	- SCRON_training_pred.txt
* Depeding on stage type and datatype (test/training), it creates ROC Curves comparing CRON and SCRON e.g. ROC_compare_purchase_training_CRON-SCRON.pdf


#####################################
Generate_Summ_Stats.r
#####################################
This code is responsible for creating the summary statistics mentioned in the paper, namely, Table S1, Table S2, Table S3, Figure S1.
* Input files required are:
	- EmailCampaign_MainData.csv
	- EmailSumm_Proj.csv
