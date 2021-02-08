###############################
This help file describes all the sub-routines used in the codes described in README_Codes.txt
###############################

joint_subroutine_random_mod.r
This subroutine contains updated versions of calculation of MCMC steps, a modified version of set of codes in joint_subroutine_random.r. It is used to run learning and prediction codes of all four models.

spikeslab_param_sig_faster.r
This subroutine produces the indicators of regression components which are to be rejected through spike and slab application. This subroutine is particularly written for the code when all 67000 random components are updated jointly.

spikeslab_random_joint.r
This subroutine is used for producing MCMC update for the fixed and random effect coefficients in the part A of the spike and slab application. It is used in the code SCORN_learning_partA.r

joint_subroutine_random.r
This subroutine contains codes for running MCMC steps of the learning and prediction codes. Since most of the these codes were modified later, this subroutine is now
mostly unused.

check_random_new2_code.r 
This subroutine contains codes for updating gamma (random components). Can be avoided in this application.