
source('/media/w/WKMACLTMP/R code/02_setupObjects_2_other objects.r')


#---------------------------------------------------------------------------------------------------------------------------------
# to run HCR simulations  (simple runs)  -------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

source('/media/w/WKMACLTMP/R code/03_runMSE.r')
source('/media/w/WKMACLTMP/R code/03.0_runMSE test.r')


#---------------------------------------------------------------------------------------------------------------------------------
# to run HCR simulations  (screen through range of F and B values)-------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# base case
source('/media/w/WKMACLTMP/R code/HCR testing/03.2_run sim batch2.0.r')
source('/media/w/WKMACLTMP/R code/HCR testing/03.2_run sim batch2.2.r')
source('/media/w/WKMACLTMP/R code/HCR testing/03.2_run sim batch2.4.r')
source('/media/w/WKMACLTMP/R code/HCR testing/03.2_run sim batch2.6.r')
source('/media/w/WKMACLTMP/R code/HCR testing/03.2_run sim batch2.8.r')
source('/media/w/WKMACLTMP/R code/HCR testing/03.2_run sim batch3.0.r')
source('/media/w/WKMACLTMP/R code/HCR testing/03.2_run sim batch3.2.r')



source('/media/w/WKMACLTMP/R code/08_outputDiagnostics batch.r')


# base case + reversible biol
source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch2.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch2.2.r')
source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch2.4.r')
source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch2.6.r')
source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch2.8.r')
source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch3.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch3.2.r')

source('/media/w/WKMACLTMP/R code/Sensitivity perm/03.2_run sim batch.r')

source('/media/w/WKMACLTMP/R code/Sensitivity perm/08_outputDiagnostics batch.r')


# base case + sensitivity recruitment
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/03.2_run sim batch2.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/03.2_run sim batch2.2.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/03.2_run sim batch2.4.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/03.2_run sim batch2.6.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/03.2_run sim batch2.8.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/03.2_run sim batch3.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/03.2_run sim batch3.2.r')




# base case + using a reference selection pattern
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/03.2_run sim batch2.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/03.2_run sim batch2.2.r')
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/03.2_run sim batch2.4.r')
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/03.2_run sim batch2.6.r')
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/03.2_run sim batch2.8.r')
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/03.2_run sim batch3.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/03.2_run sim batch3.2.r')

source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/08_outputDiagnostics batch.r')

# recovery from Blim
source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/03.2_run sim batch2.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/03.2_run sim batch2.2.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/03.2_run sim batch2.4.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/03.2_run sim batch2.6.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/03.2_run sim batch2.8.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/03.2_run sim batch3.0.r')
source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/03.2_run sim batch3.2.r')

source('/media/w/WKMACLTMP/R code/Sensitivity recovery from Blim/08_outputDiagnostics batch.r')

# sensitivity banking borrowing
source('/media/w/WKMACLTMP/R code/Sensitivity BB/03.2_run sim batch2.2.r')
source('/media/w/WKMACLTMP/R code/Sensitivity BB/08_outputDiagnostics batch.r')


# sensitivity var constraint
source('/media/w/WKMACLTMP/R code/Sensitivity Varconst/03.2_run sim batch2.2.r')
source('/media/w/WKMACLTMP/R code/Sensitivity Varconst/08_outputDiagnostics batch.r')

# sensitivity STF
source('/media/w/WKMACLTMP/R code/Sensitivity STF/03.2_run sim batch2.2.r')



#---------------------------------------------------------------------------------------------------------------------------------
# for MSY runs -------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# base case
source('/media/w/WKMACLTMP/R code/MSY runs/09_equilibrium perceived F0_10.r')
source('/media/w/WKMACLTMP/R code/MSY runs/09_equilibrium perceived F11_20.r')
source('/media/w/WKMACLTMP/R code/MSY runs/09_equilibrium perceived F21_30.r')
source('/media/w/WKMACLTMP/R code/MSY runs/09_equilibrium perceived F31_40.r')
source('/media/w/WKMACLTMP/R code/MSY runs/09_equilibrium perceived F41_50.r')



source('/media/w/WKMACLTMP/R code/09_2_extract results_equilibrium perceived F.r')

# base case + reversible biol
source('/media/w/WKMACLTMP/R code/Sensitivity perm/09_equilibrium perceived F.r')

# base case - Ricker SR model only
source('/media/w/WKMACLTMP/R code/Sensitivity recruitment/09_equilibrium perceived F Ricker only.r')

# base case + using a reference selection pattern
source('/media/w/WKMACLTMP/R code/Sensitivity ref selection pattern/09_equilibrium perceived F.r')