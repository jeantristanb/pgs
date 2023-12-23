Rscript ../script_prs_cmdl.r --data  sim_1.phe --pheno  Pheno --covariables  Sex,HTN --best_prs sim_1_sd3.best,sim_1_sd4.best --head_bestprs sd3,sd4 --quantile_number 5 --out comparison_sd
Rscript ../script_prs_cmdl.r --data  sim_1.phe --pheno  Pheno --covariables  Sex --best_prs sim_1_sd3.best,sim_1_sd4.best --head_bestprs sd3,sd4 --quantile_number 5 --out comparison_sd_htn --pop_header HTN
Rscript ../script_prs_cmdl.r --data  sim_1.phe --pheno  Pheno --covariables  Sex --best_prs sim_1_sd3.best --head_bestprs sd3 --quantile_number 5 --out comparison_sd_htn_sd3 --pop_header HTN

