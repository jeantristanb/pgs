# compute and plot result of PRSice-2
## input :
* take data, pheno and covariable for raw data, take for 
* take input for PRSCise-2, `.best` result
## algorithm :
 * variable transformation
  * transform pheno (`--tr_var`) 
  * used glm with covariable  and extracted residuals 
  * transorm covariable
 * compared PRS from PRSCise-2 and Residuals for each input of PRSCise-2 :
  * variance explained of models as : lm(Pheno residual~PRS) and plot
  * plot of residuals by PRS quantile : PRS by quantile, and comparaison of quantile with phenotype see `--quantile_number`, comparison with quantile 1
  
## output :
 * plot of variance
 * plot of decile 
 * file contained information relative to results

## argument :
* `--data` : type="character", default=NULL, help="phenotype file, each column separated by a tabulation or space"
* `--best_prs`  type="character", default=NULL, help="1 file or more from best result PRScise-2 separed by a comma"
* `--head_bestprs` : type="character", default=NULL, help="header name to give at each prs file, must be in same order than `--best_prs`"
* `--pheno` , type="character", default=NULL, help= "henotype header of data"
* `--covariables`, type="character", default=NULL, help="list of covariable separated by a comma, must be present in data"
* `--plotqt_type`, type="character", default="Coefficient", help="plot Coefficient, OR or mean between different first quantile of PRS and other quantile, Effect/OR : effect or OR are result of LM"
* `--quantile_number`: type="integer", default=10, help="number of quantile to plot"
* `--tr_var`: type="character", default="nullfct", help="transforms variables to perform lm and residualisations"
* `--tr_res`: type="character", default="nullfct", help : "transform residuals after glm"
* `--prs_header` :  type="character", default="PRS", help="PRS header in PRS file"
* `--out`, type="character", default="save_", help="output file name 
* `--pop_header` type="character", default="All", help = 'split analyse by one phenotype, by default none'
 
## example 
see [see folder : data_test](./data_test/)
