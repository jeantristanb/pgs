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
* `--type_out`, type="character", default="Coefficient", help="used Coefficient of lm or OR to showed increased"
* `--quantile_number`: type="integer", default=10, help="number of quantile to plot"
* `--tr_var`: type="character", default="nullfct", help="transforms variables to perform lm and residualisations"
* `--tr_res`, type="character", default="nullfct", help : "transform residuals after glm"
* `--out`, type="character", default="save_", help="output file name 
 
```

Rscript pgs/PRSice-2/plotresult/script_prs_cmdl.r --data GWAS_pheno.txt --pheno  DBP_statusAHM321 --covariables age,Age2,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10 --best_prs file1.best,file2.best,file3.best --head_bestprs African,Europea,Trans --quantile_number 5 --out DBP
```
![distribution by quanile](./imagesexample/DBP_dist_quant.pdf)
![distribution quantitative](./imagesexample/DBP_varexplained.pdf)
