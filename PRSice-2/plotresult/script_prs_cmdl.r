rm(list=ls())

## function to defined path
thisPath <- function() {
  stub <- function() {}
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    scriptPath<-normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    scriptPath <-dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    scriptPath <-dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
  return(scriptPath)
}

check_andloadpackage<-function(packages){
        #packages <- c("devtools","optparse", "shinydashboard", "nVennR", "nVennR", "rsvg", "ggrepel", "shiny", "DT", "Cairo", "RColorBrewer", "Cairo", 'grid', 'lattice', "shinyjs")
        # Install packages not yet installed
        installed_packages <- packages %in% rownames(installed.packages())
        if (any(installed_packages == FALSE)) {
                cat("package need in ubunu :\n * R : apt-get install r-dev \n* apt-get install librsvg2-dev libcurl4-openssl-dev libxml2-dev libssl-dev\n")
                packagetoinstall=packages[!installed_packages]
                if('devtools' %in% packagetoinstall){
                        install.packages('devtools')
                }
                installed_packages <- packages %in% rownames(installed.packages())
                packagetoinstall=packages[!installed_packages]
                install.packages(packagetoinstall[!(packagetoinstall %in% c('rsvg','devtools'))])
                if('rsvg' %in% packagetoinstall){
                        devtools::install_github('https://github.com/jeroen/rsvg')
                }
        }
        invisible(lapply(packages, library, character.only = TRUE))
}
check_andloadpackage(c('optparse', 'psychometric', 'ggplot2'))

option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="phenotype file, each column separated by a tabulation or space", metavar="character"),
  make_option(c("--best_prs"), type="character", default=NULL,
              help="1 pr more than 1 files from PRScise separed by a comma", metavar="character"),
  make_option(c("--head_bestprs"), type="character", default=NULL,
              help="header name to give at each prs file, must be in same order than --best_prs", metavar="character"),
  make_option(c("--pheno"), type="character", default=NULL,
              help="phenotype header of data", metavar="character"),
  make_option(c("--covariables"), type="character", default=NULL,
              help="list of covariable separated by a comma, must be present in data", metavar="character"),
  make_option(c("--type_out"), type="character", default="Coefficient",
              help="used Coefficient of lm or mean to showed increased", metavar="character"),
  make_option(c("--quantile_number"), type="integer", default=10,
              help="number of quantile to plot", metavar="character"),
  make_option(c("--tr_var"), type="character", default="nullfct",
              help="transforms variables to perform lm and residualisations", metavar="character"),
  make_option(c("--tr_res"), type="character", default="nullfct",
              help="transform residuals after glm", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="save_",
              help="output file name [default= %default]", metavar="character")
);

script.dir=thisPath()
fileutils<-paste(script.dir, '/utils.r',sep='')
if(!file.exists(fileutils)){
 cat('utils.r contains function not found\nexit') 	
 q('no', 5)
}
#library('psychometric')
#library(ggplot2)
source(fileutils)


#listfile_prs<-paste(c("GWAS321_PRS_Sep22/1.SBP/Multi_PAGE/SBP_statusAHM321__PRSout.best","GWAS321_PRS_Sep22/1.SBP/EUR_Evan/SBP_statusAHM321_EUR_Evan_PRSout.best"),collapse=',')
#prsice2_header<-paste(c('Page','Evangelou'), collapse=',')


#listfile_prs<-paste(c("DBP/AfriPRS_MetaGU_RE2/DBP_statusAHM321_AfriPRS_MetaGU_RE2_PRSout.best","DBP/EUR_Evan/DBP_statusAHM321_EUR_Evan_PRSout.best", "DBP/Multi_PAGE/DBP_statusAHM321__PRSout.best"), collapse=',')
#prsice2_header<-paste(c('African','European', 'Multi'), collapse=',')
#pheno<-"DBP_AHM"
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

listfile_prs<-opt[['best_prs']]
prsice2_header<-opt[['head_bestprs']]
pheno<-opt[['pheno']]
filedata<-opt[['data']]
covariable<-opt[['covariables']];covariable_tmp<-strsplit(covariable, split=',')[[1]]
type_coeff<-opt[['type_out']]
NQuant<-opt[['quantile_number']]
out=opt[['out']]
fcttr1<-opt[['tr_var']]
fcttr2<-opt[['tr_res']]


pheno_tr<-paste(pheno,'res',sep='_')
data<-read.table(filedata, header=T)
if(!(pheno %in% names(data))){
cat('pheno ',pheno,' not in header of ', filedata, '\n')
cat('list data ', names(data),'\n')
cat('--- exit ----')
q('no', 3)
}
if(!all(covariable_tmp %in% names(data))){
cat('one or more covariable in not in data file \n')
cat('covariable with issue :',paste(covariable_tmp[!(covariable_tmp %in% names(data))],collapse=','), '\n')
q('no', 3)
}

#covariable<-paste("age","Age2","sex","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",sep=',')

## readding 
StuctPrscice=extract_prsicefile(listfile_prs, prsice2_header, filedata, pheno, covariable)
DataIPrscice<-extract_prscice_info(StuctPrscice$listefile,StuctPrscice$studies)


## add PRS
nullfct<-function(x)return(x)

# transforamtion of you 
dataset<-StuctPrscice$studies


debugvar<-T
var_struct<-'All'
lmres<-computed_lm(StuctPrscice, pheno,  StuctPrscice$covar,StuctPrscice$studies, get(fcttr1), get(fcttr2), NQuant, var_struct)
plot_prscice2(lmres$reslm,StuctPrscice$studies, type_coeff, 'All')
ggsave(paste(out,'_dist_quant.pdf',sep=''))
d<-do_lmprsice(lmres$data_tr,lmres$pheno_tr, dataset, var_struct)
plot_varexplained(d, StuctPrscice$studies, var_struct)
ggsave(paste(out,'_varexplained.pdf',sep=''))
write.csv(d, file=paste(out,'_varexplained.csv',sep=''))









