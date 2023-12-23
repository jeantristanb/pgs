source('../utils.r')
outext='sim_1'
if(!file.exists(paste(outext,'.phe',sep=''))){
Pheno=rnorm(10000)
HTN=sample(c(0,1), length(Pheno),prob=c(0.45,0.55),replace=T)
Sex=sample(c(0,1), length(Pheno),prob=c(0.7,0.3),replace=T)
Pheno=(Pheno-min(Pheno))/(max(Pheno)-min(Pheno))+1
Pheno=Pheno-HTN*0.2+ rnorm(length(Pheno),sd=0.1)*HTN 
Pheno=Pheno+Sex*0.2
Pheno=Pheno-min(Pheno)+1
PRS=Pheno+ rnorm(length(Pheno),sd=0.3)
PRS2=Pheno+ rnorm(length(Pheno),sd=0.4)
cor(Pheno,PRS)
FID=paste("IID",1:length(PRS),sep='')
dffinal=data.frame(FID=FID,IID=FID,Pheno=Pheno,PRS=PRS,PRS2=PRS2,Sex=Sex,HTN=HTN)
getlm(dffinal,'Pheno','PRS',NULL)
getlm(dffinal,'Pheno','PRS',NULL, cov=c('HTN','Sex'))
getlm(dffinal,'Pheno','PRS2',NULL, cov=c('HTN','Sex'))
write.table(dffinal[,c('FID','IID', 'Pheno','Sex','HTN')], sep='\t', row.names=F, col.names=T, quote=F, file=paste(outext,'.phe',sep=''))
write.table(dffinal[,c('FID','IID', 'PRS')], sep='\t', row.names=F, col.names=T, quote=F, file=paste(outext,'_sd3.best',sep=''))
dfbest<-dffinal[,c('FID','IID', 'PRS2')];names(dfbest)<-c('FID','IID', 'PRS')
write.table(dfbest, sep='\t', row.names=F, col.names=T, quote=F, file=paste(outext,'_sd4.best',sep=''))
}




