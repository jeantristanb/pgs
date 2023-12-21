require(psychometric)

plotvenn<-function(Data){
 listres=list('Admixture'=Data$FID[!is.na(Data$NoAdm) & !Data$NoAdm], 'Relatdness'=Data$FID[!is.na(Data$noRelatded) & !Data$noRelatded], 'Neigbour'=Data$FID[!is.na(Data$PassNeighbour) & !Data$PassNeighbour], 'Pcs'=Data$FID[!is.na(Data$PassPcsFilt) & !Data$PassPcsFilt], 'EigenStrat'=Data$FID[!is.na(Data$EigPcsFilter) & !Data$EigPcsFilter])
plotVenn(listres,showPlot=T, systemShow=F,outFile='tmpvendiag.svg')
}


#data.m<-DataInd;pheno<-PhenoPlot;NQuant<-10;cov<-Cov
get_quantile2<-function(x, NQuant){
 #cat('NA', length(x[is.na(x)]), '\n')
 tmp<-trunc(rank(x)/length(x)*NQuant)+1
 tmp[tmp>NQuant]<-NQuant
 #cat('Na2',paste(table(tmp,useNA='always'),collapse=','),'\n')
 tmp<-factor(tmp)
 levels(tmp)<-sort(unique(as.numeric(as.character(tmp))))
 tmp
}

get_quantile_prsice<- function(x, num.quant, quant.ref=NA){
    quant <- as.numeric(cut(x,
                            breaks = unique(quantile(
                                x, probs = seq(0, 1, 1 / num.quant)
                            )),
                            include.lowest = T))
    if(is.na(quant.ref) | is.null(quant.ref)){
        quant.ref <- ceiling(num.quant / 2)
    }
    quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
    return(quant)
}




get_infoquant<-function(data.m,pheno,Cov,NQuant){
 data.m<-na.omit(data.m[,c('FID','IID',pheno,Cov,'PRS')])
 data.m$quant<-get_quantile2(data.m$PRS, NQuant)
 AggrQuantMean<-aggregate(as.formula(paste(pheno,"~quant")),data=data.m, mean)
 AggrQuantSd<-aggregate(as.formula(paste(pheno,"~quant")),data=data.m, sd)
 AggrQuantLength<-aggregate(as.formula(paste(pheno,"~quant")),data=data.m, length)
 meansd<-merge(merge(AggrQuantMean,AggrQuantSd,by='quant'), AggrQuantLength,by='quant')
 names(meansd)<-c('quant','mean','sd','N')
 meansd$mean_l<-meansd$mean - 1.96*meansd$sd/sqrt(meansd$N)
 meansd$mean_u<-meansd$mean + 1.96*meansd$sd/sqrt(meansd$N)
 family <- gaussian
 reg <- summary(glm(as.formula(paste(pheno,"~ .")), family, data = na.omit(data.m[, c(pheno, 'quant', Cov)])))
 coef.quantiles <- c(0,reg$coefficients[2:nrow(reg$coefficients), 1])
 ci <- c(0,(1.96 * reg$coefficients[2:nrow(reg$coefficients), 2]))
 ci.quantiles.u <- coef.quantiles + ci
 ci.quantiles.l <- coef.quantiles - ci
 QNam<-names(coef.quantiles)
 QNam[grep('quant', QNam)]<-gsub('quant', '', QNam[grep('quant', QNam)])
 QNam[1]<-"1"
 CoefQuant<-data.frame(quant=QNam, Coef=coef.quantiles, ci.u=ci.quantiles.u, ci.l=ci.quantiles.l, OR=exp(coef.quantiles), OR_U=exp(ci.quantiles.u), OR_L=exp(ci.quantiles.l))
 merge(meansd,CoefQuant, by='quant')
}

plot_quantile<-function(quantile_resume, fill, type="coeff", Pop="All"){
  theme_sam <- theme_classic()+theme(axis.title=element_text(face="bold", size=18),
                              axis.text=element_text(size=14),
                              legend.title=element_text(face="bold", size=18),
                              legend.text=element_text(size=14),
                              axis.text.x=element_text(angle=45, hjust=1)
                              )
  if(type=="Coefficient"){y = "Coef";ymin = "ci.u";ymax = "ci.l";yhor=0}
  if(type=='Or'){
      y = "OR"
      ymin = "OR_L"
      ymax = "OR_U"
      yhor=1
  }
 if(type=='Mean'){
      y = "mean"
      ymin = "mean_l"
      ymax = "mean_u"
 }
 #[1] quant   mean    sd      Coef    ci.u    ci.l    studies
 if(Pop=='All')quantiles.plot <-     ggplot(quantile_resume, aes_string(x = "quant",y = y,ymin = ymin,ymax = ymax, fill=fill,color=fill))+ geom_point(width=0.2,position=position_dodge(width=0.5))+geom_errorbar(width = 0.2, position=position_dodge(width=0.5)) +theme_sam
 else {
	names(quantile_resume)[names(quantile_resume)=="Pop"]<-Pop
	if(length(unique(quantile_resume[,fill]))>1)quantiles.plot <-     ggplot(quantile_resume, aes_string(x = "quant",y = y,color=fill,fill=fill))+ geom_errorbar(position=position_dodge2(width=0.5),  aes_string(ymin = ymin,ymax = ymax), width=0.5)+geom_point(position=position_dodge2(width=0.5), aes_string(color=fill)) #+theme_sam 
	else  quantiles.plot <-     ggplot(quantile_resume, aes_string(x = "quant",y = y,ymin = ymin,ymax = ymax, fill=Pop,color=Pop))+ geom_point(position=position_dodge(width=0.5),aes_string(color=Pop, shape=Pop))+geom_errorbar(position=position_dodge(width=0.5)) +theme_sam
 }
 if(type %in% c('Or','Coefficient'))quantiles.plot<-quantiles.plot + geom_hline(yintercept=yhor, linetype="dashed", color = "red")
 return(quantiles.plot)
}

RankNormTr<-function(a)qnorm((rank(a)-0.5)/(length(a)-2*0.5+1))
nulltr=function(x)x
#Var<-var<-"og_friedewald_qc";fid<-"fid";data<-DataFilt

addresiduals<-function(data, var ,covar, varnewnam,fcttr1=nulltr, fcttr2=nulltr){
 data2<-data[,c(var, covar)]
 data2[,var]=fcttr1(data2[,var])
 data$num_tmp215<-paste('ind',1:nrow(data),sep='')
 rownames(data2)<-data$num_tmp215
 lmres<-lm(paste(var,"~",  paste(covar, collapse='+')) ,data=data2)
 tmpdf<-data.frame(names(lmres$residuals), fcttr2(lmres$residuals))
 names(tmpdf)<-c('num_tmp215',varnewnam)
 tmpdata<-merge(data,tmpdf,all=T,by="num_tmp215")
 tmpdata
}

getlm<-function(data1,head1,head2,head, cov=NULL){
if(debugvar)cat('getlm', head1, head2, head, cov,'\n')
 pval_format<-function(x){
  if(x>10**-2)return(as.character(round(x, 2)))
  else return(sprintf('%.1e', x))
 } 
 data1<-na.omit(data1[,c(head1,head2,cov)])
 lmres<-lm(as.formula(paste(head1,'~',paste(c(head2,cov), collapse="+"))), data=data1)
 tnlmres<-anova(lmres)
 tmpglm<-summary(lmres)
 estimate<-tmpglm$coefficients[2,1]
 pval<-pf(tmpglm$fstatistic[1L],tmpglm$fstatistic[2L], tmpglm$fstatistic[3L], lower.tail = FALSE)
 infoI<-CI.Rsq(tmpglm$r.squared, tmpglm$df[2],tmpglm$df[1])
 tmp<-data.frame(n=length(lmres$residuals), r2.adj_perc=round(tmpglm$r.squared*100,4),estimate=round(estimate,1),pval=pval_format(pval),lower=round(infoI[3]*100,1),upper=round(infoI[4]*100,1))
 if(!is.null(head))names(tmp)<-paste(names(tmp),head,sep='.')
 if(length(cov)>0){ 
 covvar2<-rownames(tmpglm$coefficients)[3:nrow(tmpglm$coefficients)]
 tmp2<-matrix(c(tmpglm$coefficients[3:nrow(tmpglm$coefficients),1], tmpglm$coefficients[3:nrow(tmpglm$coefficients),2], tmpglm$coefficients[3:nrow(tmpglm$coefficients),3]), ,nrow=1)
 colnames(tmp2)<- paste(rep(c('beta','se','pval'), each=length(Cov)),covvar2,sep="_")
 tmp<-cbind(tmp,tmp2)
 print(tmp)
}
tmp
}



extract_prsicefile<-function(listfileI, listheadI, filepheno, var, covar, info_ind=NULL, fcttr1=nulltr, fcttr2=nulltr){
 listfile=strsplit(listfileI, split=',')[[1]]
 listhead=strsplit(listheadI, split=',')[[1]]
 if(length(listfile)!=length(listhead)){
  cat('error list head != list file')
  cat('list of best prscie2 must be with a comma', listfile,'\n')
  cat('name of header of each file must be separated with a comma', listfile,'\n')
  quit('no', 5)
 }
 if(!is.null(covar))covar=strsplit(covar, split=',')[[1]]
 if(!is.null(info_ind))info_ind=strsplit(info_ind, split=',')[[1]]
 pheno=strsplit(var, split=',')[[1]]
 DataPheno<-read.table(filepheno, header=T)
 DataPheno<-DataPheno[,c(names(DataPheno)[1:2],pheno, covar, info_ind)]
 DataF<-DataPheno
 DataF2<-DataPheno
 Cmt<-1
 for(File in listfile){
    Data<-read.table(File, header=T)[,c(1,2,4)]
    names(Data)[3]<-listhead[Cmt]
    DataF<-merge(DataF,Data, by=c(1,2),all=T)
    Cmt<-Cmt+1
 }
 for(ind in info_ind)DataF[,ind]<-as.character(DataF[,ind])
 return(list(data=DataF, pheno=pheno, covar=covar, studies=listhead, listefile=listfile, info_ind=info_ind))
}


computed_lm_all<-function(DataAll, pheno,covar,liststudies, fcttr1, fcttr2, NQuant, residuals=T){
#[1] "FID"                    "IID"                    "MW"
#[4] "eGFR_sdfilter8_pcs5_15"
#[1] "FID"                          "IID"
#[3] "eGFR_sdfilter8_pcs5_15"       "eGFR_sdfilter8_pcs5_15_trans"  \
 listpheno<-c(pheno,covar)
 listpheno_not<-listpheno[!(listpheno %in% names(DataAll))]
 if(length(listpheno_not)>0){
   print(names(DataAll))
   print(listpheno_not)
   q(2)
 }

 
 if(debugvar)cat("computed_lm_all : end\n")
 varnewnam<-paste(pheno,'_trans',sep='')
 if(!is.null(covar))DataAll<-addresiduals(DataAll, pheno,covar, varnewnam,fcttr1, fcttr2)
 else {
	DataAll[,varnewnam]<-fcttr1(DataAll[,pheno])
 }
 if(residuals==T){
    covar=c()
    pheno=varnewnam
 }
 Cmt<-1
 for(studies in liststudies){
	   Data<-DataAll[,c(names(DataAll)[1:2],studies,pheno,covar)]
	   names(Data)[1:3]<-c('FID','IID', 'PRS')
	 if(debugvar)cat("computed_lm_all ",studies,": begin get_infoquant\n")
	   reslm<-get_infoquant(Data,pheno,covar,NQuant)
	 if(debugvar)cat("computed_lm_all ",studies,": end get_infoquant\n")
	   reslm$studies<-studies
           if(Cmt==1)allResLM<-reslm  
	   else allResLM<-rbind(allResLM,reslm)
	 if(debugvar)cat("computed_lm_all ",studies,": rbind\n")
           Cmt<-Cmt+1
 }
 if(debugvar)cat("computed_lm_all : end\n")
 return(list(reslm=allResLM, data_tr=DataAll, pheno_tr=varnewnam))
}

computed_lm<-function(Struct, pheno,covar,liststudies, fcttr1, fcttr2, NQuant,Pop){
 if(debugvar)print('begin computed_lm')
 DataAll<-Struct$data
 #print(Pop)
 #print(Pop=='All')
 if(is.null(Pop) || Pop=='All'){
	 print('not sub structure')
	 lmres=computed_lm_all(DataAll, pheno,covar,liststudies, fcttr1, fcttr2, NQuant)
	 lmres$pop<-"All"
	 return(list(reslm=lmres$reslm, data_tr=lmres$data_tr, pheno_tr=lmres$pheno_tr, pop=lmres$pop))
 }else {
	 Cmt<-1
	 DataAll<-DataAll[!is.na(DataAll[,Pop]),]
	 TmpRes<-computed_lm_all(DataAll, pheno,covar,liststudies, fcttr1, fcttr2, NQuant) 
	 AllResLM<-TmpRes$reslm;AllResLM$Pop<-"All"
	 DataTr<-TmpRes$data_tr;DataTr$Pop<-"All"
	 listpheno=TmpRes$pheno_tr
	 for(var in unique(DataAll[,Pop])){
	   if(debugvar)cat('computed_lm begin var ',var, '\n')
	   TmpRes<-computed_lm_all(DataAll[DataAll[,Pop]==var,], pheno,covar,liststudies, fcttr1, fcttr2, NQuant) 
	   PopResLM<-TmpRes$reslm;PopResLM$Pop<-var;AllResLM<-rbind(AllResLM,PopResLM)
	   PopTr<-TmpRes$data_tr;PopTr$Pop<-var;DataTr<-rbind(DataTr,PopTr)
	   if(debugvar)cat('computed_lm end',var, '\n')
	   Cmt=Cmt+1
	}
	if(debugvar)cat('end computed_lm', Pop, '\n')
	return(list(reslm=AllResLM, data_tr=DataTr, pheno_tr=listpheno, pop=Pop))
 }
}
subsample_lm<-function(str_lm, list_pop){
   if(is.null(list_pop))return(str_lm)
   if(length(list_pop)==0)return(str_lm)
   list(reslm=str_lm$reslm[str_lm$reslm$Pop %in% list_pop,], data_tr=str_lm$data_tr[str_lm$data_tr$Pop %in% list_pop ,], pheno_tr=str_lm$pheno_tr, pop=str_lm$pop)
}

plot_prscice2<-function(dglm,liststudies, type, Pop){
	cat('studies',dglm$studies,'\n')
	#print('list studies ',liststudies)
	plot_quantile(dglm[dglm$studies %in% liststudies,], 'studies', type, Pop)
}

plot_density<-function(data, pheno, Pop){
 if(Pop=='All')p<-ggplot(data, aes_string(x=pheno)) +  geom_density()
 else {
         data[,Pop]<-as.factor(data[,Pop])
	 p<-ggplot(data, aes_string(x=pheno, color=Pop, fill=Pop)) +  geom_density()
	}
 return(p) 
}

do_lmprsice_sub<-function(dglm, pheno, listfile){
 Cmt<-1
  if(debugvar==T){
	  print('begin do_lmprsice_sub')
 print(listfile)
 print(pheno)
 print(head(dglm))
  }
 for(studies in listfile){
      if(debugvar)cat('begin : getlm ',studies, '\n')
      dlm<-getlm(dglm,pheno,studies,NULL)
      if(debugvar)cat('end : getlm ',studies, '\n')
      if(Cmt==1)dlmf<-dlm
      else dlmf<-rbind(dlmf,dlm)
      if(debugvar)cat('end : getlm rbind ',studies, '\n')
      Cmt<-Cmt+1
 }
 dlmf$studies<-listfile
 return(dlmf)
}

do_lmprsice<-function(dglm, pheno, listfile, Pop){ 
  if(debugvar==T){
	  print('begin do_lmprsice')
  print(listfile)
  }
  if(Pop=='All'){
	  dlmall<-do_lmprsice_sub(dglm, pheno, listfile)
	  print(dlmall)
	  dlmall$Pop<-'All'
  }else {
    dglm2<-dglm[!is.na(dglm[,"Pop"]),]
    Cmt<-1
    for(var in unique(as.character(dglm2[,"Pop"]))){
       if(debugvar)cat('begin : do_lmprsice',var, '\n')
       dlmallsub<-do_lmprsice_sub(dglm2[dglm2[,"Pop"]==var,],pheno,listfile) 
       if(debugvar)cat('end do_lmprsice',var, '\n')
       dlmallsub$Pop<-var
       if(Cmt==1)dlmall<-dlmallsub else dlmall<-rbind(dlmall,dlmallsub)
       Cmt<-Cmt+1
    }
  }
 return(dlmall)
}


plot_varexplained<-function(dlmf, listfile, pop){
	dlmf<-dlmf[dlmf$studies %in% listfile,]
	if(pop=='All')p<-ggplot(dlmf, aes_string(x="studies", y="r2.adj_perc", fill="studies")) +geom_bar(stat="identity")+theme_minimal() +xlab('') + ylab('% r2')
	else {
		p<-ggplot(dlmf, aes_string(x="studies", y="r2.adj_perc", fill="Pop")) +geom_bar(stat="identity",position = "dodge")+theme_minimal() +  +xlab('') + ylab('% r2')

	}
        return(p)
}

extract_prscice_info<-function(listfile, listhead){
	Cmt<-1
	for(file in listfile){
	 filelog<-gsub('.best$', '.log', file)
          datalog<-readLines(filelog)
	  nbvarincluded<-strsplit(grep(" variant(s) included ", datalog,fixed=T,value=T),split=" ")[[1]][1]
	 filesummary<-gsub('.best$', '.summary', file)
          data<-read.table(filesummary, header=T, sep='\t')
	  data$nbvarincluded<-nbvarincluded
	  data$studies<-listhead[Cmt]
	  if(Cmt==1)dataF<-data
	  else dataF<-rbind(dataF,data)
          Cmt<-Cmt+1
	}
  dataF[,-1]
}

