library(GECKO)

n1in<-round(mean(sumstat_1$N))
n2in<-round(mean(sumstat_2$N))
nsin<-0
Weightin = T #is always set for GECKO to improve the efficiency
Fix_Vein = T #(if two studies have non-overlapped samples, otherwise false)
Test = T #Test for significance or not
###Need to specify the number of individuals within each study: n1,n2, and number of the overlapping individuals in the two studies:nsin
###if the two samples are from the separate studies, nsin = 0, and Fix_Vein = 1
###if the two samples are from the same study, nsin need to be specified

Result<-GECKO_R(sumstat_1,sumstat_2,n1in,n2in,nsin,ldscore,Weightin,Fix_Vein,Test)
