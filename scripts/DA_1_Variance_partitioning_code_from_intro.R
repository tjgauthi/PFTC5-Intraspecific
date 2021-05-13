#Code from our introduction
m.SLA<-lmer(log(SLA)~1+(1|Species/Individual/Leaves)+(1|Site/Plot),data=d,na.action=na.omit)
variances.SLA<-c(unlist(lapply(VarCorr(m.SLA),diag)),attr(VarCorr(m.LMA),”sc”)^2)
#raw variances
Variances.SLA
#%Variance 
(var.comp.SLA<-variances.SLA/sum(variances.SLA))
# this is the same as: print(VarCorr(m.LMA),comp=”Variance”)
