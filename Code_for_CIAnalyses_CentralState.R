
##------------------------------------------------------------------------##
## This script reproduces the analyses done to Central State variables    ##
## Full-subset analysis for each coral group with 1992-2017 data-set      ##
## Model selection, variable importance scores and variance contribution  ##
## for the selected predictors in the Best Models                         ##
##------------------------------------------------------------------------##


library(dplyr)
library(RCurl)
require(car)
require(doBy)
require(gplots)
require(RColorBrewer)
require(gridExtra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(mgcv)
library(ggpubr)
library(FSSgam)

##----------------------------------------------
#Load data-set for analyses of State variables
##----------------------------------------------
CIdata=read.csv("CIdata_full.csv") 


##----------------------------------------------
#Model Central State variables
##----------------------------------------------
Central= CIdata %>% filter(NewRegion=="central") %>%
  as.data.frame() %>% droplevels()

##Check Correlations among potential predictors
cpreds=c("DHW_exp","COTS_exp","Cycl_exp",#All acute and cumulative metrics
         "CIacute","TCI3","TCI5", #
         "TCI3_COTS","TCI5_COTS",
         "TCI3_cycl","TCI5_cycl", 
         "TCI3_dhw", "TCI5_dhw",
         "Sal","Mud","DINe","Chla","PAR","DINriv", "Rivinf","sedWQI", ##Environmental data  
         "Dist.coast","Dist.barr","Rel.Dist.shelf",#Spatial 
         "time","LONGITUDE","LATITUDE") 

##Correlation coefficients
C_corpred<-round(cor(Central[,cpreds],use="complete.obs"),2)
##Plot and correlation values
cc=corrplot(C_corpred, method="number",type = "upper", tl.cex = 0.6,tl.col="black",number.cex=0.55,title="Central reefs")## 
#C.Cor <- rcorr(as.matrix(Central[,cpreds],use="complete.obs"))
#C.Cor$P## p-values

#Supplementary Figure S1A
jpeg(file = "CorrCentral.jpeg")
corrplot(C_corpred, method="number",type = "upper", tl.cex = 0.6,tl.col="black",number.cex=0.55,title="Central reefs",mar = c(0,0,1,0), number.digits = 2)
dev.off()


## model cover as proportions with beta-distribution
# transform fraction with 0 and 1 to <1 and >0
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

## Transform and select state variables for modeling based on correlations
Centdat=Central %>%
  mutate(sqrt.DHW=sqrt(DHW_exp),
         sqrt.COTS=sqrt(COTS_exp),
         sqrt.Cycl=sqrt(Cycl_exp),
         sqrt.DINe=sqrt(DINe),
         sqrt.MUD=sqrt(Mud),
         log.sal=log(Sal+1),
         log.DINriv=log(DINriv+1), 
         mPAR=(PAR^3),
         ACBX=transform01(ACBX/100),
         ACTO=transform01(ACTO/100),
         CBRN=transform01(CBRN/100),
         MSE=transform01(MSE/100),
         Total=transform01(COVER.HC/100))%>%
  dplyr::select(c(LATITUDE,zone,FULLREEF_ID,REEF,REEF_SITE_NO,time,sqrt.DHW,sqrt.COTS,sqrt.Cycl,CIacute,TCI5,
                  COTSpa,TCI5_COTS,TCI5_cycl,TCI3_dhw,TCI5_dhw,log.sal,sqrt.MUD,sqrt.DINe,Chla,mPAR,log.DINriv,
                  ACTO,MSE,ACBX,CBRN,Total))%>%
  droplevels()


#----------------------------------------------------------------------
##Full-subset analyses (Run on HPC separately for each coral group) 
#----------------------------------------------------------------------

#predictors
Centrans.preds=c("sqrt.DHW","sqrt.COTS","sqrt.Cycl","CIacute","TCI5","TCI5_COTS","TCI5_cycl","TCI3_dhw","TCI5_dhw","log.sal","sqrt.MUD","sqrt.DINe","Chla","mPAR","log.DINriv")
factor.vars= c("zone","COTSpa")

###############
## Total hard-coral cover 
Model.Total=gam(Total~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                  s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Centdat)

CTotal_model.set=generate.model.set(use.dat=Centdat,
                                    max.predictors=5, 
                                    test.fit=Model.Total, 
                                    k=3, 
                                    pred.vars.cont=Centrans.preds,
                                    pred.vars.fact = factor.vars,
                                    factor.smooth.interactions=F,
                                    smooth.smooth.interactions=F,
                                    cov.cutoff=0.28,
                                    null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

CTotal.list=fit.model.set(CTotal_model.set)## FITS mod
##extract Variable Importance Scores
TOTAL.var.imp=CTotal.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CTotal=CTotal.list$mod.data.out
mod.table.CTotal=mod.table.CTotal[order(mod.table.CTotal$AICc),]
CTotal.less.2AICc=mod.table.CTotal[which(mod.table.CTotal$delta.AICc<2),]


###############
## Acropora tabular
Model.ACTO=gam(ACTO~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                 s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Centdat)

CActo_model.set=generate.model.set(use.dat=Centdat,
                                   max.predictors=5, 
                                   test.fit=Model.ACTO, 
                                   k=3, 
                                   pred.vars.cont=Centrans.preds,
                                   pred.vars.fact = factor.vars,
                                   factor.smooth.interactions=F,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

CActo.list=fit.model.set(CActo_model.set)## FITS mod
##extract Variable Importance Scores
ACTO.var.imp=CActo.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CActo=CActo.list$mod.data.out
mod.table.CActo=mod.table.CActo[order(mod.table.CActo$AICc),]
CActo.less.2AICc=mod.table.CActo[which(mod.table.CActo$delta.AICc<2),]


###############
##Model of Acropora branching
Model.ACBX=gam(ACBX~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                 s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Centdat)

CAcbr_model.set=generate.model.set(use.dat=Centdat,
                                   max.predictors=5, 
                                   test.fit=Model.ACBR, 
                                   k=3, 
                                   pred.vars.cont=Centrans.preds,
                                   pred.vars.fact = factor.vars,
                                   factor.smooth.interactions=F,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

CAcbr.list=fit.model.set(CAcbr_model.set)## FITS mod
##extract Variable Importance Scores
ACBR.var.imp=CAcbr.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CAcbr=CAcbr.list$mod.data.out
mod.table.CAcbr=mod.table.CAcbr[order(mod.table.CAcbr$AICc),]
CAcbr.less.2AICc=mod.table.CAcbr[which(mod.table.CAcbr$delta.AICc<2),]

###############
## Other non-branching corals
Model.CBRN=gam(CBRN~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                 s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Centdat)

Ccbrn_model.set=generate.model.set(use.dat=Centdat,
                                   max.predictors=5, 
                                   test.fit=Model.CBRN, 
                                   k=3, 
                                   pred.vars.cont=Centrans.preds,
                                   pred.vars.fact = factor.vars,
                                   factor.smooth.interactions=F,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

Ccbrn.list=fit.model.set(Ccbrn_model.set)## FITS mod
##extract Variable Importance Scores
CBRN.var.imp=Ccbrn.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.Ccbrn=Ccbrn.list$mod.data.out
mod.table.Ccbrn=mod.table.Ccbrn[order(mod.table.Ccbrn$AICc),]
Ccbrn.less.2AICc=mod.table.Ccbrn[which(mod.table.Ccbrn$delta.AICc<2),]


###############
## Massive-submassive-encrusting
Model.MSE=gam(MSE~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Centdat)

CMse_model.set=generate.model.set(use.dat=Centdat,
                                  max.predictors=5, 
                                  test.fit=Model.MSE, 
                                  k=3, 
                                  pred.vars.cont=Centrans.preds,
                                  pred.vars.fact = factor.vars,
                                  factor.smooth.interactions=F,
                                  smooth.smooth.interactions=F,
                                  cov.cutoff=0.28,
                                  null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

CMse.list=fit.model.set(CMse_model.set)## FITS mod
##extract Variable Importance Scores
MSE.var.imp=CMse.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.CMse=CMse.list$mod.data.out
mod.table.CMse=mod.table.CMse[order(mod.table.CMse$AICc),]
CMse.less.2AICc=mod.table.CMse[which(mod.table.CMse$delta.AICc<2),]


#-----------------------------------
##Gather data from Models for Variable Importance Scores
#-----------------------------------

ALLcorals.vis<- do.call("cbind",mget(ls(pattern = "var.imp*")))
CSF_vis=as.data.frame(ALLcorals.vis)
names(CSF_vis)<- substr(colnames(CSF_vis), 1, 4)
CSF_vis$predictor <- rownames(CSF_vis)

CSF_vis$predictor=dplyr::recode(CSF_vis$predictor,
                                      sqrt.COTS="COTS density",
                                      sqrt.Cycl="Cyclone exposure",
                                      sqrt.DHW="DHW",
                                      COTSpa="COTS outbreaks",
                                      CIacute="Acute events",
                                      TCI5_COTS="5yr-outbreaks",
                                      TCI5_cycl="5yr-storms",
                                      TCI5_dhw="5yr-bleaching risk",
                                      TCI3_dhw="3yr-bleaching risk",
                                      TCI5="5yr-acute events",
                                      sqrt.DINe="Environ DIN",
                                      mPAR="PAR",
                                      log.DINriv="River DIN",
                                      log.sal="Salinity Index",
                                      sqrt.MUD="Sediments")


rownames(CSF_vis)=CSF_vis$predictor
##Order columns and remove predictor column
CSF_vis=CSF_vis[c(2,1,3,4,5)]
##Convert to matrix
CSF_vis <- data.matrix(CSF_vis)
##Transpose for heatmap
CSF_vis<-t(CSF_vis)
predscsf=c("COTS density","Cyclone exposure","DHW","COTS outbreaks","Acute events","3yr-bleaching risk",
           "5yr-outbreaks","5yr-storms","5yr-bleaching risk","5yr-acute events","Chla","Environ DIN","PAR","River DIN",
           "Salinity Index","Sediments","zone") 

##Figure S4.A
dev.new()
# 
ggsave(file="VarImp_CSF.jpeg",width = 9, height = 6, units = c("in"),dpi = 300,
       heatmap.2(CSF_vis[,predscsf],notecex=0.3,  dendrogram ="none",
                 col=colorRampPalette(c("white","lavender", "lightsteelblue1", "cornflowerblue", "slateblue4"))(15),
                 main="Central state\n (1992-2017)",
                 trace="none",key.title = "",keysize=0.5,
                 notecol="black",key=T,
                 sepcolor = "black",margins=c(11,11), lhei=c(2,6),lwid=c(2,6),cexRow = 2,cexCol = 2,Rowv=FALSE,Colv=FALSE))



#-----------------------------------
##Gather data for  Best Models
#-----------------------------------

CSF<- do.call("rbind",mget(ls(pattern = "2AICc*")))
CSFtable=CSF[c(1,3,9,11,5,7)]
CSFtable$Group <- substr(rownames(CSFtable), 1, 3)

##Table S2.Central
CSFtable%<>%as.data.frame(row.names = NULL)%>%
  rename(predictors=modname)%>%
  mutate(Model="CSF")%>%
  arrange(match(Group, c("Total","ACT", "ACB", "CBR","MSE", desc(delta.AICc))))%>%
  dplyr::select(8,7,everything())%>%
  glimpse()



#-----------------------------------------------------------------------
##Re-fit Best mod and Estimate Variance Contribution to select variables 
#-----------------------------------------------------------------------

CSFtable$predictors[CSFtable$Group=="Total"]
#[1] log.DINriv+sqrt.Cycl+sqrt.DHW+TCI5_COTS+TCI5_cycl

mod=gam(Total~
          s(log.DINriv, k = 3, bs = "cr") +  
          s(sqrt.Cycl, k = 3, bs = "cr") + 
          s(sqrt.DHW, k = 3, bs = "cr") +
          s(TCI5_COTS, k = 3, bs = "cr") +
          s(TCI5_cycl, k = 3, bs = "cr") +
          s(time, bs = "cr",k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO,bs = "re") + 
          s(LATITUDE, k = 5), family = betar(link = "logit"), data=Centdat,method="REML")

### remove each predictor at a time
mp_dinr=update(mod, ~.-s(log.DINriv, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_cyc=update(mod, ~.-s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_dhw=update(mod, ~.-s(sqrt.DHW, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_cots=update(mod, ~.-s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-4])
mp_cyc5=update(mod, ~.-s(TCI5_cycl, k = 3, bs = "cr"),sp=mod$sp[-5])


listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))
listmods=listmods[c(1:5)]
names(listmods) <- c('mcots', 'mcyl', 'mcyc5', 'mdhw', 'mdin')


#function to estimate Explained variance per term 
devexp<-function(x,y){
  ((summary(y)$dev.expl)-(summary(x)$dev.expl))/(summary(y)$dev.expl)
}

pred_dev= as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:5){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[6,1]="total"
pred_dev[6,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]

rm(list=setdiff(ls(), c("CSFtable","Centdat","mod","devexp")))


##----------------------------------
##Get Model predictions and plot
##----------------------------------


## pred Cycl
tdcycl <- expand.grid(sqrt.Cycl =seq(min(Centdat$sqrt.Cycl),max(Centdat$sqrt.Cycl),length.out = 50),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      sqrt.DHW=median(mod$model$sqrt.DHW),
                      log.DINriv=median(mod$model$log.DINriv),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits1 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl = tdcycl%>%data.frame(fits1)%>%
  group_by(sqrt.Cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
  upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

## pred Cycl5
tdcycl5 <- expand.grid(TCI5_cycl=seq(min(Centdat$TCI5_cycl),max(Centdat$TCI5_cycl),length.out = 50),
                       TCI5_COTS=median(mod$model$TCI5_COTS),
                       sqrt.Cycl=median(mod$model$sqrt.Cycl),
                       sqrt.DHW=median(mod$model$sqrt.DHW),
                       log.DINriv=median(mod$model$log.DINriv),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits2 <- predict.gam(mod, newdata=tdcycl5, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl5 = tdcycl5%>%data.frame(fits2)%>%
  group_by(TCI5_cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred Cots5
tdcots5 <- expand.grid(TCI5_COTS=seq(min(Centdat$TCI5_COTS),max(Centdat$TCI5_COTS),length.out = 50),
                       TCI5_cycl=median(mod$model$TCI5_cycl),
                       sqrt.Cycl=median(mod$model$sqrt.Cycl),
                       sqrt.DHW=median(mod$model$sqrt.DHW),
                       log.DINriv=median(mod$model$log.DINriv),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits3 <- predict.gam(mod, newdata=tdcots5, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcots5 = tdcots5%>%data.frame(fits3)%>%
  group_by(TCI5_COTS)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred dhw
tdhw <- expand.grid(sqrt.DHW=seq(min(Centdat$sqrt.DHW),max(Centdat$sqrt.DHW),length.out = 50),
                    TCI5_cycl=median(mod$model$TCI5_cycl),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    TCI5_COTS=median(mod$model$TCI5_COTS),
                    log.DINriv=median(mod$model$log.DINriv),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits4 <- predict.gam(mod, newdata=tdhw, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pdhw = tdhw%>%data.frame(fits4)%>%
  group_by(sqrt.DHW)%>%
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred din
tdin <- expand.grid(log.DINriv=seq(min(Centdat$log.DINriv),max(Centdat$log.DINriv),length.out = 50),
                    TCI5_cycl=median(mod$model$TCI5_cycl),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    TCI5_COTS=median(mod$model$TCI5_COTS),
                    sqrt.DHW=median(mod$model$sqrt.DHW),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits5 <- predict.gam(mod, newdata=tdin, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pdin = tdin%>%data.frame(fits5)%>%
  group_by(log.DINriv)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##Backtransform response for plotting
backtransform.est <- function(x, n) {
  y <- (x * n - 0.5) / (n - 1)
  return(y)
}


listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pcots5', 'pcyl', 'pcycl5', 'pdhw', 'pdin')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat))*100),
                       lower=(backtransform.est(lower,nrow(Centdat))*100),
                       upper=(backtransform.est(upper,nrow(Centdat))*100))
}

##PLOT PREDICTED VALUES (Figure 3)

png(file="Total_csf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

##--COTS
pcots5=listmods[[1]]
plot(pcots5$TCI5_COTS, pcots5$response, type="n", lwd=3, ylim=(c(0,25)),
     main="",ylab="% Total coral cover",xlab = "5yr-Outbreaks",cex.axis=2,cex.lab=2)
polygon(c(pcots5$TCI5_COTS, rev(pcots5$TCI5_COTS)), 
        c(pcots5$lower,rev(pcots5$upper)), col="magenta",
        border=NA)
lines(pcots5$TCI5_COTS, pcots5$response,  lwd=1)
quantile(Centdat$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)


##--Cyclones
pcyl=listmods[[2]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,25)),
     main="",ylab="",xlab = "Cyclone exposure",cex.axis=2,cex.lab=2)
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Centdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=0, col="black", lty=2)

##STORMS
pcyl5=listmods[[3]]
plot(pcyl5$TCI5_cycl, pcyl5$response, type="n", lwd=3, ylim=(c(0,25)),
     main="",ylab="",xlab = "5yr-Storms",cex.axis=2,cex.lab=2)
polygon(c(pcyl5$TCI5_cycl, rev(pcyl5$TCI5_cycl)), 
        c(pcyl5$lower,rev(pcyl5$upper)), col="darkorange",
        border=NA)
lines(pcyl5$TCI5_cycl, pcyl5$response,  lwd=1)
quantile(Centdat$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)

##DHW
pdhw=listmods[[4]]
plot(pdhw$sqrt.DHW^2, pdhw$response, type="n", lwd=3, ylim=(c(0,25)),
     main="",ylab="",xlab = "DHW",cex.axis=2,cex.lab=2)
polygon(c(pdhw$sqrt.DHW^2, rev(pdhw$sqrt.DHW^2)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$sqrt.DHW^2, pdhw$response,  lwd=1)
quantile(Centdat$sqrt.DHW^2,probs=c(0.50,0.95))
abline(v=0.94, col="black", lty=2)
abline(v=5.29, col="black", lty=2)

rm(list=setdiff(ls(), c("CSFtable","Centdat","backtransform.est","devexp")))

##----------------------------------
###ACROPORA TABULAR-Best Mod
##----------------------------------
CSFtable$predictors[CSFtable$Group=="ACT"]
#[1] log.DINriv+sqrt.Cycl+sqrt.DHW+TCI5_COTS+TCI5_cycl

mod=gam(ACTO~
          s(log.DINriv,k=3, bs='cr') +  
          s(sqrt.Cycl,k=3, bs='cr') + 
          s(sqrt.DHW,k=3, bs='cr') +
          s(TCI5_COTS,k=3, bs='cr') +
          s(TCI5_cycl,k=3, bs='cr') +
          s(time, k=5, bs='cr') +s(FULLREEF_ID, bs='re') + s(REEF_SITE_NO, bs='re') + 
          s(LATITUDE, k=5), family = betar(link = 'logit'), 
        data=Centdat,method='REML')


##Explained variance per term
mp_dinr=update(mod, ~.-s(log.DINriv, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_cyc=update(mod, ~.-s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_dhw=update(mod, ~.-s(sqrt.DHW, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_cots=update(mod, ~.-s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-4])
mp_cyc5=update(mod, ~.-s(TCI5_cycl, k = 3, bs = "cr"),sp=mod$sp[-5])


listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:5)]
names(listmods) <- c('mcots', 'mcyl', 'mcyc5', 'mdhw', 'mdin')

pred_dev= as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:5){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[6,1]="total"
pred_dev[6,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)

pred_dev[]



##Model Predictions
## pred Cycl
tdcycl <- expand.grid(sqrt.Cycl=seq(min(Centdat$sqrt.Cycl),max(Centdat$sqrt.Cycl),length.out = 50),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      log.DINriv=median(mod$model$log.DINriv),
                      sqrt.DHW=median(mod$model$sqrt.DHW),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits1 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl = tdcycl%>%data.frame(fits1)%>%
  group_by(sqrt.Cycl)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


tdcots <- expand.grid(TCI5_COTS=seq(min(Centdat$TCI5_COTS),max(Centdat$TCI5_COTS),length.out = 50),
                      sqrt.Cycl=median(mod$model$sqrt.Cycl),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      log.DINriv=median(mod$model$log.DINriv),
                      sqrt.DHW=median(mod$model$sqrt.DHW),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits2 <- predict.gam(mod, newdata=tdcots , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))


pcots = tdcots%>%data.frame(fits2)%>%
  group_by(TCI5_COTS)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred 
tdin <- expand.grid(log.DINriv=seq(min(Centdat$log.DINriv),max(Centdat$log.DINriv),length.out = 50),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    TCI5_COTS=median(mod$model$TCI5_COTS),
                    TCI5_cycl=median(mod$model$TCI5_cycl),
                    sqrt.DHW=median(mod$model$sqrt.DHW),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits3 <- predict.gam(mod, newdata=tdin , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))


pdin = tdin%>%data.frame(fits3)%>%
  group_by(log.DINriv)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred DHW
tdhw <- expand.grid(sqrt.DHW=seq(min(Centdat$sqrt.DHW),max(Centdat$sqrt.DHW),length.out = 50),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    TCI5_COTS=median(mod$model$TCI5_COTS),
                    TCI5_cycl=median(mod$model$TCI5_cycl),
                    log.DINriv=median(mod$model$log.DINriv),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits4 <- predict.gam(mod, newdata=tdhw , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)'))

pdhw = tdhw%>%data.frame(fits4)%>%
  group_by(sqrt.DHW)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred Cycl5
tdcycl2 <- expand.grid(TCI5_cycl=seq(min(Centdat$TCI5_cycl),max(Centdat$TCI5_cycl),length.out = 50),
                       TCI5_COTS=median(mod$model$TCI5_COTS),
                       sqrt.Cycl=median(mod$model$sqrt.Cycl),
                       log.DINriv=median(mod$model$log.DINriv),
                       sqrt.DHW=median(mod$model$sqrt.DHW),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits5 <- predict.gam(mod, newdata=tdcycl2, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl2 = tdcycl2%>%data.frame(fits5)%>%
  group_by(TCI5_cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##Plot mod
listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## list of all model fits with names
names(listmods) <- c('pcots', 'pcyl','pcyl2','pdhw', 'pdin')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat))*100),
                       lower=(backtransform.est(lower,nrow(Centdat))*100),
                       upper=(backtransform.est(upper,nrow(Centdat))*100))
}


png(file="acto_csf_pred2.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

#---CoTS
pcots=listmods[[1]]
plot(pcots$TCI5_COTS, pcots$response, type="n", lwd=3, 
     ylim=(c(0,7)), cex.axis=2,cex.lab=2,
     main="",ylab="% cover Acropora tabular",xlab = "5yr Outbreaks")
polygon(c(pcots$TCI5_COTS, rev(pcots$TCI5_COTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$TCI5_COTS, pcots$response,  lwd=1)
quantile(Centdat$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)


#---Cyclones
pcyl=listmods[[2]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,7)), cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "Cyclone exposure")
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Centdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=0, col="black", lty=2)

##--storms
pcyl2=listmods[[3]]
plot(pcyl2$TCI5_cycl, pcyl2$response, type="n", lwd=3, ylim=(c(0,7)),
     main="",ylab="",xlab = "5yr-Storms",cex.axis=2,cex.lab=2)
polygon(c(pcyl2$TCI5_cycl, rev(pcyl2$TCI5_cycl)), 
        c(pcyl2$lower,rev(pcyl2$upper)), col="darkorange",
        border=NA)
lines(pcyl2$TCI5_cycl, pcyl2$response,  lwd=1)
quantile(Centdat$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)


##--DHW
pdhw=listmods[[4]]
plot(pdhw$sqrt.DHW^2, pdhw$response, type="n", lwd=3, ylim=(c(0,7)), cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "DHW")
polygon(c(pdhw$sqrt.DHW^2, rev(pdhw$sqrt.DHW^2)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$sqrt.DHW^2, pdhw$response,  lwd=1)
quantile(Centdat$sqrt.DHW^2,probs=c(0.50,0.95))
abline(v=0.94, col="black", lty=2)
abline(v=5.29, col="black", lty=2)
dev.off()

rm(list=setdiff(ls(), c("CSFtable","Centdat","backtransform.est","devexp")))

##----------------------------------
##ACROPORA BRANCHING-- Best Mod
##----------------------------------
CSFtable$predictors[CSFtable$Group=="ACB"]
##sqrt.Cycl+TCI5_COTS+TCI5_cycl+TCI5_dhw (most parsimonious)

mod=gam(ACBX ~ s(sqrt.Cycl, k = 3, bs = "cr") + s(TCI5_COTS, k = 3, bs = "cr")+ s(TCI5_cycl, k = 3, bs = "cr")+ 
          s(TCI5_dhw,k=3, bs= "cr")+ s(time, bs = "cr", k = 5)+s(FULLREEF_ID, bs = "re")+
          s(REEF_SITE_NO, bs = "re")+ s(LATITUDE, k = 5),family=betar(link="logit"),  data=Centdat, method="REML")


mp_cycl= update(mod, ~. - s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_cots= update(mod, ~. - s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_cycl5= update(mod, ~. - s(TCI5_cycl, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_dhw5= update(mod, ~. - s(TCI5_dhw, k = 3, bs = "cr"),sp=mod$sp[-4])


listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:4)]
names(listmods) <- c('mcots', 'mcyl', 'mcyc5', 'mdhw')

pred_dev= as.data.frame(matrix(nrow = 4, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:4){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[5,1]="total"
pred_dev[5,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]


##Model Predictions

tdcycl <- expand.grid(sqrt.Cycl =seq(min(Centdat$sqrt.Cycl),max(Centdat$sqrt.Cycl),length.out = 50),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      TCI5_dhw=median(mod$model$TCI5_dhw),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits1 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl = tdcycl%>%data.frame(fits1)%>%
  group_by(sqrt.Cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



tdcots <- expand.grid(TCI5_COTS=seq(min(Centdat$TCI5_COTS),max(Centdat$TCI5_COTS),length.out = 50),
                      sqrt.Cycl=median(mod$model$sqrt.Cycl),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      TCI5_dhw=median(mod$model$TCI5_dhw),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits2 <- predict.gam(mod, newdata=tdcots , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcots = tdcots%>%data.frame(fits2)%>%
  group_by(TCI5_COTS)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



tdcycl <- expand.grid(TCI5_cycl=seq(min(Centdat$TCI5_cycl),max(Centdat$TCI5_cycl),length.out = 50),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      sqrt.Cycl=median(mod$model$sqrt.Cycl),
                      TCI5_dhw=median(mod$model$TCI5_dhw),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits3 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl2 = tdcycl%>%data.frame(fits3)%>%
  group_by(TCI5_cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


tdhw <- expand.grid(TCI5_dhw=seq(min(Centdat$TCI5_dhw),max(Centdat$TCI5_dhw),length.out = 50),
                    TCI5_COTS=median(mod$model$TCI5_COTS),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    TCI5_cycl=median(mod$model$TCI5_cycl),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits4 <- predict.gam(mod, newdata=tdhw , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)'))


pdhw = tdhw%>%data.frame(fits4)%>%
  group_by(TCI5_dhw)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##Plot Mods (Fifure 3)
listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))##
names(listmods) <- c('pcots','pcyl','pdhw','pdin')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat))*100),
                       lower=(backtransform.est(lower,nrow(Centdat))*100),
                       upper=(backtransform.est(upper,nrow(Centdat))*100))
}


png(file="acbx_csf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))


##--COTS
pcots=listmods[[1]]
plot(pcots$TCI5_COTS, pcots$response, type="n", lwd=3, ylim=(c(0,4)),cex.axis=2,cex.lab=2,
     main="",ylab="% cover Acropora branching",xlab = "5yr-Outbreaks")
polygon(c(pcots$TCI5_COTS, rev(pcots$TCI5_COTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$TCI5_COTS, pcots$response,  lwd=1)
quantile(Centdat$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)


##pcyl
pcyl=listmods[[2]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,4)),cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "Cyclone exposure")
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Centdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=0, col="black", lty=2)


##pcyl2
pcyl2=listmods[[3]]
plot(pcyl2$TCI5_cycl, pcyl2$response, type="n", lwd=3, ylim=(c(0,4)), cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "5y-Cyclone events")
# confidence 
polygon(c(pcyl2$TCI5_cycl, rev(pcyl2$TCI5_cycl)), 
        c(pcyl2$lower,rev(pcyl2$upper)), col="darkorange",
        border=NA)
lines(pcyl2$TCI5_cycl, pcyl2$response,  lwd=1)
quantile(Centdat$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)

######
pdhw=listmods[[4]]
plot(pdhw$TCI5_dhw, pdhw$response, type="n", lwd=3, ylim=(c(0,4)),cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "5yr-Bleaching-risk")
polygon(c(pdhw$TCI5_dhw, rev(pdhw$TCI5_dhw)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$TCI5_dhw, pdhw$response,  lwd=1)
quantile(Centdat$TCI5_dhw,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=1, col="black", lty=2)
dev.off()

rm(list=setdiff(ls(), c("CSFtable","Centdat","backtransform.est","devexp")))


##----------------------------------
##--OTHER BRANCHING-- Best Mod
##----------------------------------

CSFtable$predictors[CSFtable$Group=="CBR"]
##[1] log.DINriv+sqrt.Cycl+sqrt.DHW+TCI5_COTS+TCI5_cycl


mod=gam(CBRN ~ s(log.DINriv, k = 3, bs = "cr") + s(sqrt.Cycl, k = 3,bs = "cr") + s(sqrt.DHW, k = 3, bs = "cr")+ 
                s(TCI5_COTS, k = 3, bs = "cr") + s(TCI5_cycl, k = 3, bs = "cr") + 
                s(time, bs = "cr", k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO, bs = "re") + 
                s(LATITUDE, k = 5),family = betar(link = "logit"),  data=Centdat, method="REML")


mp_dine= update(mod, ~. - s(log.DINriv, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_cycl= update(mod, ~. - s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_dhw= update(mod, ~. - s(sqrt.DHW, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_cots5= update(mod, ~. - s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-4])
mp_cycl5= update(mod, ~. - s(TCI5_cycl, k = 3, bs = "cr"),sp=mod$sp[-5])



listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:5)]
names(listmods) <- c('mcots', 'mcyl', 'mcyc5', 'mdhw','mdine')

pred_dev= as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:5){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[6,1]="total"
pred_dev[6,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]


##Model Predictions

tdtc5cy <- expand.grid(TCI5_cycl=seq(min(Centdat$TCI5_cycl),max(Centdat$TCI5_cycl),length.out = 50),
                       sqrt.Cycl =mean(mod$model$sqrt.Cycl),
                       sqrt.DHW=mean(mod$model$sqrt.DHW),
                       log.DINriv=mean(mod$model$log.DINriv),
                       TCI5_COTS=mean(mod$model$TCI5_COTS),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE),
                       LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits <- predict.gam(mod, newdata=tdtc5cy, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))#

pcyc5 = tdtc5cy%>%data.frame(fits)%>%
  group_by(TCI5_cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



## pred Cycl
tdcycl <- expand.grid(sqrt.Cycl =seq(min(Centdat$sqrt.Cycl),max(Centdat$sqrt.Cycl),length.out = 50),
                      TCI5_cycl =mean(mod$model$TCI5_cycl),
                      sqrt.DHW=mean(mod$model$sqrt.DHW),
                      log.DINriv=mean(mod$model$log.DINriv),
                      TCI5_COTS=mean(mod$model$TCI5_COTS),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE),
                      LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()


fits <- predict.gam(mod, newdata=tdcycl, type='response',se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))#

pcyc= tdcycl%>%data.frame(fits)%>%
  group_by(sqrt.Cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

####
tdtc5co <- expand.grid(TCI5_COTS=seq(min(Centdat$TCI5_COTS),max(Centdat$TCI5_COTS),length.out = 50),
                       sqrt.Cycl =mean(mod$model$sqrt.Cycl),
                       sqrt.DHW=mean(mod$model$sqrt.DHW),
                       log.DINriv=mean(mod$model$log.DINriv),
                       TCI5_cycl=mean(mod$model$TCI5_cycl),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE),
                       LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()


fits <- predict.gam(mod, newdata=tdtc5co, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))#

pcots5= tdtc5co%>%data.frame(fits)%>%
  group_by(TCI5_COTS)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##DHW
tddhw <- expand.grid(sqrt.DHW=seq(min(Centdat$sqrt.DHW),max(Centdat$sqrt.DHW),length.out = 50),
                     sqrt.Cycl =mean(mod$model$sqrt.Cycl),
                     TCI5_COTS=mean(mod$model$TCI5_COTS),
                     log.DINriv=mean(mod$model$log.DINriv),
                     TCI5_cycl=mean(mod$model$TCI5_cycl),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE),
                     LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(mod, newdata=tddhw, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))#

pdhw= tddhw%>%data.frame(fits)%>%
  group_by(sqrt.DHW)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred DINr
tdin<- expand.grid(log.DINriv =seq(min(Centdat$log.DINriv),max(Centdat$log.DINriv),length.out = 50),
                      sqrt.Cycl =mean(mod$model$sqrt.Cycl),
                      TCI5_COTS=mean(mod$model$TCI5_COTS),
                      sqrt.DHW=mean(mod$model$sqrt.DHW),
                      TCI5_cycl=mean(mod$model$TCI5_cycl),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE),
                      LONGITUDE=mean(mod$model$LONGITUDE))%>%
  distinct()%>%
  glimpse()

fits <- predict.gam(Cbrnmod, newdata=tdin, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))#

pdinr= tdin%>%data.frame(fits)%>%
  group_by(log.DINriv)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## list of all model fits with names
names(listmods) <- c('pcots','pcyl','pcyl2','pdhw','pdin')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat))*100),
                       lower=(backtransform.est(lower,nrow(Centdat))*100),
                       upper=(backtransform.est(upper,nrow(Centdat))*100))
}

png(file="cbrn_csf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

##Plot Figure 3
##--Cyclones
pcyl=listmods[[2]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,8)),cex.axis=2,cex.lab=2,
     main="",ylab="% cover other branching",
     xlab = "Cyclone exposure")
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Centdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=0, col="black", lty=2)

##--Storms
pcyl2=listmods[[3]]
plot(pcyl2$TCI5_cycl, pcyl2$response, type="n", lwd=3, ylim=(c(0,8)),cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "5y- Storms")
polygon(c(pcyl2$TCI5_cycl, rev(pcyl2$TCI5_cycl)), 
        c(pcyl2$lower,rev(pcyl2$upper)), col="darkorange",
        border=NA)
lines(pcyl2$TCI5_cycl, pcyl2$response,  lwd=1)
quantile(Centdat$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)

##--DHW
pdhw=listmods[[4]]
plot(pdhw$sqrt.DHW^2, pdhw$response, type="n", lwd=3, ylim=(c(0,8)),cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "DHW")
polygon(c(pdhw$sqrt.DHW^2, rev(pdhw$sqrt.DHW^2)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$sqrt.DHW^2, pdhw$response,  lwd=1)
quantile(Centdat$sqrt.DHW^2,probs=c(0.50,0.95))
abline(v=0.94, col="black", lty=2)
abline(v=5.29, col="black", lty=2)
dev.off()


rm(list=setdiff(ls(), c("CSFtable","Centdat","backtransform.est","devexp")))

##--------------------------------
###-- MSE-- BestMod
##--------------------------------
CSFtable$predictors[CSFtable$Group=="MSE"]
##[1] sqrt.Cycl+sqrt.DHW+sqrt.DINe+TCI5_COTS+TCI5_cycl


mod=gam(MSE~s(sqrt.Cycl, k=3, bs='cr') + s(sqrt.DHW, k=3, bs='cr') + s(sqrt.DINe, k=3, bs='cr') + s(TCI5_COTS, k=3, bs='cr') +  
          s(TCI5_cycl, k=3, bs='cr') +s(time, k=5, bs='cr') +s(FULLREEF_ID, bs='re') + 
          s(REEF_SITE_NO, bs='re') + s(LATITUDE, k=5), family = betar(link = 'logit'),data=Centdat,method='REML')



mp_cycl= update(mod, ~. - s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_dhw= update(mod, ~. - s(sqrt.DHW, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_din= update(mod, ~. - s(sqrt.DINe, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_cots5= update(mod, ~. - s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-4])
mp_cycl5= update(mod, ~. - s(TCI5_cycl, k = 3, bs = "cr"),sp=mod$sp[-5])



listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:5)]
names(listmods) <- c('mcots', 'mcyl', 'mcyc5', 'mdhw', 'mdin')


#function to estimate Explained variance per term 
devexp<-function(x,y){
  ((summary(y)$dev.expl)-(summary(x)$dev.expl))/(summary(y)$dev.expl)
}

pred_dev= as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:5){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[6,1]="total"
pred_dev[6,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]

rm(list=setdiff(ls(), c("CSFtable","Centdat","mod","devexp")))



## pred Cycl
tdcycl <- expand.grid(sqrt.Cycl =seq(min(Centdat$sqrt.Cycl),max(Centdat$sqrt.Cycl),length.out = 50),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      sqrt.DINe=median(mod$model$sqrt.DINe),
                      sqrt.DHW=median(mod$model$sqrt.DHW),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl = tdcycl%>%data.frame(fits)%>%
  group_by(sqrt.Cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred 
tdcots <- expand.grid(TCI5_COTS=seq(min(Centdat$TCI5_COTS),max(Centdat$TCI5_COTS),length.out = 50),
                      sqrt.Cycl=median(mod$model$sqrt.Cycl),
                      sqrt.DINe=median(mod$model$sqrt.DINe),
                      sqrt.DHW=median(mod$model$sqrt.DHW),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()


## Exclude reef sites for prediction to reduce CI 
fits2 <- predict.gam(mod, newdata=tdcots , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcots = tdcots%>%data.frame(fits2)%>%
  group_by(TCI5_COTS)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



## pred 
tdin <- expand.grid(sqrt.DINe=seq(min(Centdat$sqrt.DINe),max(Centdat$sqrt.DINe),length.out = 50),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    TCI5_COTS=median(mod$model$TCI5_COTS),
                    sqrt.DHW=median(mod$model$sqrt.DHW),
                    TCI5_cycl=median(mod$model$TCI5_cycl),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits3 <- predict.gam(mod, newdata=tdin , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))


pdin = tdin%>%data.frame(fits3)%>%
  group_by(sqrt.DINe)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



tdhw <- expand.grid(sqrt.DHW=seq(min(Centdat$sqrt.DHW),max(Centdat$sqrt.DHW),length.out = 50),
                    TCI5_COTS=median(mod$model$TCI5_COTS),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    TCI5_cycl=median(mod$model$TCI5_cycl),
                    sqrt.DINe=median(mod$model$sqrt.DINe),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits4 <- predict.gam(mod, newdata=tdhw , type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)'))


pdhw = tdhw%>%data.frame(fits4)%>%
  group_by(sqrt.DHW)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred Cycl
tdcycl2 <- expand.grid(TCI5_cycl=seq(min(Centdat$TCI5_cycl),max(Centdat$TCI5_cycl),length.out = 50),
                       TCI5_COTS=median(mod$model$TCI5_COTS),
                       sqrt.Cycl=median(mod$model$sqrt.Cycl),
                       sqrt.DINe=median(mod$model$sqrt.DINe),
                       sqrt.DHW=median(mod$model$sqrt.DHW),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits5 <- predict.gam(mod, newdata=tdcycl2, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pcyl2 = tdcycl2%>%data.frame(fits5)%>%
  group_by(TCI5_cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## list of all model fits with names
names(listmods) <- c('pcots','pcyl','pcyl2','pdhw','pdin')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Centdat))*100),
                       lower=(backtransform.est(lower,nrow(Centdat))*100),
                       upper=(backtransform.est(upper,nrow(Centdat))*100))
}

##PLOT Figure 3
png(file="mse_csf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

##--COTS
pcots=listmods[[1]]

plot(pcots$TCI5_COTS, pcots$response, type="n", lwd=3, ylim=(c(0,10)),cex.axis=2,cex.lab=2,
     main="",ylab="% cover Massive/Encrusting",xlab = "5yr-Outbreaks")
polygon(c(pcots$TCI5_COTS, rev(pcots$TCI5_COTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$TCI5_COTS, pcots$response,  lwd=1)
quantile(Centdat$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)


##-Cyclones
pcyl=listmods[[2]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,10)),cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "Cyclone exposure")
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Centdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=0, col="black", lty=2)

##Storms
pcyl2=listmods[[3]]
plot(pcyl2$TCI5_cycl, pcyl2$response, type="n", lwd=3, ylim=(c(0,10)),cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "5y- Storm")
polygon(c(pcyl2$TCI5_cycl, rev(pcyl2$TCI5_cycl)), 
        c(pcyl2$lower,rev(pcyl2$upper)), col="darkorange",
        border=NA)
lines(pcyl2$TCI5_cycl, pcyl2$response,  lwd=1)
quantile(Centdat$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)


##--DHW
pdhw=listmods[[4]]
plot(pdhw$sqrt.DHW^2, pdhw$response, type="n", lwd=3, ylim=(c(0,10)),cex.axis=2,cex.lab=2,
     main="",ylab="",xlab = "DHW")
polygon(c(pdhw$sqrt.DHW^2, rev(pdhw$sqrt.DHW^2)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$sqrt.DHW^2, pdhw$response,  lwd=1)
quantile(Centdat$sqrt.DHW^2,probs=c(0.50,0.95))
abline(v=0.94, col="black", lty=2)
abline(v=5.29, col="black", lty=2)
dev.off()
