
##------------------------------------------------------------------------##
## This script reproduces the analyses done to Southern State variables   ##
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
#Model Southern State variables
##----------------------------------------------
South= CIdata %>% filter(NewRegion!="central") %>%
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
S_corpred<-round(cor(South[,cpreds],use="complete.obs"),2)

##Plot and correlation values
sc=corrplot(S_corpred, method="number",type = "upper", tl.cex = 0.6,tl.col="black",number.cex=0.55,title="South reefs")## 

#Supplementary Figure S1B
jpeg(file = "CorrSouth.jpeg")
corrplot(S_corpred, method="number",type = "upper", tl.cex = 0.6,tl.col="black",number.cex=0.55,title="South reefs",mar = c(0,0,1,0), number.digits = 2)
dev.off()

###############################################
## based on correlations keep these predictors:
###############################################
South.preds=c("DHW_exp","COTS_exp","Cycl_exp",#Acute
              "CIacute","TCI5","TCI5_COTS","TCI5_cycl","TCI5_dhw",#Acute cumulative
              "Sal","Mud","DINe","Chla","PAR","DINriv",## Environ
              "Dist.coast","time","LONGITUDE","LATITUDE") # ALL continuous predictors.


## model cover as proportions with beta-distribution
# transform fraction with 0 and 1 to <1 and >0
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}


## Transform and select state variables for modeling based on correlations
Southdat=South %>%
  mutate(sqrt.DHW=sqrt(DHW_exp),
         sqrt.COTS=sqrt(COTS_exp),
         sqrt.Cycl=sqrt(Cycl_exp),
         sqrt.DINe=sqrt(DINe),
         sqrt.MUD=sqrt(Mud),
         sqrt.sal=sqrt(Sal)%>%
         log.DINriv=log(DINriv+1), 
         mPAR=(PAR^3),
         ACBX=transform01(ACBX/100),
         ACTO=transform01(ACTO/100),
         CBRN=transform01(CBRN/100),
         MSE=transform01(MSE/100),
         Total=transform01(COVER.HC/100))%>%
  dplyr::select(c(LATITUDE,zone,FULLREEF_ID,REEF,REEF_SITE_NO,time,sqrt.DHW,sqrt.COTS,sqrt.Cycl,CIacute,TCI5,
                  COTSpa,TCI5_COTS,TCI5_cycl,TCI5_dhw,sqrt.sal,sqrt.MUD,sqrt.DINe,Chla,mPAR,log.DINriv,Dist.coast,
                  ACTO,MSE,ACBX,CBRN,Total))%>%
  droplevels()


#----------------------------------------------------------------------
##Full-subset analyses (Run on HPC separately for each coral group) 
#----------------------------------------------------------------------

## updated predictor list (only continues-not null variables)
south.trans.preds=c("sqrt.DHW","sqrt.COTS","sqrt.Cycl","CIacute","TCI5","TCI5_COTS",
                    "TCI5_cycl","TCI5_dhw","sqrt.sal","sqrt.MUD","sqrt.DINe","Chla",
                    "mPAR","log.DINriv","Dist.coast")

factor.vars= c("zone","COTSpa")

###############
## Total hard-coral cover 
Model.Total=gam(Total~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                  s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Southdat)

STotal_model.set=generate.model.set(use.dat=Southdat,
                                    max.predictors=5, 
                                    test.fit=Model.Total, 
                                    k=3, 
                                    pred.vars.cont=south.trans.preds,
                                    pred.vars.fact = factor.vars,
                                    factor.smooth.interactions=F,
                                    smooth.smooth.interactions=F,
                                    cov.cutoff=0.28,
                                    null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

STotal.list=fit.model.set(STotal_model.set)## FITS mod
##extract Variable Importance Scores
TOTAL.var.imp=STotal.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.STotal=STotal.list$mod.data.out
mod.table.STotal=mod.table.STotal[order(mod.table.STotal$AICc),]
STotal.less.2AICc=mod.table.STotal[which(mod.table.STotal$delta.AICc<2),]


###############
## Acropora tabular
Model.ACTO=gam(ACTO~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                 s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Southdat)

SActo_model.set=generate.model.set(use.dat=Southdat,
                                   max.predictors=5, 
                                   test.fit=Model.ACTO, 
                                   k=3, 
                                   pred.vars.cont=south.trans.preds,
                                   pred.vars.fact = factor.vars,
                                   factor.smooth.interactions=F,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

SActo.list=fit.model.set(SActo_model.set)## FITS mod
##extract Variable Importance Scores
ACTO.var.imp=SActo.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.SActo=SActo.list$mod.data.out
mod.table.SActo=mod.table.SActo[order(mod.table.SActo$AICc),]
SActo.less.2AICc=mod.table.SActo[which(mod.table.SActo$delta.AICc<2),]


###############
##Model of Acropora branching
Model.ACBX=gam(ACBX~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                 s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Southdat)

SAcbr_model.set=generate.model.set(use.dat=Southdat,
                                   max.predictors=5, 
                                   test.fit=Model.ACBR, 
                                   k=3, 
                                   pred.vars.cont=south.trans.preds,
                                   pred.vars.fact = factor.vars,
                                   factor.smooth.interactions=F,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

SAcbr.list=fit.model.set(SAcbr_model.set)## FITS mod
##extract Variable Importance Scores
ACBR.var.imp=SAcbr.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.SAcbr=SAcbr.list$mod.data.out
mod.table.SAcbr=mod.table.SAcbr[order(mod.table.SAcbr$AICc),]
SAcbr.less.2AICc=mod.table.SAcbr[which(mod.table.SAcbr$delta.AICc<2),]

###############
## Other non-branching corals
Model.CBRN=gam(CBRN~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                 s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Southdat)

Scbrn_model.set=generate.model.set(use.dat=Southdat,
                                   max.predictors=5, 
                                   test.fit=Model.CBRN, 
                                   k=3, 
                                   pred.vars.cont=south.trans.preds,
                                   pred.vars.fact = factor.vars,
                                   factor.smooth.interactions=F,
                                   smooth.smooth.interactions=F,
                                   cov.cutoff=0.28,
                                   null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

Scbrn.list=fit.model.set(Scbrn_model.set)## FITS mod
##extract Variable Importance Scores
CBRN.var.imp=Scbrn.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.Scbrn=Scbrn.list$mod.data.out
mod.table.Scbrn=mod.table.Scbrn[order(mod.table.Scbrn$AICc),]
Scbrn.less.2AICc=mod.table.Scbrn[which(mod.table.Scbrn$delta.AICc<2),]


###############
## Massive-submassive-encrusting
Model.MSE=gam(MSE~s(sqrt.DHW)+ s(time,bs='cr',k=5)+s(LATITUDE,k=5)+
                s(FULLREEF_ID,bs='re')+ s(REEF_SITE_NO,bs='re'),family=betar(link="logit"),data=Southdat)

SMse_model.set=generate.model.set(use.dat=Southdat,
                                  max.predictors=5, 
                                  test.fit=Model.MSE, 
                                  k=3, 
                                  pred.vars.cont=south.trans.preds,
                                  pred.vars.fact = factor.vars,
                                  factor.smooth.interactions=F,
                                  smooth.smooth.interactions=F,
                                  cov.cutoff=0.28,
                                  null.terms="s(time,bs='cr',k=5)+s(LATITUDE,k=5)+s(FULLREEF_ID,bs='re')+s(REEF_SITE_NO,bs='re')")

SMse.list=fit.model.set(SMse_model.set)## FITS mod
##extract Variable Importance Scores
MSE.var.imp=SMse.list$variable.importance$aic$variable.weights.raw
# Model selection table and best models
mod.table.SMse=SMse.list$mod.data.out
mod.table.SMse=mod.table.SMse[order(mod.table.SMse$AICc),]
SMse.less.2AICc=mod.table.SMse[which(mod.table.SMse$delta.AICc<2),]



#-----------------------------------
##Gather data from Models for Variable Importance Scores
#-----------------------------------

ALLcorals.vis<- do.call("cbind",mget(ls(pattern = "var.imp*")))
SSF_vis=as.data.frame(ALLcorals.vis)
names(SSF_vis)<- substr(colnames(SSF_vis), 1, 4)
SSF_vis$predictor <- rownames(SSF_vis)


##Join all
ALLcorals.vis$predictor=dplyr::recode(SSF_vis$predictor,
                                      sqrt.COTS="COTS density",
                                      sqrt.Cycl="Cyclone exposure",
                                      sqrt.DHW="DHW",
                                      COTS="COTS outbreaks",
                                      TCI5_COTS="5yr-outbreaks",
                                      TCI5_cycl="5yr-storms",
                                      TCI5_dhw="5yr-bleaching risk",
                                      TCI5="5yr-acute events",
                                      CIacute="Acute events",
                                      sqrt.DINe="Environ DIN",
                                      mPAR="PAR",
                                      log.DINriv="River DIN",
                                      sqrt.sal="Salinity Index",
                                      sqrt.MUD="Sediments")

rownames(SSF_vis)=SSF_vis$predictor
##Order columns and remove predictor column
SSF_vis=SSF_vis[c(2,1,3,4,5)]
##Convert to matrix
SSF_vis <- data.matrix(SSF_vis)
##Transpose for heatmap
SSF_vis<-t(SSF_vis)
orderssf=c("COTS density","Cyclone exposure","DHW","COTS outbreaks","Acute events","5yr-outbreaks","5yr-storms",
           "5yr-bleaching risk","5yr-acute events","Chla","Environ DIN","PAR","River DIN","Salinity Index",
           "Sediments","zone","Dist.coast")


##Figure S4.B
dev.new()

ggsave(file="VarImp_SSF.jpeg",width = 9, height = 6, units = c("in"),dpi = 300,
       heatmap.2(SSF_VarImp[,orderssf],notecex=0.3,  dendrogram ="none",
                 col=colorRampPalette(c("white","lavender", "lightsteelblue1", "cornflowerblue", "slateblue4"))(15),
                 main="South state\n (1992-2017)",
                 trace="none",key.title = "",keysize=0.5,
                 notecol="black",key=T,
                 sepcolor = "black",margins=c(11,11), lhei=c(2,6),lwid=c(2,6),cexRow = 2,cexCol = 2,
                 Rowv=FALSE,Colv=FALSE))


#-----------------------------------
##Gather data for  Best Models
#-----------------------------------

SSF<- do.call("rbind",mget(ls(pattern = "2AICc*")))
SSFtable=SSF[c(1,3,9,11,5,7)]
SSFtable$Group <- substr(rownames(SSFtable), 1, 3)

##Table S2.South
SSFtable%<>%as.data.frame(row.names = NULL)%>%
  rename(predictors=modname)%>%
  mutate(Model="SSF")%>%
  arrange(match(Group, c("Total","ACT","ACB", "CBR","MSE", desc(delta.AICc))))%>%
  dplyr::select(8,7,everything())
##Select Only the Most parsimonious Models
SSFtable%<>% slice (34, 1, 3, 17, 19)


#-------------------------------------------------------------------------------
##Re-fit Best mod and Estimate Variance Contribution to select variables 
#-------------------------------------------------------------------------------

SSFtable$predictors[SSFtable$Group=="Total"]
#[1] sqrt.Cycl+TCI5+TCI5_dhw

mod=gam(Total~
          s(sqrt.Cycl, k = 3, bs = "cr") + 
          s(TCI5, k = 3, bs = "cr") +
          s(TCI5_dhw, k = 3, bs = "cr") +
          s(time, bs = "cr",k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO,bs = "re") + 
          s(LATITUDE, k = 5), family = betar(link = "logit"), data=Southdat,method="REML")

### remove each predictor at a time
mp_cyc=update(mod, ~.-s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_tc5=update(mod, ~.-s(TCI5, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_dhw=update(mod, ~.-s(TCI5_dhw, k = 3, bs = "cr"),sp=mod$sp[-4])


listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:3)]
names(listmods) <- c('mcyl', 'mdhw', 'mtc5')


#function to estimate Explained variance per term 
devexp<-function(x,y){
  ((summary(y)$dev.expl)-(summary(x)$dev.expl))/(summary(y)$dev.expl)
}

pred_dev= as.data.frame(matrix(nrow = 3, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:3){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[4,1]="total"
pred_dev[4,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]

rm(list=setdiff(ls(), c("SSFtable","Southdat","mod","devexp")))


##----------------------------------
##Get Model predictions and plot
##----------------------------------


## pred Cycl
tdcycl <- expand.grid(sqrt.Cycl =seq(min(Southdat$sqrt.Cycl),max(Southdat$sqrt.Cycl),length.out = 50),
                      TCI5=median(mod$model$TCI5),
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

## pred TC5
tdc5 <- expand.grid(TCI5=seq(min(Southdat$TCI5),max(Southdat$TCI5),length.out = 50),
                    TCI5_dhw=median(mod$model$TCI5_dhw),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits2 <- predict.gam(mod, newdata=tdc5, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pci5 = tdc5%>%data.frame(fits2)%>%
  group_by(TCI5)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


## pred dhw
tdhw <- expand.grid(TCI5_dhw=seq(min(Southdat$TCI5_dhw),max(Southdat$TCI5_dhw),length.out = 50),
                    TCI5=median(mod$model$TCI5),
                    sqrt.Cycl=median(mod$model$sqrt.Cycl),
                    time=mean(mod$model$time),
                    REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                    FULLREEF_ID=(mod$model$FULLREEF_ID),
                    LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits3 <- predict.gam(mod, newdata=tdhw, type='response', se.fit=T, exclude= c('s(REEF_SITE_NO)','s(FULLREEF_ID)'))

pdhw = tdhw%>%data.frame(fits3)%>%
  group_by(TCI5_dhw)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##----------------------------------
## plot (Figure 3)
##----------------------------------


##Backtransform response for plotting
backtransform.est <- function(x, n) {
  y <- (x * n - 0.5) / (n - 1)
  return(y)
}


listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## 
names(listmods) <- c('pci5','pcyl','pdhw')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Southdat))*100),
                       lower=(backtransform.est(lower,nrow(Southdat))*100),
                       upper=(backtransform.est(upper,nrow(Southdat))*100))
}



png(file="total_ssf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

pci5=listmods[[1]]
plot(pci5$TCI5, pci5$response, type="n", lwd=3, ylim=(c(0,35)), xlim = (c(0,6)),
     main="",ylab="% Total cover",xlab = "No. of acute 5yr period",cex.axis=2,cex.lab=2)
polygon(c(pci5$TCI5, rev(pci5$TCI5)), 
        c(pci5$lower,rev(pci5$upper)), col="grey",
        border=NA)
lines(pci5$TCI5, pci5$response,  lwd=1)
quantile(Southdat$TCI5,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=5, col="black", lty=2)


pcyl=listmods[[2]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,35)),
     main="",ylab="% Total coral cover",xlab = "Cyclone exposure",cex.axis=2,cex.lab=2)
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Southdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=3, col="black", lty=2)


pdhw=listmods[[3]]
plot(pdhw$TCI5_dhw, pdhw$response, type="n", lwd=3, ylim=(c(0,60)),
     main="",ylab="% Total coral cover",xlab = "5yr Bleaching-risk",cex.axis=2,cex.lab=2)
polygon(c(pdhw$TCI5_dhw, rev(pdhw$TCI5_dhw)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)

lines(pdhw$TCI5_dhw, pdhw$response,  lwd=1)
quantile(Southdat$TCI5_dhw,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)
dev.off()



##----------------------------------
###ACROPORA TABULAR-Best Mod
##----------------------------------

SSFtable$predictors[SSFtable$Group=="ACT"]
#[1] sqrt.Cycl+sqrt.MUD+sqrt.sal+TCI5+TCI5_dhw


mod=gam(ACTO ~ s(sqrt.Cycl, k = 3, bs = "cr")+ s(sqrt.MUD, k = 3, bs = "cr")+s(sqrt.sal, k = 3, bs = "cr")
        +s(TCI5, k = 3, bs = "cr")+s(TCI5_dhw, k = 3, bs = "cr")+s(time, bs = "cr", k = 5)+
          s(FULLREEF_ID, bs = "re")+s(REEF_SITE_NO, bs = "re")+ s(LATITUDE, k = 5),
               family=betar(link="logit"),  data=Southdat, method="REML")

### remove each predictor at a time
mp_cycl= update(mod, ~. - s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_mud= update(mod, ~. - s(sqrt.MUD, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_sal= update(mod, ~. - s(sqrt.sal, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_tc5= update(mod, ~. - s(TCI5, k = 3, bs = "cr"),sp=mod$sp[-4])
mp_dhw5= update(mod, ~. - s(TCI5_dhw, k = 3, bs = "cr"),sp=mod$sp[-5])



listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:5)]
names(listmods) <- c('mcyc', 'mdhw5', 'mmud', 'msal', 'mtc5')

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
tdcycl <- expand.grid(sqrt.Cycl =seq(min(Southdat$sqrt.Cycl),max(Southdat$sqrt.Cycl),length.out = 50),
                      sqrt.MUD=median(mod$model$sqrt.MUD),
                      sqrt.sal=median(mod$model$sqrt.sal),
                      TCI5_dhw=median(mod$model$TCI5_dhw),
                      TCI5=median(mod$model$TCI5),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits1 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyc = tdcycl%>%data.frame(fits1)%>%
  group_by(sqrt.Cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

############
tdmud <- expand.grid(sqrt.MUD=seq(min(Southdat$sqrt.MUD),max(Southdat$sqrt.MUD),length.out = 50),
                     sqrt.Cycl=median(mod$model$sqrt.Cycl),
                     sqrt.sal=median(mod$model$sqrt.sal),
                     TCI5_dhw=median(mod$model$TCI5_dhw),
                     TCI5=median(mod$model$TCI5),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits2 <- predict.gam(mod, newdata=tdmud , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pmud = tdmud %>%data.frame(fits2)%>%
  group_by(sqrt.MUD)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

############
tdsal <- expand.grid(sqrt.sal=seq(min(Southdat$sqrt.sal),max(Southdat$sqrt.sal),length.out = 50),
                     sqrt.Cycl=median(mod$model$sqrt.Cycl),
                     sqrt.MUD=median(mod$model$sqrt.MUD),
                     TCI5_dhw=median(mod$model$TCI5_dhw),
                     TCI5=median(mod$model$TCI5),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits3 <- predict.gam(mod, newdata=tdsal , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

psal = tdsal %>%data.frame(fits3)%>%
  group_by(sqrt.sal)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



############
tc5 <- expand.grid(TCI5=seq(min(Southdat$TCI5),max(Southdat$TCI5),length.out = 50),
                   sqrt.Cycl=median(mod$model$sqrt.Cycl),
                   sqrt.MUD=median(mod$model$sqrt.MUD),
                   TCI5_dhw=median(mod$model$TCI5_dhw),
                   sqrt.sal=median(mod$model$sqrt.sal),
                   time=mean(mod$model$time),
                   REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                   FULLREEF_ID=(mod$model$FULLREEF_ID),
                   LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits4 <- predict.gam(mod, newdata=tc5 , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

ptc5 = tc5 %>%data.frame(fits4)%>%
  group_by(TCI5)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


############
tcdhw <- expand.grid(TCI5_dhw=seq(min(Southdat$TCI5_dhw),max(Southdat$TCI5_dhw),length.out = 50),
                     sqrt.Cycl=median(mod$model$sqrt.Cycl),
                     sqrt.MUD=median(mod$model$sqrt.MUD),
                     TCI5=median(mod$model$TCI5),
                     sqrt.sal=median(mod$model$sqrt.sal),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits4 <- predict.gam(mod, newdata=tcdhw , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pdhw = tcdhw %>%data.frame(fits4)%>%
  group_by(TCI5_dhw)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


##Plot mod
listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## list of all model fits with names
names(listmods) <- c('pcyc', 'pdhw','pmud','psal','ptc5', )

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Southdat))*100),
                       lower=(backtransform.est(lower,nrow(Southdat))*100),
                       upper=(backtransform.est(upper,nrow(Southdat))*100))
}


##PLOT
png(file="acto_ssf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

ptc5=listmods[[5]]
plot(ptc5$TCI5, ptc5$response, type="n", lwd=3, ylim=(c(0,8)), xlim = (c(0,6)),
     main="",ylab="Acropora Tabular",xlab = "5yr-acute events",cex.axis=2,cex.lab=2)
# confidence 
polygon(c(ptc5$TCI5, rev(ptc5$TCI5)), 
        c(ptc5$lower,rev(ptc5$upper)), col="grey",
        border=NA)
lines(ptc5$TCI5, ptc5$response,  lwd=1)
quantile(Southdat$TCI5,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=5, col="black", lty=2)


pcyc=listmods[[1]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,6)),
     main="",ylab="",
     xlab = "Cyclone exposure",cex.axis=2,cex.lab=2)
# confidence 
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Southdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=3, col="black", lty=2)


pdhw=listmods[[2]]
plot(pdhw$TCI5_dhw, pdhw$response, type="n", lwd=3, #ylim=(c(0,8)),
     main="",ylab="",xlab = "Bleaching-risk events",cex.axis=2,cex.lab=2)
# confidence 
polygon(c(pdhw$TCI5_dhw, rev(pdhw$TCI5_dhw)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$TCI5_dhw, pdhw$response,  lwd=1)
quantile(Southdat$TCI5_dhw,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=2, col="black", lty=2)
dev.off()


rm(list=setdiff(ls(), c("SSFtable","Southdat","backtransform.est","devexp")))


##----------------------------------
##ACROPORA BRANCHING-- Best Mod
##----------------------------------

SSFtable$predictors[SSFtable$Group=="ACB"]
#[1] sqrt.Cycl+sqrt.DHW+TCI5_COTS+TCI5_cycl


mod=gam(ACBX~ s(sqrt.Cycl, k = 3, bs = "cr") + 
          s(sqrt.DHW, k = 3, bs = "cr") +
          s(TCI5_COTS, k = 3, bs = "cr") +
          s(TCI5_cycl, k = 3, bs = "cr") +
          s(time, bs = "cr",k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO,bs = "re") + 
          s(LATITUDE, k = 5), family = betar(link = "logit"), 
        data=Southdat,method="REML")


### remove each predictor at a time
mp_cycl= update(mod, ~. - s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_dhw= update(mod, ~. - s(sqrt.DHW, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_cots5= update(mod, ~. - s(TCI5_COTS, k = 3,bs="cr"),sp=mod$sp[-3])
mp_cyc5= update(mod, ~. - s(TCI5_cycl, k = 3, bs="cr"),sp=mod$sp[-4])



listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:4)]
names(listmods) <- c('mcots', 'mcycl', 'mcyc5', 'mdhw')

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

## pred Cycl
tdcycl <- expand.grid(sqrt.Cycl =seq(min(Southdat$sqrt.Cycl),max(Southdat$sqrt.Cycl),length.out = 50),
                      sqrt.DHW=median(mod$model$sqrt.DHW),
                      TCI5_COTS=median(mod$model$TCI5_COTS),
                      TCI5_cycl=median(mod$model$TCI5_cycl),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits1 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcycl = tdcycl%>%data.frame(fits1)%>%
  group_by(sqrt.Cycl)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

############
tdtcots <- expand.grid(TCI5_COTS=seq(min(Southdat$TCI5_COTS),max(Southdat$TCI5_COTS),length.out = 40),
                       sqrt.DHW=median(mod$model$sqrt.DHW),
                       sqrt.Cycl=median(mod$model$sqrt.Cycl),
                       TCI5_cycl=median(mod$model$TCI5_cycl),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits2 <- predict.gam(mod, newdata=tdtcots, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcots= tdtcots %>%data.frame(fits2)%>%
  group_by(TCI5_COTS)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


############
tdtc5cy <- expand.grid(TCI5_cycl=seq(min(Southdat$TCI5_cycl),max(Southdat$TCI5_cycl),length.out = 40),
                       sqrt.DHW=median(mod$model$sqrt.DHW),
                       sqrt.Cycl=median(mod$model$sqrt.Cycl),
                       TCI5_COTS=median(mod$model$TCI5_COTS),
                       time=mean(mod$model$time),
                       REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                       FULLREEF_ID=(mod$model$FULLREEF_ID),
                       LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits3 <- predict.gam(mod, newdata=tdtc5cy, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyc5= tdtc5cy%>%data.frame(fits3)%>%
  group_by(TCI5_cycl)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


######
tddhw <- expand.grid(sqrt.DHW=seq(min(Southdat$sqrt.DHW),max(Southdat$sqrt.DHW),length.out = 40),
                     TCI5_cycl=median(mod$model$TCI5_cycl),
                     sqrt.Cycl=median(mod$model$sqrt.Cycl),
                     TCI5_COTS=median(mod$model$TCI5_COTS),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

## Exclude reef sites for prediction to reduce CI 
fits4 <- predict.gam(mod, newdata=tddhw, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pdhw = tddhw%>%data.frame(fits4)%>%
  group_by(sqrt.DHW)%>% #only change here
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()

##Plot Mods (Figure 3)
listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))##
names(listmods) <- c('pcots','pcycl','pcyc5','pdhw')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Southdat))*100),
                       lower=(backtransform.est(lower,nrow(Southdat))*100),
                       upper=(backtransform.est(upper,nrow(Southdat))*100))
}


png(file="acbx_csf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))


##PLOT
png(file="acbx_ssf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))


pcycl=listmods[[2]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,6)),
     main="",ylab="%Acropora branching",
     xlab = "Hours of exposure",cex.axis=2,cex.lab=2)
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Southdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=3, col="black", lty=2)


pcyc5=listmods[[2]]
plot(pcyc5$TCI5_cycl, pcyc5$response, type="n", lwd=3, ylim=(c(0,6)),
     main="",ylab="",xlab = "storms 5yr period",cex.axis=2,cex.lab=2)
# confidence 
polygon(c(pcyc5$TCI5_cycl, rev(pcyc5$TCI5_cycl)), 
        c(pcyc5$lower,rev(pcyc5$upper)), col="darkorange",
        border=NA)
lines(pcyc5$TCI5_cycl, pcyc5$response,  lwd=1)
quantile(Southdat$TCI5_cycl,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=3, col="black", lty=2)


pdhw=listmods[[4]]
plot(pdhw$sqrt.DHW^2, pdhw$response, type="n", lwd=3, ylim=(c(0,6)),
     main="",ylab="",xlab = "DHW",cex.axis=2,cex.lab=2)
polygon(c(pdhw$sqrt.DHW^2, rev(pdhw$sqrt.DHW^2)), 
        c(pdhw$lower,rev(pdhw$upper)), col="red",
        border=NA)
lines(pdhw$sqrt.DHW^2, pdhw$response,  lwd=1)
quantile(Southdat$sqrt.DHW^2,probs=c(0.50,0.95))
abline(v=0.16, col="black", lty=2)
abline(v=3.66, col="black", lty=2)


pcots=listmods[[1]]
plot(pcots$TCI5_COTS, pcots$response, type="n", lwd=3, ylim=(c(0,6)),
     main="",ylab="",xlab = "outbreak 5yr period",cex.axis=2,cex.lab=2)
polygon(c(pcots$TCI5_COTS, rev(pcots$TCI5_COTS)), 
        c(pcots$lower,rev(pcots$upper)), col="magenta",
        border=NA)
lines(pcots$TCI5_COTS, pcots$response,  lwd=1)
quantile(Southdat$TCI5_COTS,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=4, col="black", lty=2)
dev.off()


rm(list=setdiff(ls(), c("SSFtable","Southdat","backtransform.est","devexp")))


##----------------------------------
##--OTHER BRANCHING-- Best Mod
##----------------------------------

SSFtable$predictors[SSFtable$Group=="CBR"]
##[1] sqrt.Cycl+sqrt.DHW+TCI5

mod=gam(CBRN~s(sqrt.Cycl, k = 3, bs = "cr") + 
           s(sqrt.DHW, k = 3, bs = "cr") +
           s(TCI5, k = 3, bs = "cr") +
           s(time, bs = "cr",k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO,bs = "re") + 
           s(LATITUDE, k = 5), family = betar(link = "logit"), 
         data=Southdat,method="REML")


mp_cycl= update(mod, ~. - s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_dhw= update(mod, ~. - s(sqrt.DHW, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_tc5= update(mod, ~. - s(TCI5, k = 3,bs="cr"),sp=mod$sp[-3])


listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:3)]
names(listmods) <- c('mcyl','mdhw','mtc5')

pred_dev= as.data.frame(matrix(nrow = 3, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:3){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[4,1]="total"
pred_dev[4,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]


##Model Predictions

## pred Cycl
tdcycl <- expand.grid(sqrt.Cycl =seq(min(Southdat$sqrt.Cycl),max(Southdat$sqrt.Cycl),length.out = 50),
                      sqrt.DHW=median(mod$model$sqrt.DHW),
                      TCI5=median(mod$model$TCI5),
                      time=mean(mod$model$time),
                      REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                      FULLREEF_ID=(mod$model$FULLREEF_ID),
                      LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

fits1 <- predict.gam(mod, newdata=tdcycl, type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pcyl = tdcycl%>%data.frame(fits1)%>%
  group_by(sqrt.Cycl)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



tddhw <- expand.grid(sqrt.DHW=seq(min(Southdat$sqrt.DHW),max(Southdat$sqrt.DHW),length.out = 50),
                     sqrt.Cycl=median(mod$model$sqrt.Cycl),
                     TCI5=median(mod$model$TCI5),
                     time=mean(mod$model$time),
                     REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                     FULLREEF_ID=(mod$model$FULLREEF_ID),
                     LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()

 
fits2 <- predict.gam(mod, newdata=tddhw , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

pdhw = tddhw %>%data.frame(fits2)%>%
  group_by(sqrt.DHW)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()



tc5 <- expand.grid(TCI5=seq(min(Southdat$TCI5),max(Southdat$TCI5),length.out = 50),
                   sqrt.Cycl=median(mod$model$sqrt.Cycl),
                   sqrt.DHW=median(mod$model$sqrt.DHW),
                   time=mean(mod$model$time),
                   REEF_SITE_NO=(mod$model$REEF_SITE_NO),
                   FULLREEF_ID=(mod$model$FULLREEF_ID),
                   LATITUDE=mean(mod$model$LATITUDE))%>%
  distinct()%>%
  glimpse()


fits3 <- predict.gam(mod, newdata=tc5 , type='response', se.fit=T, exclude= c("s(REEF_SITE_NO)","s(FULLREEF_ID)"))

ptc5 = tc5 %>%data.frame(fits3)%>%
  group_by(TCI5)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  ungroup()


listmods=setNames(lapply(ls(pattern="\\p"), get), paste(ls(pattern="\\p")))## list of all model fits with names
names(listmods) <- c('pcyl','pdhw','ptc5')

for (i in names(listmods)){
  listmods[[i]]=mutate(listmods[[i]], response=(backtransform.est(response,nrow(Southdat))*100),
                       lower=(backtransform.est(lower,nrow(Southdat))*100),
                       upper=(backtransform.est(upper,nrow(Southdat))*100))
}

##PLOT
png(file="cbrn_ssf_pred.png",
    width=900, height=250,res =100)
par(mfrow=c(1,4))

ptc5=listmods[[3]]
plot(ptc5$TCI5, ptc5$response, type="n", lwd=3, ylim=(c(0,8)),
     main="",ylab="% other branching",xlab = "5yr-acute events",cex.axis=2,cex.lab=2)
polygon(c(ptc5$TCI5, rev(ptc5$TCI5)), 
        c(ptc5$lower,rev(ptc5$upper)), col="grey",
        border=NA)
lines(ptc5$TCI5, ptc5$response,  lwd=1)
quantile(Southdat$TCI5,probs=c(0.50,0.95))
abline(v=1, col="black", lty=2)
abline(v=5, col="black", lty=2)


pcyl=listmods[[1]]
plot(pcyl$sqrt.Cycl^2, pcyl$response, type="n", lwd=3, ylim=(c(0,8)),
     main="",ylab="",
     xlab = "Cyclone exposure",cex.axis=2,cex.lab=2)
# confidence 
polygon(c(pcyl$sqrt.Cycl^2, rev(pcyl$sqrt.Cycl^2)), 
        c(pcyl$lower,rev(pcyl$upper)), col="darkorange",
        border=NA)
lines(pcyl$sqrt.Cycl^2, pcyl$response,  lwd=1)
quantile(Southdat$sqrt.Cycl^2,probs=c(0.50,0.95))
abline(v=0, col="black", lty=2)
abline(v=3, col="black", lty=2)
dev.off()


rm(list=setdiff(ls(), c("CSFtable","Southdat","backtransform.est","devexp")))

##--------------------------------
###-- MSE-- BestMod
##--------------------------------
SSFtable$predictors[SSFtable$Group=="MSE"]
##[1] sqrt.Cycl+sqrt.sal+TCI5_COTS+TCI5_dhw


mod=gam(MSE~   s(sqrt.Cycl, k = 3, bs = "cr") + 
          s(sqrt.sal, k = 3, bs = "cr") +
          s(TCI5_COTS, k = 3, bs = "cr") +  
          s(TCI5_dhw, k = 3, bs = "cr") + 
          s(time, bs = "cr",k = 5) + s(FULLREEF_ID, bs = "re") + s(REEF_SITE_NO,bs = "re") + 
          s(LATITUDE, k = 5), family = betar(link = "logit"), 
        data=Southdat,method="REML")

mp_cycl= update(mod, ~. - s(sqrt.Cycl, k = 3, bs = "cr"),sp=mod$sp[-1])
mp_sal= update(mod, ~. - s(sqrt.sal, k = 3, bs = "cr"),sp=mod$sp[-2])
mp_tc5= update(mod, ~. - s(TCI5_COTS, k = 3, bs = "cr"),sp=mod$sp[-3])
mp_dhw= update(mod, ~. - s(TCI5_dhw, k = 3, bs = "cr"),sp=mod$sp[-4])


listmods=setNames(lapply(ls(pattern="\\mp_"), get), paste(ls(pattern="\\mp_")))## list of all model fits with names
listmods=listmods[c(1:4)]
names(listmods) <- c('mcyl', 'mdhw', 'msal', 'mtc5')



pred_dev= as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(pred_dev)=c("pred","devex")
for (i in 1:4){
  dev=devexp(listmods[[i]],mod)*100
  pred_dev[i, 1] <- names(listmods)[i]
  pred_dev[i, 2] <- dev
}
pred_dev[5,1]="total"
pred_dev[5,2]=pred_dev%>%summarize_if(is.numeric, sum, na.rm=TRUE)
pred_dev[]##MOd mostly explained by Null vars



