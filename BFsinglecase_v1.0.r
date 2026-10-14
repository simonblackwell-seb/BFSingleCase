# Analyses for Study: "Using Bayesian Sequential Single Case Methodology 
# to Personalize Low-intensity Psychological Interventions for Depression: 
# Initial Pilot Work"

# Registration: https://clinicaltrials.gov/ct2/show/NCT04779437
# Study data/materials: https://osf.io/qcwfx/
#
# Contact: Simon Blackwell, simon.blackwell@uni-goettingen.de
#
# v1.0 13.10.26

#### load libraries ####
library(readxl)
library(xlsx)
library(devEMF)
library(ggplot2)
library(ggpubr)

# The BayesSingleSub package was version 0.6.3, commit 3f44cf5ec218d0c5da061b0291c2e13e950fdf8e
# from https://github.com/richarddmorey/BayesSingleSub
# If you can't work out a way to install it successfully, get in contact and I can help!

library(BayesSingleSub)

#### load data ####
mydata <- read_excel("singlecaseseries2023-10-12.xlsx")


#### add variable labels ####

mydata$gender <- factor(mydata$gender,levels = c("1","2"),labels = c("Male","Female"))
mydata$germanspeaker <- factor(mydata$germanspeaker,levels = c("1","2"),labels = c("Yes","No"))
mydata$education <- factor(mydata$education, levels = c("1","2","3","4","5","6"), labels = c("None","School-leaving certificate","Vocational college certificate","Bachelors degree","Masters degree","PhD"))
mydata$employment <- factor(mydata$employment,levels = c("1","2","3","4","5"), labels=c("Unemployed/Seeking employment","Student/trainee","Self-employed","Employed","Retired"))
mydata$currenttx <- factor(mydata$currenttx,levels = c("1","2"),labels=c("No","Yes"))
mydata$pasttx <- factor(mydata$pasttx,levels = c("1","2"),labels=c("No","Yes"))
mydata$everdiag <- factor(mydata$everdiag,levels = c("1","2"),labels=c("No","Yes"))
mydata$relstatus<-factor(mydata$relstatus,levels=c("1","2","3"),labels=c("Single","Single, in a stable relationship","Married/civil partnership/cohabiting"))
mydata$income<-factor(mydata$income,levels=c("1","2","3","4","5","6"),labels=c("up to 899€","900-2099€","2100-3099€","3100-4099€","4100-5099€","5100€ or more"))
mydata$nationality<-factor(mydata$nationality,levels=c("1","3"),labels=c("German","Polish"))



#### Reliabilities and Rescoring ####


#' Function to produce nicely formatted alpha values
#' @param name Name of the variable (for output display)
#' @param titems A vector of the names of the individual columns for the items
#' @param tdata The dataframe in which these items can be found
palpha<-function(name="?",titems,tdata=mydata){
  a<-psych::alpha(tdata[titems])
  capture.output(aci<-psych::alpha.ci(a$total$raw_alpha,n.obs = length(a$scores),n.var=a$nvar,digits=3),file='NUL')
  cat(name,": ",round(a$total$raw_alpha,2), " [",round(aci$lower.ci,2),",",round(aci$upper.ci,2),"]",sep="")
}

# Note that the QIDS needs to be rescored as the scoring algorithm in the online platform was incorrect

mydata<-transform(mydata,QIDSst1r_0=apply(mydata[c("QIDSst1_0","QIDSst2_0","QIDSst3_0","QIDSst4_0")],1,FUN=max))
mydata$QIDSst6r_0<-ifelse(mydata$QIDSst6_0>3,mydata$QIDSst6_0-5,mydata$QIDSst6_0)
mydata$QIDSst8r_0<-ifelse(mydata$QIDSst7_0>3,mydata$QIDSst7_0-4,mydata$QIDSst7_0)
mydata<-transform(mydata,QIDSst7r_0=apply(mydata[c("QIDSst6r_0","QIDSst8r_0")],1,FUN=max))
mydata<-transform(mydata,QIDSst14r_0=apply(mydata[c("QIDSst13_0","QIDSst14_0")],1,FUN=max))
mydata<-transform(mydata,QIDSsttotal_0=9*apply(mydata[c("QIDSst1r_0","QIDSst5_0","QIDSst7r_0","QIDSst8_0","QIDSst9_0","QIDSst10_0","QIDSst11_0","QIDSst12_0","QIDSst14r_0")],1,FUN=mean,na.rm=T))
palpha("QIDS baseline",c("QIDSst1r_0","QIDSst5_0","QIDSst7r_0","QIDSst8_0","QIDSst9_0","QIDSst10_0","QIDSst11_0","QIDSst12_0","QIDSst14r_0"),mydata)

mydata<-transform(mydata,QIDSst1r_14=apply(mydata[c("QIDSst1_14","QIDSst2_14","QIDSst3_14","QIDSst4_14")],1,FUN=max))
mydata$QIDSst6r_14<-ifelse(mydata$QIDSst6_14>3,mydata$QIDSst6_14-5,mydata$QIDSst6_14)
mydata$QIDSst8r_14<-ifelse(mydata$QIDSst7_14>3,mydata$QIDSst7_14-4,mydata$QIDSst7_14)
mydata<-transform(mydata,QIDSst7r_14=apply(mydata[c("QIDSst6r_14","QIDSst8r_14")],1,FUN=max))
mydata<-transform(mydata,QIDSst14r_14=apply(mydata[c("QIDSst13_14","QIDSst14_14")],1,FUN=max))
mydata<-transform(mydata,QIDSsttotal_14=9*apply(mydata[c("QIDSst1r_14","QIDSst5_14","QIDSst7r_14","QIDSst8_14","QIDSst9_14","QIDSst10_14","QIDSst11_14","QIDSst12_14","QIDSst14r_14")],1,FUN=mean,na.rm=T))
palpha("QIDS baseline",c("QIDSst1r_14","QIDSst5_14","QIDSst7r_14","QIDSst8_14","QIDSst9_14","QIDSst10_14","QIDSst11_14","QIDSst12_14","QIDSst14r_14"),mydata)

mydata<-transform(mydata,QIDSst1r_28=apply(mydata[c("QIDSst1_28","QIDSst2_28","QIDSst3_28","QIDSst4_28")],1,FUN=max))
mydata$QIDSst6r_28<-ifelse(mydata$QIDSst6_28>3,mydata$QIDSst6_28-5,mydata$QIDSst6_28)
mydata$QIDSst8r_28<-ifelse(mydata$QIDSst7_28>3,mydata$QIDSst7_28-4,mydata$QIDSst7_28)
mydata<-transform(mydata,QIDSst7r_28=apply(mydata[c("QIDSst6r_28","QIDSst8r_28")],1,FUN=max))
mydata<-transform(mydata,QIDSst14r_28=apply(mydata[c("QIDSst13_28","QIDSst14_28")],1,FUN=max))
mydata<-transform(mydata,QIDSsttotal_28=9*apply(mydata[c("QIDSst1r_28","QIDSst5_28","QIDSst7r_28","QIDSst8_28","QIDSst9_28","QIDSst10_28","QIDSst11_28","QIDSst12_28","QIDSst14r_28")],1,FUN=mean,na.rm=T))
palpha("QIDS baseline",c("QIDSst1r_28","QIDSst5_28","QIDSst7r_28","QIDSst8_28","QIDSst9_28","QIDSst10_28","QIDSst11_28","QIDSst12_28","QIDSst14r_28"),mydata)

mydata<-transform(mydata,QIDSst1r_42=apply(mydata[c("QIDSst1_42","QIDSst2_42","QIDSst3_42","QIDSst4_42")],1,FUN=max))
mydata$QIDSst6r_42<-ifelse(mydata$QIDSst6_42>3,mydata$QIDSst6_42-5,mydata$QIDSst6_42)
mydata$QIDSst8r_42<-ifelse(mydata$QIDSst7_42>3,mydata$QIDSst7_42-4,mydata$QIDSst7_42)
mydata<-transform(mydata,QIDSst7r_42=apply(mydata[c("QIDSst6r_42","QIDSst8r_42")],1,FUN=max))
mydata<-transform(mydata,QIDSst14r_42=apply(mydata[c("QIDSst13_42","QIDSst14_42")],1,FUN=max))
mydata<-transform(mydata,QIDSsttotal_42=9*apply(mydata[c("QIDSst1r_42","QIDSst5_42","QIDSst7r_42","QIDSst8_42","QIDSst9_42","QIDSst10_42","QIDSst11_42","QIDSst12_42","QIDSst14r_42")],1,FUN=mean,na.rm=T))
palpha("QIDS baseline",c("QIDSst1r_42","QIDSst5_42","QIDSst7r_42","QIDSst8_42","QIDSst9_42","QIDSst10_42","QIDSst11_42","QIDSst12_42","QIDSst14r_42"),mydata)

#Correct DARS item scoring from 1 to 5 to 0 to 4
mydata$DARStotal_0<-as.numeric(mydata$DARStotal_0)-17 #correct the 1 to 5 to 0 to 4 scoring
mydata$DARStotal_14<-as.numeric(mydata$DARStotal_14)-17 #correct the 1 to 5 to 0 to 4 scoring
mydata$DARStotal_28<-as.numeric(mydata$DARStotal_28)-17 #correct the 1 to 5 to 0 to 4 scoring
mydata$DARStotal_42<-as.numeric(mydata$DARStotal_42)-17 #correct the 1 to 5 to 0 to 4 scoring

mydata$GAD7total_0<-as.numeric(mydata$GAD7total_0)
RCI <- -6.49

mydata$QIDS_0A <- mydata$QIDSsttotal_14-mydata$QIDSsttotal_0
mydata$QIDS_AB <- mydata$QIDSsttotal_28-mydata$QIDSsttotal_14
mydata$QIDS_BC <- mydata$QIDSsttotal_42-mydata$QIDSsttotal_28
mydata$QIDS_0C <- mydata$QIDSsttotal_42-mydata$QIDSsttotal_0
mydata$RC_0A <- ifelse(mydata$QIDS_0A < RCI,TRUE,FALSE)
mydata$RC_AB <- ifelse(mydata$QIDS_AB < RCI,TRUE,FALSE)
mydata$RC_BC <- ifelse(mydata$QIDS_BC < RCI,TRUE,FALSE)
mydata$RC_0C <- ifelse(mydata$QIDS_0C < RCI,TRUE,FALSE)

#### Basic descriptives ####

vlist1<-c("age","gender","nationality","education","employment","relstatus","income","currenttx","pasttx","everdiag","QIDSsttotal_0","DARStotal_0","GAD7total_0")

#' Outputs Table 1 as dataframe
#' Continuous variables are presented as M(SD), categorical as N(%)
#' @param variablelist Vector of names of variables to be included in table
makeTable1<-function(variablelist)
{
  allrdata<-mydata
  #allrdata$arm<-"x"
  tHeaders<-c("",unique(allrdata$arm))
  narms<-length(unique(allrdata$arm))
  allRows<-vector()
  for (i in 1:narms)
  {
    allRows<-c(allRows,paste("n = ",table(allrdata$arm)[i][[1]],sep=""))
    
  }
  allRows<-rbind(c(" ",allRows),c(" ",rep("M (SD) / N (%)",narms)))
  
  for (i in 1:length(variablelist))
  {
    #if (variablelist[i]%in%labeldata$variable){
    #  thisRow<-labeldata[labeldata$variable==variablelist[i],]$label
    #}
    #else
    #{
      thisRow<-variablelist[i]
    #}
    
    #if numeric/integer work out mean and SD
    if (class(allrdata[variablelist[i]][[1]])=="integer"|class(allrdata[variablelist[i]][[1]])=="numeric"){
      for (j in 1:narms){
        thisRow<-c(thisRow,paste(round(mean(allrdata[allrdata$arm==unique(allrdata$arm)[j],variablelist[i]],na.rm=TRUE),2)," (",round(sd(allrdata[allrdata$arm==unique(allrdata$arm)[j],variablelist[i]],na.rm=TRUE),2),")",sep=""))
        
      }
      allRows<-rbind(allRows,thisRow)
    }else#if not numeric/integer do N(%)
    {
      thisRow<-c(thisRow,rep(" ",narms))
      allRows<-rbind(allRows,thisRow)
      for (k in 1:nlevels(allrdata[variablelist[i]][[1]])){
        theselevels<-levels(allrdata[variablelist[i]][[1]])
        thisRow<-theselevels[[k]]
        ltable<-table(allrdata[c(variablelist[i],"arm")])
        for (j in 1:narms){
          #tryCatch(
          #{
          thisRow<-c(thisRow,paste(ltable[theselevels[k],j]," (",round((prop.table(ltable[,j])[[theselevels[k]]]*100),2),"%)",sep=""))
          #},
          #error=function(cond) {
          #    thisRow<-c(thisRow,paste(0," (",0,"%)",sep=""))
          #})
          
          
        }
        allRows<-rbind(allRows,thisRow)
      }
      
    }
    
  }
  
  table1pt1<-data.frame(allRows)
  colnames(table1pt1)<-tHeaders
  return(table1pt1)
  
}#end of makeTable1<-function(variablelist)


#output table as excel file if wanted (to insert into paper)
#write.xlsx(makeTable1(vlist1),"Table1_13.10.26.xlsx",row.names = FALSE)
table(mydata$gender)
psych::describe(mydata$age)

#### Single Case Data ####

##### Daily QIDS #####

###### Load Data ######

qidsdata <- read_excel("dailyQIDS2023-10-12.xlsx")

#get unique participant ids and put in ascending order
pids <- unique (mydata$id)
pids <- pids[order(as.numeric(pids))]

#cycle through participants and produce daily depression analyses and graphs

# list for storing all participants daily data in, for later
alldailyqids <-list()
alldailyBFs <-list()


#prepare autocorrelation (lag 1) variables and max BF values
mydata$acfA <- NA
mydata$acfB <- NA
mydata$acfC <- NA
mydata$maxBFmAB <- NA
mydata$maxBFsAB <- NA
mydata$maxBFmBC <- NA
mydata$maxBFsBC <- NA
mydata$minBFmAB <- NA
mydata$minBFsAB <- NA
mydata$minBFmBC <- NA
mydata$minBFsBC <- NA
mydata$meanQIDSA <- NA
mydata$meanQIDSB <- NA
mydata$meanQIDSC <- NA
mydata$slopeQIDSA <- NA
mydata$slopeQIDSB <- NA
mydata$slopeQIDSC <- NA

allGraphs <- list()
set.seed(2021)

for (ii in 1:length(pids)){
  
  tid <- pids[ii]
  #find when each phase ends for this individual
  endA <- mydata[mydata$id == tid,]$QIDSstwhendone_14  
  endB <- mydata[mydata$id == tid,]$QIDSstwhendone_28  
  endC <- mydata[mydata$id == tid,]$QIDSstwhendone_42  
  startB <- endA + 1
  startC <- endB + 1
  
  #Prepare daily score sequences (making sure NAs are there for missing days)

  dailyqids <- data.frame(c(seq(0,as.numeric(endA)),seq((as.numeric(endA) + 1),as.numeric(endB)),seq((as.numeric(endB) + 1),as.numeric(endC))))
  colnames(dailyqids) <- "day"
  
  allqids <- data.frame(cbind(qidsdata[qidsdata$userid == tid,]$whendone,qidsdata[qidsdata$userid == tid,]$total))
  dailyqids <- merge(dailyqids,allqids,by.x = "day",by.y = "X1",all.x = TRUE)  
  colnames(dailyqids) <- c("day","score")
  dailyqids$score <- as.numeric(dailyqids$score)
  dailyqids$day <- as.numeric(dailyqids$day)
  
  #calculate autocorrelations and add to main dataframe
  
  mydata[mydata$id == tid,]$acfA <- (acf(dailyqids[dailyqids$day < endA+ 1,]$score,plot=FALSE,na.action = na.exclude))$acf[2]
  mydata[mydata$id == tid,]$acfB <- (acf(dailyqids[dailyqids$day > endA & dailyqids$day < endB +1,]$score,plot=FALSE,na.action = na.exclude))$acf[2]
  mydata[mydata$id == tid,]$acfC <- (acf(dailyqids[dailyqids$day > endB & dailyqids$day < endC +1,]$score,plot=FALSE,na.action = na.exclude))$acf[2]

  #sequential analyses
  
  qidsA <- dailyqids[dailyqids$day < endA+ 1,]$score
  qidsB <- dailyqids[dailyqids$day > endA & dailyqids$day < endB +1,]$score
  qidsC <- dailyqids[dailyqids$day > endB & dailyqids$day < endC +1,]$score
  
  #add to list
  alldailyqids[[pids[ii]]]$qidsA <- qidsA
  alldailyqids[[pids[ii]]]$qidsB <- qidsB
  alldailyqids[[pids[ii]]]$qidsC <- qidsC
  
  BFsSlope<-vector()
  BFsMean<-vector()
  for (i in 2:length(qidsB)){
    
      output<-trendtest.Gibbs.AR(qidsA,qidsB[1:i],iterations=5000,leftSidedInt = T,leftSidedSlp = T,r.scaleInt = 1,r.scaleSlp = 1,return.chains = T,return.onesided = T,progress=F)
      BFsSlope<-c(BFsSlope,exp(output$logbfOnesided)[[2]])
      output2<-BayesSingleSub::ttest.Gibbs.AR(qidsA,qidsB[1:i],iterations=5000,return.chains=TRUE,leftSided = TRUE,r.scale=1,return.onesided = TRUE)
      BFsMean<-c(BFsMean,exp(output2$logbfOnesided))
    
  }
  
  #find most extreme BF for mean and slope
  maxBFmAB <- max(BFsMean)
  minBFmAB <- min(BFsMean)

  maxBFsAB <- max(BFsSlope[!is.na(BFsSlope)])
  minBFsAB <- min(BFsSlope[!is.na(BFsSlope)])

  mydata[mydata$id == tid,]$maxBFmAB <- maxBFmAB
  mydata[mydata$id == tid,]$maxBFsAB <- maxBFsAB
  mydata[mydata$id == tid,]$minBFmAB <- minBFmAB
  mydata[mydata$id == tid,]$minBFsAB <- minBFsAB
  mydata[mydata$id == tid,]$meanQIDSA <- mean(qidsA,na.rm = TRUE)
  days<-seq(0,length(qidsA)-1)
  mydata[mydata$id == tid,]$slopeQIDSA <- (lm(qidsA~days))$coefficients[2][[1]]
  
  mydata[mydata$id == tid,]$meanQIDSB <- mean(qidsB,na.rm = TRUE)
  days<-seq(0,length(qidsB)-1)
  mydata[mydata$id == tid,]$slopeQIDSB <- (lm(qidsB~days))$coefficients[2][[1]]
  
  BFsSlope2<-vector()
  BFsMean2<-vector()
  for (i in 2:length(qidsC)){
      output<-trendtest.Gibbs.AR(qidsB,qidsC[1:i],leftSidedInt = T,leftSidedSlp = T,r.scaleInt = 1,r.scaleSlp = 1,return.chains = T,return.onesided = T,progress=F)
      BFsSlope2<-c(BFsSlope2,exp(output$logbfOnesided)[[2]])
      output2<-BayesSingleSub::ttest.Gibbs.AR(qidsB,qidsC[1:i],iterations=5000,return.chains=TRUE,leftSided = T,r.scale=1,return.onesided = TRUE)
      BFsMean2<-c(BFsMean2,exp(output2$logbfOnesided))
    
  }
  #find most extreme BF for mean and slope
  maxBFmBC <- max(BFsMean2)
  minBFmBC <- min(BFsMean2)

  
  maxBFsBC <- max(BFsSlope2[!is.na(BFsSlope2)])
  minBFsBC <- min(BFsSlope2[!is.na(BFsSlope2)])

  mydata[mydata$id == tid,]$maxBFmBC <- maxBFmBC
  mydata[mydata$id == tid,]$maxBFsBC <- maxBFsBC
  mydata[mydata$id == tid,]$minBFmBC <- minBFmBC
  mydata[mydata$id == tid,]$minBFsBC <- minBFsBC
  
  mydata[mydata$id == tid,]$meanQIDSC <- mean(qidsC,na.rm = TRUE)
  days<-seq(0,length(qidsC)-1)
  mydata[mydata$id == tid,]$slopeQIDSC <- (lm(qidsC~days))$coefficients[2][[1]]
  
  

  #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  colourkey<-c("#000000","#0072B2","#009E73")
  allQIDS<-c(qidsA,qidsB,qidsC)
  dayssincestart<-seq(0,length(allQIDS)-1)
  phase<-c(rep("A",length(qidsA)),rep("B",length(qidsB)),rep("C",length(qidsC)))
  allQIDSd<-data.frame(cbind(allQIDS,dayssincestart,phase))
  colnames(allQIDSd) <- c("QIDS",colnames(allQIDSd)[2:3])
  allQIDSd$QIDS<-as.integer(as.character(allQIDSd$QIDS))
  allQIDSd$dayssincestart<-as.numeric(as.character(allQIDSd$dayssincestart))
  (qidsplot<-ggplot(allQIDSd,aes(x=dayssincestart,y=QIDS,colour=phase))+geom_point(size=3)
    +geom_line(size=0.5)
    +scale_colour_manual(values=colourkey)
    +geom_smooth(method="lm",size=1,se=FALSE)
    +theme(axis.text.y = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.text.x = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.title.y = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.title.x=element_blank())
    +theme(legend.title = element_blank())+#element_text(colour="black", size=rel(1.75), face="bold"))+
      theme(legend.position = "none")
    +theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "lightgrey"), 
      axis.line = element_line(size = 0.5, colour = "black", linetype=1)
      
    )
    )
  
  allQIDSd<-cbind(allQIDSd,c(rep(NA,length(qidsA)),c(NA,BFsMean),c(NA,BFsMean2)))
  allQIDSd<-cbind(allQIDSd,c(rep(NA,length(qidsA)),c(NA,BFsSlope),c(NA,BFsSlope2)))
  colnames(allQIDSd)<-c("QIDS",colnames(allQIDSd)[2:3],"BFsMean","BFsSlope")
  (BFmplot<-ggplot(allQIDSd,aes(x=dayssincestart,y=BFsMean,colour=phase))+geom_point(size=3)
    +geom_line(size=0.5)
    +scale_colour_manual(values=colourkey)
    +theme(axis.text.y = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.text.x = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.title.y = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.title.x=element_blank())
    +theme(legend.title = element_blank())+#element_text(colour="black", size=rel(1.75), face="bold"))+
      theme(legend.position = "none")
    +theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "lightgrey"), 
      axis.line = element_line(size = 0.5, colour = "black", linetype=1)
      
    )
    +labs(y = "BF Mean"))
  
  
  (BFslplot<-ggplot(allQIDSd,aes(x=dayssincestart,y=BFsSlope,colour=phase))+geom_point(size=3)
    +geom_line(size=0.5)
    +scale_colour_manual(values=colourkey)
    +theme(axis.text.y = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.text.x = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.title.y = element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(axis.title.x=element_text(size = rel(2.5),colour="black",face="bold"))
    +theme(legend.title = element_blank())+#element_text(colour="black", size=rel(1.75), face="bold"))+
      theme(legend.position = "none")
    +theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "lightgrey"), 
      axis.line = element_line(size = 0.5, colour = "black", linetype=1)
      
    )
    +labs(y = "BF Slope", x="Days since Start"))
  
  
  Ngraphs<-ggarrange(qidsplot,BFmplot,BFslplot,heights = c(1,1,1),ncol = 1, nrow = 3,align = "v")
  allGraphs[[ii]] <- Ngraphs
  alldailyBFs[[ii]]<-allQIDSd
  # uncomment to save graphs
  #emf(file=paste0("SCgraph",tid,".emf"),width=10,height=12)
  #print(Ngraphs)
  #dev.off()
  
}
# uncomment to display graphs:
# allGraphs




##### Simulate continuations #####

###### Functions ######

####Continuation sim 5####

#'
#' @param thispre The data points from the baseline phase (numerical vector)
#' @param nB The number of data points in before 'final' should be reached
#' @param istep interval between final qids score in simulations
#' @param extra The number of extra data points after nB
#' @param nsims The number of simulations
contSCsim5 <- function(thispre,nB,istep,extra=0,nsims=1000){
  allSims<-list()
  thissd<-sd(thispre,na.rm=TRUE)
  fit <-arima(thispre,order=c(1,0,0),xreg=1:length(thispre))
  startmean<-mean(thispre,na.rm=TRUE)
  #set start mean and scale sd to desired value
  preforreg<-data.frame(cbind(thispre,seq(1,length(thispre))))
  mod<-lm(thispre~V2,data=preforreg)
  preendpoint<-mod$coefficients[[1]]+mod$coefficients[[2]]*length(pre)
  nochangeqids<-mod$coefficients[[1]]+mod$coefficients[[2]]*(length(pre)+nB)
  allendqids<-seq(nochangeqids,0,by=-istep)
  slopes<-(allendqids-preendpoint) / (nB-1)
  #calculate pre standard deviation in terms of deviation from regression line
  ss<-mod$residuals*mod$residuals
  spre<-sqrt(sum(ss)/(length(mod$residuals)-1))
  for (k in 1:length(allendqids)){
    thisSim<-list()
    #set slope for intervention phase
    slope<-slopes[[k]]
    postendpoint<-preendpoint+slope*nB
    
    for (i in 1:nsims){
      
      
      ar.sim<-arima.sim(model=fit$model,n=nB+extra,sd=spre)
      post<-ar.sim+preendpoint
      for(j in 2:length(post)){
        post[j]<-post[j]+(slope)*(j-1)
      }
      #make data more realistic by restricting to whole numbers and making sure no values < 0
      
      post<-round(post)
      post[post < 0] <- 0
      thisSim[[i]] <- list()
      thisSim[[i]]$pre <- thispre
      thisSim[[i]]$post <- post
      thisSim[[i]]$final<-postendpoint
      thisSim[[i]]$nB<-nB
      thisSim[[i]]$extra<-extra
      thisSim[[i]]$allendqids<-allendqids
    } #end of for (i in 1:nsims){
    allSims[[k]]<-thisSim
  }#end of for (k in 1:length(allendqids)){
  return(allSims)
} #end of contSCsim5


checkSims2 <- function(sim,index=0,pnumber=999){
  thissim <- sim
  thissample <- sample(length(thissim),1)
  nB<-thissim[[1]]$nB
  cat("\n Sample ",thissample, ":\n\n")
  cat("Mean phase A = ",mean(thissim[[thissample]]$pre,na.rm=TRUE)," \n")
  cat("SD phase A = ",sd(thissim[[thissample]]$pre,na.rm=TRUE)," \n")
  cat("Mean phase B = ",mean(thissim[[thissample]]$post,na.rm=TRUE)," \n")
  cat("SD phase B = ",sd(thissim[[thissample]]$post,na.rm=TRUE)," \n")
  cat("See Plots for plot (gap shows break between phases)\n\n")

  #From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  colourkey<-c("#000000","#0072B2","#009E73")
  allQIDS<-c(thissim[[thissample]]$pre,thissim[[thissample]]$post)
  dayssincestart<-seq(0,length(allQIDS)-1)
  Phase<-c(rep("A",length(thissim[[thissample]]$pre)),rep("Initial Continuation",nB),rep("Further Continuation",length(thissim[[thissample]]$post)-nB))
  allQIDSd<-data.frame(cbind(allQIDS,dayssincestart,Phase))
  colnames(allQIDSd) <- c("QIDS",colnames(allQIDSd)[2:3])
  allQIDSd$QIDS<-as.integer(as.character(allQIDSd$QIDS))
  allQIDSd$dayssincestart<-as.numeric(as.character(allQIDSd$dayssincestart))
  (qidsplot<-ggplot(allQIDSd,aes(x=dayssincestart,y=QIDS,colour=Phase))+geom_point(size=3)
    +geom_line(size=0.5)
    +scale_colour_manual(values=colourkey)
    +geom_smooth(method="lm",size=1,se=FALSE)
    +theme(axis.text.y = element_text(size = rel(1.5),colour="black",face="bold"))
    +theme(axis.text.x = element_text(size = rel(1.5),colour="black",face="bold"))
    +theme(axis.title.y = element_text(size = rel(1.5),colour="black",face="bold"))
    +theme(axis.title.x=element_blank())
    +theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "lightgrey"), 
      axis.line = element_line(size = 0.5, colour = "black", linetype=1))
    +geom_vline(xintercept=nB+length(thissim[[thissample]]$pre),linetype="dashed")
    +labs(title=paste0("Sample continuation for p",ifelse(pnumber==999,"?",pnumber),", endqids = ",ifelse(index==0,"?",round(thissim[[1]]$allendqids[[index]],2)))))
  print(qidsplot)
  
  allmeans<-vector()
  
  for (j in 1:(length(thissim[[1]]$pre))){
    thismean<-vector()
    for (i in 1:length(thissim)){
      thismean<-c(thismean,thissim[[i]]$pre[j])
      
    }
    allmeans<-c(allmeans,mean(thismean,na.rm=TRUE))
  }
  for (j in 1:length(thissim[[1]]$post)){
    thismean<-vector()
    for (i in 1:length(thissim)){
      thismean<-c(thismean,thissim[[i]]$post[j])
      
    }
    allmeans<-c(allmeans,mean(thismean,na.rm=TRUE))
  }
  
  
  allQIDS<-allmeans
  dayssincestart<-seq(0,length(allQIDS)-1)
  Phase<-c(rep("A",length(thissim[[thissample]]$pre)),rep("Initial Continuation",nB),rep("Further Continuation",length(thissim[[thissample]]$post)-nB))
  allQIDSd<-data.frame(cbind(allQIDS,dayssincestart,Phase))
  colnames(allQIDSd) <- c("QIDS",colnames(allQIDSd)[2:3])
  allQIDSd$QIDS<-as.integer(as.character(allQIDSd$QIDS))
  allQIDSd$dayssincestart<-as.numeric(as.character(allQIDSd$dayssincestart))
  (qidsplot<-ggplot(allQIDSd,aes(x=dayssincestart,y=QIDS,colour=Phase))+geom_point(size=3)
    +geom_line(size=0.5)
    +scale_colour_manual(values=colourkey)
    +geom_smooth(method="lm",size=1,se=FALSE)
    +theme(axis.text.y = element_text(size = rel(1.5),colour="black",face="bold"))
    +theme(axis.text.x = element_text(size = rel(1.5),colour="black",face="bold"))
    +theme(axis.title.y = element_text(size = rel(1.5),colour="black",face="bold"))
    +theme(axis.title.x=element_blank())
    +theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "lightgrey"), 
      axis.line = element_line(size = 0.5, colour = "black", linetype=1))
    +geom_vline(xintercept=nB+length(thissim[[thissample]]$pre),linetype="dashed")
    +labs(title=paste0("Averaged continuation for p",ifelse(pnumber==999,"?",pnumber),", endqids = ",ifelse(index==0,"?",round(thissim[[1]]$allendqids[[index]],2)))))
  print(qidsplot)
  
}#end of checkSims2()

BFseqsim2 <- function(sim,rscaleS=1,rscaleM=1){
  sims<-sim
  simResults<-list()
  simResults$sim<-sim
  simResults$rscaleM<-rscaleM
  simResults$rscaleS<-rscaleS
  simResults$results<-list()
  nsims<-length(sims)
  
  minlength <- 3
  maxlength <- length(sims[[1]]$post)
  
  
  cat("\n Sim number: ")
  for (j in 1:nsims)
  {
    simResults$results[[j]]<-list()
    ns<-vector()
    BFSlopes<-vector()
    BFMeans <- vector()
    if (j%%100==0)
    {
      cat(j,",",sep="");
    }
    pre <-sims[[j]]$pre
    post <-sims[[j]]$post
    
    for (i in minlength:maxlength){
      
      capture.output(output<-trendtest.Gibbs.AR(pre,post[1:i],leftSidedInt = T,leftSidedSlp = T,r.scaleInt = rscaleS,r.scaleSlp = rscaleS,return.chains = T,return.onesided = T,progress=F))
      BFSlopes <- c(BFSlopes,(exp(output$logbfOnesided))[[2]])
      capture.output(output2<-BayesSingleSub::ttest.Gibbs.AR(pre,post[1:i],return.chains=TRUE,leftSided = T,r.scale=rscaleM,return.onesided = TRUE))
      BFMeans <- c(BFMeans,exp(output2$logbfOnesided))
      ns <- c(ns,i)
      
    }
    
    
    
    simResults$results[[j]]$BFMeans<-BFMeans
    simResults$results[[j]]$BFSlopes<-BFSlopes
    simResults$results[[j]]$ns<-ns
  }
  return(simResults)
}#end of BFseqsim2()


#' @param simresult A list of simulated sequential Bayesian analyses produced by BFseqsim2()
#' @param nmin Nmin
#' @param nmax Nmax
#' @param bffailM BF fail for the mean
#' @param bffailS BF fail for the slope
#' @param failrule can be either "either", "slope" or "mean"
#' @param bfsuccessM BF success for the mean
#' @param bfsuccessS BF succcess for the slope
#' @param successrule can be "either", "slope", or "mean"
#' @returns returns a dataframe
makesimTablev2<-function(simresult,nmin,nmax,bffailM,bffailS,failrule="either",bfsuccessM,bfsuccessS,successrule="either",rNames=c("NA")){
  
  failhead<-c(paste0(seq(nmin,nmax),"F"),paste0(seq(nmin,nmax),"S"))
  allrows<-vector()
  for (k in 1:length(simresult)){
    thisrow<-vector()
    sims<-simresult[[k]]$results
    
    
    
    fails<-rep(0,length(failhead)/2)
    successes<-rep(0,length(failhead)/2)
    nsims<-length(sims)
    
    for (j in 1:nsims){
      slopefailed<-FALSE
      meanfailed<-FALSE
      slopesuccess<-FALSE
      meansuccess<-FALSE
      for (i in nmin:nmax){
        thisi <- which(sims[[j]]$ns==i)
        #get BFs at this sample size
        thisslopeBF<-sims[[j]]$BFSlopes[thisi]
        thismeanBF <-sims[[j]]$BFMeans[thisi]
        
        if (thisslopeBF > bffailS){
          slopefailed<-TRUE
        }
        if (thismeanBF > bffailM){
          meanfailed<-TRUE
        }
        if (thisslopeBF < bfsuccessS){
          slopesuccess<-TRUE
        }
        if (thismeanBF < bfsuccessM){
          meansuccess<-TRUE
        }
        
        #check for failure
        if (failrule=="either"){
          if (slopefailed|meanfailed){
            fails[i+1-nmin]<-fails[i+1-nmin]+1
            break
          }
        }else if (failrule=="slope"){
          if (slopefailed){
            fails[i+1-nmin]<-fails[i+1-nmin]+1
            break
          }
        }else if (failrule=="mean"){
          if (meanfailed){
            fails[i+1-nmin]<-fails[i+1-nmin]+1
            break
          }
        }else if (failrule=="both"){
          if (meanfailed&slopefailed){
            fails[i+1-nmin]<-fails[i+1-nmin]+1
            break
          }
        }
        
        #check for success
        if (successrule=="either"){
          if (slopesuccess|meansuccess){
            successes[i+1-nmin]<-successes[i+1-nmin]+1
            break
          }
        }else if (successrule=="slope"){
          if (slopesuccess){
            successes[i+1-nmin]<-successes[i+1-nmin]+1
            break
          }
        }else if (successrule=="mean"){
          if (meansuccess){
            successes[i+1-nmin]<-successes[i+1-nmin]+1
            break
          }
        }else if (successrule=="both"){
          if (meansuccess&slopesuccess){
            successes[i+1-nmin]<-successes[i+1-nmin]+1
            break
          }
        }
        
        
      }
      
    }
    for (ii in 2:length(fails)) {
      fails[[ii]]<-fails[[ii]]+fails[[ii-1]]
    }
    for (ii in 2:length(successes)) {
      successes[[ii]]<-successes[[ii]]+successes[[ii-1]]
    }
    thisrow<-c(fails,successes)
    thisrow<-thisrow*100/nsims
    thisrow<-c(thisrow,100-(thisrow[length(fails)]+thisrow[length(successes)+length(fails)]))
    allrows<-rbind(allrows,thisrow)
  }
  
  allrows<-data.frame(allrows)
  colnames(allrows)<-c(failhead,"Nmax")
  if (length(rNames)>1){
    rownames(allrows)<-rNames
  }
  return(allrows)
  
}#end of makesimTablev2()

####all participants####

contSims<-list()
contResults<-list()
set.seed(2021)

for (ii in 1:length(pids)) {
  cat("participant: ",ii,"\n\n")
  thisqidsA<-alldailyqids[[ii]]$qidsA
  contSims[[ii]]<-contSCsim5(thispre=thisqidsA,nB=14,istep=2,extra=14,nsims=1000)
  contResults[[ii]]<-list()
  for (i in 1:length(contSims[[ii]])){
    cat("simulation: ",i,"\n\n")
    contResults[[ii]][[i]]<-BFseqsim2(contSims[[ii]][[i]],rscaleS = 1, rscaleM = 1)
  }
}

# save if wanted 
# save(contResults,file="contresultsall_13.10.26.RData")


#check how simulated data looks
for (j in 1:length(pids)){
  for (i in 1:length(contSims[[j]])){
    checkSims2(contSims[[j]][[i]],index=i,pnumber=j)
  }
}


# Show results for range of parameters. Here I just show the final set of parameters found in the paper.
# Prior to this a process of trial and error was used testing out different Nmins, Nmaxs, BF boundaries etc
# Until these specific parameters were arrived at:

makesimTablev2(contResults[[1]],6,21,3,12,"either",1/3,1/3,"either",round(contResults[[1]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[2]],10,28,3,8,"either",1/10,1/10,"either",round(contResults[[2]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[3]],4,14,4,9,"either",1/3,1/3,"either",round(contResults[[3]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[3]],15,28,8,12,"either",1/3,1/3,"either",round(contResults[[3]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[4]],4,28,8,12,"either",1/9,1/9,"either",round(contResults[[4]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[5]],6,28,8,16,"either",1/3,1/3,"either",round(contResults[[5]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[5]],4,28,8,12,"either",1/3,1/3,"either",round(contResults[[5]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[6]],9,14,3,3,"either",1/8,1/8,"either",round(contResults[[6]][[1]]$sim[[1]]$allendqids,2))
makesimTablev2(contResults[[7]],3,14,3,4,"slope",1/3,1/3,"either",round(contResults[[7]][[1]]$sim[[1]]$allendqids,2))


#####Get data for table 4#####

for (i in 1:length(pids)){
  cat("\nParticipant: ",i,"\n\n")
  thispre<-contResults[[i]][[1]]$sim[[1]]$pre
  fit <-arima(thispre,order=c(1,0,0),xreg=1:length(thispre))
  cat("AR1: ",fit$coef[[1]],"\n")
  preforreg<-data.frame(cbind(thispre,seq(1,length(thispre))))
  mod<-lm(thispre~V2,data=preforreg)
  ss<-mod$residuals*mod$residuals
  spre<-sqrt(sum(ss)/(length(mod$residuals)-1))
  cat("SDslope: ",spre,"\n")
  cat("Slope: ",mod$coefficients[[2]],"\n")
  preendpoint<-mod$coefficients[[1]]+mod$coefficients[[2]]*length(pre)
  nochangeqids<-mod$coefficients[[1]]+mod$coefficients[[2]]*(length(pre)+contResults[[i]][[1]]$sim[[1]]$nB)
  cat("Pre endpoint: ",round(preendpoint,2),"\n")
  cat("H0 end QIDS (change): ",round(nochangeqids,2)," (",round((preendpoint-nochangeqids),2),")","\n")
  slopes<-(contResults[[i]][[1]]$sim[[1]]$allendqids-preendpoint) / (contResults[[i]][[1]]$sim[[1]]$nB-1)
  #calculate pre standard deviation in terms of deviation from regression line
  cat("End QIDS: \n")
  cat("EndQIDS\tDecrease\tSlope\n")
  for (j in 1:length(slopes)){
    cat(round(contResults[[i]][[1]]$sim[[1]]$allendqids[j],2),"\t",round(preendpoint-contResults[[i]][[1]]$sim[[1]]$allendqids[j],2),"\t",round(slopes[j],2),"\n")
  }
  cat("BFs:\n")
  print(alldailyBFs[i])
}





#### Compare weekly and daily QIDS ####
rows <- c(" "," ","QIDS","","Daily QIDS"," ","Most extreme BF"," ")
rows <- rbind(rows,c(" "," ","Score","Reliable change?","Mean","Slope","Mean","Slope"))
rows <-rbind(rows,c("Participant number"," "," "," "," "," "," "," "))
pids <- unique(mydata$id)
pids<-pids[order(as.numeric(pids))]
for (i in 1:length(pids)){
  thisdata <- mydata[mydata$id == pids[i],]
  label1 <- "Post-CBM"
  label2 <- "Post-CCT"
  if (thisdata$arm == "CCTfirst"){
    label1 <- label2
    label2 <- "Post-CBM"
  }
  thisrow <- c(as.character(i)," "," "," "," "," "," "," ")
  rows <- rbind(rows,thisrow)
  thisrow <- c(" ","Pre-Baseline",thisdata$QIDSsttotal_0," "," "," "," "," ")
  rows <-rbind(rows,thisrow)
  thisrow <- c(" ","Post-Baseline",thisdata$QIDSsttotal_14,thisdata$RC_0A,round(thisdata$meanQIDSA,2),round(thisdata$slopeQIDSA,2)," "," ")
  rows <-rbind(rows,thisrow)
  tmaxmBF<-ifelse(thisdata$maxBFmAB>(1/thisdata$minBFmAB),thisdata$maxBFmAB,thisdata$minBFmAB)
  tmaxsBF<-ifelse(thisdata$maxBFsAB>(1/thisdata$minBFsAB),thisdata$maxBFsAB,thisdata$minBFsAB)
  thisrow <- c(" ",label1,thisdata$QIDSsttotal_28,thisdata$RC_AB,round(thisdata$meanQIDSB,2),round(thisdata$slopeQIDSB,2),round(tmaxmBF,2),round(tmaxsBF,2))
  rows <-rbind(rows,thisrow)
  tmaxmBF<-ifelse(thisdata$maxBFmBC>(1/thisdata$minBFmBC),thisdata$maxBFmBC,thisdata$minBFmBC)
  tmaxsBF<-ifelse(thisdata$maxBFsBC>(1/thisdata$minBFsBC),thisdata$maxBFsBC,thisdata$minBFsBC)
  thisrow <- c(" ",label2,thisdata$QIDSsttotal_42,thisdata$RC_BC,round(thisdata$meanQIDSC,2),round(thisdata$slopeQIDSC,2),round(tmaxmBF,2),round(tmaxsBF,2))
  rows <-rbind(rows,thisrow)
  
}
table3<-as.data.frame(rows)
#output table as excel file if wanted (to insert into paper)
#write.xlsx(table3,"Table3_13.10.26.xlsx",row.names = FALSE)


##### All outcome measures supplement table####

measures<-c("QIDSst","GAD7","DARS","PMH","PANASPweek")
tpids<-c(mydata[mydata$arm=="CBMfirst",]$id[order(as.numeric(mydata[mydata$arm=="CBMfirst",]$id))],mydata[mydata$arm=="CCTfirst",]$id[order(as.numeric(mydata[mydata$arm=="CCTfirst",]$id))])
ends<-c("0","14","28","42")
rows<-c("","",tpids)
for (i in 1:length(measures)){
  rows<-rbind(rows,c(measures[i],"",rep("",length(tpids))))
  for (j in 1:4){
    tvar<-paste0(measures[i],"total_",ends[j])
    thisrow<-c("",as.character(j))
    for (k in 1:length(tpids)){
      thisrow<-c(thisrow,mydata[mydata$id==tpids[k],tvar])
      
    }
    rows<-rbind(rows,thisrow)
  }
}
rows<-data.frame(rows)
View(rows)

#output as table if needed
#write.xlsx(rows,"TableS1_13.10.26.xlsx",col.names = FALSE,row.names = FALSE)

#RCIS
#QIDS = 6.49
#GAD7 = 5.84
#PMH = use patient sample from https://link.springer.com/article/10.1186/s40359-016-0111-x
# SD = 0.66, alpha = 0.91, RCI = 0.548828 x 9 = 4.94
#PANAS-18, using baseline data from Westermann et al., alpha = 0.926, SD = 12.777, RCI = 9.63
#DARS from Rizvi et al., study 3 (depressed sample), alpha = 0.96, SD = 13.95 , RCI = 7.73
# pooling SDS from Rizvi TRD and non-TRD samples to get SD:
# https://www.statology.org/pooled-standard-deviation-in-r/
#pooled <- sqrt(((37-1)*15.5^2 + (15-1)*8.8^2) / (37+15-2))

