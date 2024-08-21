lapply(c("tidyverse","ggthemes","ggpubr","lme4","nlme","piecewiseSEM"),require,character.only=T) #load packages
setwd() #insert directory where data files are saved

#Figure 2 Irrigation, Precipitation vs Year
t2L=read.csv("PRECIP_IRRIGATION_AF_YEAR_MORE2.csv")#Long format
t2L$Time="Full Year (Dec-Nov)";
t2L$Time[str_detect(t2L$Metric,"1")]="Mosquito season (Apr-Nov)"
t2L$IrrPre="Irrig";
t2L$IrrPre[str_detect(t2L$Metric,"PRECIP")]="Precip"
t2Lm=aggregate(Value~Metric,data=t2L,FUN=mean)
t2Lm$Time=rep(unique(t2L$Time),2)
t2Lm$IrrPre=rep(unique(t2L$IrrPre),each=2)
lb1=data.frame(YEAR=1997.9,Value=c(2.9,1.2),Let=c("A","B"),IrrPre=NA,Time=unique(t2L$Time))
Fig2=ggplot(t2L,aes(x=YEAR,y=Value,col=IrrPre))+theme_few()+
  geom_line()+geom_text(data=lb1,aes(label=Let),col="black")+
  geom_pointrange(aes(ymin=Value-SE,ymax=Value+SE),size=.5)+
  geom_hline(data=t2Lm,aes(yintercept=Value,col=IrrPre),linetype="dashed")+
  facet_wrap(~Time,nrow=2,scales="free_y")+
  scale_color_manual(values=c("red","blue"),name="",labels=c("Irrigation","Precipitation"))+
  labs(x="Year", y = "Irrigation or Precipitation (Acre-feet)")+
  theme(text=element_text(size=20),legend.position="right")+
  scale_x_continuous(expand=c(0.01,0.15),breaks=1998:2010);Fig2
ggsave("Figure 2.pdf",Fig2,height=9,width=11)

t1=read.csv("PRECIP_IRRIGATION_AF_YEAR_MORE.csv")#Wide format
mean(t1$AW_ACRE1)/(mean(t1$PRECIP_ACRE1)+mean(t1$AW_ACRE1))*100
#31% 
mean(t1$AW_ACRE1)/mean(t1$PRECIP_ACRE1)*100
sd(t1$AW_ACRE)/sqrt(length(t1$AW_ACRE))

#Figure 3 TAR abundance vs Year and Irrigation
t3z=read.csv("LAT_LONG_YEAR_TAR_COUNTS_WITHLANDCOVER.csv")#load in data
t3za=aggregate(cbind(TAR_COUNT,PIP_QNX_COUNT,ERY_COUNT,
                     DAU_CENTER_LAT,DAU_CENTER_LONG)~DAU_YEAR+DAU_CODE,data=t3z,FUN=mean) #ave by DAU-YR
t3za2=aggregate(cbind(TAR_COUNT,PIP_QNX_COUNT,ERY_COUNT,
                      DAU_CENTER_LAT,DAU_CENTER_LONG)~DAU_CODE,data=t3za,FUN=mean) #ave by DAU
t3a2=aggregate(cbind(APPLIED_WATER_ACRE,PRECIP_AF_ACRE,MEAN_TEMPC,DEV_TOTAL_PERCENT_DAU,
                     OPENWATER_PERCENT_DAU,WETLAND_PERCENT_DAU,RICE_PERCENT_DAU,AG2010,
                     DAU_CENTER_LAT,DAU_CENTER_LONG)~DAU_CODE,data=t3z,FUN=mean) #Ave by DAU
t3a=merge(t3za2,t3a2) #merge TAR & other columns
t3a$APPLIED_WATER_ACREp.001=t3a$APPLIED_WATER_ACRE+.001


t3e = t3z[t3z$DAU_CODE%in%c(89,100,255),] %>% group_by(DAU_CODE,YEAR) %>% 
  summarize(SE=sd(TAR_COUNT)/(n()^.5),TAR_COUNT=mean(TAR_COUNT) )

f3a=ggplot(t3e[t3e$YEAR>2002,],aes(x=YEAR,y=TAR_COUNT,color=as.factor(DAU_CODE)))+theme_few()+
  geom_pointrange(aes(ymin=TAR_COUNT-SE,ymax=TAR_COUNT+SE),size=.5)+
  geom_line()+labs(y=expression(italic("C. tarsalis")~'abundance'),x=NULL)+
  scale_y_continuous(trans="log10",expand = c(0,0),limits=c(.2,500),labels = scales::comma,
                     breaks = c(.2,.5,1,2,5,10,20,50,100,200,500))+
  scale_x_continuous(expand=c(.02,.02),breaks = 2003:2010)+
  scale_color_manual(name = "Example DAUs:",values=c("lightgreen","darkgreen","black"),
                     labels=c("Low irrig.","Med. irrig.","High irrig."))+
  theme(legend.position="top")+
  theme(plot.margin=unit(c(25,5.5,5.5,5,5),"pt"))+#unit(c(t,r,b,l),"pt"), default: all 5.5pt
  theme(axis.text.x=element_blank())+
  annotate(x=2003,y=400,geom="text",label="A",size=5)+
  theme(text=element_text(size=20));f3a

#Fitted model in Fig 3b
Fs3=gls(log10(TAR_COUNT)~log10(APPLIED_WATER_ACREp.001),
        correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
        data=t3a, method = "ML");summary(Fs3)
rsquared(Fs3)
c(10^coef(Fs3)[2],10^(3*coef(Fs3)[2]))

nd3=data.frame(APPLIED_WATER_ACREp.001=c(.001*2^(0:10),3.2))
nd3$pTAR=predict(Fs3,newdata=nd3)
nd3$SE=predict(Fs3,newdata=nd3,se.fit=T)$se.fit

#Figure 3b
f3b=ggplot(t3a, aes(y=TAR_COUNT, x=APPLIED_WATER_ACREp.001))+theme_few()+ 
  geom_point(data=t3a[t3a$DAU_CODE%in%c(89,100,255),],
             col=c("darkgreen","black","lightgreen"),size=6)+
  geom_point(aes(col=MEAN_TEMPC),size=3)+
  scale_colour_gradient2(low="blue",mid="lightgray",high="red",midpoint=18)+
  theme(legend.position="top")+
  geom_line(data=nd3,aes(x=APPLIED_WATER_ACREp.001,y=10^(pTAR)))+
  geom_ribbon(data=nd3,aes(x=APPLIED_WATER_ACREp.001,y=10^(pTAR),
                           ymin=10^(pTAR-SE),ymax=10^(pTAR+SE)),alpha=.1)+
  labs(x=NULL,y=NULL,color="Temp, C")+
  scale_x_continuous(trans="log10",expand = c(0.02,0.0),limits=c(.001,3.2),
                     breaks = c(.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2))+
  scale_y_continuous(trans="log10",expand = c(0,0),limits=c(.02,650),labels = scales::comma,
                     breaks = c(.02,.05,.1,.2,.5,1,2,5,10,20,50,100,200,500))+
  annotate(x=.001,y=450,geom="text",label="B",size=5)+
  theme(axis.text.x=element_blank())+
  theme(text=element_text(size=20));f3b

#Univariate models
Ts2a=gls(log10(TAR_COUNT)~log10(APPLIED_WATER_ACREp.001),
        correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
        data=t3a, method = "ML");summary(Ts2a)
Ts2b=gls(log10(TAR_COUNT)~log10(PRECIP_AF_ACRE),
        correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
        data=t3a, method = "ML");summary(Ts2b)
Ts2c=gls(log10(TAR_COUNT)~MEAN_TEMPC,
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts2c)
Ts2d=gls(log10(TAR_COUNT)~log10(WETLAND_PERCENT_DAU),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts2d)
Ts2e=gls(log10(TAR_COUNT)~log10(DEV_TOTAL_PERCENT_DAU),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts2e)
Ts2f=gls(log10(TAR_COUNT)~log10(OPENWATER_PERCENT_DAU),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts2f)

#extract info from TAR models above
Ts2=data.frame(Predictor=NA,Intercept=NA,Slope=NA,SE=NA,P_value=NA,R2=NA)
for (i in 1:6) {
  modeli=paste0("Ts2",letters[i])
  Ts2[i,]=c(rownames(summary(get(paste0("Ts2",letters[i])))$tTable)[2],
            coef(get(modeli))["(Intercept)"],
            coef(get(modeli))[2],
            summary(get(modeli))$tTable[2,"Std.Error"],
            summary(get(modeli))$tTable[2,"p-value"],
            rsquared(get(modeli))$R.squared)  }
write.csv(Ts2,"Table S2 TAR.csv")

#Table S3
Ts3=gls(log10(TAR_COUNT)~log10(APPLIED_WATER_ACRE+0.001)+log10(PRECIP_AF_ACRE)+
          MEAN_TEMPC+log10(WETLAND_PERCENT_DAU)+
          log10(DEV_TOTAL_PERCENT_DAU)+
          log10(OPENWATER_PERCENT_DAU),
        correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
        data=t3a, method = "ML");summary(Ts3)
rsquared(Ts3)
shapiro.test(residuals(Ts3))

#Examine quadratic for Temp
t3a$MEAN_TEMPC2=t3a$MEAN_TEMPC^2
Ts3q=gls(log10(TAR_COUNT)~log10(APPLIED_WATER_ACRE+0.001)+log10(PRECIP_AF_ACRE)+
           (MEAN_TEMPC)+MEAN_TEMPC2+log10(WETLAND_PERCENT_DAU)+
           log10(DEV_TOTAL_PERCENT_DAU)+
           log10(OPENWATER_PERCENT_DAU),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), data=t3a, method = "ML");summary(Ts3q)
AIC(Ts3,Ts3q) #model w/out Temp2 better

#Univariate PIP_QNX_COUNT
Ts4a=gls(log10(PIP_QNX_COUNT+.1)~log10(APPLIED_WATER_ACREp.001),
        correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
        data=t3a, method = "ML");summary(Ts4a)
rsquared(Ts4a)
c(10^coef(Ts4a)[2],10^(3*coef(Ts4a)[2]))

Ts4b=gls(log10(PIP_QNX_COUNT+.1)~log10(PRECIP_AF_ACRE),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts4b)

Ts4c=gls(log10(PIP_QNX_COUNT+.1)~MEAN_TEMPC,
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts4c)

Ts4d=gls(log10(PIP_QNX_COUNT+.1)~log10(WETLAND_PERCENT_DAU),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts4d)

Ts4e=gls(log10(PIP_QNX_COUNT+.1)~log10(DEV_TOTAL_PERCENT_DAU),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts4e)

Ts4f=gls(log10(PIP_QNX_COUNT+.1)~log10(OPENWATER_PERCENT_DAU),
         correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
         data=t3a, method = "ML");summary(Ts4f)

#extract info from PIP/QNX models above
Ts4=data.frame(Predictor=NA,Intercept=NA,Slope=NA,SE=NA,P_value=NA,R2=NA)
for (i in 1:6) {
  modeli=paste0("Ts4",letters[i])
  Ts4[i,]=c(rownames(summary(get(paste0("Ts4",letters[i])))$tTable)[2],
            coef(get(modeli))["(Intercept)"],
            coef(get(modeli))[2],
            summary(get(modeli))$tTable[2,"Std.Error"],
            summary(get(modeli))$tTable[2,"p-value"],
            rsquared(get(modeli))$R.squared)  }
write.csv(Ts4,"Table S4 PIP_QNX.csv")

t3ex = t3z[t3z$DAU_CODE%in%c(92,100,186),] %>% group_by(DAU_CODE,YEAR) %>% 
  summarize(SE=sd(PIP_QNX_COUNT)/(n()^.5),PIP_QNX_COUNT=mean(PIP_QNX_COUNT) )

f3c=ggplot(t3ex[t3ex$YEAR>2002,],aes(x=YEAR,y=PIP_QNX_COUNT,color=as.factor(DAU_CODE)))+theme_few()+
  geom_pointrange(aes(ymin=PIP_QNX_COUNT-SE,ymax=PIP_QNX_COUNT+SE),size=.5)+
  geom_line()+labs(y=expression(italic("C. pipiens")~'complex abundance'),x="Year")+
  scale_y_continuous(trans="log10",expand = c(0,0),limits=c(.2,500),labels = scales::comma,
                     breaks = c(.2,.5,1,2,5,10,20,50,100,200,500))+
  scale_x_continuous(expand=c(.02,.02),breaks = 2003:2010)+
  scale_color_manual(name = "Example DAUs:",values=c("lightgreen","darkgreen","black"),
                     labels=c("Low irrigation","Medium irrigation","High irrigation"))+
  theme(legend.position="none")+
  annotate(x=2003,y=400,geom="text",label="C",size=5)+
  theme(text=element_text(size=20));f3c

nd3x=data.frame(APPLIED_WATER_ACREp.001=c(.001*2^(0:10),3.2))
nd3x$pPIP=predict(Ts4a,newdata=nd3x)
nd3x$SE=predict(Ts4a,newdata=nd3x,se.fit=T)$se.fit

t3a[,c("DAU_CODE","APPLIED_WATER_ACREp.001","PIP_QNX_COUNT")]
f3d=ggplot(t3a, aes(y=PIP_QNX_COUNT, x=APPLIED_WATER_ACREp.001))+theme_few()+ 
  geom_point(data=t3a[t3a$DAU_CODE%in%c(92,100,186),],
             col=c("darkgreen","black","lightgreen"),size=6)+
  geom_point(aes(col=MEAN_TEMPC),size=3)+
  scale_colour_gradient2(low="blue",mid="lightgray",high="red",midpoint=18)+
  theme(legend.position="none")+
  geom_line(data=nd3x,aes(x=APPLIED_WATER_ACREp.001,y=10^(pPIP)))+
  geom_ribbon(data=nd3x,aes(x=APPLIED_WATER_ACREp.001,y=10^(pPIP),
                            ymin=10^(pPIP-SE),ymax=10^(pPIP+SE)),alpha=.1)+
  labs(x="Irrigation, Acre-feet",y=NULL,color="Temp, C")+
  scale_x_continuous(trans="log10",expand = c(0.02,0.0),limits=c(.001,3.3),
                     breaks = c(.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2))+
  scale_y_continuous(trans="log10",expand = c(0,0),limits=c(.02,650),labels = scales::comma,
                     breaks = c(.02,.05,.1,.2,.5,1,2,5,10,20,50,100,200,500))+
  theme(text=element_text(size=20))+
  annotate(x=.001,y=450,geom="text",label="D",size=5);f3d

Fig3=ggarrange(f3a,f3b,f3c,f3d,ncol=2,nrow=2);Fig3#Figure 3--------------
ggsave("Figure 3.pdf",Fig3,height=12,width=16)

#Table S5
Ts5=gls(log10(PIP_QNX_COUNT+.1)~log10(APPLIED_WATER_ACRE+0.001)+log10(PRECIP_AF_ACRE)+
          (MEAN_TEMPC)+log10(WETLAND_PERCENT_DAU)+
          log10(DEV_TOTAL_PERCENT_DAU)+
          log10(OPENWATER_PERCENT_DAU),
        correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
        data=t3a, method = "ML");summary(Ts5)
rsquared(Ts5)

#Table S6
Ts7=gls(log10(ERY_COUNT+.1)~log10(APPLIED_WATER_ACRE+0.001)+log10(PRECIP_AF_ACRE)+
          (MEAN_TEMPC)+log10(WETLAND_PERCENT_DAU)+
          log10(DEV_TOTAL_PERCENT_DAU)+
          log10(OPENWATER_PERCENT_DAU),
        correlation = corExp(1, form = ~ DAU_CENTER_LAT+DAU_CENTER_LONG, nugget=F), 
        data=t3a, method = "ML");summary(Ts7)
rsquared(Ts7)


#Figure 4 TAR CV
t4B=read.csv("temp_mosq_try.csv",head=TRUE,sep=",")#load in human index seasonal data 
#Subset to two sites: 1 each DAU_Year 256_2010, 257_2010
t4_2=t4B[t4B$LAT_LONG_YEAR%in%c("35.5775_-119.3475_2010","35.4325_-118.887222_2010"),]
t4_2$HighLow="High";t4_2$HighLow[t4_2$LAT_LONG_YEAR=="35.5775_-119.3475_2010"]="Low"
f4a=ggplot(t4_2,aes(x=WEEK,y=TAR_COUNT,col=HighLow))+theme_few()+
  geom_line()+geom_point(size=3)+lims(y=c(0,30))+
  labs(x="Week",y=expression(paste(italic("C. tarsalis "),"abundance")))+
  scale_color_manual(values=c("red","blue"),name="Example DAUs:",
                     labels=c("Low irrigation","High irrigation"))+
  scale_x_continuous(expand = c(0.02,0.0),limits=c(14,46))+
  annotate(x=14,y=29,geom="text",label="A",size=5)+
  theme(text=element_text(size=20),legend.position="top");f4a

s1=read.csv("CO2_MOSQ_RAW_TRAPS_USED.csv",head=TRUE,sep=",");#head(s)
DLLT=read.csv("DAU_LAT_LONG_TRAPS.csv")
s1=merge(s1,DLLT)
s1$PIP_QNX_COUNT=s1$PIP_COUNT+s1$QNX_COUNT
s2=aggregate(cbind(TAR_COUNT,PIP_QNX_COUNT,ERY_COUNT)~LAT_LONG_YEAR+DAU_CODE+YEAR,
             data=s1,FUN=mean) #ave by DAU-YR
s2a=s1 %>% group_by(LAT_LONG_YEAR,DAU_CODE) %>% summarize(TAR_SD=sd(TAR_COUNT),
                                                          PIP_QNX_SD=sd(PIP_QNX_COUNT),
                                                          ERY_SD=sd(ERY_COUNT))
s2=merge(s2,s2a)
s2$TAR_CV=100*s2$TAR_SD/s2$TAR_COUNT
s2$PIP_QNX_CV=100*s2$PIP_QNX_SD/s2$PIP_QNX_COUNT
s2$ERY_CV=100*s2$ERY_SD/s2$ERY_COUNT
s2=s2[s2$TAR_COUNT!=0,]

s3a=aggregate(cbind(TAR_COUNT,PIP_QNX_COUNT,ERY_COUNT,TAR_CV)~DAU_CODE+YEAR,data=s2,FUN=mean)#avg across DAU-YEAR
s3=aggregate(cbind(TAR_COUNT,PIP_QNX_COUNT,ERY_COUNT,TAR_CV)~DAU_CODE,data=s3a,FUN=mean)#avg across DAU

s3b=read.csv("DAU_ACROSS_YEARS4_ALL.csv")
s3=merge(s3,s3b)
s3$IRRIGATION_AFp.001=s3$IRRIGATION_AF+.001

Fm4=gls(log10(TAR_CV)~log10(IRRIGATION_AFp.001),
        correlation = corExp(1, form = ~ LAT_CENTER + LONG_CENTER, nugget=F), 
        data=s3, method = "ML");summary(Fm4)
shapiro.test(residuals(Fm4))
rsquared(Fm4)
c(10^coef(Fm4)[2],10^(3*coef(Fm4)[2]))


nd4=data.frame(IRRIGATION_AFp.001=c(.001*2^(0:10),2.76))
nd4$pCV=predict(Fm4,newdata=nd4)
nd4$SE=predict(Fm4,newdata=nd4,se.fit=T)$se.fit

f4b=ggplot(s3,aes(x=IRRIGATION_AFp.001,y=TAR_CV))+theme_few()+
  geom_point(size=2.5)+
  geom_point(data=s3[s3$DAU_CODE%in%c(256,257),],col=c("blue","red"),size=4)+
  geom_line(data=nd4,aes(x=IRRIGATION_AFp.001,y=10^(pCV)))+
  geom_ribbon(data=nd4,aes(x=IRRIGATION_AFp.001,y=10^(pCV),
                           ymin=10^(pCV-SE),ymax=10^(pCV+SE)),alpha=.1)+
  labs(x="Irrigation (Acre-feet)",
    y=expression(paste("CV ", italic("C. tarsalis "),"abundance")))+
  scale_x_continuous(trans="log10",expand = c(0.02,0.0),limits=c(.001,2.76),
                     breaks = c(.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2))+
  scale_y_continuous(trans="log10",expand = c(0,0),limits=c(65,450),labels = scales::comma,
                     breaks = c(75,100,200,300,400))+
  annotate(x=.001,y=410,geom="text",label="B",size=5)+
  theme(text=element_text(size=20));f4b
Fig4=ggarrange(f4a,f4b,nrow=2);Fig4
ggsave("Figure 4.pdf",Fig4,height=9,width=10)

#Table S5 TAR CV
Ts5=gls(log10(TAR_CV)~log10(IRRIGATION_AFp.001)+
           log10(PRECIP_AF)+
           TEMPC_AVG+
           log10(WETLAND_PERCENT)+
           log10(DEV_TOTAL)+
           log10(OPEN_WATER_PERCENT),
         correlation = corExp(1, form = ~ LAT_CENTER + LONG_CENTER, nugget=F), 
         data=s3, method = "ML");summary(Ts5)
shapiro.test(residuals(Ts5))
rsquared(Ts5)

#Table S6 PIP_QNX CV
s3a_PIP_QNX=aggregate(cbind(PIP_QNX_CV)~DAU_CODE+YEAR,data=s2,FUN=mean)#avg across DAU-YEAR
s3_PIP_QNX=aggregate(cbind(PIP_QNX_CV)~DAU_CODE,data=s3a_PIP_QNX,FUN=mean)#avg across DAU
s3_PIP_QNX=merge(s3_PIP_QNX,s3b)
s3_PIP_QNX$IRRIGATION_AFp.001=s3_PIP_QNX$IRRIGATION_AF+.001

Ts6=gls(log10(PIP_QNX_CV)~log10(IRRIGATION_AFp.001)+
          log10(PRECIP_AF)+
          TEMPC_AVG+
          log10(WETLAND_PERCENT)+
          log10(DEV_TOTAL)+
          log10(OPEN_WATER_PERCENT),
        correlation = corExp(1, form = ~ LAT_CENTER + LONG_CENTER, nugget=F), 
        data=s3_PIP_QNX, method = "ML");summary(Ts6)
shapiro.test(residuals(Ts6))
rsquared(Ts6)

#Table S7 ERY CV
s3a_ERY=aggregate(cbind(ERY_CV)~DAU_CODE+YEAR,data=s2,FUN=mean)#avg across DAU-YEAR
s3_ERY=aggregate(cbind(ERY_CV)~DAU_CODE,data=s3a_ERY,FUN=mean)#avg across DAU
s3_ERY=merge(s3_ERY,s3b)
s3_ERY$IRRIGATION_AFp.001=s3_ERY$IRRIGATION_AF+.001

Ts7=gls(log10(ERY_CV)~log10(IRRIGATION_AFp.001)+
          log10(PRECIP_AF)+
          TEMPC_AVG+
          log10(WETLAND_PERCENT)+
          log10(DEV_TOTAL)+
          log10(OPEN_WATER_PERCENT),
        correlation = corExp(1, form = ~ LAT_CENTER + LONG_CENTER, nugget=F), 
        data=s3_ERY, method = "ML");summary(Ts7)
shapiro.test(residuals(Ts7))
rsquared(Ts7)

#Figure 5---------------------------------------
t5a=read.csv("FIPS_MOSQ_WNV_CASES.csv",head=TRUE,sep=",")#load in human index seasonal data 
t5=t5a[t5a$OVERLAP>75&t5a$AG_2010>0,] #subsets to 52 counties with some overlap with DAU.
t5$APPLIED_WATER_ACREp.001=t5$APPLIED_WATER_ACRE+.001

Fm5=gls(log10(INCID_2004_2010+0.1)~log10(APPLIED_WATER_ACREp.001),
            correlation = corExp(1, form = ~ LAT+LONG, nugget=F), 
        data=t5, method = "ML");summary(Fm5)
rsquared(Fm5)

nd5=data.frame(APPLIED_WATER_ACREp.001=c(.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2))
nd5$pWNV=predict(Fm5,newdata=nd5)
nd5$SE=predict(Fm5,newdata=nd5,se.fit=T)$se.fit

Fig5=ggplot(t5, aes(x=APPLIED_WATER_ACREp.001,y=INCID_2004_2010))+theme_few()+
  geom_point(aes(col=MEAN_TEMPC),size=5)+
  scale_colour_gradient2(low="blue",mid="lightgray",high="red",midpoint=10)+
  geom_line(data=nd5,aes(x=APPLIED_WATER_ACREp.001,y=10^(pWNV)))+
  geom_ribbon(data=nd5,aes(x=APPLIED_WATER_ACREp.001,y=10^(pWNV),
                           ymin=10^(pWNV-SE),ymax=10^(pWNV+SE)),alpha=.1)+
  labs(x="Irrigation, Acre-feet",y="WNV human incidence",col="Temp.")+
  scale_x_continuous(expand=c(0,0),trans="log10",limits=c(.0008,2.2),
                     breaks = c(.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2))+
  scale_y_continuous(expand=c(0,0),trans="log10",limits=c(.1,210),
                     breaks = c(.01,.02,.05,.1,.2,.5,1,2,5,10,20,50,100,200))+
  theme(text=element_text(size=20));Fig5
ggsave("Figure 5.pdf",Fig5,height=9,width=12)

#Table S11
FmTs11=gls(log10(INCID_2004_2010+0.1)~log10(APPLIED_WATER_ACREp.001)+
          log10(PRECIP_AF_ACRE)+
          MEAN_TEMPC+log10(WETLAND_COVER)+
          log10(DEV_TOTAL_PERCENT_DAU)+
          log10(OPENWATER_PERCENT_DAU),
          correlation = corExp(1, form = ~ LAT+LONG, nugget=F), data=t5, method = "ML");summary(FmTs11)
rsquared(FmTs11)

FmTs10a=gls(log10(INCID_2004_2010+0.1)~log10(APPLIED_WATER_ACREp.001), 
            correlation = corExp(1, form = ~ LAT+LONG, nugget=F), 
            data=t5, method = "ML");summary(FmTs10a)

FmTs10b=gls(log10(INCID_2004_2010+0.1)~log10(PRECIP_AF_ACRE),
            correlation = corExp(1, form = ~ LAT+LONG, nugget=F), 
            data=t5, method = "ML");summary(FmTs10b)

FmTs10c=gls(log10(INCID_2004_2010+0.1)~log10(DEV_TOTAL_PERCENT_DAU),
            correlation = corExp(1, form = ~ LAT+LONG, nugget=F), 
            data=t5, method = "ML");summary(FmTs10c)

FmTs10d=gls(log10(INCID_2004_2010+0.1)~MEAN_TEMPC,
            correlation = corExp(1, form = ~ LAT+LONG, nugget=F), 
            data=t5, method = "ML");summary(FmTs10d)

FmTs10e=gls(log10(INCID_2004_2010+0.1)~log10(WETLAND_COVER),
            correlation = corExp(1, form = ~ LAT+LONG, nugget=F), 
            data=t5, method = "ML");summary(FmTs10e)

FmTs10f=gls(log10(INCID_2004_2010+0.1)~log10(OPENWATER_PERCENT_DAU),
            correlation = corExp(1, form = ~ LAT+LONG, nugget=F), 
            data=t5, method = "ML");summary(FmTs10f)

#extract info from models above
TS10=data.frame(Predictor=NA,Intercept=NA,Slope=NA,SE=NA,P_value=NA,R2=NA)
for (i in 1:6) {
  modeli=paste0("FmTs10",letters[i])
  TS10[i,]=c(rownames(summary(get(paste0("FmTs10",letters[i])))$tTable)[2],
            coef(get(modeli))["(Intercept)"],
            coef(get(modeli))[2],
            summary(get(modeli))$tTable[2,"Std.Error"],
            summary(get(modeli))$tTable[2,"p-value"],
            rsquared(get(modeli))$R.squared)  }
write.csv(TS10,"Table S10.csv")

