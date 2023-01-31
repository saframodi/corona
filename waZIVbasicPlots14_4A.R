source('/work/smodi/scripts/covid/func.R')

####################################################################################
#       panel 3
fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP//',pattern='\\.silentMM.IGGA.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.silentMM.IGGA.csv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
silentGA=loadFiles(fs)

fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.silentMM.IGDM.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.silentMM.IGDM.csv$',full.names = T))
silentDM=loadFiles(fs)

fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.MM.IGGA.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.MM.IGGA.csv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
allGA=loadFiles(fs)

fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.MM.IGDM.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.MM.IGDM.csv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
allDM=loadFiles(fs)


fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.tsv$',full.names = T)
fs2=list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
fs2=fs2[!grepl('/P',fs2)]
f=c(fs,fs2)
f=f[!grepl('AP',f)&!grepl('CT',f)]
s=as.integer(factor(ifelse(grepl('100',f),'control',ifelse(grepl('sev',f),'severe','mild/moderate'))
                    ,levels = c('control','mild/moderate','severe')))




pdf('/work/smodi/OurCovid/pdfWOziv/REPNEWplot34old.pdf',width=7.5,height=9)
par(mfrow=c(3,2))
par(cex=0.7,mar=c(3,5,3,2),mgp=c(3.5,1,0))
j=1;ti=c('All mutations IGA/G', 'Silent mutations IGA/G')
r=NULL
for(iso in list(allGA[1:79,],silentGA[1:79,])){
  k=1
  y=list()
  for(i in c('NNANN','NNTNN','NNCNN','NNGNN')){
    y[[k]]=iso[s==1,i]*256
    y[[k+1]]=iso[s==2,i]*256
    y[[k+2]]=iso[s==3,i]*256
    k=k+3
  }
  groups=paste(rep(c('A','T','C','G'),each=3),1:3)
  group=c();value=c();for(i in 1:12){
    group=c(group,rep(groups[i],length(y[[i]])))
    value=c(value,y[[i]])
  }
  q=pairwise.wilcox.test(value,group,p.adjust.method = 'BH')
  q$p.value[is.na(q$p.value)]=1
  q=data.frame(q$p.value);r=NULL
  for(i in c('A','T','C','G')){
    z=q[startsWith(rownames(q),i),startsWith(colnames(q),i)]
    for(a in rownames(z))for(b in colnames(z))if(z[a,b]<=0.05)r=rbind(r,data.frame(
      a=a,b=b,pv=z[a,b] ,j=j   ))
  }
  par(mar=c(3,5,2,2))
  boxplot(y,col=rep(c(3,'darksalmon',2),4),outline=F
          ,xaxt='n',ylab='Mutability',xlab='',las=1,main=ti[j])
  axis(side = 1,at = c(2,5,8,11),labels=c("A",'T','C','G'),lwd.ticks = T)
  if(j==1){
    
    mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
    legend('topright',legend = c('control','mild','severe'),fill=c(3,'darksalmon',2),box.col = NA        )
    
  }
  if(j==2){
    mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
  }
  for(i in 1:nrow(r)){print(r[i,])
    if(substr(r$a[i],1,1)=='A')d=0
    if(substr(r$a[i],1,1)=='T')d=3
    if(substr(r$a[i],1,1)=='C')d=6
    if(substr(r$a[i],1,1)=='G')d=9
    print(c(d+as.integer(substr(r$a[i],nchar(r$a[i]),100)),
            d+as.integer(substr(r$b[i],nchar(r$b[i]),100))))
    lines(c(d+as.integer(substr(r$a[i],nchar(r$a[i]),100)),
           d+as.integer(substr(r$b[i],nchar(r$b[i]),100))),
         rep(0.275+(as.integer(substr(r$a[i],nchar(r$a[i]),100))+
                     as.integer(substr(r$b[i],nchar(r$b[i]),100))-2)*0.015,2))
    text((d+as.integer(substr(r$a[i],nchar(r$a[i]),100))+
          d+as.integer(substr(r$b[i],nchar(r$b[i]),100)))/2,
         0.2775+(as.integer(substr(r$a[i],nchar(r$a[i]),100))+
                as.integer(substr(r$b[i],nchar(r$b[i]),100))-2)*0.015,ifelse(
                  r$pv[i]>0.01,'*',ifelse(r$pv[i]>0.001,'**','***')))
  }
  j=j+1
}

r=function(x,dy=0.0015){
  q=pairwise.wilcox.test(c(rowMeans(x[s==1,colnames(x)%in%WRC_GYW]),
                         rowMeans(x[s==2,colnames(x)%in%WRC_GYW]),
                         rowMeans(x[s==3,colnames(x)%in%WRC_GYW]),
                         rowMeans(x[s==1,colnames(x)%in%WA_TW]),
                         rowMeans(x[s==2,colnames(x)%in%WA_TW]),
                         rowMeans(x[s==3,colnames(x)%in%WA_TW]),
                         rowMeans(x[s==1,colnames(x)%in%SYC_GRS]),
                         rowMeans(x[s==2,colnames(x)%in%SYC_GRS]),
                         rowMeans(x[s==3,colnames(x)%in%SYC_GRS]),
                         rowMeans(x[s==1,colnames(x)%in%REST]),
                         rowMeans(x[s==2,colnames(x)%in%REST]),
                         rowMeans(x[s==3,colnames(x)%in%REST])),
                       paste0(rep(c(rep(1,sum(s==1)),rep(2,sum(s==2)),rep(3,sum(s==3))),4),
                             c(rep('WRC',nrow(x)),rep('WA',nrow(x)),
                                     rep('SYC',nrow(x)),rep('Neu',nrow(x)) )),p.adjust.method = 'BH')
  q=data.frame(q$p.value);q[is.na(q)]=1;rownames(q)=paste0('X',rownames(q))
  r=NULL;for(i in c('WRC','WA','SYC','Neu')){
    y=q[endsWith(rownames(q),i),endsWith(colnames(q),i)]
    for(i in rownames(y))for(j in colnames(y))if(y[i,j]<0.05)r=rbind(r,
          data.frame(group=substr(i,3,6),a=as.integer(substr(i,2,2)),b=as.integer(substr(j,2,2))
            ,pv=y[i,j]))
  }  
  r$sig=ifelse(r$pv>0.01,'*',ifelse(r$pv>0.001,'**','***'))
  for(i in 1:nrow(r)){
    if(r$group[i]=='WRC'){dy=0.0018}else dy=0.0015
    if(r$group[i]=='WRC'){dd=0}
    if(r$group[i]=='WA'){dd=3}
    if(r$group[i]=='SYC'){dd=6}
    if(r$group[i]=='Neu'){dd=9}
    lines(dd+c(r$a[i],r$b[i]),(r$a[i]+r$b[i]-2)*0.0001+c(dy,dy))
    text(dd+(r$a[i]+r$b[i])/2,dy+0.00005+(r$a[i]+r$b[i]-2)*0.0001,r$sig[i])
  }
}

x=allDM
boxplot(rowMeans(x[s==1,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==2,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==3,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==1,colnames(x)%in%WA_TW]),
        rowMeans(x[s==2,colnames(x)%in%WA_TW]),
        rowMeans(x[s==3,colnames(x)%in%WA_TW]),
        rowMeans(x[s==1,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==2,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==3,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==1,colnames(x)%in%REST]),
        rowMeans(x[s==2,colnames(x)%in%REST]),
        rowMeans(x[s==3,colnames(x)%in%REST]),
        col=rep(c(3,'darksalmon',2),4),xaxt='n',outline=F,ylab='Mean mutability',xlab='Hotspot',
        ylim=c(0,0.002),las=1,main='All mutations IGD/M')
axis(side = 1,at = c(2,5,8,11),labels=c("WRC/GYW","WA/TW","SYC/GRS","Neutral"),lwd.ticks = T)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
r(x)
x=silentDM
boxplot(rowMeans(x[s==1,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==2,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==3,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==1,colnames(x)%in%WA_TW]),
        rowMeans(x[s==2,colnames(x)%in%WA_TW]),
        rowMeans(x[s==3,colnames(x)%in%WA_TW]),
        rowMeans(x[s==1,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==2,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==3,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==1,colnames(x)%in%REST]),
        rowMeans(x[s==2,colnames(x)%in%REST]),
        rowMeans(x[s==3,colnames(x)%in%REST]),
        col=rep(c(3,'darksalmon',2),4),xaxt='n',outline=F,ylab='Mean mutability',xlab='Hotspot',
        ylim=c(0,0.0025),las=1,main='Silent mutations IGD/M')
axis(side = 1,at = c(2,5,8,11),labels=c("WRC/GYW","WA/TW","SYC/GRS","Netural"),lwd.ticks = T)
r(x)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)

x=allGA[1:79,]
boxplot(rowMeans(x[s==1,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==2,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==3,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==1,colnames(x)%in%WA_TW]),
        rowMeans(x[s==2,colnames(x)%in%WA_TW]),
        rowMeans(x[s==3,colnames(x)%in%WA_TW]),
        rowMeans(x[s==1,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==2,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==3,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==1,colnames(x)%in%REST]),
        rowMeans(x[s==2,colnames(x)%in%REST]),
        rowMeans(x[s==3,colnames(x)%in%REST]),
        col=rep(c(3,'darksalmon',2),4),xaxt='n',outline=F,ylab='Mean mutability',xlab='Hotspot',
        ylim=c(0,0.002),las=1,main='All mutations IGA/G')
axis(side = 1,at = c(2,5,8,11),labels=c("WRC/GYW","WA/TW","SYC/GRS","Neutral"),lwd.ticks = T)
mtext('E', side = 3, line = 0.5, adj = 0, cex = 1.1)
r(x)
x=silentGA[1:79,]
boxplot(rowMeans(x[s==1,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==2,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==3,colnames(x)%in%WRC_GYW]),
        rowMeans(x[s==1,colnames(x)%in%WA_TW]),
        rowMeans(x[s==2,colnames(x)%in%WA_TW]),
        rowMeans(x[s==3,colnames(x)%in%WA_TW]),
        rowMeans(x[s==1,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==2,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==3,colnames(x)%in%SYC_GRS]),
        rowMeans(x[s==1,colnames(x)%in%REST]),
        rowMeans(x[s==2,colnames(x)%in%REST]),
        rowMeans(x[s==3,colnames(x)%in%REST]),
        col=rep(c(3,'darksalmon',2),4),xaxt='n',outline=F,ylab='Mean mutability',xlab='Hotspot',
        ylim=c(0,0.0025),las=1,main='Silent mutations IGA/G')
axis(side = 1,at = c(2,5,8,11),labels=c("WRC/GYW","WA/TW","SYC/GRS","Netural"),lwd.ticks = T)
mtext('F', side = 3, line = 0.5, adj = 0, cex = 1.1)
r(x)
dev.off()

for(i in list())


w=c();t=c()
for(i in 1:4)for(j in 1:4){
  if(i==1){x=allDM}else if(i==2){x=allGA}else if(i==3){x=silentDM}else if(i==4)x=silentGA
  if(j==1){y=WRC_GYW}else if(j==2){y=WA_TW}else if(j==3){y=SYC_GRS}else if (j==4)y=REST
  w=c(w,wilcox.test(rowMeans(x[s==1,colnames(x)%in%y]),
                    rowMeans(x[s!=1,colnames(x)%in%y]))$p.value)
  t=c(t,t.test(rowMeans(x[s==1,colnames(x)%in%y]),
                    rowMeans(x[s!=1,colnames(x)%in%y]))$p.value)
  
}

##########################################################################
#     Panel 5
dir='/work/smodi/OurCovid/28_10_2021_long/'
fs=list.files(dir,pattern='\\.MM.csv$',full.names = T)
lightAll=loadFiles(fs[grepl('LIGHT',fs)])
lightAll$stage=ifelse(grepl('LIGHTCov',fs[grepl('LIGHT',fs)]),'case','control')
lightAll$X=NULL;lightAll$id=NULL
lightrs=lightAll



t=lightrs[,grepl('tar',colnames(lightrs))]
t$stage=lightrs$stage
pall=traintest(t[1:60,],pv=1,iterations = 50)


dir='/work/smodi/OurCovid/28_10_2021_long/'
fs=list.files(dir,pattern='\\.silentMM.csv$',full.names = T)
lightAll=loadFiles(fs[grepl('LIGHT',fs)])
lightAll$stage=ifelse(grepl('LIGHTCov',fs[grepl('LIGHT',fs)]),'case','control')
lightAll$X=NULL;lightAll$id=NULL
t=lightAll[,grepl('tar',colnames(lightAll))]
t$stage=lightrs$stage
p=traintest(t[1:60,],pv=1,iterations = 50)

pdf('/work/smodi/OurCovid/pdfWOziv/NEWplotLight.pdf',width=7.5,height=3)
par(cex=0.7)
par(mfrow=c(1,2))
boxplot(pall,las=2,ylim=c(0,1))
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(p,las=2,ylim=c(0,1))
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()

t1=lightAll[,grepl('tar',colnames(lightAll))]
t=allGA[,grepl('sub',colnames(allGA))]
t=t[1:60,]
t$id=1:60;t1$id=1:60
hl=merge(t,t1,by='id');hl$id=NULL
hl=hl[,colSums(hl)>0]
hl$stage=c(rep('control',28),rep('case',32))
p=traintest(hl,0.1,iterations = 50)
boxplot(p)


fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP//',pattern='\\.silentMM.IGGA.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.silentMM.IGGA.csv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
silentGA=loadFiles(fs)

fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.silentMM.IGDM.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.silentMM.IGDM.csv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
silentDM=loadFiles(fs)

fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.MM.IGGA.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.MM.IGGA.csv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
allGA=loadFiles(fs)
allGA$id=fs
fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.MM.IGDM.csv$',full.names = T)
fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.MM.IGDM.csv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
allDM=loadFiles(fs)


fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='\\.tsv$',full.names = T)
fs2=list.files('/work/smodi/OurCovid/28_10_2021_longREP/more/',pattern='\\.tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
fs2=fs2[!grepl('/P',fs2)]
fs2=fs2[!grepl('AP',fs2)&!grepl('CT',fs2)]
f=c(fs,fs2)
s=as.integer(factor(ifelse(grepl('100',f),'control',ifelse(grepl('sev',f),'severe','mild/moderate'))
                    ,levels = c('control','mild/moderate','severe')))

z=allDM[1:79,]
z$stage=ifelse(s==1,'move',ifelse(s==2,'control','case'))
z=z[s!=1,]
z$X=NULL;z$id=NULL
z=z[,colnames(z)=='stage'|substr(colnames(z),1,5)%in%SYC_GRS]
con=z[z$stage=='control',]
z=z[z$stage!='control',]
p=list();f=list()
for(i in 1:20){
  q=rbind(z[z$stage=='case',],con[order(runif(nrow(con)))[1:12],])
  #q=rbind(q,x2)
  q2=leaveOneOut(q,pv=0.1)
  p[[i]]=data.frame(F1=F1_Score(q$stage,q2,positive='case'),
                    acc=Accuracy(q$stage,q2),
                    sen=Sensitivity(q$stage,q2,positive='case'),
                    spe=Specificity(q$stage,q2,positive='case'))
  f[[i]]=extractFeatures(q,1)
  boxplot(rbindlist(p),main=i)
}
q=f
boxplot(rbindlist(p))
if(run=='save'){
  write.csv(data.frame(rbindlist(p)),'/work/smodi/OurCovid/REPseverityPrediction.csv')

}
t=allGA[,grepl('sub',colnames(allGA))]
f=allGA$id
s=as.integer(factor(ifelse(grepl('100',f)|grepl('/P',f),'control',ifelse(grepl('sev',f),'severe','mild/moderate'))
                    ,levels = c('control','mild/moderate','severe')))
t$stage=c(ifelse(s==2|s==3,'case','control'))
substitution=traintest(t[1:60,],pv=0.001,iterations = 50)
boxplot(substitution)
pv=0.001
x=t[1:60,]
k=rep(T,ncol(x))
if(pv<1){
  k=c()
  for(i in 1:ncol(x))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
    k[is.na(k)]=1
    k=k<pv
}
r=train(stage~.,x[,k],method='glmnet',tuneLength=5)
val=predict(r,t[61:nrow(t),]);valTrue=t[61:nrow(t),'stage']
sum(val==valTrue)

if(run=='save')write.csv(substitution,'/work/smodi/OurCovid/REPsubstitutionPrediction.csv')
t=silentGA[,!grepl('N',colnames(silentGA))&!grepl('sub',colnames(silentGA))]
t$stage=c(ifelse(s==1,'control','case'),rep('control',nrow(t)-length(s)))
t$X=NULL;t$id=NULL
silent=traintest(t[1:60,],pv=0.001,iterations = 50)
boxplot(silent)

pv=0.001
x=t[1:60,]
k=rep(T,ncol(x))
if(pv<1){
  k=c()
  for(i in 1:ncol(x))
    if(colnames(x)[i]=='stage'){
      k=c(k,0)
    }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
    k[is.na(k)]=1
    k=k<pv
}
r=train(stage~.,x[,k],method='glmnet',tuneLength=5)
valSilent=predict(r,t[61:nrow(t),]);valSilentTrue=t[61:nrow(t),'stage']
sum(valSilent==valSilentTrue)
if(run=='save')write.csv(silent,'/work/smodi/OurCovid/REPsilentPrediction.csv')
f=q
f=data.frame(rbindlist(q))
if(run=='save')write.csv(f,'/work/smodi/OurCovid/REPFofPrediction.csv')



write.csv(rbind(data.frame(accuracy=Accuracy(val,valTrue),sensitivity=Sensitivity(valTrue,val,positive='case'),
           specificity=Specificity(valTrue,val,positive='case'),F1=F1_Score(valTrue,
                                          val,positive='case'),group='validation'),
data.frame(accuracy=Accuracy(valSilent,valSilentTrue),
           sensitivity=Sensitivity(valSilentTrue,valSilent,positive='case'),
           specificity=Specificity(valSilentTrue,valSilent,positive='case'),
           F1=F1_Score(valSilentTrue,valSilent,positive='case'),group='validationSilent')),
  '/work/smodi/OurCovid/28_10_2021_long/more/REPvalidation.csv')



pdf('/work/smodi/OurCovid/pdfWOziv/REPNEWplot5validation.pdf',width=7.5,height=3)
par(cex=0.7)
par(mfrow=c(1,2))
barplot(c(F1_Score(valTrue,val,positive='case'),
        Accuracy(val,valTrue),Sensitivity(valTrue,val,positive='case'),
        Specificity(valTrue,val,positive='case')),
        names.arg=c('F1 score','accuracy','sensitivity','specificity'),ylim=c(0,1),
        las=2)
box()
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
barplot(c(F1_Score(valSilentTrue,valSilent,positive='case'),
          Accuracy(valSilent,valSilentTrue),Sensitivity(valSilentTrue,valSilent,positive='case'),
          Specificity(valSilentTrue,valSilent,positive='case')),
        names.arg=c('F1 score','accuracy','sensitivity','specificity'),ylim=c(0,1),
        las=2)
box()
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()

substitution=read.csv('/work/smodi/OurCovid/REPsubstitutionPrediction.csv')
silent=read.csv('/work/smodi/OurCovid/REPsilentPrediction.csv')
p=read.csv('/work/smodi/OurCovid/REPseverityPrediction.csv')
f=read.csv('/work/smodi/OurCovid/REPFofPrediction.csv')
library(plyr)
z=allDM
fs=list.files('/work/smodi/OurCovid/28_10_2021_long/',pattern='\\.tsv$',full.names = T)
fs2=list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='\\.tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
fs2=fs2[!grepl('/P',fs2)]
fs2=fs2[!grepl('AP',fs2)&!grepl('CT',fs2)]
f=c(fs,fs2)
s=as.integer(factor(ifelse(grepl('100',f),'control',ifelse(grepl('sev',f),'severe','mild/moderate'))
                    ,levels = c('control','mild/moderate','severe')))
s=s[1:79]

z$stage=s
z=z[s!=1,]
z$X=NULL;z$id=NULL
z=z[,colnames(z)=='stage'|substr(colnames(z),1,5)%in%SYC_GRS]
con=z[z$stage==2,]
z=z[z$stage!=2,]


f=read.csv('/work/smodi/OurCovid/FofPrediction.csv')
t=data.frame(table(f$gene))
t=t[t$Freq>=8,]
t=t[order(t$Freq,decreasing = T),]
q=list();j=1
zz=count(f[f$gene%in%t$Var1,],'gene','s1')
zz=zz[zz$freq>0,];zz=zz[order(zz$freq,decreasing = T),]
for(i in zz$gene){
  q[[j]]=con[,i];q[[j+1]]=z[,i];j=j+2
}
boxplot(q)
substitution$X=NULL;silent$X=NULL;p$X=NULL
pdf('/work/smodi/OurCovid/pdfWOziv/NEWplot5.pdf',width=7.5,height=3)
par(cex=0.7)
par(mar=c(5,3,3,1))
layout(t(c(1,2,3,4,4)))
boxplot(substitution,las=2,ylim=c(0,1),las=2)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(silent,names=names(substitution),las=2,ylim=c(0,1),las=2)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot((p),names=names(substitution),las=2,ylim=c(0,1),las=2)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
par(mar=c(5,5,3,1),mgp=c(4,0.5,0))
boxplot(q,col=rep(c('darksalmon',2),7),ylim=c(0,0.00121),las=1,
        xaxt='n',outline=F,ylab='Mutability frequency')
axis(side = 1,at = 1.5+2*(0:7),labels=zz$gene,lwd.ticks = T,las=2)
legend('topleft',legend = c('mild','severe'),fill=c('darksalmon',2),box.col = NA)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()




