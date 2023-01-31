source('/work/smodi/scripts/covid/func.R')
pv=function(p){
  value=c();group=c();gene=c();
  j=1;jj=1
  for(i in p){
    value=c(value,i)
    group=c(group,rep(paste0('X',j),length(i)))
    gene=c(gene,rep(gsub(' ',';',gsub('-','_',names(p)[jj])),length(i)))
    j=j+1;jj=jj+1;if(j==4)j=1
  }
  z=data.frame(value=value,group=group,gene=gene,class=paste0(gene,group))
  z=z[!is.na(z$gene),]
  q=pairwise.wilcox.test(value,z$class,p.adjust.method = 'BH')
  q$p.value[is.na(q$p.value)]=1
  q=data.frame(q$p.value)
  rownames(q)=gsub('-','.',gsub(' ','.',rownames(q)))
  colnames(q)=gsub('\\.',';',colnames(q))
  r=NULL
  for(i in unique(z$gene)){
    i=gsub('-','.',gsub(' ','.',i))
    y=q[substr(rownames(q),1,nchar(rownames(q))-2)==i,
        substr(colnames(q),1,nchar(colnames(q))-2)==i]
    for(a in rownames(y))for(b in colnames(y))if(y[a,b]<=0.05)r=rbind(r,data.frame(
      a=a,b=b,pv=y[a,b]    ))
  }
  r$gene=gsub(';',' ',gsub('_','-',substr(r$a,1,nchar(r$a)-2)))
  r$a=as.integer(substr(r$a,nchar(r$a),nchar(r$a)))
  r$b=as.integer(substr(r$b,nchar(r$b),nchar(r$b)))
  r$sig=ifelse(r$pv>0.05,'',ifelse(r$pv>0.01,'*',ifelse(r$pv>0.001,'**','***')))
  return(r)
}
print(pv(d));print(pv(J));print(pv(vj));print(pv(vjjl));print(pv(v))

########################################################################################
# 		Create basic data


f=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
f=f[grepl('/IGH',f)]
y=as.integer(as.factor(ifelse(grepl('100',f),'control',ifelse(grepl('seve',f),'severe','mild/moderate'))))
f=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
f=f[grepl('/IGL',f)]



p=mclapply(f,function(i){
  r=data.frame(id=i)
  rep=read.delim(i)
  l=sort(rep[!is.na(rep$junction_aa_length),]$junction_aa_length)
  r=merge(r,data.frame(jl10=l[0.1*length(l)],
                       jl50=median(l),jl90=l[0.9*length(l)]))
  l=sort(rep[!is.na(rep$v_identity),]$v_identity)
  r=merge(r,data.frame(vi10=l[0.1*length(l)],
                       vi50=median(l),vi90=l[0.9*length(l)]))
  
},mc.cores=50)
library(data.table)
p=data.frame(rbindlist(p))
library(openxlsx)
z=p
z$id=NULL
z[,4]=(1-z[,4]);z[,5]=(1-z[,5]);z[,6]=(1-z[,6])
basics=z



p=mclapply(f,function(i){
  id=gsub('\\.fasta.tab','',gsub('/work/smodi/crohn/igblast//','',i))
  rep=read.delim(i)
  rep$v=getGene(rep$v_call)
  rep$j=getGene(rep$j_call)
  rep$VF=getFamily(rep$v_call)
  rep$vj=paste(getGene(rep$v_call),getGene(rep$j_call))
  r=data.frame(table(rep$vj))
  r$Freq=r$Freq/sum(r$Freq)
  r$id=id
  return(r)
},mc.cores=50)
p=data.frame(rbindlist(p))
m=merge(levels(as.factor(p$id)),levels(as.factor(p$Var1)))
m$p=paste(m$x,m$y)
colnames(m)=c('id','Var1','p')
p$p=paste(p$id,p$Var1)
m=merge(m,p,by='p',all=T)
m$p=NULL;m$Var1.y=NULL;m$id.y=NULL
colnames(m)=c('id','vj','freq')
m$freq=ifelse(is.na(m$freq),0,m$freq)
z=m
library(plyr)
l=count(z,'vj','freq')
l=l[order(l$freq,decreasing = T),]
l=l[1:50,]
z=z[z$vj%in%l$vj,]
id=data.frame(id=f,stage=y)
z=merge(z,id,by='id')
p=list()
j=1
for(i in unique(z$vj)){
  p[[j*3]]=z[z$stage==3&z$vj==i,]$freq
  p[[j*3-1]]=z[z$stage==2&z$vj==i,]$freq
  p[[j*3-2]]=z[z$stage==1&z$vj==i,]$freq
  j=j+1
}
names(p)=rep(unique(z$vj),each=3)
vj=p

p=mclapply(f,function(i){
  id=gsub('\\.fasta.tab','',gsub('/work/smodi/crohn/igblast//','',i))
  rep=read.delim(i)
  rep$v=getGene(rep$v_call)
  rep$j=getGene(rep$j_call)
  rep$VF=getFamily(rep$v_call)
  rep$vj=(getGene(rep$v_call))
  r=data.frame(table(rep$vj))
  r$Freq=r$Freq/sum(r$Freq)
  r$id=id
  return(r)
},mc.cores=50)
p=data.frame(rbindlist(p))
m=merge(levels(as.factor(p$id)),levels(as.factor(p$Var1)))
m$p=paste(m$x,m$y)
colnames(m)=c('id','Var1','p')
p$p=paste(p$id,p$Var1)
m=merge(m,p,by='p',all=T)
m$p=NULL;m$Var1.y=NULL;m$id.y=NULL
colnames(m)=c('id','vj','freq')
m$freq=ifelse(is.na(m$freq),0,m$freq)
z=m
z=merge(z,id,by='id')
library(plyr)
l=count(z,'vj','freq')
l=l[order(l$freq,decreasing = T),]
l=l[1:50,]
z=z[z$vj%in%l$vj,]
p=list()
j=1
for(i in unique(z$vj)){
  p[[j*3]]=z[z$stage==3&z$vj==i,]$freq
  p[[j*3-1]]=z[z$stage==2&z$vj==i,]$freq
  p[[j*3-2]]=z[z$stage==1&z$vj==i,]$freq
  j=j+1
}
names(p)=rep(unique(z$vj),each=3)
v=p
###############################################################################
#		Plot S1

pdf('/work/smodi/OurCovid/pdfWOziv/NEWplotS1IGL.pdf',width=7.5,height=7.5)
layout(rbind(c(1,2),c(3,3),c(4,4)))
par(mar=c(4.5,4.5,2,1))
z=basics
boxplot(z$jl10[y==1],z$jl10[y==2],z$jl10[y==3],
        z$jl50[y==1],z$jl50[y==2],z$jl50[y==3],
        z$jl90[y==1],z$jl90[y==2],z$jl90[y==3],
        col=rep(c(3,'darksalmon',2),4),xaxt='n',ylim=c(10,15),
        outline=F,ylab='Junction AA length',xlab='Precentile',las=1)
axis(side = 1,at = c(2,5,8),labels=c("10","50","90"),lwd.ticks = T)
legend('topleft',legend = c(
  'control','mild','severe'),
  fill=c(3,'darksalmon',2),box.col = NA)
group=paste(rep(c(rep(1,sum(y==1)),rep(2,sum(y==2)),rep(3,sum(y==3))),3)
            ,c(rep(10,length(y)),rep(50,length(y)),rep(90,length(y))))
r=pairwise.wilcox.test(c(z$jl10[y==1],z$jl10[y==2],z$jl10[y==3],
                         z$jl50[y==1],z$jl50[y==2],z$jl50[y==3],
                         z$jl90[y==1],z$jl90[y==2],z$jl90[y==3]),group
                       ,p.adjust.method = 'BH')
lines(c(7,8),c(14.5,14.5))
lines(c(4,6),c(14,14))
text(c(7.5,5),c(14.7,14.2),c('*','*'))

mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(z$vi90[y==1],z$vi90[y==2],z$vi90[y==3],
        z$vi50[y==1],z$vi50[y==2],z$vi50[y==3],
        z$vi10[y==1],z$vi10[y==2],z$vi10[y==3],
        col=rep(c(3,'darksalmon',2),4),xaxt='n',
        outline=F,ylab='V distance from germline',xlab='Precentile',las=1)
axis(side = 1,at = c(2,5,8),labels=c("10","50","90"),lwd.ticks = T)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
pairwise.wilcox.test(c(z$vi90[y==1],z$vi90[y==2],z$vi90[y==3],
                       z$vi50[y==1],z$vi50[y==2],z$vi50[y==3],
                       z$vi10[y==1],z$vi10[y==2],z$vi10[y==3]),
                     paste(rep(c(rep(1,sum(y==1)),rep(2,sum(y==2)),rep(3,sum(y==3))),3)
                           ,c(rep(10,length(y)),rep(50,length(y)),rep(90,length(y)))),p.adjust.method = 'BH')


p=v
par(mar=c(4,4.5,2.5,1))
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),37),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:36)*3.5,2+c(0:36)*3.5,c(0:36)*3.5+3)),ylim=c(0,0.35))
axis(side = 1,at = 3.5*c(0:36)+2,labels=gsub('IGL','',names(p)[3*1:37-1]),lwd.ticks = T,las=2)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
r=pv(p);dy=0.3;ry=0.015
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3;ddy=(r$a[i]+r$b[i]-2)*ry
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy+ddy,dy+ddy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+ddy+ry/2,r$sig[i])
}

p=vj
par(mar=c(5.5,4.5,1,1))
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),50),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:49)*3.5,2+c(0:49)*3.5,c(0:49)*3.5+3)),ylim=c(0,0.2))
axis(side = 1,at = 3.5*c(0:49)+2,labels=gsub('IGL','',names(p)[3*1:50-1]),lwd.ticks = T,las=2)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
r=pv(p);dy=0.15;ry=0.01
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3;ddy=(r$a[i]+r$b[i]-2)*ry
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy+ddy,dy+ddy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+ddy+ry/2,r$sig[i])
}

dev.off()

####################################################################################
source('/work/smodi/scripts/covid/func.R')
########################################################################################
# 		Create basic data


f=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
f=f[grepl('/IGH',f)]
y=as.integer(as.factor(ifelse(grepl('100',f),'control',ifelse(grepl('seve',f),'severe','mild/moderate'))))
f=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
f=f[grepl('/IGK',f)]



p=mclapply(f,function(i){
  r=data.frame(id=i)
  rep=read.delim(i)
  l=sort(rep[!is.na(rep$junction_aa_length),]$junction_aa_length)
  r=merge(r,data.frame(jl10=l[0.1*length(l)],
                       jl50=median(l),jl90=l[0.9*length(l)]))
  l=sort(rep[!is.na(rep$v_identity),]$v_identity)
  r=merge(r,data.frame(vi10=l[0.1*length(l)],
                       vi50=median(l),vi90=l[0.9*length(l)]))
  
},mc.cores=50)
library(data.table)
p=data.frame(rbindlist(p))
library(openxlsx)
z=p
z$id=NULL
z[,4]=(1-z[,4]);z[,5]=(1-z[,5]);z[,6]=(1-z[,6])
basics=z



p=mclapply(f,function(i){
  id=gsub('\\.fasta.tab','',gsub('/work/smodi/crohn/igblast//','',i))
  rep=read.delim(i)
  rep$v=getGene(rep$v_call)
  rep$j=getGene(rep$j_call)
  rep$VF=getFamily(rep$v_call)
  rep$vj=paste(getGene(rep$v_call),getGene(rep$j_call))
  r=data.frame(table(rep$vj))
  r$Freq=r$Freq/sum(r$Freq)
  r$id=id
  return(r)
},mc.cores=50)
p=data.frame(rbindlist(p))
m=merge(levels(as.factor(p$id)),levels(as.factor(p$Var1)))
m$p=paste(m$x,m$y)
colnames(m)=c('id','Var1','p')
p$p=paste(p$id,p$Var1)
m=merge(m,p,by='p',all=T)
m$p=NULL;m$Var1.y=NULL;m$id.y=NULL
colnames(m)=c('id','vj','freq')
m$freq=ifelse(is.na(m$freq),0,m$freq)
z=m
library(plyr)
l=count(z,'vj','freq')
l=l[order(l$freq,decreasing = T),]
l=l[1:50,]
z=z[z$vj%in%l$vj,]
id=data.frame(id=f,stage=y)
z=merge(z,id,by='id')
p=list()
j=1
for(i in unique(z$vj)){
  p[[j*3]]=z[z$stage==3&z$vj==i,]$freq
  p[[j*3-1]]=z[z$stage==2&z$vj==i,]$freq
  p[[j*3-2]]=z[z$stage==1&z$vj==i,]$freq
  j=j+1
}
names(p)=rep(unique(z$vj),each=3)
vj=p

p=mclapply(f,function(i){
  id=gsub('\\.fasta.tab','',gsub('/work/smodi/crohn/igblast//','',i))
  rep=read.delim(i)
  rep$v=getGene(rep$v_call)
  rep$j=getGene(rep$j_call)
  rep$VF=getFamily(rep$v_call)
  rep$vj=(getGene(rep$v_call))
  r=data.frame(table(rep$vj))
  r$Freq=r$Freq/sum(r$Freq)
  r$id=id
  return(r)
},mc.cores=50)
p=data.frame(rbindlist(p))
m=merge(levels(as.factor(p$id)),levels(as.factor(p$Var1)))
m$p=paste(m$x,m$y)
colnames(m)=c('id','Var1','p')
p$p=paste(p$id,p$Var1)
m=merge(m,p,by='p',all=T)
m$p=NULL;m$Var1.y=NULL;m$id.y=NULL
colnames(m)=c('id','vj','freq')
m$freq=ifelse(is.na(m$freq),0,m$freq)
z=m
z=merge(z,id,by='id')
library(plyr)
l=count(z,'vj','freq')
l=l[order(l$freq,decreasing = T),]
l=l[1:50,]
z=z[z$vj%in%l$vj,]
p=list()
j=1
for(i in unique(z$vj)){
  p[[j*3]]=z[z$stage==3&z$vj==i,]$freq
  p[[j*3-1]]=z[z$stage==2&z$vj==i,]$freq
  p[[j*3-2]]=z[z$stage==1&z$vj==i,]$freq
  j=j+1
}
names(p)=rep(unique(z$vj),each=3)
v=p
###############################################################################
#		Plot S1

pdf('/work/smodi/OurCovid/pdfWOziv/NEWplotS1IGK.pdf',width=7.5,height=7.5)

layout(rbind(c(1,2),c(3,3),c(4,4)))
par(mar=c(5,4.5,2,1))
z=basics
boxplot(z$jl10[y==1],z$jl10[y==2],z$jl10[y==3],
        z$jl50[y==1],z$jl50[y==2],z$jl50[y==3],
        z$jl90[y==1],z$jl90[y==2],z$jl90[y==3],
        col=rep(c(3,'darksalmon',2),4),xaxt='n',ylim=c(10,13),
        outline=F,ylab='Junction AA length',xlab='Precentile',las=1)
axis(side = 1,at = c(2,5,8),labels=c("10","50","90"),lwd.ticks = T)
legend('topleft',legend = c(
  'control','mild','severe'),
  fill=c(3,'darksalmon',2),box.col = NA)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)

group=paste(rep(c(rep(1,sum(y==1)),rep(2,sum(y==2)),rep(3,sum(y==3))),3)
            ,c(rep(10,length(y)),rep(50,length(y)),rep(90,length(y))))
r=pairwise.wilcox.test(c(z$jl10[y==1],z$jl10[y==2],z$jl10[y==3],
                         z$jl50[y==1],z$jl50[y==2],z$jl50[y==3],
                         z$jl90[y==1],z$jl90[y==2],z$jl90[y==3]),group
                       ,p.adjust.method = 'BH')
lines(c(8,9),c(12.5,12.5))
text(c(8.5),c(12.7),c('*'))


boxplot(z$vi90[y==1],z$vi90[y==2],z$vi90[y==3],
        z$vi50[y==1],z$vi50[y==2],z$vi50[y==3],
        z$vi10[y==1],z$vi10[y==2],z$vi10[y==3],
        col=rep(c(3,'darksalmon',2),4),xaxt='n',
        outline=F,ylab='V distance from germline',xlab='Precentile',las=1)
axis(side = 1,at = c(2,5,8),labels=c("10","50","90"),lwd.ticks = T)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)

pairwise.wilcox.test(c(z$vi90[y==1],z$vi90[y==2],z$vi90[y==3],
                       z$vi50[y==1],z$vi50[y==2],z$vi50[y==3],
                       z$vi10[y==1],z$vi10[y==2],z$vi10[y==3]),
                     paste(rep(c(rep(1,sum(y==1)),rep(2,sum(y==2)),rep(3,sum(y==3))),3)
                           ,c(rep(10,length(y)),rep(50,length(y)),rep(90,length(y)))),p.adjust.method = 'BH')


p=v
par(mar=c(5,4.5,1,1))
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),35),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:34)*3.5,2+c(0:34)*3.5,c(0:34)*3.5+3)))
axis(side = 1,at = 3.5*c(0:34)+2,labels=gsub('IGK','',names(p)[2*1:35-1]),lwd.ticks = T,las=2)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
r=pv(p);dy=0.2;ry=0.015
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3;ddy=(r$a[i]+r$b[i]-2)*ry
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy+ddy,dy+ddy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+ddy+ry/2,r$sig[i])
}



par(mar=c(5.5,4.5,1,1))
p=vj
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),50),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:49)*3.5,2+c(0:49)*3.5,c(0:49)*3.5+3)))
axis(side = 1,at = 3.5*c(0:49)+2,labels=gsub('IGK','',names(p)[2*1:50-1]),lwd.ticks = T,las=2)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)


r=pv(p);dy=0.08;ry=0.005
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3;ddy=(r$a[i]+r$b[i]-2)*ry
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy+ddy,dy+ddy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+ddy+ry/2,r$sig[i])
}

dev.off()


