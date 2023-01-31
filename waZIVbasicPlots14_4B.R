source('/work/smodi/scripts/covid/func.R')
if(run=='create'){
  createMutationLoad(dir='/work/smodi/OurCovid/28_10_2021_long/')
  createMutationLoad(dir='/work/smodi/OurCovid/28_10_2021_long/more/')
  
  
  f=list.files('/work/smodi/OurCovid/punmock/',pattern='tsv$',full.names = T)
  mclapply(f,function(i){
    rep=read.delim(i)
    write.csv(MM(rep,'rs'),paste0(i,'.MM.csv'))
    write.csv(MM(rep,'s'),paste0(i,'.silentMM.csv'))
    
  },mc.cores=55)
  x=loadFiles(list.files('/work/smodi/OurCovid/punmock/',pattern='\\.MM.csv$',full.names = T))
  f=list.files('/work/smodi/OurCovid/punmock/',pattern='\\.MM.csv$',full.names = T)
  s=ifelse(grepl('Sev',f),3,ifelse(grepl('Mild',f),1,2))
  punmock=x;punmocks=s
}
########################################################################################
# 		Create basic data

fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
fs2=(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
fs2=fs2[!grepl('/P',fs2)&!grepl('AP',fs2)&!grepl('CT',fs2)]

f=c(fs,fs2)
y=as.integer(as.factor(ifelse(grepl('100',f),'control',ifelse(grepl('seve',f),'severe','mild/moderate'))))


p=mclapply(f,function(i){
  id=gsub('\\.fasta.tab','',gsub('/work/smodi/crohn/igblast//','',i))
  rep=read.delim(i)
  rep$v=getGene(rep$v_call)
  rep$j=getGene(rep$j_call)
  rep$VF=getFamily(rep$v_call)
  rep$vj=paste(getGene(rep$d_call))
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
z$stage=as.integer(as.factor(ifelse(grepl('100',z$id),'control',ifelse(grepl('sev',z$id),'severe','mild/moderate'))))
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
d=p;d[[1]]=NULL;d[[1]]=NULL;d[[1]]=NULL
D=d

p=mclapply(f,function(i){
  id=gsub('\\.fasta.tab','',gsub('/work/smodi/crohn/igblast//','',i))
  rep=read.delim(i)
  rep$v=getGene(rep$v_call)
  rep$j=getGene(rep$j_call)
  rep$VF=getFamily(rep$v_call)
  rep$vj=paste(getGene(rep$j_call))
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
z$stage=as.integer(as.factor(ifelse(grepl('100',z$id),'control',ifelse(grepl('sev',z$id),'severe','mild/moderate'))))
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
J=p



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
library(alakazam)
library(plyr)
p=mclapply(f,function(i){
  a=read.delim(i)
  aa=data.frame(clone_id=a[!is.na(a$clone_id),]$clone_id)
  e=alphaDiversity(a,max_n=2500)
  e@diversity$group=gsub('/work/smodi/crohn/igblast//','',
                         gsub('\\.fasta_clone-pass.tsv','',i))
  return(e@diversity)
},mc.cores=50)
library(data.table)
p2=rbindlist(p)
p=data.frame(p2)
q=p[p$q%in%c(0,1,2,3,4),]
q$stage=ifelse(grepl('100',q$group),'control',ifelse(grepl('seve',q$group),'severe','mild/moderate'))
q$stage=as.factor(q$stage)
q$stage=as.integer(q$stage)
Q=q



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
z$stage=as.integer(as.factor(ifelse(grepl('100',z$id),'control',ifelse(grepl('sev',z$id),'severe','mild/moderate'))))
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
z$stage=as.integer(as.factor(ifelse(grepl('100',z$id),'control',ifelse(grepl('sev',z$id),'severe','mild/moderate'))))
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




p=mclapply(f,function(i){
  id=gsub('\\.fasta.tab','',gsub('/work/smodi/crohn/igblast//','',i))
  rep=read.delim(i)
  rep$v=getGene(rep$v_call)
  rep$j=getGene(rep$j_call)
  rep$VF=getFamily(rep$v_call)
  rep=rep[rep$cdr3_aa!='',]
  rep$vj=paste(getGene(rep$v_call),getGene(rep$j_call),nchar(rep$cdr3_aa))
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
z$stage=as.integer(as.factor(ifelse(grepl('100',z$id),'control',ifelse(grepl('sev',z$id),'severe','mild/moderate'))))
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
vjjl=p
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
  q=pairwise.wilcox.test(value,z$class,p.adjust.method = 'BH')
  q$p.value[is.na(q$p.value)]=1
  q=data.frame(q$p.value)
  rownames(q)=gsub('-','.',gsub(' ','.',rownames(q)))
  colnames(q)=gsub('\\.',';',colnames(q))
  r=NULL
  for(i in unique(gene)){
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


pdf('/work/smodi/OurCovid/pdfWOziv/NEWplotS1new.pdf',width=7.5,height=5.5)
layout(rbind(c(1,1,1,1,2),c(3,3,3,3,3),c(4,4,4,4,4)))
par(mar=c(5.5,4,2,1))
d=D
p=d
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),25),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:24)*3.5,2+c(0:24)*3.5,c(0:24)*3.5+3)))
axis(side = 1,at = 3.5*c(0:24)+2,labels=gsub('IGH','',names(p)[3*1:25-1]),lwd.ticks = T,las=2)
r=pv(p);dy=0.15
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy,dy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+0.01,r$sig[i])
}
legend('topright',legend = c(
  'control','mild','severe'),
  fill=c(3,'darksalmon',2),box.col = NA)

mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)

p=J
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),6),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:5)*3.5,2+c(0:5)*3.5,c(0:5)*3.5+3)))
axis(side = 1,at = 3.5*c(0:5)+2,labels=gsub('IGH','',names(p)[3*1:6]),lwd.ticks = T,las=2)
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
r=pv(p);dy=0.5
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy,dy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+0.01,r$sig[i])
}

p=vj
par(mar=c(5.5,4,1,1))
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),50),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:49)*3.5,2+c(0:49)*3.5,c(0:49)*3.5+3)))
axis(side = 1,at = 3.5*c(0:49)+2,labels=gsub('IGH','',names(p)[3*1:50-1]),lwd.ticks = T,las=2)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
r=pv(p);dy=0.07;ry=0.0075
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3;ddy=(r$a[i]+r$b[i]-2)*ry
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy+ddy,dy+ddy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+ddy+ry/2,r$sig[i])
}

p=vjjl
par(mar=c(7,4,1,1))
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),50),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:49)*3.5,2+c(0:49)*3.5,c(0:49)*3.5+3)))
axis(side = 1,at = 3.5*c(0:49)+2,labels=gsub('IGH','',names(p)[3*1:50-1]),lwd.ticks = T,las=2)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
r=pv(p);dy=0.0135;ry=0.001
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3;ddy=(r$a[i]+r$b[i]-2)*ry
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy+ddy,dy+ddy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+ddy+ry/2,r$sig[i])
}
dev.off()


pdf('/work/smodi/OurCovid/pdfWOziv/NEWplot1new.pdf',width=7.5,height=8)

layout(rbind(c(1,1,2,2,3,3),c(4,4,4,4,4,4),c(5,5,5,6,6,6)))
par(mar=c(4,4,2,1))
z=basics
boxplot(z$jl10[y==1],z$jl10[y==2],z$jl10[y==3],
        z$jl50[y==1],z$jl50[y==2],z$jl50[y==3],
        z$jl90[y==1],z$jl90[y==2],z$jl90[y==3],
        col=rep(c(3,'darksalmon',2),4),xaxt='n',
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
lines(c(1,2),c(0.05,0.05))
lines(c(1,3),c(0.08,0.08))
text(c(1.5,2),c(0.06,0.09),c('**','**'))
q=Q
boxplot(q[q$q==0&q$stage==1,]$d,
        q[q$q==0&q$stage==2,]$d,
        q[q$q==0&q$stage==3,]$d,
        q[q$q==1&q$stage==1,]$d,
        q[q$q==1&q$stage==2,]$d,
        q[q$q==1&q$stage==3,]$d,
        q[q$q==2&q$stage==1,]$d,
        q[q$q==2&q$stage==2,]$d,
        q[q$q==2&q$stage==3,]$d,
        q[q$q==3&q$stage==1,]$d,
        q[q$q==3&q$stage==2,]$d,
        q[q$q==3&q$stage==3,]$d,
        q[q$q==4&q$stage==1,]$d,
        q[q$q==4&q$stage==2,]$d,
        q[q$q==4&q$stage==3,]$d,
        col=rep(c(3,'darksalmon',2),4),xaxt='n',ylab='Diversity',xlab='q'  ,las=1 ,ylim=c(0,2400),
        at=c(1,2,3,4.5,5.5,6.5,8,9,10,11.5,12.5,13.5,15,16,17))
axis(side = 1,at = c(2,5.5,9,12.5,16),labels=c(0,1,2,3,4),lwd.ticks = T)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
pairwise.wilcox.test(c(q[q$q==0&q$stage==1,]$d,
                     q[q$q==0&q$stage==2,]$d,
                     q[q$q==0&q$stage==3,]$d,
                     q[q$q==1&q$stage==1,]$d,
                     q[q$q==1&q$stage==2,]$d,
                     q[q$q==1&q$stage==3,]$d,
                     q[q$q==2&q$stage==1,]$d,
                     q[q$q==2&q$stage==2,]$d,
                     q[q$q==2&q$stage==3,]$d,
                     q[q$q==3&q$stage==1,]$d,
                     q[q$q==3&q$stage==2,]$d,
                     q[q$q==3&q$stage==3,]$d,
                     q[q$q==4&q$stage==1,]$d,
                     q[q$q==4&q$stage==2,]$d,
                     q[q$q==4&q$stage==3,]$d),
                     paste(rep(c(rep(1,28),rep(2,39),rep(3,12)),4),rep(0:4,each=79)),p.adjust.method = 'BH')
lines(c(4.5,6.5),c(2150,2150));text(5.5,2300,'*')
lines(c(8,10),c(2000,2000));text(9,2150,'*')
lines(c(11.5,13.5),c(1800,1800));text(12.5,1950,'*')
lines(c(11.5,12.5),c(1500,1500));text(12,1650,'*')
lines(c(15,16),c(1500,1500));text(15.5,1650,'*')
lines(c(15,17),c(1800,1800));text(16,1950,'*')


par(mar=c(5,4,1,1))

p=v
boxplot(p,xaxt='n',col=rep(c(3,'darksalmon',2),50),ylab='Frequency',xlab='',las=1
        ,outline=F,at=sort(c(1+c(0:49)*3.5,2+c(0:49)*3.5,c(0:49)*3.5+3)))
axis(side = 1,at = 3.5*c(0:49)+2,labels=gsub('IGH','',names(p)[3*1:50-1]),lwd.ticks = T,las=2)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
r=pv(p);dy=0.135;ry=0.02
for(i in 1:nrow(r)){
  dd=(min((1:length(p))[names(p)==r$gene[i]])-1)/3;ddy=(r$a[i]+r$b[i]-2)*ry
  lines(c(3.5*dd+r$a[i],3.5*dd+r$b[i]),c(dy+ddy,dy+ddy))
  text((3.5*dd+r$a[i]+3.5*dd+r$b[i])/2,dy+ddy+ry/2,r$sig[i])
}
fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
fs2=(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
fs2=fs2[!grepl('/P',fs2)&!grepl('AP',fs2)&!grepl('CT',fs2)]
f=c(fs,fs2)
x=data.frame(rbindlist(mclapply(f,function(i){
  rep=read.delim(i)
  return(data.frame(IGD=sum(grepl('IGD',rep$c_call))/nrow(rep),
                    IGM=sum(grepl('IGM',rep$c_call))/nrow(rep),
                    IGG=sum(grepl('IGG',rep$c_call))/nrow(rep),
                    IGA=sum(grepl('IGA',rep$c_call))/nrow(rep)                    ))
},mc.cores=45)))
x$stage=ifelse(grepl('100',f),'control',ifelse(grepl('seve',f),'severe','mild/moderate'))
x$stage=factor(x$stage)
iso=x


x=iso
y=as.integer(x$stage)

z=rbind(read.csv('/work/smodi/OurCovid/28_10_2021_long/allMutationLoadConst4.csv'),
        read.csv('/work/smodi/OurCovid/28_10_2021_long/more/allMutationLoadConst4.csv'))
z=z[!grepl('/P',z$X)&!grepl('more//Cov_5H',z$X)&!grepl('more//Cov_6H',z$X)
    &!grepl('AP',z$X)&!grepl('CT',z$X),]
fs=list.files('/work/smodi/OurCovid/28_10_2021_long/',pattern='\\.tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
fs2=list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='\\.tsv$',full.names = T)
fs2=fs2[!grepl('/P',fs2)&!grepl('AP',fs2)&!grepl('CT',fs2)]
f=c(fs,fs2)
z$stage=ifelse(grepl('100',f),'control',ifelse(grepl('sev',f),'severe','mild/moderate'))


boxplot(x$IGD[y==1],x$IGD[y==2],x$IGD[y==3],
        x$IGM[y==1],x$IGM[y==2],x$IGM[y==3],
        x$IGG[y==1],x$IGG[y==2],x$IGG[y==3],
        x$IGA[y==1],x$IGA[y==2],x$IGA[y==3],
        col=rep(c(3,'darksalmon',2),4),xaxt='n',
        outline=F,ylab='Frequency',xlab='Isotype',las=1)
ig=c('IGD','IGM','IGG','IGA')

value=c(x$IGD[y==1],x$IGD[y==2],x$IGD[y==3],
  x$IGM[y==1],x$IGM[y==2],x$IGM[y==3],
  x$IGG[y==1],x$IGG[y==2],x$IGG[y==3],
  x$IGA[y==1],x$IGA[y==2],x$IGA[y==3])
group=paste(rep(c('IGD','IGM','IGG','IGA'),each=length(y)),rep(c(rep(1,sum(y==1)),rep(2,sum(y==2)),rep(3,sum(y==3))),4)
)
q=pairwise.wilcox.test(value,group,method='BH')
q$p.value[is.na(q$p.value)]=1
q=data.frame(q$p.value);r=NULL
for(i in c('IGD','IGM','IGG','IGA')){
  z=q[startsWith(rownames(q),i),startsWith(colnames(q),i)]
  for(a in rownames(z))for(b in colnames(z))if(z[a,b]<=0.05)r=rbind(r,data.frame(
    a=a,b=b,pv=z[a,b] ,j=j   ))
}
for(i in 1:nrow(r)){
  if(substr(r$a[i],1,3)=='IGD'){d=0;yy=0.2}
  if(substr(r$a[i],1,3)=='IGM')d=3
  if(substr(r$a[i],1,3)=='IGG'){d=6;yy=0.6}
  if(substr(r$a[i],1,3)=='IGA')d=9
  print(c(d+as.integer(substr(r$a[i],nchar(r$a[i]),100)),
          d+as.integer(substr(r$b[i],nchar(r$b[i]),100))))
  lines(c(d+as.integer(substr(r$a[i],nchar(r$a[i]),100)),
          d+as.integer(substr(r$b[i],nchar(r$b[i]),100))),
        rep(yy,2))
  text((d+as.integer(substr(r$a[i],nchar(r$a[i]),100))+
         d+as.integer(substr(r$b[i],nchar(r$b[i]),100)))/2,yy+0.05,ifelse(r$pv[i]<0.01,'**','*'))
}
axis(side = 1,at = c(2,5,8,11),labels=c("IGD","IGM","IGG","IGA"),lwd.ticks = T)
#legend('topleft',legend = c(
#  'control','mild','severe'),
#  fill=c(3,'darksalmon',2),box.col = NA)
mtext('E', side = 3, line = 0.5, adj = 0, cex = 1.1)
z=rbind(read.csv('/work/smodi/OurCovid/28_10_2021_long/allMutationLoadConst4.csv'),
        read.csv('/work/smodi/OurCovid/28_10_2021_long/more/allMutationLoadConst4.csv'))
z=z[!grepl('/P',z$X)&!grepl('more//Cov_5H',z$X)&!grepl('more//Cov_6H',z$X)
    &!grepl('AP',z$X)&!grepl('CT',z$X),]
fs=list.files('/work/smodi/OurCovid/28_10_2021_long/',pattern='\\.tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
fs2=list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='\\.tsv$',full.names = T)
fs2=fs2[!grepl('/P',fs2)&!grepl('AP',fs2)&!grepl('CT',fs2)]
f=c(fs,fs2)
z$stage=ifelse(grepl('100',f),'control',ifelse(grepl('sev',f),'severe','mild/moderate'))

boxplot(z[z$stage=='control','s_IGD'],
        z [z$stage=='mild/moderate','s_IGD'],z[z$stage=='severe', 's_IGD'],
        z[z$stage=='control','s_IGM'],
        z [z$stage=='mild/moderate','s_IGM'],z[z$stage=='severe', 's_IGM'],
        z[z$stage=='control','s_IGG'],
        z [z$stage=='mild/moderate','s_IGG'],z[z$stage=='severe', 's_IGG'],
        z[z$stage=='control','s_IGA'],
        z [z$stage=='mild/moderate','s_IGA'],z[z$stage=='severe', 's_IGA'],col=rep(
          c(3,'darksalmon',2),4),xaxt='n',outline=F,
        ylab='Silent V gene distance from germline',xlab='Isotype',ylim=c(0,0.05),las=1)
axis(side = 1,at = c(2,5,8,11),labels=c("IGD","IGM","IGG","IGA"),lwd.ticks = T)


value=c(z[z$stage=='control','s_IGD'],
        z [z$stage=='mild/moderate','s_IGD'],z[z$stage=='severe', 's_IGD'],
        z[z$stage=='control','s_IGM'],
        z [z$stage=='mild/moderate','s_IGM'],z[z$stage=='severe', 's_IGM'],
        z[z$stage=='control','s_IGG'],
        z [z$stage=='mild/moderate','s_IGG'],z[z$stage=='severe', 's_IGG'],
        z[z$stage=='control','s_IGA'],
        z [z$stage=='mild/moderate','s_IGA'],z[z$stage=='severe', 's_IGA'])
group=paste(rep(c('IGD','IGM','IGG','IGA'),each=length(y)),rep(c(rep(1,sum(y==1)),rep(2,sum(y==2)),rep(3,sum(y==3))),4)
)
q=pairwise.wilcox.test(value,group,method='BH')
q$p.value[is.na(q$p.value)]=1
q=data.frame(q$p.value);r=NULL
for(i in c('IGD','IGM','IGG','IGA')){
  z=q[startsWith(rownames(q),i),startsWith(colnames(q),i)]
  for(a in rownames(z))for(b in colnames(z))if(z[a,b]<=0.05)r=rbind(r,data.frame(
    a=a,b=b,pv=z[a,b] ,j=j   ))
}
for(i in 1:nrow(r)){
  if(substr(r$a[i],1,3)=='IGD'){d=0;yy=0.02}
  if(substr(r$a[i],1,3)=='IGM')d=3
  if(substr(r$a[i],1,3)=='IGG'){d=6;yy=0.04}
  if(substr(r$a[i],1,3)=='IGA')d=9
  print(c(d+as.integer(substr(r$a[i],nchar(r$a[i]),100)),
          d+as.integer(substr(r$b[i],nchar(r$b[i]),100))))
  lines(c(d+as.integer(substr(r$a[i],nchar(r$a[i]),100)),
          d+as.integer(substr(r$b[i],nchar(r$b[i]),100))),
        rep(yy,2))
  text((d+as.integer(substr(r$a[i],nchar(r$a[i]),100))+
         d+as.integer(substr(r$b[i],nchar(r$b[i]),100)))/2,yy+0.005,'***')
}

mtext('F', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()

######################################################################################
source('/work/smodi/scripts/covid/func.R')
fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
#fs=c(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
target=data.frame(rbindlist(mclapply(fs,function(i){
  x=read.delim(i)
  z=AAtarget(x)
  
  return(z)
},mc.cores=50  ),fill = T))
target[is.na(target)]=0
target$stage=as.factor(ifelse(grepl('100',fs)|grepl('/P',fs),'control','case'))
z=target
p2=traintest(z,pv=0.001,iterations = 100);boxplot(p2)
r=train(stage~.,z,method='glmnet')
z3=cov;z3$sequence_alignment_aa=cov$VH.or.VHH
for(i in 1:103){
  for(a in .AA)
    z3[,paste0(a,i)]=ifelse(substr(z3$sequence_alignment_aa,i,i)==a,1,0)
}

fs=c(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
target2=data.frame(rbindlist(mclapply(fs,function(i){
  x=read.delim(i)
  z=AAtarget(x)
  
  return(z)
},mc.cores=50  ),fill = T))
target2[is.na(target2)]=0
target2$stage=as.factor(ifelse(grepl('100',fs)|grepl('/P',fs),'control','case'))
p=predict(r,target2)


p=predict(r,z3);sum(p=='case')/length(p)
x=extractFeatures(z,0.001)
x$a=factor(substr(x$gene,1,1))
x$p=as.integer(substr(x$gene,2,5))
cov=read.csv('~/covAbdab310122.csv',encoding='UTF8')
m=read.csv('~/cov.table.aligned.csv',header=F)
cov$V1=1:nrow(cov)
mm=merge(m,cov,by='V1')
cov=mm
cov$VH.or.VHH=cov$V2
for(i in 1:nrow(x)){
  a=as.character(x$a[i]);p=x$p[i]
  cov[paste0(a,p)]=ifelse(substr(cov$VH.or.VHH,p,p)==a,1,0)
}

s=rep(0,nrow(cov))
for(i in 1:nrow(x))s=s+x$s1[i]*cov[,x$gene[i]]
summary(s)
bneu=c();for(i in c(0,-10000,-20000,-30000,-40000,-50000,-75000,-100000,-150000,-200000,-300000)){
  a=data.frame(table(cov[s<(i),]$Neutralising.Vs));a$Freq=a$Freq/sum(a$Freq)
  bneu=c(bneu,sum(a[grepl('SARS-CoV2',a$Var1),]$Freq)/sum(a$Freq))
};bneu
fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)&!grepl('AP',fs)&!grepl('CT',fs)]
m=mclapply(fs,function(i){
  k=i
  rep=read.delim(i)
  rep$seq=gsub('\\-','N',gsub('\\.','N',rep$sequence_alignment))
  rep$seq=as.character(translate(DNAStringSet(rep$seq),if.fuzzy.codon = 'X'))
  for(i in 1:nrow(x)){
    ;a=as.character(x$a[i]);p=x$p[i]
    rep[paste0(a,p)]=ifelse(substr(rep$seq,p,p)==a,1,0)
  }
  s=rep(0,nrow(rep))
  for(i in 1:nrow(x))s=s+x$s1[i]*rep[,x$gene[i]]
  return(s)
},mc.cores=50)
w=c();for(i in m)w=c(w,i)
k=c();for(i in m)k=c(k,sum(i<0)/length(i))
sum(s<0)/length(s)
for(i in m)print(sum(i<0)/length(i));
boxplot(k[1:28],k[29:60])
m=mclapply(fs[1:28],function(i){
  k=i
  rep=read.delim(i)
  rep$seq=gsub('\\-','N',gsub('\\.','N',rep$sequence_alignment))
  rep$seq=as.character(translate(DNAStringSet(rep$seq),if.fuzzy.codon = 'X'))
  for(i in 1:nrow(x)){
    v=as.character(x$v[i]);a=as.character(x$a[i]);p=x$p[i]
    rep[paste0(a,p)]=ifelse(substr(rep$seq,p,p)==a,1,0)
  }
  s=rep(0,nrow(rep))
  for(i in 1:nrow(x))s=s+x$s1[i]*rep[,x$gene[i]]
  h2=hist(s,main=k,breaks=1000,plot=F)
  h2$counts=cumsum(h2$counts)
  return(data.frame(breaks=h2$breaks[2:length(h2$breaks)],counts=h2$counts))
},mc.cores=50)

plot((-h$breaks[2:length(h$breaks)]),h$counts/max(h$counts),pch='.',
     xlab='Ab anti covid19 score',ylab='Cumulative fraction',las=1)
j=1
for(i in m){
  lines((-i$breaks),i$counts/max(i$counts),col=ifelse(grepl('Spike',fs[j]),2,4))
  j=j+1
}
lines((-h$breaks[2:length(h$breaks)]),(h$counts/max(h$counts)),col=2)
legend('topright',legend = c('CoV-abdab Abs','Repertoires'),
       fill=c(2,3),box.col = NA)

######################################################################################
#         targeting AAA
source('/work/smodi/scripts/covid/func.R')
fs=list.files('/work/jenkins/align_n_annotate/130/',pattern='tsv$',full.names = T)
#
#target=data.frame(rbindlist(
q=mclapply(fs,function(i){
  x=read.delim(i)
  z=AAtarget(x[startsWith(x$v_call,'IGHV1'),])
  colnames(z)=paste0('IGHV1',colnames(z))
  for(i in paste0('IGHV',c(2:5))){
    y=AAtarget(x[startsWith(x$v_call,i),])
    colnames(y)=paste0(i,colnames(y))
    z=merge(z,y)
  }
  return(z)
},mc.cores=25  )


#,fill = T))
target[is.na(target)]=0
target$stage=0
z0=target
fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
#fs=c(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
target=data.frame(rbindlist(mclapply(fs,function(i){
  x=read.delim(i)
  z=AAtarget(x[startsWith(x$v_call,'IGHV1'),])
  colnames(z)=paste0('IGHV1',colnames(z))
  for(i in paste0('IGHV',c(2:5))){
    y=AAtarget(x[startsWith(x$v_call,i),])
    colnames(y)=paste0(i,colnames(y))
    z=merge(z,y)
  }
  return(z)
},mc.cores=50  ),fill = T))
target[is.na(target)]=0
target$stage=as.factor(ifelse(grepl('100',fs)|grepl('/P',fs),'control','case'))
z=target




y=rbind(z,z0)



#####################################################################


source('/work/smodi/scripts/covid/func.R')

fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
#fs=c(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
target=data.frame(rbindlist(mclapply(fs,function(i){
  x=read.delim(i)
  z=AAtarget(x[startsWith(x$v_call,'IGHV1'),])
  colnames(z)=paste0('IGHV1',colnames(z))
  for(i in paste0('IGHV',c(2:5))){
    y=AAtarget(x[startsWith(x$v_call,i),])
    colnames(y)=paste0(i,colnames(y))
    z=merge(z,y)
  }
  return(z)
},mc.cores=50  ),fill = T))
target[is.na(target)]=0
target$stage=as.factor(ifelse(grepl('100',fs)|grepl('/P',fs),'control','case'))
z=target
p2=traintest(z,pv=0.001,iterations=50);boxplot(p2)
x=extractFeatures(z,0.001)
r=train(stage~.,z[,colnames(z)%in%c('stage',x$gene)],method='glmnet',tuneLength=5)
fs=c(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
target=data.frame(rbindlist(mclapply(fs,function(i){
  x=read.delim(i)
  z=AAtarget(x[startsWith(x$v_call,'IGHV1'),])
  colnames(z)=paste0('IGHV1',colnames(z))
  for(i in paste0('IGHV',c(2:5))){
    y=AAtarget(x[startsWith(x$v_call,i),])
    colnames(y)=paste0(i,colnames(y))
    z=merge(z,y)
  }
  return(z)
},mc.cores=50  ),fill = T))
target[is.na(target)]=0
target$stage=as.factor(ifelse(grepl('100',fs)|grepl('/P',fs),'control','case'))



x$v=factor(substr(x$gene,1,5))
x$a=factor(substr(x$gene,6,6))
x$p=as.integer(substr(x$gene,7,9))
library(openxlsx)
a=openxlsx::read.xlsx('/work/smodi/scripts/covid/2022_08_03 SARS-CoV-2 Abs table for Gur.xlsx',sheet=2)
cov=a
cov$VH.or.VHH=cov$heavy_chain_seq_imgt
cov$Heavy.V.Gene=cov$v_gene_heavy
cov=read.csv('~/Downloads/CoV-AbDab_011121.csv',encoding='UTF8')
cov=read.csv('~/Downloads/CoV-AbDab_',encoding='UTF8')

#cov=read.csv('~/covmeital.csv')
cov=read.csv('~/covAbdab310122.csv',encoding='UTF8')
m=read.csv('~/cov.table.aligned.csv',header=F)
cov$V1=1:nrow(cov)
mm=merge(m,cov,by='V1')
cov=mm
#cov=cov[grepl('uman',cov$Origin),]
#cov=cov[!grepl('SARS-CoV2',cov$Doesn.t.Bind.to),]
cov$VH.or.VHH=cov$V2
for(i in 1:nrow(x)){
  v=as.character(x$v[i]);a=as.character(x$a[i]);p=x$p[i]
  cov[paste0(v,a,p)]=ifelse(startsWith(cov$Heavy.V.Gene,v)&
                              substr(cov$VH.or.VHH,p,p)==a,1,0)
}

s=rep(0,nrow(cov))
for(i in 1:nrow(x))s=s+x$s1[i]*cov[,x$gene[i]]
h=hist(s,breaks=10*-3000:1000)
h$counts=cumsum(h$counts)
covn=cov[grepl('SARS-CoV2',cov$Neutralising.Vs),]
sn=s[grepl('SARS-CoV2',cov$Neutralising.Vs)]
hn=hist(sn,breaks=10*-3000:1000)
hn$counts=cumsum(hn$counts)

plot((-h$breaks[2:length(h$breaks)]),(h$counts/max(h$counts)),pch='.')
lines((-h$breaks[2:length(h$breaks)]),(h$counts/max(h$counts)))
lines((-hn$breaks[2:length(hn$breaks)]),(hn$counts/max(hn$counts)),col='red')


plot((-h$breaks[2:length(h$breaks)]),(h$counts),pch='.')
lines((-h$breaks[2:length(h$breaks)]),(h$counts))
lines((-hn$breaks[2:length(hn$breaks)]),(hn$counts),col='red')
plot((-hn$breaks[2:length(hn$breaks)]),(h$counts/hn$counts))
lines((-hn$breaks[2:length(hn$breaks)]),(h$counts/hn$counts),col='blue')

neu=c()
for(i in c(100,50,0,-50,-100,-200,-300,-400,-1000)){
  a=data.frame(table(cov[s<(i),]$Neutralising.Vs));a$Freq=a$Freq/sum(a$Freq)
  neu=c(neu,sum(a[grepl('SARS-CoV2',a$Var1),]$Freq)/sum(a$Freq))
};neu
dev.off()



bneu=c()
for(i in c(0,-10000,-20000,-30000,-40000,-50000)){
  a=data.frame(table(cov[s<(i),]$Neutralising.Vs));a$Freq=a$Freq/sum(a$Freq)
  bneu=c(bneu,sum(a[grepl('SARS-CoV2',a$Var1),]$Freq)/sum(a$Freq))
}

nneu=c()
for(i in c(0,-1500,-2500,-5000,-20000)){
  a=data.frame(table(cov[s<(i),]$Neutralising.Vs));a$Freq=a$Freq/sum(a$Freq)
  nneu=c(nneu,sum(a[!grepl('SARS-CoV2',a$Var1),]$Freq)/sum(a$Freq))
}
fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)&!grepl('AP',fs)&!grepl('CT',fs)]
m=mclapply(fs,function(i){
  k=i
  rep=read.delim(i)
  
  rep$seq=gsub('\\-','N',gsub('\\.','N',rep$sequence_alignment))
  rep$seq=as.character(translate(DNAStringSet(rep$seq),if.fuzzy.codon = 'X'))
  #rep$seq=gsub('X','',rep$seq)
  for(i in 1:nrow(x)){
    v=as.character(x$v[i]);a=as.character(x$a[i]);p=x$p[i]
    rep[paste0(v,a,p)]=ifelse(startsWith(rep$v_call,v)&
                                substr(rep$seq,p,p)==a,1,0)
  }
  s=rep(0,nrow(rep))
  for(i in 1:nrow(x))s=s+x$s1[i]*rep[,x$gene[i]]
  return(s)
},mc.cores=50)
w=c();for(i in m)w=c(w,i)
k=c();for(i in m)k=c(k,sum(i<0)/length(i))
sum(s<0)/length(s)
for(i in m)print(sum(i<0)/length(i));
boxplot(k[1:28],k[29:60])
pdf('/work/smodi/OurCovid/pdfWOziv/NEWplot2.additional.pdf',width=7.5,height=5)
par(mfrow=c(1,2),mgp=c(3,0.5,0))
boxplot(-w,-s,outline=F,names=c('All repertoires','CoV Abdab'),las=1,
        ylab='Score',ylim=c(0,2000))
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
barplot(neu,names.arg = c(0,5000,10000,15000),ylim=c(0,1),
        xlab='Score',ylab='Fraction neutralizing',las=1,axis.lty=1);box()
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()
m=mclapply(fs[1:28],function(i){
  k=i
  rep=read.delim(i)
  rep$seq=gsub('\\-','N',gsub('\\.','N',rep$sequence_alignment))
  rep$seq=as.character(translate(DNAStringSet(rep$seq),if.fuzzy.codon = 'X'))
  for(i in 1:nrow(x)){
    v=as.character(x$v[i]);a=as.character(x$a[i]);p=x$p[i]
    rep[paste0(v,a,p)]=ifelse(startsWith(rep$v_call,v)&
                                substr(rep$seq,p,p)==a,1,0)
  }
  s=rep(0,nrow(rep))
  for(i in 1:nrow(x))s=s+x$s1[i]*rep[,x$gene[i]]
  h2=hist(s,main=k,breaks=1000,plot=F)
  h2$counts=cumsum(h2$counts)
  return(data.frame(breaks=h2$breaks[2:length(h2$breaks)],counts=h2$counts))
},mc.cores=50)

plot((-h$breaks[2:length(h$breaks)]),h$counts/max(h$counts),pch='.',
     xlab='Ab anti covid19 score',ylab='Cumulative fraction',las=1)
j=1
for(i in m){
  lines((-i$breaks),i$counts/max(i$counts),col=ifelse(grepl('Spike',fs[j]),2,4))
  j=j+1
}
lines((-h$breaks[2:length(h$breaks)]),(h$counts/max(h$counts)),col=2)
legend('topright',legend = c('CoV-abdab Abs','Repertoires'),
       fill=c(2,3),box.col = NA)

p=predict(r,target[1:31,])

####################################################################################
#	plot 2

pdf('/work/smodi/OurCovid/pdfWOziv/NEWplot2.pdf',width=7.5,height=5.5)
par(cex=0.7)
layout(rbind(c(1,2,3),c(4,4,4)))
boxplot(p2,ylim=c(0,1),las=2)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
p=p[1:31];target=target[1:31,]
barplot(c(F1_Score(target$stage,p,positive='case'),
          sum(p==target$stage)/nrow(target),
          Sensitivity(target$stage,p,positive='case'),
          Specificity(target$stage,p,positive='case')),
        names.arg=c('F1 score','accuracy','sensitivity','specificity'),ylim=c(0,1),axis.lty=1,
        las=2)
box()
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
plot((-h$breaks[2:length(h$breaks)]),h$counts/max(h$counts),pch='.',
     xlab='Ab anti covid19 score',ylab='Cumulative fraction',las=1)
for(i in m)
  lines((-i$breaks),i$counts/max(i$counts),col=3)
lines((-h$breaks[2:length(h$breaks)]),(h$counts/max(h$counts)),col=2)
legend('topright',legend = c('CoV-Abdab Abs','repertoires'),
       fill=c(2,3),box.col = NA)
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)

x$pos=log10(abs(x$s1))*(-1)*x$s1/abs(x$s1)
x=x[order(x$pos,decreasing = T),]
barplot(x$pos,names.arg = paste0(substr(x$gene,4,5),' ',x$a,x$p),axis.lty=1,
        col=ifelse(x$pos>0,2,3),las=2,ylim=c(-5.5,5.5),ylab='Log10 coefficient')
box()
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()



pdf('/work/smodi/OurCovid/pdfWOziv/NEWplot2correct.pdf',width=7.5,height=5.5)
par(cex=0.7)
layout(rbind(c(1,1,2,2,3,3,4,4),rep(5,8)))
boxplot(p2,ylim=c(0,1),las=2)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
target=target[1:31,]
barplot(c(F1_Score(target$stage,p,positive='case'),
          sum(p==target$stage)/nrow(target),
          Sensitivity(target$stage,p,positive='case'),
          Specificity(target$stage,p,positive='case')),
        names.arg=c('F1 score','accuracy','sensitivity','specificity'),ylim=c(0,1),axis.lty=1,
        las=2)
box()
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(k[1:28],k[29:60],sum(s<0)/length(s),outline=F,names=c('control',
      'covid-19','CoV-Abdab'),las=2,
        ylab='Fraction predicted as SARS-CoV2 Abs')
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)
plot((-h$breaks[2:length(h$breaks)]),(h$counts),pch='.',xlim=c(-2500,5000),
     ylab='Cumulative count',xlab='Score')
lines((-h$breaks[2:length(h$breaks)]),(h$counts))
lines((-hn$breaks[2:length(hn$breaks)]),(hn$counts),col='red')
legend('topright',legend=c('All','SARS-CoV2\nneutralizing'),fill=c(1,2),bty='n',cex = 0.7)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
box()
x$pos=log10(abs(x$s1))*(-1)*x$s1/abs(x$s1)
x=x[order(x$pos,decreasing = T),]
barplot(x$pos,names.arg = paste0(substr(x$gene,4,5),' ',x$a,x$p),axis.lty=1,cex.names=0.65,
        col=ifelse(x$pos>0,2,3),las=2,ylim=c(-5.5,5.5),ylab='Log10 coefficient')
box()
mtext('E', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()


q=c();for(i in -100:100*100)q=c(q,sum(grepl('SARS-CoV2',cov[s>i,]$Neutralising.Vs))/sum(s>i))

dev.off();plot(-100:100*100,q,pch='.');lines(-100:100*100,q)






plot(x$p,x$pos,pch='')
text(x$p,x$pos,labels=x$a,col=ifelse(x$pos>0,2,3))
library(seqLogo)
#write.csv(x,'~/coefOfAAfreq.csv')
x=read.csv('~/coefOfAAfreq.csv')
k=list()
for(i in 1:103){
  y=data.frame(aa=x[x$p==i,]$a,val=x[x$p==i,]$s1)
  rownames(y)=y$aa;y$aa=NULL
  if(nrow(y)==0){k[[i]]=data.frame(X=0)
  }else   k[[i]]=data.frame(t(y))
}
k=data.table::rbindlist(k,fill=T)
k=data.frame(k)
k[is.na(k)]=0


pdf('/work/smodi/OurCovid/pdfWOziv/NEWplot2new.pdf',width=7.5,height=8)
par(cex=0.7)
layout(rbind(c(1,1,2,2,3),c(4,4,4,5,5),c(6,6,6,6,6)))
boxplot(p2,ylim=c(0,1),las=2)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
target=target[1:31,]
barplot(c(F1_Score(target$stage,p,positive='case'),
          sum(p==target$stage)/nrow(target),
          Sensitivity(target$stage,p,positive='case'),
          Specificity(target$stage,p,positive='case')),
        names.arg=c('F1 score','accuracy','sensitivity','specificity'),ylim=c(0,1),axis.lty=1,
        las=2)
box()
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
boxplot(-w,-s,outline=F,names=c('All repertoires','CoV Abdab'),las=2,
        ylab='Score',ylim=c(0,2000))
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)

plot((-h$breaks[2:length(h$breaks)]),h$counts/max(h$counts),pch='.',
     xlab='Ab anti covid19 score',ylab='Cumulative fraction',las=1)
for(i in m)
  lines((-i$breaks),i$counts/max(i$counts),col=3)
lines((-h$breaks[2:length(h$breaks)]),(h$counts/max(h$counts)),col=2)
legend('topright',legend = c('CoV-Abdab Abs','repertoires'),
       fill=c(2,3),box.col = NA)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
barplot(neu,names.arg = c(0,25000,50000,100000),ylim=c(0,1),
        xlab='Score',ylab='Fraction neutralizing',las=1,axis.lty=1);box()
mtext('E', side = 3, line = 0.5, adj = 0, cex = 1.1)

x$pos=log10(abs(x$s1))*(-1)*x$s1/abs(x$s1)
x=x[order(x$pos,decreasing = T),]
barplot(x$pos,names.arg = paste0(substr(x$gene,4,5),' ',x$a,x$p),axis.lty=1,
        col=ifelse(x$pos>0,2,3),las=2,ylim=c(-5.5,5.5),ylab='Log10 coefficient')
box()
mtext('F', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()

x=x3
x$v=factor(substr(x$gene,1,5))
x$a=factor(substr(x$gene,6,6))
x$p=as.integer(substr(x$gene,7,9))
x$log=log10(abs(x$s1))
x$sign=ifelse(x$s1<0,1,-1)
div='';for(m in sort(unique(x$p))){
  y=x[x$p==m,]
  jplus=200;jminus=200;
  for(j in 1:nrow(y)){
    div=paste0(div,'<p style="',
             'position:absolute;left:',15*m,'px;top:',
             ifelse(y$sign[j]==1,(jplus-15*y$log[j]),jminus),'px;font-size:15px;',
             'transform: scale(1, ',y$log[j],');">',y$a[j],'</p>' )
    if(y$sign[j]==1){
      jplus=jplus-15*y$log[j]
    }else jminus=jminus+15*y$log[j]
  }
};div

  
#####################################################################################
x=openxlsx::read.xlsx('/home/bcrlab/smodi/nicwulab-SARS-CoV-2_Abs-05c663b/code/CoV_Encoder/data/RBD-HA_new_full_CDR.xlsx')
x=x[!is.na(x$VH_nuc),]
write.table(paste0('>',1:nrow(x),'\n',x$VH_nuc),'~/nuc/covnuc.fasta',quote=F,row.names=F,col.names=F)
rep=read.delim('~/nuc/covnuc.fasta.tab')
