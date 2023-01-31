source('/work/smodi/scripts/covid/func.R')
f=list.files('/work/smodi/covid/igblast',pattern='tab$',full.names=T)
s=read.csv('/work/smodi/covid/igblast/all/Plasmablast.csv')
library(msa)
zz=mclapply(1:14,function(i){
  x=read.delim(f[i]) 
  #y=s[s$id==i,]
  y=s
  x=x[startsWith(x$v_call,'IGH')&x$sequence_id%in%y$cell_id,]
  if(nrow(x)==0)return(NULL)
  return(data.frame(v=x$v_call,j=x$j_call,seq=x$v_sequence_alignment,germ=x$v_germline_alignment,
                    all=x$sequence_alignment_aa,cdr3=x$cdr3_aa,
                    dup=ifelse(is.null(x$duplicate_count),rep(1,nrow(x))
                               ,x$duplicate_count),id=rep(i,nrow(x))))
},mc.cores=50)
l=data.frame(rbindlist(zz))
l$v=alakazam::getGene(l$v);l$j=alakazam::getGene(l$j);
l$p=paste(l$v,l$j,l$cdr3)
plasmablast=l
zz=mclapply(1:14,function(i){
  x=read.delim(f[i]) 
  #y=s[s$id==i,]
  y=s
  x=x[startsWith(x$v_call,'IGH')&!x$sequence_id%in%y$cell_id,]
  if(nrow(x)==0)return(NULL)
  return(data.frame(v=x$v_call,j=x$j_call,seq=x$v_sequence_alignment,germ=x$v_germline_alignment,
                    all=x$sequence_alignment_aa,cdr3=x$cdr3_aa,
                    dup=ifelse(is.null(x$duplicate_count),rep(1,nrow(x))
                               ,x$duplicate_count),id=rep(i,nrow(x))))
},mc.cores=50)
l=data.frame(rbindlist(zz))
l$v=alakazam::getGene(l$v);l$j=alakazam::getGene(l$j);
l$p=paste(l$v,l$j,l$cdr3)
all=l

library(shazam)
modelPlasma <- createTargetingModel(plasmablast,sequenceColumn = 'seq',germlineColumn = 'germ',
                              model='rs', vCallColumn="v")
modelAll <- createTargetingModel(all,sequenceColumn = 'seq',germlineColumn = 'germ',
                                    model='rs', vCallColumn="v")
source('/work/smodi/scripts/covid/func.R')
p=data.frame(t(data.frame(modelPlasma@mutability)))
a=data.frame(t(data.frame(modelAll@mutability)))
#write.csv(a,'~/allSCmutability.csv')
#write.csv(p,'~/plasmablastmutability.csv')



a=read.csv('~/allSCmutability.csv')
p=read.csv('~/plasmablastmutability.csv')
a$X=NULL;p$x=NULL;a[is.na(a)]=0;p[is.na(p)]=0
b=boxplot(rowMeans(a[,colnames(a)%in%WRC_GYW]),
        rowMeans(p[,colnames(p)%in%WRC_GYW]),
        rowMeans(a[,colnames(a)%in%WA_TW]),
        rowMeans(p[,colnames(p)%in%WA_TW]),
        rowMeans(a[,colnames(a)%in%SYC_GRS]),
        rowMeans(p[,colnames(p)%in%SYC_GRS]),
        rowMeans(a[,colnames(a)%in%REST]),
        rowMeans(p[,colnames(p)%in%REST]),
        col=rep(c(3,2),4),xaxt='n',outline=F,ylab='Mean mutability',xlab='Hotspot',
        ylim=c(0,0.002),las=1,main='')

#write.csv(a,'~/allSCmutability.csv')
#write.csv(p,'~/plasmablastmutability.csv')
#############################################################################################
loadOur = function(stage='covid'){
  f=list.files('/work/smodi/OurCovid/28_10_2021_long/',pattern='tsv$',full.names=T)
  f2=list.files('/work/smodi/OurCovid/28_10_2021_long/more',pattern='tsv$',full.names=T)
  if(stage=='ourPlasma'){
    f=list.files('/work/smodi/covid/igblast',pattern='tab$',full.names=T)
    s=read.csv('/work/smodi/covid/igblast/all/Plasmablast.csv')
    library(msa)
    zz=mclapply(1:14,function(i){
      x=read.delim(f[i]) 
      #y=s[s$id==i,]
      y=s
      x=x[startsWith(x$v_call,'IGH')&x$sequence_id%in%y$cell_id,]
      if(nrow(x)==0)return(NULL)
      return(data.frame(v=x$v_call,j=x$j_call,all=x$sequence_alignment_aa,cdr3=x$cdr3_aa,
                        dup=ifelse(is.null(x$duplicate_count),rep(1,nrow(x))
                                   ,x$duplicate_count),id=rep(i,nrow(x))))
    },mc.cores=50)
    l=data.frame(rbindlist(zz))
    l$v=alakazam::getGene(l$v);l$j=alakazam::getGene(l$j);
    l$p=paste(l$v,l$j,l$cdr3)
    return(l)
    
  }
  if(stage=='covid'){
    f2=f2[grepl('Cov_',f2)]
    f=f[grepl('Co',f)&!grepl('LIGHT',f)]
    f=f[!grepl('100',f)]
  }else if(stage=='control'){
    f2=f2[!grepl('Cov_',f2)&!grepl('CT',f2)]
    f=f[grepl('100',f)&!grepl('LIGHT',f)]
  }else if(stage=='sc'){
    f=list.files('/work/smodi/covid/igblast',pattern='tab$',full.names=T)
    f2=c()
  }else if (stage=='natalie'){
    l=read.delim('~/Natalie/all.tsv')
    l=l[startsWith(l$v_call,'IGH'),]
    l=(data.frame(v=l$v_call,j=l$j_call,cdr3=substr(
      l$junction_aa,2,nchar(l$junction_aa)-1),
      dup=ifelse(is.null(l$duplicate_count),rep(1,nrow(l))
                 ,l$duplicate_count),id=l$patient_alias))
    l$v=alakazam::getGene(l$v);l$j=alakazam::getGene(l$j);
    l$p=paste(l$v,l$j,l$cdr3)
    return(l)
  }
  f=c(f,f2)
  library(parallel);library(data.table)
  l=data.frame(rbindlist(mclapply(f,function(i){
    x=read.delim(i)
    x=x[nchar(as.character(x$cdr3_aa))>5&x$productive,]
    return(data.frame(v=x$v_call,j=x$j_call,all=x$sequence_alignment_aa,cdr3=x$cdr3_aa,
                      dup=ifelse(is.null(x$duplicate_count),rep(1,nrow(x))
                                 ,x$duplicate_count),id=rep(i,nrow(x))))
  },mc.cores=50)))
  l$v=alakazam::getGene(l$v);l$j=alakazam::getGene(l$j);
  l$p=paste(l$v,l$j,l$cdr3)
  return(l)
}

clonesByVJCDR3=function(l,thresh=0.25,method='hamminig'){
  r=NULL
  for(i in 1:nrow(cov)){
    z=stringdist::stringdist(cov$CDRH3[i],l[l$v==cov$v[i]&l$j==cov$j[i]&
                                              nchar(l$cdr3)==nchar(cov$CDRH3[i]),]$cdr3,method='hamming')
    if(sum(z<=(thresh*nchar(cov$CDRH3[i])))>0){
      print(paste(length(z),sum(z<=(thresh*nchar(cov$CDRH3[i])))))
      x=l[l$v==cov$v[i]&l$j==cov$j[i]&
            nchar(l$cdr3)==nchar(cov$CDRH3[i]),]
      x=x[z<=(thresh*nchar(cov$CDRH3[i])),]
      r=rbind(r,data.frame(index=i,id=cov$p[i],dupsum=sum(x$dup),reps=length(unique(x$id)),
                           dup=paste0(x$dup,collapse=','),
                           dists=paste0(z[z<=(thresh*nchar(cov$CDRH3[i]))],   collapse=','),
                           seqids=paste0(rownames(x),collapse=','),
                           repids=paste0((x$id),collapse=',')))
    }
  }
  return(r)
}

convert=function(rep){
  index=unlist(lapply(1:nrow(rep),function(i){return(rep(rep$index[i],
                                                         length(unlist(strsplit(rep$dup[i],',')))))}))
  id=unlist(lapply(1:nrow(rep),function(i){return(rep(rep$id[i],
                                                      length(unlist(strsplit(rep$dup[i],',')))))}))
  
  y=data.frame(index=index,id=id,
               seqids=unlist(strsplit(paste0(rep$seqids,collapse=','),',')),
               dup=as.integer(unlist(strsplit(paste0(rep$dup,collapse=','),','))),
               repids=unlist(strsplit(paste0(rep$repids,collapse=','),',')),
               dist=unlist(strsplit(paste0(rep$dists,collapse=','),',')))
  return(y)
}

library(parallel);library(data.table)
cov=read.csv('/work/smodi/scripts/covidClones/CoV-AbDab_310122.csv')
cov$v=alakazam::getGene(cov$Heavy.V.Gene)
cov$j=alakazam::getGene(cov$Heavy.J.Gene)
cov$p=paste(cov$v,cov$j,cov$CDRH3)
cov=cov[!duplicated(cov$p),]

plas=loadOur('ourPlasma')
covid=loadOur()
control=loadOur('control')
sc=loadOur('sc')
rcovid=clonesByVJCDR3(covid)
rcon=clonesByVJCDR3(control)
rsc=clonesByVJCDR3(sc)
rplas=clonesByVJCDR3(plas)
con=convert(rcon)
cvd=convert(rcovid)
con$stage='control'
cvd$stage='covid'
all=rbind(con,cvd)
allmeta=rbind(plyr::count(control,'id','dup'),
              plyr::count(covid,'id','dup'))
#write.csv(all,'~/reverse.all.csv')
#write.csv(allmeta,'~/reverse.allmeta.csv')

csc$stage='single cell'
cpl$stage='plasma'
all=rbind(cpl,csc)
allmeta=rbind(plyr::count(sc,'id','dup'),
              plyr::count(plas,'id','dup'))
#write.csv(all,'~/reverseCS.all.csv')
#write.csv(allmeta,'~/reverseCS.allmeta.csv')

d=list.files('~/covid19/tcr/cellranger-4.0.0/modi/yyw/')
v=d[endsWith(d,'vdj')]
e=d[endsWith(d,'exp')]
e1=paste0('~/covid19/tcr/cellranger-4.0.0/modi/yyw/',e)
d=list.files('/work/yyw/weizmann/')
d=d[startsWith(d,'r2')]
v=d[endsWith(d,'vdj')]
e=d[endsWith(d,'exp')]
e2=paste0('/work/yyw/weizmann/',e,'/')
library(Seurat);library(DropletUtils);library(SingleR);library(dplyr)
i=2
h=MonacoImmuneData()

library(parallel)
x=mclapply(1:14,function(i){
  r1 <- CreateSeuratObject(Read10X(paste0(e1[i],'/outs/filtered_feature_bc_matrix//')),
                           min.cells = 3,
                           min.genes = 200,project=paste(i,1))
  r2 <- CreateSeuratObject(Read10X(paste0(e2[i],'/outs/filtered_feature_bc_matrix//')),
                           min.cells = 3,
                           min.genes = 200,project=paste(i,2))
  r=merge(x=r1,y=r2,project=as.character(i))
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = r), value = TRUE)
  percent.mito <- colSums(r[mito.genes, ]) /
    colSums(r)
  subdata <- AddMetaData(object = r, metadata = percent.mito, col.name =
                           "percent.mito")
  p <- SingleR(test = as.SingleCellExperiment(subdata), ref = h, assay.type.test=1,
               labels = h$label.fine)
  subdata <- AddMetaData(object = subdata, metadata = p$first.labels, col.name =
                           "group")
  
},mc.cores=14)
y=x[[1]]
for(i in 2:14){print(i);y=merge(y,x[[i]])}


subdata=y
subdata=subdata[,grepl('B',subdata$group)]
subdata=y[,y@meta.data[,'percent.mito']<0.05&
            y@meta.data[,'nCount_RNA']>=200&y@meta.data[,'nCount_RNA']<=2500&
            grepl('B',y$group)|y$group=='Plasmablasts']

subdata
subdata <- NormalizeData(object = subdata,
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4)
subdata <- FindVariableFeatures(object = subdata,
                                mean.function = ExpMean,
                                dispersion.function = LogVMR,
                                x.low.cutoff = 0.0125,
                                x.high.cutoff = 3,
                                y.cutoff = 0.5)
subdata <- ScaleData(object = subdata,
                     vars.to.regress = c("nUMI", "percent.mito"))
subdata <- RunPCA(object = subdata,
                  pc.genes = subdata@var.genes,
                  do.print = TRUE,
                  pcs.print = 1:5,
                  genes.print = 5)
p=RunUMAP(subdata,dims=1:10)
q=data.frame(x=p@reductions$umap@cell.embeddings[,1],
             y=p@reductions$umap@cell.embeddings[,2],group=p$group)
#write.csv(q,'~/umapnew.csv')
##################################################################################
##      plot last figure            ##############################################
##################################################################################
a=function(){
pdf('/work/smodi/OurCovid/pdfWOziv/NEWlastFig.pdf',width=7.5,height=6)
layout(rbind(c(1,1,1,2,2,3,3),c(4,4,4,4,4,5,5)))
par(mar=c(4.1,4,2,1))
par(mgp=c(3,0.8,0))
all=read.csv('~/reverse.all.csv')
allmeta=read.csv('~/reverse.allmeta.csv')
allmeta=allmeta[!grepl('P',allmeta$id)|grepl('P4',allmeta$id)|grepl('P9',allmeta$id),]
all$n=nchar(substr(gsub(' ','',substr(gsub('IGHJ','                                                '
                                           ,all$id),30,170)),2,60))
all=all[as.integer(all$dist)/all$n<0.15,]
x=merge(plyr::count(all,'repids','dup'),allmeta,by.x='repids',by.y='id',all=T)
x[is.na(x)]=0
x$freq=x$freq.x/x$freq.y
x=x[x$freq.y>2000,]
x$stage=!(grepl('/P',x$repids)|grepl('100',x$repids))
z=x
y=pROC::roc(z$stage~z$freq)
table(x$stage,x$freq>median(x$freq))
barplot(sort(x$freq)+0.0000007,col=ifelse(x[order(x$freq),]$stage,2,3),log='y',las=2
        ,ylim=c(5e-7,2e-2),ylab='Fraction close to known antibodies',xlab='sample');box()
legend('topleft',legend = c('control','covid'),fill=c(3,2),box.col = NA)
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
plot(y$specificities,y$sensitivities,xlim=c(1,-0.1),pch='.',xlab='Specificity',ylab='Sensitivity')
lines(y$specificities,y$sensitivities      )
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)

all=read.csv('~/reverseCS.all.csv')
allmeta=read.csv('~/reverseCS.allmeta.csv')
all$n=nchar(substr(gsub(' ','',substr(gsub('IGHJ','                                                '
                                           ,all$id),30,170)),2,60))
all=all[as.integer(all$dist)/all$n<0.15,]
x=merge(plyr::count(all,'repids','dup'),allmeta,by.x='repids',by.y='id',all=T)
x[is.na(x)]=0
x$freq=x$freq.x/x$freq.y
zz=x
sum(x$freq.x[15:28])/sum(x$freq.y[15:28])
x=x[x$freq.y>2000,]
sum(x$freq.x)/sum(x$freq.y)
zzz=x
x$stage=!(grepl('/P',x$repids)|grepl('100',x$repids))
table(x$stage,x$freq>median(x$freq))
barplot(sort(x$freq)+0.0000007,col=ifelse(x[order(x$freq),]$stage,2,3),log='y',las=2
        ,ylim=c(5e-7,2e-2),ylab='Fraction close to known antibodies',xlab='sample');box()
#abline(h=median(z$freq))
mtext('C', side = 3, line = 0.5, adj = 0, cex = 1.1)

p=read.csv('~/umapnew.csv')
par(mar=c(5.5,4,2,1))
par(cex=0.8)

p$gr=ifelse(grepl('B',p$group)|p$group=='Plasmablasts',p$group,'Non B cells')
plot(p$x,p$y,col=2+as.integer(as.factor(p$gr)),las=1,
     pch='+',xlab='UMAP 1',ylab='UMAP 2',ylim=c(-18,6))
legend('bottom',legend = gsub('memory','\nmemory',as.character(levels(as.factor(p$gr)))),
       fill=2+as.integer(as.factor(
  levels(as.factor(p$gr)))),box.col = NA,ncol=3,cex=0.8,)
mtext('D', side = 3, line = 0.5, adj = 0, cex = 1.1)
par(mgp=c(3,0.8,0))
barplot(c(sum(z[!z$stage,]$freq.x)/sum(z[!z$stage,]$freq.y),
          sum(z[z$stage,]$freq.x)/sum(z[z$stage,]$freq.y),
          sum(zzz$freq.x)/sum(zzz$freq.y),
          sum(zz[15:28,]$freq.x)/sum(zz[15:28,]$freq.y)),names.arg = c(
            'control','covid','   covid    \nsingle-cell','plasma'
          ),las=2,ylab='Fraction close to known antibodies',ylim=c(0,0.005),axis.lty=1)
box()
mtext('E', side = 3, line = 0.5, adj = 0, cex = 1.1)

dev.off()
};a()

pdf('/work/smodi/OurCovid/pdfWOziv/NEWlastFigSup.pdf',width=7.5,height=6)
layout(mat = matrix(c(1,1,2),nrow=1))
all=read.csv('~/reverse.all.csv')
allmeta=read.csv('~/reverse.allmeta.csv')
all$n=nchar(substr(gsub(' ','',substr(gsub('IGHJ','                                                '
                                         ,all$id),30,170)),2,60))
all=all[as.integer(all$dist)/all$n<0.15,]
all=all[!grepl('100',all$repids)&!grepl('/P',all$repids),]
barplot(table(plyr::count(all,'id','dup')$freq),xlab='Clone size',main='',las=1,
        ylab='Frequency',axis.lty=1,ylim=c(0,60));box()
mtext('A', side = 3, line = 0.5, adj = 0, cex = 1.1)
x=c();for(i in unique(all$id))x=c(x,length(unique(all[all$id==i,]$repids)))
barplot(table(x),axis.lty=1,las=1,xlab='Number of repertoires',ylim=c(0,100),ylab='Frequency');box()
mtext('B', side = 3, line = 0.5, adj = 0, cex = 1.1)
dev.off()
###############################################################

source('/work/smodi/scripts/covid/func.R')
fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
fs=fs[grepl('/IGH',fs)]
all=data.frame(data.table::rbindlist(parallel::mclapply(fs,function(f){
  x=read.delim(f)
  x=x[nchar(x$cdr3_aa)>0&nchar(x$cdr3_aa)<40,]
  x=data.frame(id=rep(f,nrow(x)),v=alakazam::getGene(x$v_call),
                           j=alakazam::getGene(x$j_call),cdr=x$cdr3_aa,dup=x$duplicate_count)
  x$len=nchar(x$cdr)
  x=x[x$cdr!='',]
  if(nrow(x)==0)return(NULL)
  x$cluster=paste(x$cdr,x$v,x$j)
  y=plyr::count(x,'cluster','dup')
  y$id=f
  y$cdr=gsub(' ','',substr(gsub(' ','                                                                  ',y$cluster),1,40))
  y$len=nchar(y$cdr)
  y$cluster=substr(y$cluster,2+y$len,100)
  y$cluster=paste(y$cluster,y$len)
  y$freq=y$freq/sum(y$freq)
  return(y)
},mc.cores=40)))

consensus=function(z){
  r=''
  for(i in 1:nchar(z$cdr[1])){
    x=data.frame(table(substr(z$cdr,i,i)))
    r=paste0(r,x[order(x$Freq,decreasing = T),]$Var1[1])
  }
  return(r)
}
create=function(z){
  z=mclapply(unique(z$cluster),function(c){
    x=z[z$cluster==c,]
    if(nrow(x)==1)return(x)
    d=stringDist(x$cdr,method='hamming')
    h=hclust(d)
    x$cluster=paste(x$cluster,cutree(h,h =x$len[1]*0.15))
    return(x)
  },mc.cores=45)
  z=data.frame(rbindlist(z))
  k=data.frame(table(z$cluster))
  k=k[k$Freq>2,]
  k$seq=unlist(mclapply(1:nrow(k),function(i){
    return(consensus(z[z$cluster==k$Var1[i],]))
  },mc.cores=45))
  z=merge(z,k[,c(1,3)],by.x='cluster',by.y='Var1')
  k$col=paste(k$seq,k$Var1)
  z$col=paste(z$seq,z$cluster)
  y=data.frame(col=k$col)
  for(i in unique(z$id))y=merge(y,plyr::count(z[z$id==i,],'col','freq'),by='col',all=T)
  
  
  rownames(y)=y$col;y$col=NULL;colnames(y)=unique(z$id)
  y[is.na(y)]=0
  return(y)
}
assign=function(test,names){
  p=mclapply(1:nrow(names),function(i){
    x=test[test$cluster==names$cluster[i],]
    if(nrow(x)==0)return(0)
    d=as.matrix(Biostrings::stringDist(c(names$cdr[i],x$cdr),method='hamming'))
    d=d[1,];d=d[2:length(d)]
    x=x[d<0.15*x$len[1],]
    return(sum(x$freq))
  },mc.cores=45)
  return(unlist(p))
  
}

res=NULL;for(j in 1:50){
i=order(runif(length(unique(all$id))))
train=all[all$id%in%unique(all$id)[i[1:50]],]
test=all[all$id%in%unique(all$id)[i[51:60]],]
print(i[1:50])
print(i[51:60])
y=create(train)
names=data.frame(name=rownames(y))
names$cdr=gsub(' ','',substring(gsub(' ','                                                         '
      ,names$name),1,40))
names$v=gsub(' ','',substr(gsub(' ','                ',substr(names$name,nchar(names$cdr)+2,100)),1,15))
names$cluster=paste(names$cdr,names$v)
names$cluster=paste(names$cluster,
     gsub(' ','',substr(gsub(' ','                ',substr(names$name,nchar(names$cluster)+2,100)),1,15)
          ),nchar(names$cdr))         
names$cluster=substr(names$cluster,nchar(names$cdr)+2,200)
y=data.frame(t(y))
y$stage=ifelse(grepl('IGHCov',rownames(y)),'covid','control')
k=c();for(i in 1:(ncol(y)-1))k=c(k,t.test(y[y$stage=='covid',i],y[y$stage!='covid',i])$p.value);k[is.na(k)]=1
r=train(stage~.,y[,c(k<0.2,T)],method='glmnet')
z=data.frame(col=names$name)
for(i in unique(test$id))z[,i]=assign(test[test$id==i,],names)
rownames(z)=colnames(y)[1:(ncol(y)-1)]
z$col=NULL
z=data.frame(t(z))
z$stage=ifelse(grepl('IGHCov',rownames(z)),'covid','control')
print(rownames(z)%in%rownames(y))
p=predict(r,z)
res=rbind(res,data.frame(i=rep(j,10),stage=z$stage,prediction=p))
write.csv(res,'~/resNewClones.csv')
n=c();for(i in 1:j)n=c(n,sum(res[res$i==i,]$stage==res[res$i==i,]$prediction)/10)
boxplot(n);print(res)
}
res=read.csv('~/resNewClones.csv')
library(MLmetrics);n=c();sp=c();se=c();F1=c();for(i in 1:max(res$i)){
  n=c(n,sum(res[res$i==i,]$stage==res[res$i==i,]$prediction)/10)
  sp=c(sp,Specificity(res[res$i==i,]$stage,res[res$i==i,]$prediction,positive = 'covid'))
  se=c(se,Sensitivity(res[res$i==i,]$stage,res[res$i==i,]$prediction,positive = 'covid'))
  F1=c(F1,F1_Score(res[res$i==i,]$stage,res[res$i==i,]$prediction,positive = 'covid'))
}
pdf('/work/smodi/OurCovid/pdfWAziv/NEWsupclones.pdf',width=6,height=4)
boxplot(list(F1,n,sp,se),names=c('F1 Score','Accuracy','Specificity','Sensitivity'),
        ylim=c(0,1),xlab,las=1);print(res)
dev.off()



#####################################################################################################
#####################################################################################################
f=list.files('/work/smodi/OurCovid/28_10_2021_long/',pattern='tsv$',full.names=T)
f2=list.files('/work/smodi/OurCovid/28_10_2021_long/more',pattern='tsv$',full.names=T)
f=f[!grepl('IGL',f)&!grepl('IGK',f)]
d=read.csv('/work/smodi/Ourcovidsamples.csv')
d$id=gsub('H_collapsed_cloned_w_nf.tsv','',gsub('/work/smodi/cancer/','',gsub(
  '/work/smodi/covid/poria/heavy/','',gsub('_cloned_w_nf.tsv','',gsub(
    'H_collapsed.tsv','',gsub('/work/smodi/covidAyelet/','',d$processed.file))))))
f=c(f,f2)
f=gsub('_cloned_w_filtered_seqs.tsv','',gsub('/work/smodi/OurCovid/28_10_2021_long//IGH','',f))
f=gsub('/work/smodi/OurCovid/28_10_2021_long/more/IGH','',f)
f=gsub('H','',gsub('_collapsed','',gsub('severe','',f)))
f=gsub('/work/smodi/OurCovid/28102021long/more/','',gsub('_','',f))
d$id=gsub(' ','',gsub('H','',gsub('\\.tsv','',gsub('_','',gsub('H_collapsed','',d$id)))))
sum(d$id%in%f)
f=f[!grepl('P4',f)&!grepl('P8',f)&!grepl('P9',f)&!grepl('P1',f)]
d=d[d$id%in%f,]
q=read.delim('/work/smodi/OurCovidpazit.csv',sep=',')
q$id=paste0('Cov',q$id)
qq=merge(q,d,by='id',all=T)
qq=qq[qq$id%in%d$id,]
write.csv(qq,'/work/smodi/OurCovidcomb.csv')
#############################################
z=read.csv('/work/smodi/OurCovidcomb.csv')
z=z[!grepl('P',z$id)&!grepl('CT',z$id),]
z$id=gsub('Cov','Cov_',z$id)
library(parallel)
add=mclapply(z$id,function(i){
  a=list.files('/home/bcrlab/yyw/practice/corona/sample_fastq/')
  a=a[startsWith(a,paste0(i,'_'))]
  r=data.frame(id=i)
  k=1
  for(j in 1:length(a)){
    r[,paste0('filename',k)]=a[j]
    system(paste0('cp ','/home/bcrlab/yyw/practice/corona/sample_fastq/',a[j],
                  ' /work/smodi/SRA/',a[j]))
    k=k+1
  }
  a=list.files('/home/bcrlab/yyw/practice/coloncancer/sample_fastq/')
  a=a[startsWith(a,paste0(i,'_'))]
  for(j in 1:length(a)){
    r[,paste0('filename',k)]=a[j]
    system(paste0('cp ','/home/bcrlab/yyw/practice/coloncancer/sample_fastq/',a[j],
                  ' /work/smodi/SRA/',a[j]))
    k=k+1
  }
  return(r)
},mc.cores=40)


add=mclapply(z$id,function(i){
  a=list.files('/home/bcrlab/yyw/practice/corona/sample_fastq/')
  a=a[startsWith(a,paste0(i,'_'))]
  r=data.frame(id=i)
  k=1
  if(length(a)>0)
  for(j in 1:length(a)){
    r[,paste0('filename',k)]=a[j]
    k=k+1
  }
  a=list.files('/home/bcrlab/yyw/practice/coloncancer/sample_fastq/')
  a=a[startsWith(a,paste0(i,'_'))]
  if(length(a)>0)
  for(j in 1:length(a)){
    r[,paste0('filename',k)]=a[j]
    k=k+1
  }
  return(r)
},mc.cores=40)
add=data.frame(data.table::rbindlist(add,fill=T))
z=merge(z,add,by='id')
write.csv(z,'~/sra.csv')


f=list.files('/work/smodi/SRASC')
n=gsub('R','',unique(substr(f,1,13)))
z=lapply(n,function(j){
  rv=data.frame(id=paste0(j,'-vdj'))
  x=f[startsWith(f,j)&grepl('vdj',f)]
  for(i in 1:length(x))rv[,paste0('filename',i)]=x[i]
  re=data.frame(id=paste0(j,'-exp'))
  x=f[startsWith(f,j)&grepl('exp',f)]
  for(i in 1:length(x))re[,paste0('filename',i)]=x[i]
  return(rbind(rv,re))            
})
z=data.frame(data.table::rbindlist(z))
write.csv(z,'~/allscbind222.csv')







library(shazam)
library(parallel)
MM = function(rep,model){
  ret=data.frame(id='')
  rep=rep[!is.na(rep$clone_id),]
  x=collapseClones(rep,method='mostCommon') 
  print('starting muatbility model')
  model <- createTargetingModel(x,sequenceColumn = 'clonal_sequence',germlineColumn = 'clonal_germline',
                                model=model, vCallColumn="v_call")
  
  x=model@substitution
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('substi',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  x=model@targeting
  for(y in 1:5){
    z=data.frame(t(x[y,]))
    colnames(z)=paste0('target',rownames(x)[y],'_',colnames(z))
    ret=merge(ret,z)
  }
  x=data.frame(model@mutability)
  z=data.frame(t(x))
  ret=merge(ret,z)
  return(ret)
}
run=''
if(run=='runMMOnMore'){
  fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/more',pattern='s.tsv$',full.names=T)
  #fs=fs[grepl('/P',fs)]
  mclapply(fs,function(i){
    rep=read.delim(i)
    rep$c_call=gsub('IGH','IG',toupper(rep$c_call))
    rep$c_call=gsub('H','',gsub('E','',rep$c_call))
    #write.csv(MM(rep,'rs'),paste0(i,'.MM.csv'))
    #write.csv(MM(rep,'s'),paste0(i,'.silentMM.csv'))
    write.csv(MM(rep[rep$c_call=='IGD'|rep$c_call=='IGM',],'rs'),paste0(i,'.MM.IGDM.csv'))
    write.csv(MM(rep[rep$c_call=='IGD'|rep$c_call=='IGM',],'s'),paste0(i,'.silentMM.IGDM.csv'))
    write.csv(MM(rep[rep$c_call=='IGA'|rep$c_call=='IGG',],'rs'),paste0(i,'.MM.IGGA.csv'))
    write.csv(MM(rep[rep$c_call=='IGA'|rep$c_call=='IGG',],'s'),paste0(i,'.silentMM.IGGA.csv'))
    
  },mc.cores=60)
  
  
  fs=list.files('/work/smodi/OurCovid/28_10_2021_longREP/',pattern='s.tsv$',full.names=T)
  mclapply(fs,function(i){
    rep=read.delim(i)
    rep$c_call=gsub('IGH','IG',toupper(rep$c_call))
    #write.csv(MM(rep,'rs'),paste0(i,'.MM.csv'))
    #write.csv(MM(rep,'s'),paste0(i,'.silentMM.csv'))
    rep$c_call=gsub('H','',gsub('E','',rep$c_call))
    write.csv(MM(rep[rep$c_call=='IGD'|rep$c_call=='IGM',],'rs'),paste0(i,'.MM.IGDM.csv'))
    write.csv(MM(rep[rep$c_call=='IGD'|rep$c_call=='IGM',],'s'),paste0(i,'.silentMM.IGDM.csv'))
    write.csv(MM(rep[rep$c_call=='IGA'|rep$c_call=='IGG',],'rs'),paste0(i,'.MM.IGGA.csv'))
    write.csv(MM(rep[rep$c_call=='IGA'|rep$c_call=='IGG',],'s'),paste0(i,'.silentMM.IGGA.csv'))
    
  },mc.cores=60)
  
}

if(run=='byVgene'){
  fs=list.files('/work/smodi/OurCovid/28_10_2021_long/additionals/',pattern='tsv$',full.names=T)
  #fs=fs[grepl('/P',fs)]
  mclapply(fs,function(i){
    rep=read.delim(i)
    rep$c_call=gsub('IGH','IG',toupper(rep$c_call))
    rep$c_call=gsub('H','',gsub('E','',rep$c_call))
    write.csv(MM(rep[startsWith(rep$v_call,'IGHV1'),],'rs'),paste0(i,'.IGHV1.csv'))
    write.csv(MM(rep[startsWith(rep$v_call,'IGHV2'),],'rs'),paste0(i,'.IGHV2.csv'))
    write.csv(MM(rep[startsWith(rep$v_call,'IGHV3'),],'rs'),paste0(i,'.IGHV3.csv'))
    write.csv(MM(rep[startsWith(rep$v_call,'IGHV4'),],'rs'),paste0(i,'.IGHV4.csv'))
    write.csv(MM(rep[startsWith(rep$v_call,'IGHV5'),],'rs'),paste0(i,'.IGHV5.csv'))
    
  },mc.cores=60)
  
  
  
}



########################################################################################
########################################################################################

source('/work/smodi/scripts/covid/func.R')


fs=list.files('/work/smodi/OurCovid/28_10_2021_long/additionals//',pattern='\\.IGHV1.csv$',full.names = T)
IGHV1=loadFiles(fs)
IGHV1$X=NULL;IGHV1$id=NULL
IGHV1$stage=factor(ifelse(grepl('100',fs),'control','case'))
V1=traintest(IGHV1[1:60,],pv=0.001,iterations = 50)

fs=list.files('/work/smodi/OurCovid/28_10_2021_long/additionals//',pattern='\\.IGHV2.csv$',full.names = T)
IGHV2=loadFiles(fs)
IGHV2$X=NULL;IGHV2$id=NULL
IGHV2$stage=factor(ifelse(grepl('100',fs),'control','case'))
V2=traintest(IGHV2[1:60,],pv=0.005,iterations = 50)

fs=list.files('/work/smodi/OurCovid/28_10_2021_long/additionals//',pattern='\\.IGHV3.csv$',full.names = T)
IGHV3=loadFiles(fs)
IGHV3$X=NULL;IGHV3$id=NULL
IGHV3$stage=factor(ifelse(grepl('100',fs),'control','case'))
V3=traintest(IGHV3[1:60,],pv=0.001,iterations = 50)

fs=list.files('/work/smodi/OurCovid/28_10_2021_long/additionals//',pattern='\\.IGHV4.csv$',full.names = T)
IGHV4=loadFiles(fs)
IGHV4$X=NULL;IGHV4$id=NULL
IGHV4$stage=factor(ifelse(grepl('100',fs),'control','case'))
V4=traintest(IGHV4[1:60,],pv=0.001,iterations = 50)

fs=list.files('/work/smodi/OurCovid/28_10_2021_long/additionals//',pattern='\\.IGHV5.csv$',full.names = T)
IGHV5=loadFiles(fs)
IGHV5$X=NULL;IGHV5$id=NULL
IGHV5$stage=factor(ifelse(grepl('100',fs),'control','case'))
V5=traintest(IGHV5[1:60,],pv=0.001,iterations = 50)


par(mfrow=c(1,5),mar=c(5.5,3.5,3,1));
boxplot(V1,las=2,ylim=c(0,1),main='IGHV1');
par(mar=c(5.5,1,3,1))
boxplot(V2,yaxt='n',las=2,ylim=c(0,1),main='IGHV2');
boxplot(V3,yaxt='n',las=2,ylim=c(0,1),main='IGHV3');
boxplot(V4,yaxt='n',las=2,ylim=c(0,1),main='IGHV4');
boxplot(V5,yaxt='n',las=2,ylim=c(0,1),main='IGHV5');
