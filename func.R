
library(MLmetrics)
library(caret)
library(parallel)
library(data.table)
library(alakazam)
library(shazam)
library(randomForest)
library(Biostrings)

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


MM = function(rep,model){
  ret=data.frame(id='')
  x=rep
  print('starting muatbility model')
  model <- createTargetingModel(x, model=model, vCallColumn="v_call")
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
  fs=list.files('/work/smodi/OurCovid/28_10_2021_long/more',pattern='s.tsv$',full.names=T)
  fs=fs[grepl('/P',fs)]
  mclapply(fs,function(i){
      rep=read.delim(i)
      rep$c_call=gsub('IGH','IG',toupper(rep$c_call))
      write.csv(MM(rep[rep$c_call=='IGA'|rep$c_call=='IGG',],'rs'),paste0(i,'.MM.IGGA.csv'))
      write.csv(MM(rep[rep$c_call=='IGA'|rep$c_call=='IGG',],'s'),paste0(i,'.silentMM.IGGA.csv'))
    
  },mc.cores=60)
  
  
}


.AA=c('A','C','D','E','F','G','H','I','W','K','L','M','N','P','Q','R','S','T','V','Y')
.target=paste0(rep(.AA,each=103),1:103)
AAtarget=function(rep){
  rep$seq=gsub('\\-','N',gsub('\\.','N',rep$sequence_alignment))
  rep$seq=as.character(translate(DNAStringSet(rep$seq),if.fuzzy.codon = 'X'))
  m=NULL
  for(i in 1:103){
    k=data.frame(table(substr(rep$seq,i,i)))
    k$Var1=paste0(k$Var1,i)
    m=rbind(m,k)
  }
  m=rbind(m,data.frame(Var1=.target,Freq=rep(0,length(.target))))
  m=m[!grepl('X',m$Var1)&!grepl('\\*',m$Var1),]
  m$Freq=m$Freq/sum(m$Freq)
  m=m[!duplicated(m$Var1),]
  rownames(m)=m$Var1
  m$Var1=NULL
  m=data.frame(t(m))
  return(m)
}
AAtarget2=function(rep){
  rep$seq=rep$sequence_alignment_aa
  m=NULL
  for(i in 1:103){
    k=data.frame(table(substr(rep$seq,i,i)))
    k$Var1=paste0(k$Var1,i)
    m=rbind(m,k)
  }
  m=rbind(m,data.frame(Var1=.target,Freq=rep(0,length(.target))))
  m=m[!grepl('X',m$Var1)&!grepl('\\*',m$Var1),]
  m$Freq=m$Freq/sum(m$Freq)
  m=m[!duplicated(m$Var1),]
  rownames(m)=m$Var1
  m$Var1=NULL
  m=data.frame(t(m))
  return(m)
}

traintest = function(data,pv=0.05,method='glmnet',fracTest=0.175,iterations=50){
  return(rbindlist(mclapply(1:iterations,function(j){
    test=data.frame(stage='case')
    while(sum(test$stage=='case')/nrow(test)<0.2|
          sum(test$stage=='case')/nrow(test)>0.8){
      i=runif(nrow(data))
      x=data[order(i)[(nrow(data)*fracTest):nrow(data)],]
      test=data[-order(i)[(nrow(data)*fracTest):nrow(data)],]
    }
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
    p=0
    if(method=='RF'){
      p=data.frame(true=test$stage,prediction=predict(randomForest(stage~.,x[,k]),test),
                   iteration=rep(j,nrow(test)))
      
    }else{
      p=data.frame(true=test$stage,prediction=
                     predict(train(stage~.,x[,k],method=method,tuneLength=5),test),
                   iteration=rep(j,nrow(test)))
    }
    return(data.frame(F1=F1_Score(test$stage,p$prediction,positive='case'),
                      accuracy=Accuracy(test$stage,p$prediction),
                      sensitivity=Sensitivity(test$stage,p$prediction,positive='case'),
                      specificity=Specificity(test$stage,p$prediction,positive='case')))
  },mc.cores=50)))
}


selfPtraintest = function(data,fracTest=0.175,iterations=50){
  return(rbindlist(mclapply(1:iterations,function(j){
    test=data.frame(stage='case')
    while(sum(test$stage=='case')/nrow(test)<0.2|
          sum(test$stage=='case')/nrow(test)>0.8){
      i=runif(nrow(data))
      x=data[order(i)[(nrow(data)*fracTest):nrow(data)],]
      test=data[-order(i)[(nrow(data)*fracTest):nrow(data)],]
    }
    k=c()
    for(i in 1:ncol(x))
      if(colnames(x)[i]=='stage'){
        k=c(k,0)
      }else k=c(k,t.test(x[x$stage!='case',i],x[x$stage=='case',i])$p.value)
    k[is.na(k)]=1
    p=list()
    o=1
    for(pv in c(0.0001,0.001,0.01,0.05,0.025,0.1)){
      g=cut(c(1:nrow(x)),5)
      a=0
      for(j in levels(g))
        if(ncol(x[g!=j,k<pv])>1)
          a=a+sum(predict(train(stage~.,x[g!=j,k<pv],method='glmnet'),
                  x[g==j,k<pv])==x$stage[g==j])
      p[[o]]=a/nrow(x)
      o=o+1
    }
    p=unlist(p)
    pv=c(0.0001,0.001,0.01,0.05,0.025,0.1)[p==max(p)][1]
    p=data.frame(true=test$stage,prediction=
                     predict(train(stage~.,x[,k<pv],method='glmnet',tuneLength=5),test),
                   iteration=rep(j,nrow(test)))
    
    return(data.frame(F1=F1_Score(test$stage,p$prediction,positive='case'),
                      accuracy=Accuracy(test$stage,p$prediction),
                      sensitivity=Sensitivity(test$stage,p$prediction,positive='case'),
                      specificity=Specificity(test$stage,p$prediction,positive='case')))
  },mc.cores=50)))
}








targetFromFiles = function(fs){
  x=mclapply(fs,function(i){
    return(AAtarget(read.delim(i)))
  },mc.cores=40)
  x=rbindlist(x,fill=T)
  x=data.frame(x)
  x[is.na(x)]=0
  return(x)
}
loadFiles = function(fs){
  x=mclapply(fs,function(i){
    read.csv(i)
  },mc.cores=40)
  x=rbindlist(x,fill=T)
  x=data.frame(x)
  x[is.na(x)]=0
  return(x)
}

leaveOneOut = function(data,pv=0.05,method='glmnet'){
  return(unlist(mclapply(1:nrow(data),function(j){
    x=data[-j,]
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
    if(method=='RF')
      return(predict(randomForest(stage~.,x[,k]),data[j,]))
    return(predict(train(stage~.,x[,k],method=method,tuneLength=5),data[j,]))
  },mc.cores=40)))
}
extractFeatures=function(x,pv=0.05,method='glmnet'){
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
  r=train(stage~.,x[,k],method=method,tuneLength=5)
  co=(coef(r$finalModel,r$bestTune$lambda))
  x=as.matrix(co)
  x=data.frame(x)
  x$gene=rownames(x)
  x=x[x$gene!='(Intercept)',]
  x=x[x[,1]!=0,]
}

N=c('A','G','C','T')
#R = A,G ; Y = C,T; M = A,C ; K = G,T; S =C,G; W = A,T; H = A,C,T ; B = C,G,T ; V = A,C,G ; D = A,G,T; N = A,C,G,T
mer=paste0(rep(N,each=256),rep(N,each=64),rep(N,each=16),rep(N,each=4),N)
WA_TW=mer[substr(mer,2,3)%in%c('AA','TA')|substr(mer,3,4)%in%c('TA','TT')]
WRC_GYW=mer[substr(mer,1,3)%in%paste0(rep(c('A','T'),each=2),c('A','G'),'C')|
              substr(mer,3,5)%in%paste('G',rep(c('C','T'),each=2),c('A','T'))]
SYC_GRS=mer[substr(mer,1,3)%in%paste0(rep(c('C','G'),each=2),c('C','T'),'C')|
              substr(mer,3,5)%in%paste('G',rep(c('A','G'),each=2),c('C','G'))]
REST=mer[!mer%in%c(WA_TW,WRC_GYW,SYC_GRS)]


createMutationLoad = function(dir='/work/smodi/OurCovid/28_10_2021/'){
  
  fs=list.files(dir,pattern='tsv$',full.names = T)
  fs=fs[!grepl('LIGHT',fs)]
  z=data.frame(rbindlist(mclapply(fs,function(i){
    rep=read.delim(i)
    rep$c_call=toupper(rep$c_call)
    rep$c_call=gsub('hIGD','IGD',gsub('hIGGE','IGG',gsub('hIGA','IGA',gsub('hIGM','IGM',rep$c_call))))
    rep=rep[rep$consensus_count>1,]
    v=observedMutations(rep, sequenceColumn="sequence_alignment",
                        germlineColumn="germline_alignment_d_mask",
                        regionDefinition=IMGT_V_BY_CODONS,
                        frequency=TRUE,                          
                        nproc=1)
    v=v[,startsWith(colnames(v),'mu_freq')|colnames(v)%in%c('c_call')]
    m=data.frame(t(colMeans(v[,2:ncol(v)])))
    for(i in c('IGD','IGM','IGG','IGA')){
      m[,paste0('s_',i)]=mean(rowMeans(v[v$c_call==i,startsWith(colnames(v),'mu_freq')&endsWith(colnames(v),'s')]))
      m[,paste0('r_',i)]=mean(rowMeans(v[v$c_call==i,startsWith(colnames(v),'mu_freq')&endsWith(colnames(v),'r')]))
    }
    for(i in c('IGD','IGM','IGG','IGA'))
      m[,i]=mean(rep[rep$c_call==i,]$v_identity)
    for(f in c(0.1,0.5,0.9))
      for(i in c('IGD','IGM','IGG','IGA')){
        y=(rowMeans(v[v$c_call==i,startsWith(colnames(v),'mu_freq')&endsWith(colnames(v),'s')]))
        m[,paste('s_',i,'_',f)]=sort(y)[length(y)*f]
        y=(rowMeans(v[v$c_call==i,startsWith(colnames(v),'mu_freq')&endsWith(colnames(v),'r')]))
        m[,paste('r_',i,'_',f)]=sort(y)[length(y)*f]
      }
    
    return(m)
  },mc.cores=40)))
  z$stage=factor(ifelse(grepl('100',fs)|grepl('/P',fs),'control','case'))
  rownames(z)=fs
  write.csv(z,paste0(dir,'/allMutationLoadConst5.csv'))
}



#########################################################################
##########################################################################
if(run=='check'){
  ######################################################################
  library(parallel);library(caret)
  fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
  fs=fs[grepl('/IGH',fs)]
  fs2=(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
  fs2=fs2[!grepl('/P',fs2)&!grepl('AP',fs2)&!grepl('CT',fs2)]
  x=data.frame(data.table::rbindlist(mclapply(48:48,function(i){
    r=read.delim(fs[i])
    r=r[nchar(r$cdr3)<=90,]
    r=r[,colnames(r)%in%c('v_sequence_alignment','v_germline_alignment','cdr3')]
    return(r)},mc.cores=40)))
  x=x[nchar(x$v_sequence_alignment)<450,]
  x=x[nchar(x$v_germline_alignment)<450,]
  y=data.frame(germ=x$v_germline_alignment)
  t=data.frame(pos=rep(1:450,each=4),base=rep(c('A','C','G','T'),450))
  z=y
  y$germ=NULL
  for(i in 1:450)for(n in c('A','C','G','T'))y[,paste(n,i)]=ifelse(
    substr(x$v_germline_alignment,i,i)==n,1,0 )
  for(i in 1:90)for(n in c('A','C','G','T'))y[,paste('cdr3',n,i)]=ifelse(
    substr(x$cdr3,i,i)==n,1,0 )
  z$germ=NULL
  y[is.na(y)]=0;z[is.na(z)]=0
  for(i in 1:450)for(n in c('A','C','G','T'))z[,paste(n,i)]=ifelse(
    substr(x$v_sequence_alignment,i,i)==n,1,0 )
  i=runif(nrow(y))
  trainY=y[i<0.8,];trainZ=z[i<0.8,]
  testY=y[i>=0.8&i<0.9,];testZ=z[i>=0.8&i<0.9,]
  
  suppressPackageStartupMessages(library(keras))
  enc_input = layer_input(shape = ncol(y))
  enc_output = enc_input %>% 
    layer_dense(units=ncol(y), activation = "tanh") %>% 
   # layer_dense(units=ncol(y), activation = "tanh") %>% 
    layer_dense(units=ncol(y), activation = "tanh") %>% 
    layer_dense(units=ncol(z), activation = "tanh")# %>% 
  
  model = keras_model(enc_input, enc_output)
  summary(model);model %>% compile(  loss=loss_mean_absolute_error, 
                                     optimizer='adam' ,run_eagerly = T )
  model %>% fit(  x = as.matrix(trainY),  y = as.matrix(trainZ),
                  validation_data=list(as.matrix(testY),as.matrix(testZ)),epochs = 500,
                  verbose = 2,batch_size=sum(i<0.8))
  a=model%>%predict(as.matrix(y[i>=0.9,]))
  a[a>0.8]=1
  a[a<0.5]=0
  zz=as.matrix(z[i>=0.9,])
  k=c();for(j in 1:nrow(zz))k=c(k,cor((a[j,]),(zz[j,])))
  yy=as.matrix(y[i>0.9,])
  k2=c();for(j in 1:nrow(zz))k2=c(k2,cor((yy[j,1:ncol(zz)]),(zz[j,])))
  boxplot(k,k2);print(summary(k));print(summary(k2))
  
  ##########################################################
  library(parallel);library(caret)
  fs=list.files('/work/smodi/OurCovid/28_10_2021_long//',pattern='tsv$',full.names = T)
  fs=fs[grepl('/IGH',fs)]
  fs2=(list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='tsv$',full.names = T))
  fs2=fs2[!grepl('/P',fs2)&!grepl('AP',fs2)&!grepl('CT',fs2)]
  x=data.frame(data.table::rbindlist(mclapply(48:48,function(i){
    r=read.delim(fs[i])
    r=r[nchar(r$cdr3)<=90,]
    r=r[,colnames(r)%in%c('sequence_alignment','germline_alignment','cdr3')]
    return(r)},mc.cores=40)))
  x=x[nchar(x$sequence_alignment)<450,]
  x=x[nchar(x$germline_alignment)<450,]
  y=data.frame(germ=x$germline_alignment)
  t=data.frame(pos=rep(1:450,each=4),base=rep(c('A','C','G','T'),450))
  z=y
  y$germ=NULL
  for(i in 1:450)for(n in c('A','C','G','T'))y[,paste(n,i)]=ifelse(
    substr(x$germline_alignment,i,i)==n,1,0 )
  for(i in 1:90)for(n in c('A','C','G','T'))y[,paste('cdr3',n,i)]=ifelse(
    substr(x$cdr3,i,i)==n,1,0 )
  z$germ=NULL
  for(i in 1:450)for(n in c('A','C','G','T'))z[,paste(n,i)]=ifelse(
    substr(x$sequence_alignment,i,i)==n,1,0 )
  
  i=runif(nrow(y))
  trainY=y[i<0.8,];trainZ=z[i<0.8,]
  testY=y[i>=0.8&i<0.9,];testZ=z[i>=0.8&i<0.9,]
  
  suppressPackageStartupMessages(library(keras))
  enc_input = layer_input(shape = ncol(y))
  enc_output = enc_input %>% 
    layer_dense(units=ncol(y), activation = "tanh") %>% 
    layer_dense(units=ncol(y), activation = "tanh") %>% 
    #layer_dense(units=ncol(y), activation = "tanh") %>% 
    layer_dense(units=ncol(z), activation = "tanh")# %>% 
  
  model = keras_model(enc_input, enc_output)
  summary(model);model %>% compile(  loss=loss_mean_absolute_error, 
                                     optimizer='adam' ,run_eagerly = T )
  model %>% fit(  x = as.matrix(trainY),  y = as.matrix(trainZ),
                  validation_data=list(as.matrix(testY),as.matrix(testZ)),epochs = 500,
                  verbose = 2,batch_size=sum(i<0.8))
  a=model%>%predict(as.matrix(y[i>=0.9,]))
  a[a>0.8]=1
  a[a<0.5]=0
  zz=as.matrix(z[i>=0.9,])
  k=c();for(j in 1:nrow(zz))k=c(k,cor((a[j,]),(zz[j,])))
  yy=as.matrix(y[i>0.9,])
  k2=c();for(j in 1:nrow(zz))k2=c(k2,cor((yy[j,1:ncol(zz)]),(zz[j,])))
  boxplot(k,k2);print(summary(k));print(summary(k2))
  
  a=model%>%predict(as.matrix(y[i<=0.9,]))
  a[a>0.8]=1
  a[a<0.5]=0
  zz=as.matrix(z[i<=0.9,])
  k=c();for(j in 1:nrow(zz))k=c(k,cor((a[j,]),(zz[j,])))
  yy=as.matrix(y[i<0.9,])
  k2=c();for(j in 1:nrow(zz))k2=c(k2,cor((yy[j,]),(zz[j,])))
  boxplot(k,k2);print(summary(k));print(summary(k2))
  
  z=allGA[1:60,];z$X=NULL
  z$id=NULL;stage=((c(rep(0,28),rep(1,32))))
  i=runif(60)
  train=as.matrix(z[order(i)[1:50],]);test=as.matrix(z[order(i)[51:60],])
  trainY=as.matrix(stage[order(i)[1:50]])
  testY=as.matrix(stage[order(i)[51:60]])
  suppressPackageStartupMessages(library(keras))
  enc_input = layer_input(shape = ncol(train))
  enc_output = enc_input %>% 
    layer_dense(units=10000, activation = "relu") %>% 
    layer_dense(units=2000, activation = "relu") %>% 
    layer_dense(units=400, activation = "relu") %>% 
    layer_dense(units=100, activation = "relu") %>% 
    layer_dense(units=25, activation = "relu")%>%
    layer_dense(units=5, activation = "relu") %>%
    layer_dense(units=1, activation = "relu")
  model = keras_model(enc_input, enc_output)
  summary(model);model %>% compile(  loss=loss_mean_absolute_error(), 
                                     optimizer='adam' ,run_eagerly = T )
  model %>% fit(  x = train,  y = trainY,epochs = 100,
                  verbose = 2,batch_size=50)
  p=model%>%predict(test)
  print(p)
  print(testY)
  
  
  
  z=allGA[1:60,];z$X=NULL
  z$id=NULL;stage=((c(rep(0,28),rep(1,32))))
  i=runif(60)
  train=as.matrix(z[order(i)[1:50],]);test=as.matrix(z[order(i)[51:60],])
  trainY=as.matrix(stage[order(i)[1:50]])
  testY=as.matrix(stage[order(i)[51:60]])
  suppressPackageStartupMessages(library(keras))
  model <- keras_model_sequential()
  model %>%
    layer_dense(units=ncol(train),activation='tanh',input_shape=ncol(train))%>%
    layer_dense(units=10000, activation = "tanh") %>% 
    layer_dense(units=2000, activation = "tanh") %>% 
    layer_dense(units=400, activation = "tanh") %>% 
    layer_dense(units=100, activation = "tanh") %>% 
    layer_dense(units=25, activation = "tanh")%>%
    layer_dense(units=5, activation = "tanh") %>%
    layer_dense(units = 2, activation = "tanh", name = "bottleneck") %>%
    layer_dense(units=5, activation = "tanh") %>% 
    layer_dense(units=25, activation = "tanh") %>% 
    layer_dense(units=100, activation = "tanh") %>% 
    layer_dense(units=400, activation = "tanh") %>% 
    layer_dense(units=2000, activation = "tanh")%>%
    layer_dense(units=10000, activation = "tanh") %>%
    layer_dense(units=ncol(train), activation = "tanh")
    
    summary(model)
    model %>% compile(  loss=loss_mean_squared_error(), optimizer='adam' ,run_eagerly = T )
    model %>% fit(  x = train,  y = train, epochs = 50,verbose = 2,batch_size=50)
    intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
    intermediate_output <- predict(intermediate_layer_model, train)
    p=predict(intermediate_layer_model,test)
    dev.off()
    plot(jitter(intermediate_output[,1]),jitter(intermediate_output[,2]), col=trainY+1)
    points(jitter(p[,1]),jitter(p[,2]),col=testY+1,pch='+')
    
  
    
    train=as.matrix(allGA[1:60,!colnames(allGA)%in%c('X','id')&nchar(colnames(allGA))==5])
    trainY=c(rep(0,28),rep(1,32))
    suppressPackageStartupMessages(library(keras))
    model <- keras_model_sequential()
    model %>%
      layer_dense(units=ncol(train),activation='tanh',input_shape=ncol(train))%>%
      layer_dense(units=10000, activation = "tanh") %>% 
      layer_dense(units=2000, activation = "tanh") %>% 
      layer_dense(units=400, activation = "tanh") %>% 
      layer_dense(units=100, activation = "tanh") %>% 
      layer_dense(units=25, activation = "tanh")%>%
      layer_dense(units=5, activation = "tanh") %>%
      layer_dense(units = 1, activation = "tanh", name = "bottleneck") %>%
      layer_dense(units=5, activation = "tanh") %>% 
      layer_dense(units=25, activation = "tanh") %>% 
      layer_dense(units=100, activation = "tanh") %>% 
      layer_dense(units=400, activation = "tanh") %>% 
      layer_dense(units=2000, activation = "tanh")%>%
      layer_dense(units=10000, activation = "tanh") %>%
      layer_dense(units=ncol(train), activation = "tanh")
    
    summary(model)
    model %>% compile(  loss=loss_mean_squared_error(), optimizer='adam' ,run_eagerly = T )
    model %>% fit(  x = train,  y = train, epochs = 30,verbose = 2,batch_size=60)
    intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
    intermediate_output <- predict(intermediate_layer_model, train)
    plot(jitter(intermediate_output[,1]),jitter(intermediate_output[,2]), col=trainY+1)
    plot((intermediate_output[,1]),(intermediate_output[,2]), col=trainY+1)
    plot(intermediate_output[,1],col=1+trainY)
    p=model%>%predict(train)
    cor(train[1,],p[1,])
    plot(train[1,],p[1,])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    source('/work/smodi/scripts/covid/func.R')
    fs=list.files('/work/smodi/OurCovid/28_10_2021_long/',pattern='\\.MM.IGGA.csv$',full.names = T)
    fs=c(fs,list.files('/work/smodi/OurCovid/28_10_2021_long/more/',pattern='\\.MM.IGGA.csv$',full.names = T))
    fs=fs[!grepl('AP',fs)&!grepl('CT',fs)]
    allGA=loadFiles(fs)
    allGA=allGA[,startsWith(nchar(colnames(allGA)),'sub')]
    res=NULL;n=c()
    for(j in 1:30){
      z=allGA[1:60,];z$X=NULL
      z$id=NULL;stage=((c(rep(0,28),rep(1,32))))
      i=runif(60)
      train=as.matrix(z[order(i)[1:50],]);test=as.matrix(z[order(i)[51:60],])
      trainY=as.matrix(stage[order(i)[1:50]])
      testY=as.matrix(stage[order(i)[51:60]])
      k=c();for(i in 1:ncol(train))k=c(k,t.test(train[trainY==0,i],train[trainY==1,i])$p.value)
      k[is.na(k)]=1;train=train[,k<0.01];test=test[,k<0.01]
      suppressPackageStartupMessages(library(keras))
      model <- keras_model_sequential()
      model %>%
        layer_dense(units=ncol(train),activation='tanh',input_shape=ncol(train))%>%
        layer_dense(units=500, activation = "tanh") %>% 
        layer_dense(units=200, activation = "tanh") %>% 
        layer_dropout(rate=0.2)%>%
        layer_dense(units=75, activation = "tanh") %>% 
        #layer_dropout(rate=0.75)%>%
        layer_dense(units=25, activation = "relu")%>%
        layer_dense(units=5, activation = "relu") %>%
        layer_dense(units = 1, activation = "relu", name = "bottleneck") 
      
      summary(model)
      model %>% compile(  loss=loss_mean_squared_error(), optimizer="adam",run_eagerly = T )
      h=0.5
      model %>% fit(  x = train,  y = trainY, epochs = 500,verbose = 2,batch_size=50)
      while(h>0.05){
        h=model %>% fit(  x = train,  y = trainY, epochs = 200,verbose = 2,batch_size=50)
        h=median(h$metrics$loss)
      }
      p=model%>%predict(test)
      n=c(n,sum((unlist(data.frame(p)[,1])>=0.5)==(unlist(data.frame(testY)[,1])==1))/10)
      print(n)
      res=rbind(res,data.frame(test=testY,pred=p))
    }
    
    
}

