#library(qgraph)
library(corrplot)
library(coloc)
library(susieR)
library(data.table)
library(data.table)
library(bigsnpr)
library(ggplot2)
library(easyGgplot2)
library(igraph)
library(RColorBrewer)
library(ggnet)
library(dplyr)
p.minim=function(x)min(pchisq(as.numeric(unlist((dataset[[i]][x,6:ncol(dataset[[i]])]^2))),df=1,lower=F))
minimiser=function(x){
  p1=unlist(susie.list[[unlist(x["d1"])]]$sets$min.p[paste0("L",unlist(x["idx1"]))])
  p2=unlist(susie.list[[unlist(x["d2"])]]$sets$min.p[paste0("L",unlist(x["idx2"]))])
  min(c(p1,p2))
  
}


locus.joyplot=function(x,window.size=10000,susie.res=susie.list){
  
  start=min(x$pos)
  final.end=max(x$pos)
  log10p=x[,6:ncol(x)]
  res.all=c()
  while(start<final.end){
    
    end=start+window.size
    idx=which(x$pos>=start & x$pos<=end)
    if(nrow(log10p[idx,])>0){
      means=apply(log10p[idx,],2,function(x)x[which.max(abs(x))])
      res.tmp=c(mean(c(start,end)),means)
      res.all=rbind(res.all,res.tmp)
    }
    start=start+window.size/2
    
  }
  library(reshape2)
  res.all=as.data.frame(res.all)
  ordine=hclust(as.dist(1-cor(as.matrix(res.all[,-1]),use ="pairwise.complete.obs")),method = "ward.D2")$order
  ordine=(names(res.all)[-1]) [ordine]
  names(res.all)[1]="pos"
  
  res.all=melt(res.all,id.vars = "pos")
  
  names(res.all)=c("pos","trait","log10p")
  res.all$trait=factor(res.all$trait,levels=ordine)
  
  ggplot(res.all,aes(x=pos,y=log10p))+geom_area(fill="blue",alpha=0.4)+
    geom_hline(yintercept = -log10(5e-8),linetype="dashed",color="red")+
    geom_hline(yintercept = log10(5e-8),linetype="dashed",color="red")+
    facet_grid(scales="free",rows = "trait")+
    theme_minimal()+theme(panel.grid = element_blank())
  
  
}

cojo.ht=function( D=datasets[[1]]
                  ,plink.bin="/project/alfredo/software/plink/1.90_20210606/plink"
                  ,gcta.bin="/project/alfredo/software/GCTA/1.94.0beta/gcta64"
                  ,bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/p01_output/ukbb_all_30000_random_unrelated_white_british"
                  ,p.tresh=1e-4)
  {
  
  require(stringi)
  random.number=stri_rand_strings(n=1, length=20, pattern = "[A-Za-z0-9]")
  
  
  write(D$snp,ncol=1,file=paste0(random.number,".snp.list"))
  system(paste0(plink.bin," --bfile ",bfile," --extract ",random.number,".snp.list --maf 0.0001 --make-bed --freqx --out ",random.number))

  
  
  freqs=fread(paste0(random.number,".frqx"))
  freqs$FreqA1=(freqs$'C(HOM A1)'*2+freqs$'C(HET)')/(2*(rowSums(freqs[,c("C(HOM A1)", "C(HET)", "C(HOM A2)")])))
  D$FREQ=freqs$FreqA1[match(D$snp,freqs$SNP)]
  idx=which(D$a1!=freqs$A1  )
  D$FREQ[idx]=1-D$FREQ[idx]
  D$se=sqrt(D$varbeta)
  
  D=D[,c("snp","a1","a0","FREQ","beta","se","pvalues","N")]
  names(D)=c("SNP" , "A1" ,  "A2"  , "freq", "b"  ,  "se" ,  "p" ,   "N")
  write.table(D,file=paste0(random.number,"_sum.txt"),row.names=F,quote=F,sep="\t")
  #step1 determine independent snps
  system(paste0(gcta.bin," --bfile ",random.number," --cojo-p ",p.tresh," --extract ",random.number,".snp.list  --cojo-file ",random.number,"_sum.txt --cojo-slct --out ",random.number,"_step1"))
  
  ind.snp=fread(paste0(random.number,"_step1.jma.cojo"))
  dataset.list=list()
  dataset.list$ind.snps=ind.snp
  dataset.list$results=list()
  if(nrow(ind.snp)>1){
    for( i in 1:nrow(ind.snp)){
      
      write(ind.snp$SNP[-i],ncol=1,file=paste0(random.number,"_independent.snp"))
      print(ind.snp$SNP[-i])
      system(paste0(gcta.bin," --bfile ",random.number,"  --extract ",random.number,".snp.list  --cojo-file ",random.number,"_sum.txt  --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
      
      step2.res=fread(paste0(random.number,"_step2.cma.cojo"),data.table = FALSE)
      dataset.list$results[[i]]=step2.res
      names(dataset.list$results)[i]=ind.snp$SNP[i]
      
    }
  }else{
    
    write(ind.snp$SNP[1],ncol=1,file=paste0(random.number,"_independent.snp"))
    system(paste0(gcta.bin," --bfile ",random.number," --cojo-p ",p.tresh," --extract ",random.number,".snp.list  --cojo-file ",random.number,"_sum.txt --cojo-slct --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
    step2.res=fread(paste0(random.number,"_step2.cma.cojo"),data.table = FALSE)
    step2.res=step2.res[,c("Chr","SNP","bp","refA","freq","b","se","p","n","freq_geno")]
    dataset.list$results[[1]]=step2.res
    names(dataset.list$results)[1]=ind.snp$SNP[1]
  }
  system(paste0("rm *",random.number,"*"))
  dataset.list
  
}



colo.cojo.ht=function(conditional.dataset1=conditional.datasets[[pairwise.list[i,1]]]
                     ,conditional.dataset2=conditional.datasets[[pairwise.list[i,2]]]
                     ,p.threshold.cond=1e-6
                     ,p.threshold.orig=5e-8
                     )
{
  
  if(length(grep("pJ",names(conditional.dataset1$ind.snps)))>0){
    
    hits.t1=conditional.dataset1$ind.snps$SNP[conditional.dataset1$ind.snps$pJ<p.threshold.cond | conditional.dataset1$ind.snps$p <p.threshold.orig]
  }else{
    hits.t1=conditional.dataset1$ind.snps$SNP
  }
  
  if(length(grep("pJ",names(conditional.dataset2$ind.snps)))>0){
    
    hits.t2=conditional.dataset2$ind.snps$SNP[conditional.dataset2$ind.snps$pJ<p.threshold.cond | conditional.dataset2$ind.snps$p<p.threshold.orig]
  }else{
    hits.t2=conditional.dataset2$ind.snps$SNP
  }
  
  
  
  
  if(length(hits.t2)>0 & length(hits.t1)>0){
  
  
  coloc.results=c()
  for(i in hits.t1){
    for(j in hits.t2){
      
      D1=conditional.dataset1$results[[i]]     
      D2=conditional.dataset2$results[[j]] 
      
      if(length(grep("bC",names(D1)))>0){
        
        D1=D1[,c("SNP","Chr","bp","bC","bC_se","n","pC","freq")]
        names(D1)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
        
      }else{
        
        D1=D1[,c("SNP","Chr","bp","b","se","n","p","freq")]
        names(D1)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
        
        
      }
      D1$type="quant"
      D1$varbeta=D1$varbeta^2
      
      if(length(grep("bC",names(D2)))>0){
        
        D2=D2[,c("SNP","Chr","bp","bC","bC_se","n","pC","freq")]
        names(D2)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
      }else{
        
        D2=D2[,c("SNP","Chr","bp","b","se","n","p","freq")]
        names(D2)=c("snp","chr","position","beta","varbeta","N","pvalues","MAF")
        
        
      }
      D2$type="quant"
      D2$varbeta=D2$varbeta^2
      D1=na.omit(D1)
      D2=na.omit(D2)
      colo.res=coloc.abf(D1,D2)
      colo.res=data.frame(t(colo.res$summary))
      colo.res$hit1=i
      colo.res$hit2=j
      coloc.results=rbind(coloc.results,colo.res)
    }
  }
  }else{
    
    coloc.results=NULL
    
    
  }
  
  coloc.results
  
}





bin2lin2=function (D, dotplot = FALSE) 
{
  if (D$type != "cc") 
    stop("type != cc")
  D1 <- D
  beta = (D1$beta/sqrt(D1$varbeta)) / sqrt (D1$N * 2*(D1$MAF)*(1-D1$MAF)) 
  D1$varbeta <- (beta/(D1$beta/sqrt(D1$varbeta)))^2
  D1$beta <-beta
  
  
  z1 <- D1$beta/sqrt(D1$varbeta)
  z <- D$beta/sqrt(D$varbeta)
  if(dotplot==TRUE){
    
    plot(z1,z)
    
  }
  
  D1$quality <- 1
  D1
}

plot.cojo.ht=function(cojo.ht.obj){
  require(ggplot2)
  library(patchwork)
  if(nrow(cojo.ht.obj$ind.snps)>1){
    
    whole.dataset=c()
    for(i in 1:nrow(cojo.ht.obj$ind.snps)){
      
      tmp=cojo.ht.obj$results[[i]]
      tmp$signal=cojo.ht.obj$ind.snps$SNP[i]
      whole.dataset=rbind(whole.dataset,tmp)
      
    }
    
    p1=ggplot(cojo.ht.obj$results[[i]],aes(x=bp,y=-log10(p)))+geom_point(alpha=0.6,size=3)+theme_minimal()+geom_point(data =cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=SNP),size=6,shape=23)
    p2=ggplot(whole.dataset,aes(x=bp,y=-log10(pC),color=signal))+facet_grid(signal~.)+geom_point(alpha=0.8,size=3)+
      theme_minimal()+ggtitle("Conditioned results")
    p3=p1/p2+ plot_layout(heights = c(1, nrow(cojo.ht.obj$ind.snps)+0.2))
    
  }else{
    
    
    p3=ggplot(cojo.ht.obj$results[[1]],aes(x=bp,y=-log10(p)))+geom_point(alpha=0.6,size=3)+theme_minimal()+geom_point(data =cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=SNP),size=6,shape=23)
    
    
    
  }
  
  (p3)
  
}


dataset.munge=function( sumstats.file="/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/cd_munged_build37.txt"
                       ,map=mappa.loc
                       ,snp.lab="SNP"
                       ,chr.lab="CHR"
                       ,pos.lab="BP"
                       ,a1.lab="A1"
                       ,a0.lab="A2"
                       ,beta.lab="BETA"
                       ,se.lab="SE"
                       ,pval.lab="P"
                       ,freq.lab="FRQ"
                       ,n.lab="N"
                       ,type="quant"
                       ,sdY=1
                       ,s=NULL
                       ){




    # Load sumstat
    if(is.character(sumstats.file)){
       dataset=fread(sumstats.file)
    }else{
      dataset=sumstats.file
    }
    
    
    
    if(  !is.null(a1.lab) & a1.lab%in%names(dataset) & !is.null(a0.lab) & a0.lab%in%names(dataset) ){
       names(dataset)[match(c(a1.lab,a0.lab),names(dataset))]=c("A1","A2")
      
   }else{
      stop("a0.lab or a1.lab have not been defined or the column is missing")
    }
    if(!is.null(beta.lab)& beta.lab%in%names(dataset)){
      names(dataset)[names(dataset)==beta.lab]="BETA"
    }else{
      stop("beta.lab has not been defined or the column is missing")
    }
    
    if(!is.null(snp.lab) & snp.lab%in%names(dataset)){
      names(dataset)[names(dataset)==snp.lab]="SNP"
    }else{
      stop("snp.lab has not been defined or the column is missing")
    }
    
    
    
    
    dataset=dataset[which(dataset$SNP %in% map$SNP),]
    
    if(!is.null(chr.lab) & chr.lab%in%names(dataset)){
       names(dataset)[names(dataset)==chr.lab]="CHR"
    }else{
      dataset$CHR=map$CHR[match(dataset$SNP,map$SNP)]
    }
    
    if(!is.null(pos.lab) & pos.lab%in%names(dataset)){
      names(dataset)[names(dataset)==pos.lab]="BP"
    }else{
      dataset$BP=map$BP[match(dataset$SNP,map$SNP)]
    }
    
    if(!is.null(freq.lab) & freq.lab%in%names(dataset)){
      names(dataset)[names(dataset)==freq.lab]="FRQ"
    }else{
      dataset$FRQ=map$MAF[match(dataset$SNP,map$SNP)]
    }
    
    if("FRQ" %in% colnames(dataset)){
        
        dataset$MAF=dataset$FRQ
        dataset$MAF[dataset$MAF>0.5]=1-dataset$MAF[dataset$MAF>0.5]
    }
    

    
    
    if(!is.null(n.lab) & n.lab%in%names(dataset)){
       names(dataset)[names(dataset)==n.lab]="N"
    }else{
      N_hat<-median(1/((2*dataset$MAF*(1-dataset$MAF))*dataset$SE^2),na.rm = T)
      dataset$N=ceiling(N_hat)
    }
    
    if(!is.null(pval.lab) & pval.lab%in%names(dataset)){
      names(dataset)[names(dataset)==pval.lab]="P"
    }else{
      
      dataset$P=pchisq((dataset$BETA/dataset$SE)^2,df=1,lower=F)
      
    }
    
    
    
    
    dataset=dataset[match(map$SNP,dataset$SNP),]
    
    flip=dataset[,c("SNP","CHR","BP","A2","A1","BETA")]
    names(flip)=c("rsid","chr","pos","a0","a1","beta")
    names(map)=c("rsid","chr","pos","maf","a1","a0")
    flip.t=snp_match(sumstats = flip,info_snp = map,join_by_pos=FALSE,strand_flip = FALSE,match.min.prop=0)
    dataset=dataset[match(flip.t$rsid,dataset$SNP),]
    dataset$A1=flip.t$a1
    dataset$A2=flip.t$a0
    dataset$BETA=flip.t$beta
    dataset$varbeta=dataset$SE^2
    dataset=dataset[,c("SNP","CHR","BP","A1","A2","BETA","varbeta","P","MAF","N")]
    dataset$type=type
    if(type=="cc"){
      
      dataset$s=s
      names(dataset)=c("snp","chr","pos","a1","a0","beta","varbeta","pvalues","MAF","N","type","s")
      
    }else if(type=="quant"){
      
      dataset$sdY=sdY
      names(dataset)=c("snp","chr","pos","a1","a0","beta","varbeta","pvalues","MAF","N","type","sdY")
      
      
    }else{
      
      stop("Type has to be either 'cc' or 'quant'")
      
    }
    
   dataset
    
    
}



pleio.table=function(conditional.datasets=conditional.datasets,loc.table=NA,plot=FALSE,plot.file=NULL, index.trait=NULL){
  
  
  ## Remove duplicates
  duplicati=table(loc.table$trait)
  if(any(duplicati>1)){
    doppi=names(duplicati)[duplicati>1]
    for(i in doppi){
      tmp=conditional.datasets[[i]]$ind.snps
      tmp=tmp[which(tmp$SNP%in%loc.table$SNP[loc.table$trait==i]),]
      min.snp=tmp$SNP[which.min(tmp$pJ)]
      loc.table=loc.table[which(loc.table$trait!=i | (loc.table$trait==i & loc.table$SNP==min.snp)),]
    }
      
  }
  
  
  if(nrow(loc.table)>1){
    for(i in 1:nrow(loc.table)){
      
      if(i ==1){
        merged.datasets=conditional.datasets[[loc.table$trait[i]]]$results[[loc.table$SNP[i]]]
        if("bC"%in% names(merged.datasets)){
          
          merged.datasets=merged.datasets[,c("SNP","bC","bC_se","pC")]
          names(merged.datasets)=c("SNP",paste(c("b","se","p"),loc.table$trait[i],sep="_"))
        }else{
          
          merged.datasets=merged.datasets[,c("SNP","b","se","p")]
          names(merged.datasets)=c("SNP",paste(c("b","se","p"),loc.table$trait[i],sep="_"))
          
          
        }
        
        
      }else{
        
        merged.datasets2=conditional.datasets[[loc.table$trait[i]]]$results[[loc.table$SNP[i]]]
        if("bC"%in% names(merged.datasets2)){
          
          merged.datasets2=merged.datasets2[,c("SNP","bC","bC_se","pC")]
          names(merged.datasets2)=c("SNP",paste(c("b","se","p"),loc.table$trait[i],sep="_"))
        }else{
          
          merged.datasets2=merged.datasets2[,c("SNP","b","se","p")]
          names(merged.datasets2)=c("SNP",paste(c("b","se","p"),loc.table$trait[i],sep="_"))
          
          
        }
        merged.datasets=merge(merged.datasets,merged.datasets2,by="SNP",all = FALSE,suffixes =c(""))
        
      }
      
    }
    
    merged.datasets=na.omit(merged.datasets)
    if(is.null(index.trait)  ){
      top.snp=as.data.frame(merged.datasets[which.min(apply(merged.datasets[,grep("p_",names(merged.datasets))],1,min)),])
    }else if (length(grep(index.trait,names(merged.datasets)))==0){
      top.snp=as.data.frame(merged.datasets[which.min(apply(merged.datasets[,grep("p_",names(merged.datasets))],1,min)),])
    }else{
      top.snp=as.data.frame(merged.datasets[which.min(merged.datasets[,paste0("p_",index.trait)]),])
    }
    
    top.snp=reshape2::melt(top.snp)
    matric=matrix(unlist(strsplit(as.character(top.snp$variable),split="_")),ncol=2,byrow = T)
    
    top.snp$variable=matric[,1]
    top.snp$trait=matric[,2]
    if(plot==TRUE){
      
      mappa=conditional.datasets[[1]]$results[[1]][,c("SNP","bp")]
      mappa=na.omit(mappa)
      mappa=merge(mappa,merged.datasets[,c(1,grep("p_",names(merged.datasets)))],by="SNP")
      per.plot=melt(mappa[,-1],id.vars = "bp")
      per.plot$variable=gsub("p_","",per.plot$variable)
      pdf(plot.file,width=14,height=7*nrow(loc.table))
      p1=ggplot(per.plot,aes(x=bp,y=-log10(value),fill=variable))+geom_point(shape=21,color="black",size=4,alpha=0.9)+
        scale_fill_manual(values=RColorBrewer::brewer.pal(n=nrow(loc.table),name = "Paired"))+
        facet_grid(variable~., scales='free')+theme_minimal()+theme(axis.line = element_line(colour="black"),panel.spacing = unit(3, "lines"))
      print(p1)
      dev.off()
    }
    
    
    
  }else{
    
    merged.datasets=conditional.datasets[[loc.table$trait[1]]]$results[[loc.table$SNP[1]]]
    if("bC"%in% names(merged.datasets)){
      
      merged.datasets=merged.datasets[,c("SNP","bC","bC_se","pC")]
      names(merged.datasets)=c("SNP","b","se","p")
    }else{
      
      merged.datasets=merged.datasets[,c("SNP","b","se","p")]
      names(merged.datasets)=c("SNP","b","se","p")
      
      
    }
    top.snp=merged.datasets[which.min(merged.datasets$p),]
    names(top.snp)=c("SNP",paste(c("b","se","p"),loc.table$trait[1],sep="_"))
    top.snp=reshape2::melt(top.snp)
    matric=matrix(unlist(strsplit(as.character(top.snp$variable),split="_")),ncol=2,byrow = T)
    top.snp$variable=matric[,1]
    top.snp$trait=matric[,2]
    
  }
  top.snp
}











