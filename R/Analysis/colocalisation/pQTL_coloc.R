source("scripts/Multi_coloc_funtions.R")
library(reshape2)
library(TwoSampleMR)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

mappa=fread("data/UKBB_30k_map.tsv")

act.r=as.numeric(args[1])




loci.table=fread("output/loci_definitions/final_locus_table.tsv")
loci.table=as.data.frame(loci.table)
loci.table=loci.table[which(!(loci.table$pan.locus%in%c(112, 113, 114, 115, 116, 117, 118, 119, 120, 121))),]
loci.table$final_locus=paste(loci.table$pan.locus,loci.table$sub_locus,sep="_")
loci=unique(loci.table$final_locus)
loci=loci[501:578]
#locus=loci[act.r]


loci.table=loci.table[loci.table$final_locus==locus,]
chr=loci.table$Chr[1]

mappa.tmp=mappa[mappa$CHR==chr,]



for(i in loci.table$pan_locus){
  print(i)
  idx=which(loci.table$pan_locus==i)
  loci.table$start[idx]=min(loci.table$start[idx])
  loci.table$end[idx]=max(loci.table$end[idx])
  
}


base.traits=c("uc","cd","psc","ra","sle","t1d","jia","asthma","derma","f1","f2","f3")



colocalization.table.all=c()
colocalization.table.H4=c()

sc.res=fread(paste0("/processing_data/shared_datasets/plasma_proteome/decode/assocs_filtered/cis_eqtls/Chr",chr,"_cis_pqtls.tsv"))




bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"
file.list=system("ls /project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/*.txt",intern=T)
labels=gsub("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/","",
            gsub("_munged_build37.txt","",file.list))
reference.map=mappa.tmp
reference.map=as.data.frame(reference.map)
t.type=c("cc","cc","cc","quant","quant","quant","cc","cc","cc","cc","cc","cc")
prev.var=c(0.2696713,0.434383,0.02902916,1,1,1,0.3593954,0.2388718,0.3269585,0.5736819,0.0377603,0.3679372)
file.table=data.frame(file=file.list,trait=labels,type=t.type,prev.var=prev.var)


print(chr)


pleiotropy.table=c()
all.mr=c()
sc.res=as.data.frame(sc.res)

print(paste("Analising Locus:",locus))

loci.table.tmp=loci.table[loci.table$final_locus==locus,]
locus.info=loci.table.tmp
start=min(loci.table.tmp$start)-100000
end=max(loci.table.tmp$end)+100000
mappa.loc=reference.map[which(reference.map$CHR==chr & reference.map$BP>=start & reference.map$BP<=end),]


eqtl.results=sc.res[which(sc.res$Pos>=start & sc.res$Pos<=end),]

if(min(eqtl.results$Pval)<=5e-8){
  
  datasets=list()
  for(i in locus.info$trait){
    datasets[[i]]=dataset.munge(  sumstats.file = file.table$file[file.table$trait==i]
                                 ,map = mappa.loc
                                 ,snp.lab = "SNP"
                                 ,chr.lab = "CHR"
                                 ,pos.lab = "BP"
                                 ,a1.lab = "A2"
                                 ,a0.lab = "A1"
                                 ,beta.lab = "BETA"
                                 ,se.lab = "SE"
                                 ,pval.lab = "P"
                                 ,n.lab="N"
                                 ,type = file.table$type[file.table$trait==i]
                                 ,sdY = file.table$prev.var[file.table$trait==i]
                                 ,s = file.table$prev.var[file.table$trait==i]
    )
    
  }
  
  
  gene.list=unique(eqtl.results$gene)
  
  
  for(i in gene.list){   
    
    
    dataset=eqtl.results[eqtl.results$gene==i,]
    dataset=dataset[which(dataset$ImpMAF>=0.01),]
    if(nrow(dataset)>0){
      dataset=as.data.frame(dataset)
      dataset$CHR=chr
      if(any(dataset$Pval<5e-8)){
        
        tmp=try(dataset.munge( sumstats.file = dataset
                               ,map = mappa.loc
                               ,snp.lab = "rsids"
                               ,chr.lab = "CHR"
                               ,pos.lab = "Pos"
                               ,a1.lab = "effectAllele"
                               ,a0.lab = "otherAllele"
                               ,beta.lab = "Beta"
                               ,se.lab = "SE"
                               ,pval.lab = "Pval"
                               ,n.lab="N"
                               ,type = "quant"
                               ,freq.lab = "ImpMAF"
                               ,sdY = 1
                               ,s = 1))
        
        if( !is(object = tmp,class2 = "try-error")){
          if(any(tmp$pvalues<5e-8) ){
            datasets[[i]]=tmp
          }
        }
      }
    }
    
  }
  
  
  
  if(length(datasets)>1 & length(which(gene.list%in%names(datasets)))){
    conditional.datasets=list()
    max.loci=1
    for(i in 1:length(datasets)){
      
      if(names(datasets)[i]%in%locus.info$trait){
        tmp=cojo.ht(D=datasets[[i]],p.tresh = 1e-4,bfile = bfile)
        tmp$ind.snps=tmp$ind.snps[tmp$ind.snps$SNP%in%locus.info$SNP[locus.info$trait==names(datasets)[i]],]
        a=names(tmp$results)[names(tmp$results)%in%locus.info$SNP[locus.info$trait==names(datasets)[i]]]
        tmp$results=tmp$results[a]
      }else{
        tmp.res=datasets[[i]]
        ind.snps=tmp.res[which.min(tmp.res$pvalues),]
        ind.snps=ind.snps[,c("chr","snp","pos","a1","MAF","beta","varbeta","pvalues","N","MAF","beta","varbeta","pvalues")]
        ind.snps$LD_r=1
        ind.snps$varbeta=sqrt(ind.snps$varbeta)
        ind.snps[,12]=sqrt(ind.snps[,12])
        names(ind.snps)=c("Chr" ,"SNP", "bp", "refA","freq","b","se",  "p" ,"n","freq_geno"   ,     "bJ"   ,  "bJ_se"     ,     "pJ"    ,  "LD_r")
        results=tmp.res[,c("chr","snp","pos","a1","MAF","beta","varbeta","pvalues","N","MAF","beta","varbeta","pvalues")]
        results$varbeta=sqrt(results$varbeta)
        results[,12]=sqrt(results[,12])
        
        names(results)=c("Chr", "SNP",  "bp" ,"refA","freq","b","se","p","n", "freq_geno","bC","bC_se","pC")
        results=list(results)
        names(results)[[1]]=ind.snps$SNP
        tmp=list(ind.snps=ind.snps,results=results)
      }
      conditional.datasets[[i]]=tmp
      names(conditional.datasets)[i]=names(datasets)[i]
      max.loci=max(max.loci,nrow(tmp$ind.snps))
      
    }
    
    
    
    
    pairwise.list=unique(expand.grid(locus.info$trait,gene.list[gene.list%in%names(datasets)],stringsAsFactors = FALSE))
    final.colocs=c()
    all.colocs=list()
    k=1
    
    for(i in 1:nrow(pairwise.list)){
      
      coloc.res=colo.cojo.ht(conditional.dataset1 = conditional.datasets[[pairwise.list[i,1]]],conditional.dataset2 = conditional.datasets[[pairwise.list[i,2]]],p.threshold.orig = 5e-8)
      if(!is.null(coloc.res)){
        
        coloc.res$t1=pairwise.list[i,1]
        coloc.res$t2=pairwise.list[i,2]
        coloc.res$locus=locus
        final.colocs=rbind(final.colocs,coloc.res)
      }
    }
    if(!is.null(final.colocs)){
      
      final.colocs.H4=final.colocs[round(final.colocs$PP.H4.abf,digits = 2)>=0.90,]
      final.colocs.H4=as.data.frame(final.colocs.H4)
      if(nrow(final.colocs.H4)>0){
        final.colocs.H4$locus=locus
        colocalization.table.H4=rbind(colocalization.table.H4,final.colocs.H4)
        tmp.loc=data.frame(trait=c(final.colocs.H4$t1,final.colocs.H4$t2),SNP=c(final.colocs.H4$hit1,final.colocs.H4$hit2),stringsAsFactors = FALSE)
        tmp.loc=unique(tmp.loc)
        coloc.genes=unique(tmp.loc$trait[which(!(tmp.loc$trait%in%base.traits))])
        p.table=c()
        for(k in coloc.genes){
          bad=coloc.genes[!(coloc.genes%in%k)]
          if(length(bad)>0){
            tmp.loc2=tmp.loc[!(tmp.loc$trait%in%bad),]
          }else{
            tmp.loc2=tmp.loc
          }
          coloc.traits=unique(unlist(colocalization.table.H4[colocalization.table.H4$t2==k,c("t1","t2")]))
          tmp.loc2=tmp.loc2[tmp.loc2$trait%in%coloc.traits,]
          p.table.tmp=pleio.table(conditional.datasets =conditional.datasets,loc.table = tmp.loc2,
                                  index.trait= unique(tmp.loc2$trait[which(!(tmp.loc2$trait%in%base.traits))]),
                                  plot=TRUE, plot.file=paste0("output/pqtls/Locus_",locus,"_colocalization_",k,"_pqtl.pdf"))
          p.table.tmp$index_trait=k
          p.table=rbind(p.table,p.table.tmp)
          out.traits=coloc.traits[!(coloc.traits%in%k)]
          all.res=acast(data = p.table.tmp,trait~variable)
          all.res=as.data.frame(all.res)
          
          for(ou in out.traits){
            
            mr.res=mr_wald_ratio(b_exp = all.res[k,"b"]
                                 ,b_out = all.res[ou,"b"]
                                 ,se_exp = all.res[k,"se"]
                                 ,se_out = all.res[ou,"se"])
            mr.res$trait=ou
            mr.res$gene=k
            mr.res$locus=locus
            all.mr=rbind(all.mr,data.frame((mr.res)))
            
            
          }
          
          
          #### MR
          
          
          
        }
        #### pleiotropy plots
        p.table$locus=locus
        pleiotropy.table=rbind(pleiotropy.table,p.table)
      }
      final.colocs$locus=loci.table.tmp$locus
      colocalization.table.all=rbind(colocalization.table.all,final.colocs)
    }
  }
}  



print(locus)


write.table(pleiotropy.table,file=paste0("output/pqtls//Locus_",locus,"_pqtl_pleiotropy.tsv"),row.names=F,quote=F,sep="\t")


write.table(colocalization.table.all,file=paste0("output/pqtls/Locus_",locus,"_pqtl_coloc.tsv"),row.names=F,quote=F,sep="\t")


write.table(all.mr,file=paste0("output/pqtls/Locus_",locus,"_pqtl_MR.tsv"),row.names=F,quote=F,sep="\t")
write(locus,file="Completed_runs_pqtl.tsv",append=TRUE)









