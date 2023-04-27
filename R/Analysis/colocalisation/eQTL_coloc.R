source("scripts/Multi_coloc_funtions.R")
library(reshape2)
library(TwoSampleMR)
args = commandArgs(trailingOnly=TRUE)

an.n=as.numeric(args[1])

ct=c("BIN"
,"BMem"
,"CD4ET"
,"CD4NC"
,"CD4SOX4"
,"CD8ET"
,"CD8NC"
,"CD8S100B"
,"DC"
,"MonoC"
,"MonoNC"
,"NK"
,"NKR"
,"Plasma")

an.list=c()
for(i in 1:22){
  
  an.list=rbind(an.list,cbind(i,ct))
  
}


chr=as.numeric(an.list[an.n,1])
cell.type=an.list[an.n,2]

#eqtl.map=fread("data/onek1k_eqtl_dataset.tsv",skip=1,select = c(3,4,7,8,9,10,11))
#eqtl.map=unique(eqtl.map)
#
#names(eqtl.map)=c("SNP","snpid","CHR","BP","A1","A2","FreqA2")
#
#write.table(eqtl.map,file="data/One1k_snp.info",row.names=F,quote=F,sep="\t")

loci.table=fread("output/loci_definitions/final_locus_table.tsv")
loci.table=as.data.frame(loci.table)
loci.table=loci.table[which(!(loci.table$pan.locus%in%c(112, 113, 114, 115, 116, 117, 118, 119, 120, 121))),]
loci.table=loci.table[loci.table$Chr==chr,]


mappa=fread("data/UKBB_30k_map.tsv")

mappa=mappa[mappa$CHR==chr,]

loci.table=loci.table[order(loci.table$start),]
loci.table$final_locus=paste(loci.table$pan.locus,loci.table$sub_locus,sep="_")

loci.table=loci.table[order(loci.table$final_locus),]

for(i in loci.table$pan_locus){
  
  idx=which(loci.table$pan_locus==i)
  loci.table$start[idx]=min(loci.table$start[idx])
  loci.table$end[idx]=max(loci.table$end[idx])
  
}


base.traits=c("uc","cd","psc","ra","sle","t1d","jia","asthma","derma","f1","f2","f3")



colocalization.table.all=c()
colocalization.table.H4=c()

sc.res=fread(paste0("data/onek_scRNA/round1/",cell.type,"_chr",chr,"_correlation_results.tsv"))

eqtl.map=fread("data/One1k_snp.info")
sc.res$SNP=eqtl.map$SNP[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$A1=eqtl.map$A1[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$A2=eqtl.map$A2[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$FREQ=eqtl.map$FreqA2[match(sc.res$snpid,eqtl.map$snpid)]
sc.res$BP=eqtl.map$BP[match(sc.res$snpid,eqtl.map$snpid)]



bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"
file.list=system("ls /project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/*.txt",intern=T)
labels=gsub("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/","",
            gsub("_munged_build37.txt","",file.list))
reference.map=mappa
reference.map=as.data.frame(reference.map)
file.list=file.list
labels=labels
t.type=c("cc","cc","cc","quant","quant","quant","cc","cc","cc","cc","cc","cc")
prev.var=c(0.2696713,0.434383,0.02902916,1,1,1,0.3593954,0.2388718,0.3269585,0.5736819,0.0377603,0.3679372)
file.table=data.frame(file=file.list,trait=labels,type=t.type,prev.var=prev.var)


print(c(cell.type,chr))


sub.loci=unique(loci.table$final_locus)
pleiotropy.table=c()
all.mr=c()
for(locus in sub.loci){
  
  
  loci.table.tmp=loci.table[loci.table$final_locus==locus,]
  locus.info=loci.table.tmp
  start=min(loci.table.tmp$start)-100000
  end=max(loci.table.tmp$end)+100000
  mappa.loc=reference.map[which(reference.map$CHR==chr & reference.map$BP>=start & reference.map$BP<=end),]
  
  
  eqtl.results=sc.res[sc.res$BP>=start & sc.res$BP<=end,]
  
  if(min(eqtl.results$p.value)<=1e-6){
 
    datasets=list()
    for(i in locus.info$trait){
    datasets[[i]]=dataset.munge( sumstats.file = file.table$file[file.table$trait==i]
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
    
    
    gene.list=unique(eqtl.results$geneid)
    
    
    for(i in gene.list){   
      
      
      dataset=eqtl.results[eqtl.results$geneid==i,]
      dataset=as.data.frame(dataset)
      dataset$CHR=chr
      dataset$SE=abs(dataset$estimate/sqrt(qchisq(dataset$p.value,df=1,lower=F)))
        
      if(any(dataset$p.value<1e-6)){
        
        tmp=dataset.munge( sumstats.file = dataset
                           ,map = mappa.loc
                           ,snp.lab = "SNP"
                           ,chr.lab = "CHR"
                           ,pos.lab = "BP"
                           ,a1.lab = "A2"
                           ,a0.lab = "A1"
                           ,beta.lab = "estimate"
                           ,se.lab = "SE"
                           ,pval.lab = "p.value"
                           ,n.lab="N"
                           ,type = "quant"
                           ,freq.lab = "FREQ"
                           ,sdY = 1
                           ,s = 1
                  )
        
        if(any(tmp$pvalues<1e-6)){
           datasets[[i]]=tmp
        }
        
        }
      
    }
    
    
    
    if(length(datasets)>1 & length(which(gene.list%in%names(datasets)))){
      conditional.datasets=list()
      max.loci=1
      for(i in 1:length(datasets)){
        
        tmp=cojo.ht(D=datasets[[i]],p.tresh = 1e-4,bfile = bfile)
        if(names(datasets)[i]%in%locus.info$trait){
          tmp$ind.snps=tmp$ind.snps[tmp$ind.snps$SNP%in%locus.info$SNP[locus.info$trait==names(datasets)[i]],]
          a=names(tmp$results)[names(tmp$results)%in%locus.info$SNP[locus.info$trait==names(datasets)[i]]]
          tmp$results=tmp$results[a]
        }
        conditional.datasets[[i]]=tmp
        names(conditional.datasets)[i]=names(datasets)[i]
        max.loci=max(max.loci,nrow(tmp$ind.snps))
        
      }
      
      for(i in locus.info$trait){
        
        conditional.datasets[[i]]$ind.snps=conditional.datasets[[i]]$ind.snps[which(conditional.datasets[[i]]$ind.snps$pJ<1e-6 | conditional.datasets[[i]]$ind.snps$p<5e-8),]
      }
      
      
     
      pairwise.list=unique(expand.grid(locus.info$trait,gene.list[gene.list%in%names(datasets)],stringsAsFactors = FALSE))
      final.colocs=c()
      all.colocs=list()
      k=1
      
      for(i in 1:nrow(pairwise.list)){
        
        coloc.res=colo.cojo.ht(conditional.dataset1 = conditional.datasets[[pairwise.list[i,1]]],conditional.dataset2 = conditional.datasets[[pairwise.list[i,2]]],p.threshold.orig = 1e-6)
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
                            plot=TRUE, plot.file=paste0("output/sc_eqtl/Locus_",locus,"_colocalization_",k,"_eqtl_",cell.type,".pdf"))
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
               mr.res$cell.type=cell.type
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

}

print(c(cell.type,chr))


pleiotropy.table$cell_type=cell.type
write.table(pleiotropy.table,file=paste0("output/sc_eqtl/Chr",chr,"_Cell_type_",cell.type,"_pleiotropy.tsv"),row.names=F,quote=F,sep="\t")


colocalization.table.all$cell_type=cell.type
write.table(colocalization.table.all,file=paste0("output/sc_eqtl/Chr",chr,"_Cell_type_",cell.type,"_coloc.tsv"),row.names=F,quote=F,sep="\t")


write.table(all.mr,file=paste0("output/sc_eqtl/Chr",chr,"_Cell_type_",cell.type,"_MR.tsv"),row.names=F,quote=F,sep="\t")

write(c(chr, cell.type),ncol=2,file="Completed_runs.tsv",append=TRUE,sep="\t")






