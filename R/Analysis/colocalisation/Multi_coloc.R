###### Multicoloc 

source("scripts/Multi_coloc_funtions.R")

loci.table=fread("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/hyprcoloc_all_against_all/loci_all_traits.txt")
loci.table=loci.table[order(loci.table$pan_locus),]
for(i in loci.table$pan_locus){
  
  idx=which(loci.table$pan_locus==i)
  loci.table$start[idx]=min(loci.table$start[idx])
  loci.table$end[idx]=max(loci.table$end[idx])
  
}

loci.table=loci.table[which(!(loci.table$pan_locus%in%c(112, 113, 114, 115, 116, 117, 118, 119, 120, 121))),]




mappa=fread("data/UKBB_30k_map.tsv")



colocalization.table.all=c()
colocalization.table.H4=c()
bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"

## GWAS info
file.list=system("ls /project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/*.txt",intern=T)
labels=gsub("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/07_colocalization/munged/","",gsub("_munged_build37.txt","",file.list))

proportions=readRDS("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_sumstat_outputs/case_controls_fraction_nic.RDS")

reference.map=mappa
t.type=c("cc","cc","cc","quant","quant","quant","cc","cc","cc","cc","cc","cc")
prev.var=c(0.2696713,0.434383,0.02902916,1,1,1,0.3593954,0.2388718,0.3269585,0.5736819,0.0377603,0.3679372)

file.table=data.frame(file=file.list,trait=labels,type=t.type,prev.var=prev.var)
col.order=c("trait","Chr",  "start", "end", "SNP" ,  "bp","refA","othA","freq","b","se","p","bJ","bJ_se","pJ","LD_r" ,"n" , "pan.locus" , "sub_locus" )
final.locus.table=c()



for( locus in unique(loci.table$pan_locus)){
  
  loci.table.tmp=loci.table[loci.table$pan_locus==locus,]
  locus.info=loci.table.tmp
  #Define genomi region
  start=min(loci.table.tmp$start)-100000
  end=max(loci.table.tmp$end)+100000
  chr=loci.table.tmp$chr[1]
  n.table=c()
  ##### load and harmonise the SNP data
  
  datasets=list()
  mappa.loc=mappa[which(mappa$CHR==chr & mappa$BP>=start & mappa$BP<=end),]
  
  
  #munge files
  for(i in unique(loci.table.tmp$trait)){
    
    datasets[[i]]=dataset.munge(sumstats.file = file.table$file[file.table$trait==i]
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
                               ,s = file.table$prev.var[file.table$trait==i])
  
    
  }  

  conditional.datasets=list()
  max.loci=1
  for(i in 1:length(datasets)){
    
    tmp=cojo.ht(D=datasets[[i]],p.tresh = 1e-4,bfile = bfile)
    conditional.datasets[[i]]=tmp
    names(conditional.datasets)[i]=names(datasets)[i]
    max.loci=max(max.loci,nrow(tmp$ind.snps))
    
  }
  
  #Plot 
  pdf(paste0("output/loci_definitions/Locus_",locus,"_Conditioned_loci.pdf"),height=3.5*max.loci,width=10)
  
  for(i in 1:length(conditional.datasets)){
    p4=plot.cojo.ht(conditional.datasets[[i]])+plot_annotation(paste("Locus",locus,names(conditional.datasets)[i]))
    print(p4)
  }
  
  dev.off()
  
  
  
  
  
  # Only if there are multiple traits.
  if(length(conditional.datasets)>1){
  
    pairwise.list=t(combn(names(datasets),2))
    final.colocs=c()

    
    ### Prepare final locus.table Part I
    final.locus.table.tmp=c()
    for(k in names(conditional.datasets)){
      
      tmp=conditional.datasets[[k]]$ind.snps
      if(nrow(tmp)>1){
        tmp=tmp[tmp$p<5e-8 | tmp$pJ<1e-6,]
      }
      
      tmp$start=loci.table.tmp$start[1]
      tmp$end=loci.table.tmp$end[1]
      tmp$trait=k
      tmp$othA=NA
      for(h in 1:nrow(tmp)){
        alleles=unlist(mappa.loc[mappa.loc$SNP==tmp$SNP[h],c("A1","A2")])
        tmp$othA[h]=alleles[!(alleles%in%tmp$refA[h])]
      }
      tmp$pan.locus=locus
      tmp$sub_locus=NA
      final.locus.table.tmp=rbind(final.locus.table.tmp,tmp)
      
    }
    
  
    ### Run colocalisation adn collect results
        
    for(i in 1:nrow(pairwise.list)){
    
      coloc.res=colo.cojo.ht(conditional.dataset1 = conditional.datasets[[pairwise.list[i,1]]]
                            ,conditional.dataset2 = conditional.datasets[[pairwise.list[i,2]]]
                            ,p.threshold.cond = 1e-6
                            ,p.threshold.orig = 5e-8 
                            )
    
      if(!is.null(coloc.res)){
        coloc.res$t1=pairwise.list[i,1]
        coloc.res$t2=pairwise.list[i,2]
        final.colocs=rbind(final.colocs,coloc.res)
      }
    
    }
    
    #Identify colocalizations
    
    final.colocs.H4=final.colocs[round(final.colocs$PP.H4.abf,digits = 2)>=0.90,]
    final.colocs.H4=as.data.frame(final.colocs.H4)
    
    
    
    
    

    
    
    
    # If any colocalizations are present 
    if(nrow(final.colocs.H4)>0){
      
      # create groups 
      
      a.graph=graph_from_data_frame(final.colocs.H4[,c("hit1","hit2")],directed=F)
  
      groups=components(a.graph)
      groups=data.frame(snp=names(groups$membership),group=groups$membership)
      
      final.colocs.H4$g1=groups$group[match(final.colocs.H4$hit1,groups$snp)]
      
      
      # Add sublocus to locus table
      
      
      
      for(k in 1:nrow(final.colocs.H4)){
        
        final.locus.table.tmp$sub_locus[final.locus.table.tmp$trait==final.colocs.H4$t1[k] &
                                final.locus.table.tmp$SNP==final.colocs.H4$hit1[k]]=final.colocs.H4$g1[k]
        final.locus.table.tmp$sub_locus[final.locus.table.tmp$trait==final.colocs.H4$t2[k] &
                                          final.locus.table.tmp$SNP==final.colocs.H4$hit2[k]]=final.colocs.H4$g1[k]
        
        
      }
      # Number non colocalising loci
      if(any(is.na(final.locus.table.tmp$sub_locus))){
        idx=which(is.na(final.locus.table.tmp$sub_locus))
        pri=max(final.locus.table.tmp$sub_locus,na.rm=T)+1
        final.locus.table.tmp$sub_locus[idx]=pri:(pri+length(idx)-1)
      }
      final.locus.table.tmp=as.data.frame(final.locus.table.tmp)[,col.order]
      
      
      
      k=1
      ### Create pleiotropy table
      pleio.all=c()
      sub.loci=sort(unique(final.locus.table.tmp$sub_locus))
      for(i in sub.loci){
        tmp=pleio.table(conditional.datasets = conditional.datasets,loc.table = final.locus.table.tmp[final.locus.table.tmp$sub_locus==i,])
        tmp$sublocus=i
        pleio.all=rbind(pleio.all,tmp)
      }
      pleio.all=unique(pleio.all)
      
      pleio.all=reshape2::dcast(pleio.all,formula = SNP+sublocus+trait~variable,fill=NA)
      pleio.all=pleio.all[order(pleio.all$sublocus), ]
      pleio.all$Z=pleio.all$b/pleio.all$se
      
      library(corrplot)
      a=dcast(pleio.all[,c("SNP","trait","Z")],SNP~trait,fill = 0)
      row.names(a)=a$SNP
      pdf(paste0("output/loci_definitions/Pleiotropy_table_Locus",locus,".pdf"),width=(dim(a)[1]/dim(a)[2])*7)
      
         corrplot(t(as.matrix(a[,-1])),is.corr = F,method = "color",addgrid.col = 'darkgrey',col=COL2('RdBu', 200),col.lim=c(max(abs(a[,-1]))*-1,max(abs(a[,-1]))))
      
      dev.off()
    
      
      
      
      
### Trasformare in funzione
      
      per.plot.data=c()
      for(i in 1:nrow(final.colocs.H4)){
        tmp=conditional.datasets[[final.colocs.H4$t1[i]]]$results[[final.colocs.H4$hit1[i]]]
        tmp$label=paste(final.colocs.H4$t1[i],final.colocs.H4$hit1[i],sep="-")
        tmp$group=final.colocs.H4$g1[i]
        if(length(grep("pC",names(tmp)))>0){
      
          tmp=tmp[,c("bp","pC","label","group")]
          names(tmp)=c("bp","p","label","group")
      
        }else{
      
        tmp=tmp[,c("bp","p","label","group")]

        }
        per.plot.data=rbind(per.plot.data,tmp)
        tmp=conditional.datasets[[final.colocs.H4$t2[i]]]$results[[final.colocs.H4$hit2[i]]]
        tmp$label=paste(final.colocs.H4$t2[i],final.colocs.H4$hit2[i],sep="-")
        tmp$group=final.colocs.H4$g1[i]
        if(length(grep("pC",names(tmp)))>0){
        
          tmp=tmp[,c("bp","pC","label","group")]
          names(tmp)=c("bp","p","label","group")
        }else{
      
          tmp=tmp[,c("bp","p","label","group")]
      
        }
        per.plot.data=rbind(per.plot.data,tmp)
    
        per.plot.data=per.plot.data[order(per.plot.data$group),]
        per.plot.data$label=factor(per.plot.data$label,levels=unique(per.plot.data$label))
  
        final.colocs.H4$locus=locus
      }
      
      # parameter file
       
      pdf(paste0("output/loci_definitions/locus_",locus,"_colocalization_plot.pdf"),width=14,height=4*length(unique(per.plot.data$label)))
      
      p1=ggplot(per.plot.data,aes(x=bp,y=-log10(p),fill=as.character(group)))+
        geom_point(shape=21,alpha=0.9,color='black',size=4)+geom_hline(yintercept = 0,size=0.5)+facet_grid(label~.,scales="free")+theme_minimal()+
        scale_fill_manual(values=brewer.pal(length(unique(per.plot.data$group)),name="Paired"))+
        scale_color_manual(values=brewer.pal(length(unique(per.plot.data$group)),name="Paired"))+
        theme(strip.text.y.right = element_text(angle = 270,size=20),axis.line.y=element_line())+ggtitle(paste("Locus",locus,"Conditional regional plot"))
      print(p1)
      dev.off()
      
      ### Function plot ends
      
      
    }else{
      
      idx=which(is.na(final.locus.table.tmp$sub_locus))
      pri=1
      final.locus.table.tmp$sub_locus[idx]=pri:(pri+length(idx)-1)
      final.locus.table.tmp=as.data.frame(final.locus.table.tmp)[,col.order]
      
    }  
    if(nrow(final.colocs)>0){
       final.colocs$locus=locus
       colocalization.table.all=rbind(colocalization.table.all,final.colocs)
    }
    if(nrow(final.colocs.H4)>0){
       colocalization.table.H4=rbind(colocalization.table.H4,final.colocs.H4)
    }
    
  }else{
    
   final.locus.table.tmp=conditional.datasets[[1]]$ind.snps
   final.locus.table.tmp$start=locus.info$start
   final.locus.table.tmp$end=locus.info$end
   final.locus.table.tmp$pan.locus=locus
   final.locus.table.tmp$sub_locus=1
   final.locus.table.tmp$freq_geno=NA
   final.locus.table.tmp$bJ=NA
   final.locus.table.tmp$bJ_se=NA
   final.locus.table.tmp$pJ=NA
   final.locus.table.tmp$LD_r=NA

   
   
   alleles=unlist(mappa.loc[mappa.loc$SNP==final.locus.table.tmp$SNP,c("A1","A2")])
   final.locus.table.tmp$othA=alleles[!(alleles%in%final.locus.table.tmp$refA)]
   final.locus.table.tmp$trait=names(datasets)[1]
   final.locus.table.tmp=as.data.frame(final.locus.table.tmp)
   final.locus.table.tmp=final.locus.table.tmp[,col.order]  
   
  }
  
  final.locus.table=rbind(final.locus.table,final.locus.table.tmp)
  print(final.locus.table)
}


write.table(final.locus.table,file="output/loci_definitions/final_locus_table.tsv",row.names=F,quote=F,sep="\t")
write.table(colocalization.table.all,file="output/loci_definitions/colocalization.table.all.tsv",row.names=F,quote=F,sep="\t")
write.table(colocalization.table.H4,file="output/loci_definitions/colocalization.table.H4.tsv",row.names=F,quote=F,sep="\t")





