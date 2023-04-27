args<-commandArgs(trailingOnly = TRUE)
counter<- as.numeric(args[1])

#in this script only one function will be run, since I requested to the function to output a log file


library(MungeSumstats) #Summary statistics here were prepared for being munged by the package MungeSumstats (https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html#overview).
library('BSgenome.Hsapiens.NCBI.GRCh38')
library('SNPlocs.Hsapiens.dbSNP144.GRCh38')
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("BSgenome.Hsapiens.1000genomes.hs37d5")



#-----------run the mune function to add rsID and convert everything to the same build -------------

#very important: the output sumstats will have the effect allele coded as A2 (a-two). . 


my_paths <-  c('/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/asthma_han-2020.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/cd_build37_delange-2017.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/derma_sliz-2021.txt',
               ###
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/jia_beta_lopezisac-2020.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/psc_ji-2016.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/ra_eu_okada-2014.txt',
               ###
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/sle_beta_bentham-2015.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/t1d_chiou-2021.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/01_prepare/ready_format/uc_build37_delange-2017.txt',
               ###
               '/project/aid_sharing/AID_sharing/outputs/rev_1/05_munge/ready_munge/f1_ready_for_munge.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/05_munge/ready_munge/f2_ready_for_munge.txt',
               '/project/aid_sharing/AID_sharing/outputs/rev_1/05_munge/ready_munge/f3_ready_for_munge.txt'
               )



#get the names of the gwas
gwas_names <- c('asthma', 'cd', 'derma',
                'jia', 'psc', 'ra',
                'sle', 't1d', 'uc',
                'f1', 'f2', 'f3')


names(my_paths) <- gwas_names #give names just for checking
#builds <- get_genome_builds(as.list(my_paths)) #check the genome build, only t1d and derma are build 38

#the same order as the files!
builds <- c('GRCH37', 'GRCH37', 'GRCH38',
            'GRCH37', 'GRCH37', 'GRCH37',
            'GRCH37', 'GRCH38', 'GRCH37',
            'GRCH37', 'GRCH37', 'GRCH37')

data('sumstatsColHeaders') #load reference file coming from the package  MungeSumstats

#run function 
format_sumstats(
  path= my_paths[counter],
  convert_ref_genome = 'GRCH37',
  ref_genome = builds[counter],
  convert_small_p = F,
  compute_z = FALSE,
  force_new_z = FALSE,
  compute_n = 0L,
  convert_n_int = F,
  analysis_trait = NULL,
  INFO_filter = 0,
  FRQ_filter = 0,
  pos_se = T,   #important for coloc
  effect_columns_nonzero = T, #important for coloc 
  N_std = 5,
  N_dropNA = F,
  rmv_chr = c("X", "Y", "MT"),
  rmv_chrPrefix = TRUE,
  on_ref_genome = TRUE,
  strand_ambig_filter = FALSE,
  allele_flip_check = TRUE,
  allele_flip_drop = TRUE,
  allele_flip_z = TRUE,
  allele_flip_frq = TRUE,
  bi_allelic_filter = TRUE,
  snp_ids_are_rs_ids = TRUE,
  remove_multi_rs_snp = T,     
  frq_is_maf = TRUE,
  sort_coordinates = TRUE,
  nThread = 2,
  save_path =  paste0('/project/aid_sharing/AID_sharing/outputs/rev_1/05_munge/munged/', names(my_paths)[counter], '_munged_build37.txt' ),
  write_vcf = FALSE,
  tabix_index = FALSE,
  return_data = FALSE,
  return_format = "data.table",
  ldsc_format = FALSE,
  log_folder_ind = FALSE,
  log_folder = paste0('/project/aid_sharing/AID_sharing/outputs/rev_1/05_munge/munged/log/',names(my_paths)[counter] ),
  log_mungesumstats_msgs = TRUE,
  imputation_ind = FALSE,
  force_new = FALSE,
  mapping_file = sumstatsColHeaders
)

