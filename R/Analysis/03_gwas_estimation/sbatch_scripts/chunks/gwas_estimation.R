#this scripts is for running a GWAS with GenomicSEM package in a fast way
#it takes the GWAS summary stats and splits into chunks so that it can be run as job array 

args<-commandArgs(trailingOnly = TRUE)
counter<-args[1]

#------------------------------------total number of chunks---------------------

# how_many_chunks <- function(summary_stats, chunk_size) {
#   require(data.table)
#   #give a number to the elements of the sum_states
#   my_index <- seq_along(1: nrow(summary_stats))
#   n <- nrow(summary_stats)
#   chunks <- rep(1:ceiling(n/chunk_size),each=chunk_size)[1:n] #the last piece [1:n] is necessary because otherwise the last chunk will become long as the chunk size and will replicate it thus exceeding the length of the sumstats
#   a <- split(my_index, chunks)
#   return(length(a))
# }

#----Function for splitting summary stats --------------------------------------

split_sum_stats <- function(summary_stats, chunk_size, which_chunk) {
  require(data.table)
  #give a number to the elements of the sum_states
  my_index <- seq_along(1: nrow(summary_stats))
  n <- nrow(summary_stats)
  chunks <- rep(1:ceiling(n/chunk_size),each=chunk_size)[1:n] #the last piece [1:n] is necessary because otherwise the last chunk will become long as the chunk size and will replicate it thus exceeding the length of the sumstats
  a <- split(my_index, chunks)
  summary_stats[ (a[[which_chunk]]) , ]
}

#-------------------------------------------------------------------------------

library(GenomicSEM)
library(data.table)

aid_sumstats <- fread('/project/aid_sharing/AID_sharing/outputs/rev_1/02_munge_ldsc_sumstats/aid_sumstats.txt', data.table = F)
ldsc_model <- readRDS('/project/aid_sharing/AID_sharing/outputs/rev_1/02_munge_ldsc_sumstats/ldsc_output_rev1.RDS')

aid_model <- 'F1 =~ NA*crohn + uc  + psc  
F2 =~ NA*jia + sle + ra+ t1d 
F3 =~ NA*asthma+ derma 

F1~~F2
F1~~F3
F2~~F3
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3

F1 ~ SNP
F2 ~ SNP
F3 ~ SNP

derma~~a*derma
a>0.001'


#dim(aid_sumstats) #5050785      24
#how_many_chunks(aid_sumstats, 10000) #506

chunk_to_use <- split_sum_stats(aid_sumstats, 10000, counter )

output <- userGWAS(covstruc = ldsc_model, 
                   SNPs = chunk_to_use, 
                   model = aid_model, 
                   sub = c("F1~SNP", "F2~SNP", 'F3~SNP'), 
                   parallel = F , 
                   cores = 1)


saveRDS(output, file.path('/project/aid_sharing/AID_sharing/outputs/rev_1/03_gwas_estimation/chunks/', paste0(counter, '_gwas_rev_1.RDS')))
















