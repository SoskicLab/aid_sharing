This folder contains the script that have been used to send the job and job arrays to run the userGWAS with genomiSEM.

A) chunks/ contains a copy of the sbatch and R.script to run the userGWAS with the model specified.
B)  missing_chunks/ contains a copy of the R.script to run the userGWAS with the model specified for the chunks that fail during the job array (due to matrix multiplication errors)


