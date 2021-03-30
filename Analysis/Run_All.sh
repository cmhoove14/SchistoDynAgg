#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M choover@berkeley.edu
#$ -m beas
# 
Rscript Analysis/0-Data-Process.R  
Rscript Analysis/11-Agg-Burden_GEE.R  
Rscript Analysis/12-kappa_change_bootstrap.R  
Rscript Analysis/21-ABC-Worm-Burden-Estimation.R 1000000  
Rscript Analysis/22-ABC-Posteriors.R 1000  
Rscript Analysis/23-ABC-Mod-Comps.R  
Rscript Analysis/24-ABC-PostPred-Checks.R  

