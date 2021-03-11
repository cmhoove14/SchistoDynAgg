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
Rscript Analysis/21-ABC-Worm-Burden-Estimation.R  
Rscript Analysis/22-ABC-Sims-PostProcess.R  
Rscript Analysis/23-ABC-PostPred-Checks.R  
