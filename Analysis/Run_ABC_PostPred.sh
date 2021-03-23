#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M choover@berkeley.edu
#$ -m beas
# 
Rscript Analysis/22-ABC-Posteriors.R 1000  
Rscript Analysis/23-ABC-Mod-Comps.R  
Rscript Analysis/24-ABC-PostPred-Checks.R  
