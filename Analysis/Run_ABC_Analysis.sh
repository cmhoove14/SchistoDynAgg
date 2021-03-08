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
Rscript Analysis/22-ABC_Sims-PostProcess.R  
Rscript Analysis/2F-ABC_Figs.R  
Rscript Analysis/2F2-ABC-Mod-Comp-BayesF.R
Rscript Analysis/2F3-ABC-Mate-Prob-Comp.R
Rscript Analysis/2FS1-ABC-Data-Gen-Results.R
Rscript Analysis/2FS2-ABC-Case3-Comparison.R
Rscript Analysis/2FS3-ABC-EstW-to-obsE-Comp.R
