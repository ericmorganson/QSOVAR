#!/usr/bin/env python

import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from astropy.io import fits
#import seaborn as sns
#import astropy.io.fits as pyfit
import sys
import math as m

failed_txt = sys.argv[1] 
failed_fits = sys.argv[2].upper()
prologue="""#!/bin/bash
#SBATCH --job-name=cluster_emcee_{0}      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail        
#SBATCH --nodes={2}                    # Run all processes on two  node   
#SBATCH --ntasks={1}                   # Run 80 single tasks            
#SBATCH --ntasks-per-node=40            # Number of ntasks per node
#SBATCH --time=04:00:00              # Time limit hrs:min:sec
#SBATCH --output=serial_%j.log     # Standard output and error log
#SBATCH --partition=caps

fig_fold="figure_cc/"
name="cc_all_and_sing_{3}_serial"
run_type="normal"
fits="/home/thrush2/caps_dir/{0}_lc.fits"
date

module load python/3"""

if __name__ == "__main__":
    count=0
    failed_txts = open(failed_txt, 'r')
    failed_lines = failed_txts.readlines()
    line_count = len(failed_lines)
    print(prologue.format(failed_fits, 40*m.ceil(line_count/40), m.ceil(line_count/40),
                          failed_fits.lower()))
    print("export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}")
    line_count = len(failed_lines)
    #print(line_count)
    for line in failed_lines:
        if count%2==0:
            start_l = line.split(' ')[-1].strip()
        if count%2==1:
            next_l = line.split('_')[-1].split(".")[0]
        
        if count>=1 and count%2==1:
            #print(start_l, next_l)
            srun_script = 'srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits {0} {1} $run_type $name $fig_fold > {2}_runlogs/cc_allandsing_{2}_{0}_{1}.txt 2>&1 &'.format(start_l, next_l, failed_fits.lower())
            print(srun_script)
        count+=1
    print('wait')
    print('date')
