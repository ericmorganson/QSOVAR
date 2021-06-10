#!/bin/bash

#SBATCH --job-name=cluster_emcee_c3_fails      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on two  node	
#SBATCH --ntasks=3                   # Run 80 single tasks		
#SBATCH --ntasks-per-node=3            # Number of ntasks per node
#SBATCH --time=1:00:00              # Time limit hrs:min:sec
#SBATCH --output=serial_%j.log     # Standard output and error log
#SBATCH --partition=caps

first_row=60000
last_row=229034
delta=$((( last_row - first_row ) / (80)))
fig_fold="figure_cc/"
name="cc_all_and_sing_c3_serial"
run_type="normal"
fits="/home/thrush2/caps_dir/C3_lc.fits"

date

module load python/3
export PYTHONPATH=/home/thrush2/scratch/mypython3:${PYTHONPATH}


srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 20178 20179 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_20178_20179.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 21751 21752 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_21751_21752.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 59074 59075 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_59074_59075.txt 2>&1 &

wait

date
