#!/bin/bash
#SBATCH --job-name=cluster_emcee_C3      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail        
#SBATCH --nodes=1                    # Run all processes on two  node   
#SBATCH --ntasks=40                   # Run 80 single tasks            
#SBATCH --ntasks-per-node=40            # Number of ntasks per node
#SBATCH --time=04:00:00              # Time limit hrs:min:sec
#SBATCH --output=serial_%j.log     # Standard output and error log
#SBATCH --partition=caps

fig_fold="figure_cc/"
name="cc_all_and_sing_c3_serial"
run_type="normal"
fits="/home/thrush2/caps_dir/C3_lc.fits"
date

module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 216331 216343 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_216331_216343.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 155919 156000 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_155919_156000.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 95245 95500 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_95245_95500.txt 2>&1 &
wait
date
