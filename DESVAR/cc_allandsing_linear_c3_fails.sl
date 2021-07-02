#!/bin/bash

#SBATCH --job-name=cluster_emcee_c3_x3_fails      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on two  node	
#SBATCH --ntasks=16                   # Run 80 single tasks		
#SBATCH --ntasks-per-node=16            # Number of ntasks per node
#SBATCH --time=8:00:00              # Time limit hrs:min:sec
#SBATCH --output=serial_%j.log     # Standard output and error log
#SBATCH --partition=caps

fig_fold="figure_cc/"
name="cc_all_and_sing_c3_serial"
nameX="cc_all_and_sing_x3_serial"
run_type="normal"
fits="/home/thrush2/caps_dir/C3_lc.fits"
fitsX="/home/thrush2/caps_dir/X3_lc.fits"

date

module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 124224 124225 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_124224_124225.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 130764 130765 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_130764_130765.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 155919 155920 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_155919_155920.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 216331 216332 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_216331_216332.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 95245 95246 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_95245_95246.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 147234 147235 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_147234_147235.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 158061 158062 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_158061_158062.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 200570 200571 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_200570_200571.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 205398 205399 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_205398_205399.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 26288 26289 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_26288_26289.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 37884 37885 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_37884_37885.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 62436 62437 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_62436_62437.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 70593 70594 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_70593_70594.txt 2>&1 &

wait

date
