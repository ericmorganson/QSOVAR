#!/bin/bash

#SBATCH --job-name=cluster_emcee_c3_x3_conts      # Job name
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

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 124225 124500 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_124225_124500.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 130765 131000 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_130765_131000.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 155920 156000 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_155920_156000.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 216332 216720 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_216332_216720.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 95246 95500 $run_type $name $fig_fold > c3_runlogs/cc_allandsing_c3_95246_95500.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 147235 147500 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_147235_147500.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 158062 158500 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_158062_158500.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 200571 200740 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_200571_200740.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 205399 205620 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_205399_205620.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 26289 26500 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_26289_26500.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 37885 38000 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_37885_38000.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 62437 62500 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_62437_62500.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsX 70594 71000 $run_type $nameX $fig_fold > x3_runlogs/cc_allandsing_x3_70594_71000.txt 2>&1 &

wait

date
