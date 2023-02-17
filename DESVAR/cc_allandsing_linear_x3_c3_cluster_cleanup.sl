#!/bin/bash

#SBATCH --job-name=cluster_emcee_x3_c3_cleanup      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on two  node	
#SBATCH --ntasks=40                   # Run 80 single tasks		
#SBATCH --ntasks-per-node=40            # Number of ntasks per node
#SBATCH --time=15:00:00              # Time limit hrs:min:sec
#SBATCH --output=serial_%j.log     # Standard output and error log
#SBATCH --partition=caps


fig_fold="figure_cc/"
namex="cc_all_and_sing_x3_serial_cluster"
run_type="normal"
fitsx="/home/thrush2/caps_dir/X3_lc.fits"
namec="cc_all_and_sing_c3_serial_cluster"
fitsc="/home/thrush2/caps_dir/C3_lc.fits"
date

module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsc 155920 156000 $run_type $namec $fig_fold > c3_runlogs_cluster/cc_allandsing_cluster_c3_155500_156000.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsc 216332 216720 $run_type $namec $fig_fold > c3_runlogs_cluster/cc_allandsing_cluster_c3_216312_216720.txt 2>&1 &
#srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsc 52000 52500 $run_type $namec $fig_fold > c3_runlogs_cluster/cc_allandsing_cluster_c3_52000_52500.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsc 95246 95500 $run_type $namec $fig_fold > c3_runlogs_cluster/cc_allandsing_cluster_c3_95000_95500.txt 2>&1 &

srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsx 147235 147528 $run_type $namex $fig_fold > x3_runlogs_cluster/cc_allandsing_cluster_x3_146784_147528.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsx 158062 158688 $run_type $namex $fig_fold > x3_runlogs_cluster/cc_allandsing_cluster_x3_157944_158688.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsx 205399 205560 $run_type $namex $fig_fold > x3_runlogs_cluster/cc_allandsing_cluster_x3_204816_205560.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsx 26289 26500 $run_type $namex $fig_fold > x3_runlogs_cluster/cc_allandsing_cluster_x3_26000_26500.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsx 37885 38000 $run_type $namex $fig_fold > x3_runlogs_cluster/cc_allandsing_cluster_x3_37500_38000.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsx 62437 62500 $run_type $namex $fig_fold > x3_runlogs_cluster/cc_allandsing_cluster_x3_62000_62500.txt 2>&1 &
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fitsx 70594 71000 $run_type $namex $fig_fold > x3_runlogs_cluster/cc_allandsing_cluster_x3_70500_71000.txt 2>&1 &

wait

date
