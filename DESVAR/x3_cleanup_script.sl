#!/bin/bash
#SBATCH --job-name=cluster_emcee_X3      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail        
#SBATCH --nodes=1                    # Run all processes on two  node   
#SBATCH --ntasks=40                   # Run 80 single tasks            
#SBATCH --ntasks-per-node=40            # Number of ntasks per node
#SBATCH --time=04:00:00              # Time limit hrs:min:sec
#SBATCH --output=serial_%j.log     # Standard output and error log
#SBATCH --partition=caps

fig_fold="figure_cc/"
name="cc_all_and_sing_x3_serial"
run_type="normal"
fits="/home/thrush2/caps_dir/X3_lc.fits"
date

module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}
srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 158199 158543 $run_type $name $fig_fold > x3_runlogs/cc_allandsing_x3_158199_158543.txt 2>&1 &
#srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 37885 38000 $run_type $name $fig_fold > x3_runlogs/cc_allandsing_x3_37885_38000.txt 2>&1 &
#srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 147235 147565 $run_type $name $fig_fold > x3_runlogs/cc_allandsing_x3_147235_147565.txt 2>&1 &
#srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 205399 205620 $run_type $name $fig_fold > x3_runlogs/cc_allandsing_x3_205399_205620.txt 2>&1 &
#srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 70594 71000 $run_type $name $fig_fold > x3_runlogs/cc_allandsing_x3_70594_71000.txt 2>&1 &
#srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 62437 62500 $run_type $name $fig_fold > x3_runlogs/cc_allandsing_x3_62437_62500.txt 2>&1 &
#srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits 26289 26500 $run_type $name $fig_fold > x3_runlogs/cc_allandsing_x3_26289_26500.txt 2>&1 &
wait
date
