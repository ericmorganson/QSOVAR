#!/bin/bash

#SBATCH --job-name=cluster_emcee_c1      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail	
#SBATCH --nodes=3                    # Run all processes on two  node	
#SBATCH --ntasks=120                   # Run 80 single tasks		
#SBATCH --ntasks-per-node=40            # Number of ntasks per node
#SBATCH --time=04:00:00              # Time limit hrs:min:sec
#SBATCH --output=serial_%j.log     # Standard output and error log
#SBATCH --partition=secondary

first_row=0
last_row=60000 #121658
delta=$((( last_row - first_row ) / (SLURM_NTASKS)))
fig_fold="figure_cc/"
name="cc_all_and_sing_c1_serial"
run_type="normal"
fits="/home/thrush2/caps_dir/C1_lc.fits"

date

module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}

for m in {1..120}; do
    next_row=$(( m * delta + first_row ))
    start_row=$(( next_row - delta ))
    if [ $m -eq $SLURM_NTASKS ] && [ $next_row -lt $last_row ]
        then
        next_row=$((last_row))
        echo "new last row is $next_row"
    fi
    echo "starting aprun for run $start_row and $next_row with step $delta"
    srun --exclusive --nodes 1 --ntasks 1 python3 -u ClusterEmcee2_linear_mu_all_and_single_fast.py $fits $start_row $next_row $run_type $name $fig_fold > c1_runlogs/cc_allandsing_c1_"$start_row"_"$next_row".txt 2>&1 &
done
wait

date
