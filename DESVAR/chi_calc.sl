#!/bin/bash

#SBATCH --job-name=chi_runs      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail        
#SBATCH --nodes=1                    # Run all processes on two  node   
#SBATCH --ntasks=16                   # Run 80 single tasks             
#SBATCH --ntasks-per-node=16            # Number of ntasks per node
#SBATCH --time=8:00:00              # Time limit hrs:min:sec
#SBATCH --output=chi_debug.log     # Standard output and error log
#SBATCH --partition=caps

date

module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}

#srun --exclusive --nodes 1 --ntasks 1 
python3 -u Chi_squared_LCs_all_band_data.py > chi2_log.txt

wait

date
