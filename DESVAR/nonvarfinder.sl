#!/bin/bash

#SBATCH --job-name=nonvar_runs      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail        
#SBATCH --nodes=6                    # Run all processes on two  node   
#SBATCH --ntasks=240                   # Run 80 single tasks             
#SBATCH --ntasks-per-node=40            # Number of ntasks per node
#SBATCH --time=12:00:00              # Time limit hrs:min:sec
#SBATCH --output=nonvar_debug.log     # Standard output and error log
#SBATCH --partition=caps

date

module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}


#"C3_lc.fits", "C2_lc.fits", "C3_lc.fits",
#"E1_lc.fits", "E2_lc.fits",
#"S1_lc.fits", "S2_lc.fits",
#"X1_lc.fits", "X2_lc.fits", "X3_lc.fits"]
srun --exclusive --nodes 1 --ntasks 1 python3 -u non_var_finder.py "S1_lc.fits" "s1" > nonvar_debug_s1.txt &
srun --exclusive --nodes 1 --ntasks 1 python3 -u non_var_finder.py "S2_lc.fits" "s2" > nonvar_debug_s2.txt &
srun --exclusive --nodes 1 --ntasks 1 python3 -u non_var_finder.py "E2_lc.fits" "e2" > nonvar_debug_e2.txt &
srun --exclusive --nodes 1 --ntasks 1 python3 -u non_var_finder.py "X1_lc.fits" "x1" > nonvar_debug_x1.txt &
srun --exclusive --nodes 1 --ntasks 1 python3 -u non_var_finder.py "X2_lc.fits" "x2" > nonvar_debug_x2.txt &
srun --exclusive --nodes 1 --ntasks 1 python3 -u non_var_finder.py "X3_lc.fits" "x3" > nonvar_debug_x3.txt &

wait

date
