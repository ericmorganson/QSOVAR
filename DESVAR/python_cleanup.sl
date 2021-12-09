#!/bin/bash
#SBATCH --job-name=cleanup_script      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=thrush2@illinois.edu    # Where to send mail        
#SBATCH --nodes=1                    # Run all processes on 1 node      
#SBATCH --ntasks=1                   # Run 1 single task           
#SBATCH --ntasks-per-node=1            # Number of ntasks per node
#SBATCH --time=0:05:00              # Time limit hrs:min:sec, 12:00:00 
#SBATCH --output=clean_%j.log     # Standard output and error log
#SBATCH --partition=caps

DIR='x3'
module load python/3
export PYTHONPATH=/home/thrush2/caps_dir/mypython3:${PYTHONPATH}

sh linux_fails_cleanup.sh ${DIR}
if [ -s /home/thrush2/QSOVAR/DESVAR/${DIR}_fail_lines.txt ]; then
    srun python3 fails_scripts.py /home/thrush2/QSOVAR/DESVAR/${DIR}_fail_lines.txt ${DIR} > ${DIR}_cleanup_script.sl
else
    echo $DIR
    echo 'Already successfully finished'
    rm /home/thrush2/QSOVAR/DESVAR/${DIR}_fail_lines.txt
    rm /home/thrush2/QSOVAR/DESVAR/${DIR}_fails.txt
fi

