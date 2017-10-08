import sys

if len(sys.argv) < 4:
    print('make_nodes.py FITS ROWSTART COLOR(first letter)')
string = """#!/bin/bash\n
# declare a name for this job to be sample_job\n
#PBS -N nschwei2_job\n
# request the queue (enter the possible names, if omitted, serial is the default)\n
#PBS -q secondary\n
# request 1 node\n
#PBS -l nodes=1:ppn=24\n
# request 4 hours of cpu time\n
#PBS -l walltime=04:00:00\n
# mail is sent to you when the job starts and when it terminates or aborts \n
#PBS -m bea\n
# specify your email address\n

#PBS -M nschwei2@illinois.edu\n

cd /home/nschwei2/\n
# run the program\n
module load python/2\n
#export PYTHONPATH=$PYTHONPATH:/usr/local/python/3.6.0/lib/python3.6/site-packages:/usr/local/python/2.7.13/mkl-numpy-scipy/lib/python2.7/site-packages:/usr/local/python/2.7.13/lib/python2.7/site-packages\n"""
def makepbs(fits,start,color,string):
    begin = int(start)
    end = begin + 240
    for core in range(24):
        string += "python /home/nschwei2/ClusterEmcee2.py " + fits + " " + str(begin) + " " + str(begin + 10) + " " + color  + "&"  + "\n"
        begin += 10
    string += "wait\n" + "exit 0;\n"
    file_name = "submit_" + start + "_" + str(end) + ".pbs"
    with open(file_name, 'w') as fout:
        fout.write(string)
    return file_name
text = ""
for num in range(10):
    start = int(sys.argv[2]) + num*240
    text += "qsub " + makepbs(sys.argv[1],str(start),sys.argv[3],string) + "\n"

with open('pbs_file_list.txt','w') as fout:
    fout.write(text)
