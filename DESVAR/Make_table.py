import glob

with open('table.txt','wb') as outfile:
    for filename in glob.glob('/home/sam/Documents/Morganson_research/QSOVAR/scratch/*object.txt'): #looks for object.txt files in scratch made from ClusterEmcee2.py
        print(filename)
        with open(filename, 'rb') as infile:
            outfile.write(infile.read())  # compiles object.txt data into a table
