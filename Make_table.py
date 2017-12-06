import glob

with open('table.txt','wb') as outfile:
    for filename in glob.glob('*object.txt'):
        print(filename)
        with open(filename, 'rb') as infile:
            outfile.write(infile.read())
