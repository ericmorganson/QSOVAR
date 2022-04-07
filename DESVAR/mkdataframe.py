#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import csv

txtpath = "/home/thrush2/QSOVAR/DESVAR/scratch_cc"
csvpath = "/home/thrush2/QSOVAR/DESVAR/"
onlyfiles = [f for f in listdir(txtpath) if isfile(join(txtpath, f))]
#csv_file = "all_and_sing_df_all.csv"
pd.set_option('display.max_rows', None)
col_allband = ["band", "fits", "object", "Tau", "V", "Num_Obs", "Mu_Bright"] 
col_allband_coord = ["band", "fits", "object", "COADD_OJBECT_ID", "RA", "Dec",
                     "Tau", "V", "Num_Obs", "Mu_Bright"]
data_frame = pd.DataFrame() 
counter = 0

#csv_rows = pd.read_csv(txtpath+'/'+csv_file, names=col_allband, skiprows=1, delimiter=",")
for f in onlyfiles:
    if "optimal" in f:
        #print("optimal")
        continue    
    #if "C3" in f or "X3" in f:
    #    pass
    #else:
    #    continue
    counter+=1
    #tmp = csv_rows.iloc[[csv_count]]
    tmp_fi = pd.read_csv(txtpath+'/'+f, names=col_allband, skiprows=1, delimiter="\t")
    #print(f)
    #print(tmp_fi)
    
    for index, tmp in tmp_fi.iterrows():   #check this to make sure it's doing what you think it is!!
        #print(tmp)
        if str(tmp['fits']) == 'nan':
            print(tmp)
            print(f)
            continue
        
        with fits.open('/home/thrush2/caps_dir/'+str(tmp['fits'])+'_lc.fits') as hdul:
            ra = hdul[1].data['RA'][int(tmp['object'])]
            dec = hdul[1].data['Dec'][int(tmp['object'])]
            co_ID = hdul[1].data['COADD_OBJECT_ID'][int(tmp['object'])]
        #TODO have if statement checking if V and tau are bad. If so, print fits & tmp_fi data
        if tmp['V']<-3 or tmp['V']>2 or tmp['Tau']<0 or tmp['Tau']>4:
            print("Oh no!  Out of bounds!")
            if "optimal" in f:
                print("optimal")
            #print(tmp_fi)
            #print(tmp)
            print(f)
            with open(txtpath+'/'+f, newline='') as csvfile:
                scratch_file = csv.reader(csvfile, delimiter='\t')
                #for row in scratch_file:
                    #print('\t'.join(row))
            #continue
        tmp_new = pd.DataFrame({"band": [tmp['band']], 
                       "fits": [tmp['fits']], 
                        "object": [tmp['object']],
                        "COADD_OBJECT_ID": [co_ID], 
                        "RA": [ra], 
                        "Dec": [dec], 
                        "Tau": [tmp['Tau']], 
                        "V": [tmp['V']],
                        "Num_Obs": [tmp['Num_Obs']], 
                        "Mu_Bright": [tmp['Mu_Bright']]
                       })
    #print(tmp_new)
    #if counter%10000 == 0:
        #print(counter)
        #print(tmp_new)
    data_frame = pd.concat([data_frame,tmp_new], axis=0)
    #if counter > 100:
    #    break
data_frame = data_frame.sort_values(by=['fits', 'object'])

print(data_frame.shape)
#print(data)

data_frame.to_csv(csvpath+"all_and_sing_df_all_RA_DEC_ID_no_optimal_new.csv")
