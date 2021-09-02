#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
import pandas as pd
from astropy.io import fits
from astropy.table import Table

txtpath = "/home/thrush2/QSOVAR/DESVAR/scratch_cc"
onlyfiles = [f for f in listdir(txtpath) if isfile(join(txtpath, f))]
#csv_file = "all_and_sing_df_all.csv"
pd.set_option('display.max_rows', None)
col_allband = ["band", "fits", "object", "Tau", "V", "Num_Obs", "Mu_Bright"] 
col_allband_coord = ["band", "fits", "object", "COADD_OJBECT_ID", "RA", "Dec",
                     "Tau", "V", "Num_Obs", "Mu_Bright"]
data = pd.DataFrame() 
counter = 0

#csv_rows = pd.read_csv(txtpath+'/'+csv_file, names=col_allband, skiprows=1, delimiter=",")
for f in onlyfiles:
    counter+=1
    #tmp = csv_rows.iloc[[csv_count]]
    tmp_fi = pd.read_csv(txtpath+'/'+f, names=col_allband, skiprows=1, delimiter="\t")
    for index, tmp in tmp_fi.iterrows():
        	
        with fits.open('/home/thrush2/caps_dir/'+tmp['fits']+'_lc.fits') as hdul:
            ra = hdul[1].data['RA'][tmp['object']]
            dec = hdul[1].data['Dec'][tmp['object']]
            co_ID = hdul[1].data['COADD_OBJECT_ID'][tmp['object']]
    
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
        if counter%10000 == 0:
            print(counter)
            print(tmp_new)
        data = pd.concat([data,tmp_new], axis=0)
data = data.sort_values(by=['fits', 'object'])

print(data.shape)
#print(data)

data.to_csv("all_and_sing_df_all_RA_DEC_ID.csv")
