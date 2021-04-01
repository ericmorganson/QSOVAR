# QSOVAR
A package for analyzing astrophysical variability (particularly for quasars)

Requiremenats:
astropy
numpy
emcee

Weep all ye who enter here.

Lol, jk, here is how you run the MCMC codes in this folder.

REQUIRED FILES:
Fits files from Eric Morganson
a_b_calc.txt (located in this directory, created with Chi_squared_LCs_all_band_data.py)


If you want to run a deep fits file (C3 or X3), it will require clustering of the data points, which will be done automatically for X3 or C3 fits files.

For all  fits files (C1, C2, or others; ask Eric Morganson for these), run the code in this manner changing the fits file, start & stop rows, name of the files, and figure directory as desired:
python -u ClusterEmcee2_linear_mu_all_and_single_fast.py ../X1_lc.fits 0 500 normal castor_all_and_sing_x1_cluster figure_castor/ > castor_x1_0_500_cluster.txt 2>&1 &

The skeleton for running this is:

python -u ClusterEmcee2_linear_mu_all_and_single_fast.py fits_file start_row end_row+1 file_name_you_want_to_give figure_dir/ > output_for_debugging.txt 2>&1 &
