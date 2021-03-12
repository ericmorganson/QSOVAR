Weep all ye who enter here.

Lol, jk, here is how you run the MCMC codes in this folder.


If you want to run a deep fits file (C3 or X3), it will require clustering of the data points, which can be run in the following way:
python -u ClusterEmcee2_linear_mu_all_and_single_fast_deep.py ../X3_lc.fits 0 500 normal castor_all_and_sing_x3_cluster figure_castor/ > castor_x3_0_500_deep_cluster.txt 2>&1 &

For all other fits files besides C3 and X3 (C1, C2, or others; ask Eric Morganson for these), run the code in this manner (mind the lack of "deep" before ".py"): 
python -u ClusterEmcee2_linear_mu_all_and_single_fast.py ../X1_lc.fits 0 500 normal castor_all_and_sing_x1_cluster figure_castor/ > castor_x1_0_500_cluster.txt 2>&1 &

The skeleton for running this is:

python -u ClusterEmcee2_linear_mu_all_and_single_fast[_deep].py fits_file start_row end_row+1 file_name_you_want_to_give figure_dir/ > output_for_debugging.txt 2>&1 &
