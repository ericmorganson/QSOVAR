#/bin/bash

python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 0 100 normal castor_linear > castor_linear_0_100.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 100 200 normal castor_linear > castor_linear_100_200.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 200 300 normal castor_linear > castor_linear_200_300.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 300 400 normal castor_linear > castor_linear_300_400.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 400 500 normal castor_linear > castor_linear_400_500.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 500 600 normal castor_linear > castor_linear_500_600.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 600 700 normal castor_linear > castor_linear_600_700.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 700 800 normal castor_linear > castor_linear_700_800.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 800 900 normal castor_linear > castor_linear_800_900.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 900 1000 normal castor_linear > castor_linear_900_1000.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 1000 1100 normal castor_linear > castor_linear_1000_1100.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 1100 1200 normal castor_linear > castor_linear_1100_1200.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 1200 1300 normal castor_linear > castor_linear_1200_1300.txt 2>&1 &
python ClusterEmcee2_linear_mu_2step.py ../C1_lc.fits 1300 1400 normal castor_linear > castor_linear_1300_1400.txt 2>&1 &
