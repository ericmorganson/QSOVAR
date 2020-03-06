#/bin/bash
for i in $(seq 1 10)
do
   j="${i}_${2}"
   echo "Run $j"
   python ClusterEmcee2_newprob_mu_all_bands.py ../../C1_lc.fits 75375 75376 $1 $j
done
