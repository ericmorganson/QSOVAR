#!/bin/bash

python ClusterEmcee2_newprob_mu.py ../C1_lc.fits 17999 18000 g normal &
python ClusterEmcee2_newprob_mu.py ../C1_lc.fits 3347 3348 g normal &
python ClusterEmcee2_newprob_mu.py ../C1_lc.fits 9004 9005 g normal &
python ClusterEmcee2_newprob_mu.py ../C1_lc.fits 47662 47663 g normal &
wait
exit;
