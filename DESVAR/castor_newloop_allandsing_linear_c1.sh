#!/bin/bash

fig_fold="figure_castor/"
name="castor_all_and_sing_c1"
run_type="normal"
fits="../C1_lc_fits"
end_row=252303
nodes=10
start_row=0
row_step=`echo "(($end_row - $start_row) / $nodes)/1" | bc` #-l #DON'T YOU FUCKING DARE CHANGE THIS
#row_step=400
next_row=`echo "($start_row + $row_step)/1" | bc` #-l

for m in $(seq 1 $nodes); do
    echo "starting aprun for run $start_row and $next_row with step $row_step"
    python ClusterEmcee2_linear_mu_all_and_single_fast.py $fits $start_row $next_row $run_type $name $fig_fold > castor_allandsing_"$start_row"_"$next_row".txt 2>&1 &
    start_row=`echo "($start_row + $row_step )/1" | bc` #-l
    next_row=`echo "($next_row + $row_step)/1" | bc` #-l
done

if [ $start_row -lt $end_row ]
then
    res=`echo "($end_row - $start_row)/1" | bc `
    next_row=`echo "($start_row + $res)/1" | bc `
    echo "starting aprun for run $start_row and $next_row with step $res"
    python ClusterEmcee2_linear_mu_all_and_single_fast.py $fits $start_row $next_row $run_type $name $fig_fold > castor_allandsing_"$start_row"_"$next_row".txt 2>&1 &
# Wait for everything to finish up.
fi
wait

echo "Finished CastorLoop script"

