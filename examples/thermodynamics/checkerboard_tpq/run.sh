#!/bin/bash

n_sites=16
ks=("Gamma.D4.A1" "Gamma.D4.A2" "Gamma.D4.B1" "Gamma.D4.B2" "Gamma.D4.E"
    "M.D4.A1" "M.D4.A2" "M.D4.B1" "M.D4.B2" "M.D4.E" "Sigma.D1.A" "Sigma.D1.B"
    "X.D2.A1" "X.D2.A2" "X.D2.B1" "X.D2.B2") 

# n_sites=20
# ks=("Gamma.C4.A" "Gamma.C4.B" "Gamma.C4.Ea" "Gamma.C4.Eb"
#     "M.C4.A" "M.C4.B" "M.C4.Ea" "M.C4.Eb" "None0.C1.A" "None1.C1.A") 
seeds=`seq 1 1 20`

n_ups=`seq 0 1 $(( n_sites / 2 ))`
J=1.00
Jd=1.00


for n_up in ${n_ups[@]}; do
for k in ${ks[@]}; do
for seed in ${seeds[@]}; do

    ./build/main $n_sites $n_up $k $J $Jd $seed

done
done
done
