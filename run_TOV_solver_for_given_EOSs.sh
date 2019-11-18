#!/bin/bash

EOSnames=(
    WFF1
    WFF2
    APR4
    SLy
    ENG
    APR3
    MPA1
    ALF2
    H4
    MS1b
    MS1
)

CentralDensities=(
    3.3e-1
    2.8e-1
    2.8e-1
    2.8e-1
    2.5e-1
    2.5e-1
    2.3e-1
    2.5e-1
    2.0e-1
    2.0e-1
    2.0e-1
)

NEOS=11

for j in {0..10}; do
    ./TOV_Solver ${EOSnames[$j]} ${CentralDensities[$j]}
done

echo Done!
