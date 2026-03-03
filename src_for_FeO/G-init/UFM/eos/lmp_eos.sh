#!/bin/bash

temperature=$1
LMP=$5
# LMP=~/software/lammps-patch_21Nov2023/src/lmp_mpi
pot_file=../pot_combine.fs
pressure=$(echo "$2 * 10000" | bc)
type=$3
x=$4
if [[ $type == "hcp" ]]; then
    ifiso="aniso"
else
ifiso="iso"
fi

state="../Fe.data"
mpirun -np 20 $LMP -var xx ${x} -var pot_file ${pot_file} -var ifiso ${ifiso} -var state ${state} -var temperature ${temperature} -var tarpressure ${pressure}  -in lmp_eos.in > output_${temperature} 

# python  result.py > eos.data
