#!/bin/bash

# LMP=$1
# atomDir=$2
# potFile=$3
# element1=$4

thermo_file=$1 
potFile=$2
atomDir=$3
LMP=$4
element1=$5
element2=$6

# srun -n 1 -c 1 --exclusive --cpu-bind=cores \
# $LMP -var atomdir $atomDir -var pot_file $potFile -var element1 $element1 -var element2 $element2 -var thermo_file thermo.json -in in.ls
mpirun -np 16 $LMP -var atomdir $atomDir -var pot_file $potFile -var element1 $element1 -var element2 $element2 -var thermo_file thermo.json -in in.ls

cp "thermo.json" ${thermo_file}