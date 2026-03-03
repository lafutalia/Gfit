#!/bin/bash

thermo_file=$1 
potFile=$2
atomDir=$3
LMP=$4
element1=$5
element2=$6
mpirun -np 16  $LMP -var atomdir $atomDir -var pot_file $potFile -var element1 $element1 -var element2 $element2 -in in.thermo

cp "thermo.json" ${thermo_file}
filename=$(basename "$thermo_file")
prefix="pot-"
suffix=".json"
result=${filename#"$prefix"}
file_directory=$(dirname "$thermo_file")
meanforce_file="${file_directory}/meanforce_${result%"$suffix"}.dat"
mf-alloy.x $atomDir/para.in force.atom >out.mf
wait
rm force.atom
cp mean-force-alloy.dat $meanforce_file