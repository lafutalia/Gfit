#!/bin/bash

LMP=$3

workdir=$PWD
temperature=$1
x=$2
pot=../pot_combine.fs

eos_dir=$workdir/../eos
echo $eos_dir
state="../Fe.data"

json_file=$eos_dir/thermo.json
vol=$(jq '.vol' $json_file)
natom=$(jq '.natom' $json_file)
k=$(jq '.k' $json_file)
lx=$(jq '.Lx' $json_file)
ly=$(jq '.Ly' $json_file)
lz=$(jq '.Lz' $json_file)
##开三次方获得晶胞边长
lattice=$(echo "scale=10; e(1/3*l($vol))" | bc -l)

mpirun -np 20 $LMP -var xx ${x} -var lx ${lx} -var ly ${ly} -var lz ${lz}  -var lattice ${lattice}  -var pot_file ${pot}  -var state ${state} -var temperature ${temperature}  -in ti-liq.in > ti_output_${temperature} 

wait

echo $natom

python integrate-liq.py ${temperature} ${natom} ${x} > TI.dat



