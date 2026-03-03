#!/bin/bash

templ=$1
LMP=$2
gr_target=$3
parain=$4
read_atom=$5
element1=$6
element2=$7

cp ${read_atom} ./read.in
pot_file=../pot_combine.fs
mpirun -np $SLURM_NTASKS $LMP -var read_atom $read_atom -var templ $templ -var pot_file $pot_file -var element1 $element1 -var element2 $element2 -in in.liquid
compare-gr-v2.py  -i gr.all -t $gr_target
cp $gr_target gr-target.dat
cp $parain para.in

