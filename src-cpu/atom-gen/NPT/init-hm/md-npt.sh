#!/bin/bash

templ=$1
LMP=$2
read_atom=$3
press=$4
element1=$5
element2=$6
x=$7
ifiso=$8
cp ${read_atom} ./read.in
pot_file=../pot_combine.fs
mpirun -np 20 $LMP  -var x $x -var ifiso $ifiso -var disp_seed $RANDOM -var read_atom $read_atom -var templ $templ -var tarpress $press -var pot_file $pot_file -var element1 $element1 -var element2 $element2 -in in-liq.npt


