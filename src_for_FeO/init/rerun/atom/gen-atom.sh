#!/bin/bash

NVT_filepath='../../atom-gen/NVT'
G_filepath='../../atom-gen/NPT'
phase1="hcp"
phase4="liq"
phase3="bcc"

Meanf_prefix_list=("T=6000_p=323-Fe4O1") 
phase_list=($phase4)
for prefix in "${Meanf_prefix_list[@]}"; do
  for phase in "${phase_list[@]}";do
  mkdir -p "Meanf-${phase}_${prefix}"
  cp ${NVT_filepath}/"${phase}_${prefix}"/dump_all.atom ${NVT_filepath}/"${phase}_${prefix}"/read.in ./"Meanf-${phase}_${prefix}"/
  cp ${NVT_filepath}/"${phase}_${prefix}"/gr-target.dat ${NVT_filepath}/"${phase}_${prefix}"/thermo.json ${NVT_filepath}/"${phase}_${prefix}"/para.in ./"Meanf-${phase}_${prefix}"/
  echo "RERUN       Meanf" >"./Meanf-${phase}_${prefix}/option.dat"
  echo "dealing with Meanf for Meanf-${phase}_${prefix}"
    done
  done


Gmix_prefix_list=("T=5500_p=323-0.04" "T=5500_p=323-0.08" "T=5500_p=323-0.12" "T=5500_p=323-0.16" "T=5500_p=323-0.20")
Gmix_prefix_ref="T=5500_p=323-0.20"
Gmix_prefix_0="T=5500_p=323-0"
phase_list=($phase4)
for prefix in "${Gmix_prefix_list[@]}"; do
  for phase in "${phase_list[@]}";do
  mkdir -p "Gmix-${phase}_${prefix}"
  cp ${G_filepath}/"${phase}_${prefix}"/dump_all.atom ./"Gmix-${phase}_${prefix}"/dump_x.atom
  cp ${G_filepath}/"${phase}_${prefix}"/thermo.json ./"Gmix-${phase}_${prefix}"/thermo_x.json
  cp ${G_filepath}/"${phase}_${Gmix_prefix_ref}"/dump_all.atom ./"Gmix-${phase}_${prefix}"/dump_y.atom
  cp ${G_filepath}/"${phase}_${Gmix_prefix_ref}"/thermo.json  ./"Gmix-${phase}_${prefix}"/thermo_y.json
  cp ${G_filepath}/"${phase}_${Gmix_prefix_0}"/thermo.json ./"Gmix-${phase}_${prefix}"/thermo_0.json
  cp ${G_filepath}/"${phase}_${Gmix_prefix_0}"/read.in ./"Gmix-${phase}_${prefix}"/read.in
    echo "dealing with Gmix for ${phase}_${prefix}"
    echo "RERUN       Gmix" >"Gmix-${phase}_${prefix}"/option.dat
    done
  done






