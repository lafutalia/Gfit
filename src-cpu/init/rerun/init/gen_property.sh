#!/bin/bash

##############################################################
#add your lammps path here
LMP=~/software/lammps-patch_21Nov2023/src/lmp_mpi
##############################################################

element1=Fe
element2=O
rerunCodeDir=../../../rerun-code/
input_atom_path=("$@")
dir_pot="./pot"
pot_file="pot-*-type-*-elements-*.fs"
prefix="pot-"
suffix=".fs"
tmp_rerun="./rerun_tmp"
work_dir=$PWD
counter=0
for data_path in "${input_atom_path[@]}"; do
  tmp_rerun=./rerun_tmp-$(basename "$data_path")
  if [ -d $tmp_rerun ]; then
    rm -r $tmp_rerun
  fi
  mkdir -p $tmp_rerun
  mkdir -p ${data_path} 
  echo "Rerun MD with ${data_path} $element1"
  atom_path="../../../atom/$(basename "$data_path")"
  option_file="../atom/$(basename "$data_path")/option.dat"
  rerun_type=$(awk '/^RERUN/ {for(i=2; i<=NF; i++) {if($i!=""){print $i; break}}}' "$option_file")
  for file in $dir_pot/$pot_file; do
    echo -n "Rerun MD with $file Fe"
    filename=$(basename "$file")
    result=${filename#"$prefix"}
    thermo_name="${prefix}${result%"$suffix"}.json"
    ##如果是rerun_type!=Meanf，还需要判断meanforce文件是否存在
    if [ -f "$data_path/$thermo_name" ] && [ "$rerun_type" != "Meanf" ]; then
    ##如果是rerun_type=Meanf，还需要判断meanforce文件是否存在
      if [ -f "$data_path/meanforce_${result%"$suffix"}.dat" ]; then
        echo "文件存在,继续计算下一个"
        continue
      fi
    else
      if [ "$rerun_type" = "Meanf" ] && [  -f "$data_path/meanforce_${result%"$suffix"}.dat" ] ; then
          echo "文件存在,继续计算下一个"
          continue
      fi
      echo "文件不存在"
      dir_tmp_rerun=$tmp_rerun/"thermo_${result%"$suffix"}"
      if [ -d $dir_tmp_rerun ]; then
        rm -rf $dir_tmp_rerun
      fi
      cp -r $rerunCodeDir/$rerun_type $dir_tmp_rerun
      
      cd $dir_tmp_rerun ||exit
      thermo_file="../../$data_path/$thermo_name"
      ./rerun.sh $thermo_file $work_dir/$file $atom_path $LMP $element1 $element2 &
      cd ../..
      counter=$((counter + 1))
    fi
    
    if [[ $counter -eq 4 ]]; then
      # 等待当前文件处理完毕
      wait

      # 重置计数器
      counter=0
    fi


  done
done
wait
echo "done"



