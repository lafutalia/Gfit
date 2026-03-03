#!/bin/bash

args=("$@")

#thermojson_dir="analysis/json_dirs"    ##表示从该文件夹获取thermo文件
pair_json_file="pair-*-type-*-elements-*.json"
pair_dir="./pair_json"                   ##表示从该目录获取pair文件
database="./database"
if [ -d $database/ ]; then
rm -r $database
fi
mkdir $database
# 从命令行参数中获取所有参数

for pairfile in $pair_dir/$pair_json_file; do
  subdirs=()
  type_list=()
  filename=$(basename "$pairfile")

  prefix="${filename#pair}"           ## 去掉前缀
  echo $filename
  thermo_filename="pot${prefix}"

  for arg in "${args[@]}"; do
    file_path=${arg}/$thermo_filename
    type_list+=("$arg")

    subdirs+=("$file_path")
  done

  out_file=$database/$thermo_filename

  cut_2.py  -p $pairfile -t ${subdirs[*]} -o $out_file -n ${type_list[*]}
done


