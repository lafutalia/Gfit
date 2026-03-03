#!/bin/bash

LMP=$1
element1="Fe"
element2="O"

echo "gen atom ensemble for liquid gr"



types=("liq")
temps=(6000)
presss=(3230000)
x=("Fe4O1")

counter=0
i=0
for xx in ${x[*]};do
  temp=${temps[i]}
  press=${presss[i]}
  type=${types[i]}
  phase_name=${type}_T=${temp}_p=$(($press/10000))-${x[i]}
  gr_target=../target-gr/$phase_name.gr
  read_atom=../init-atom/$phase_name.atom
  para=../target-gr/para.in
  rm -rf $phase_name
  cp -r init-liquid-gr $phase_name
  cd $phase_name
  ./md.sh $temp $LMP $gr_target $para $read_atom $element1 $element2  &
  counter=$((counter + 1))
  cd ../
  if [[ $counter -eq 4 ]]; then
      # 等待当前文件处理完毕
      wait
      # 重置计数器
      counter=0
  fi
  i=$(($i + 1))
done

wait