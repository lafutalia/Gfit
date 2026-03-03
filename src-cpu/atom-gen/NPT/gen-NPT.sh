#!/bin/bash

LMP=$1
element1="Fe"
element2="O"




echo "gen atom ensemble for liquid hm "

counter=0
types=("liq"  "liq" "liq" "liq" "liq" "liq" )
ifisos=("iso" "iso" "iso" "iso" "iso" "iso")
temp=5500
press=3230000
x=(0 0.04 0.08 0.12 0.16 0.20)
i=0
for xx in ${x[*]};do
  type=${types[i]}
  ifiso=${ifisos[i]}
  phase_name=${type}_T=${temp}_p=$(($press/10000))-$xx
  read_atom=../init-atom/$type.atom
  rm -rf $phase_name
  cp -r init-hm $phase_name
  cd $phase_name
  ./md-npt.sh $temp $LMP $read_atom $press $element1 $element2 $xx $ifiso  &
  counter=$((counter + 1))
  cd ../
  if [[ $counter -eq 3 ]]; then
      # 等待当前文件处理完毕
      wait
      # 重置计数器
      counter=0
  fi
  i=$(($i + 1))
done



wait