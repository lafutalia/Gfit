#!/bin/bash
#SBATCH -o TI.out
#SBATCH -e TI.err
#SBATCH -J es
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1

x_list=(0 0.04 0.08 0.12 0.16 0.20)
# x_list=(0)
counter=0
LMP=$1
##每次处理三个任务
for x in "${x_list[@]}"; do
    T=5500
    p=323
    type=liq
     ./get_G.pbs $T $p $type $x $LMP &
     counter=$((counter+1))
    if [ $counter -eq 3 ]; then
        counter=0
        wait
    fi
done

python Gmix.py


