#!/bin/bash

#!/bin/bash
workdir=$(pwd)
x_list=(0 0.04 0.08 0.12 0.16 0.20)
for x in "${x_list[@]}"; do
    filepath=$workdir/"T=5500_p=323_liq_x=${x}"/TI
    cd $filepath
    cp $workdir/UFM/TI/integrate-liq.py .
    cp $workdir/UFM/G.ipynb ..
    python integrate-liq.py 5500 32000 $x > TI.dat
    cd ..
    ipython G.ipynb 
    cd $workdir
done
