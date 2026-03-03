#!/usr/bin/env bash

##############################################################
#add your lammps path here
LMP=lmp_mpi
##############################################################

workDir=$PWD  
echo "fitting is in $workDir" 
CMDGenPath="atom-gen"
modeTm=0      
if  [   ${modeTm}  -eq 1    ];
then
   ensembleList=(
    Gmix-liq_T=5500_p=323-0.04 Gmix-liq_T=5500_p=323-0.08
   Gmix-liq_T=5500_p=323-0.12 Gmix-liq_T=5500_p=323-0.16

   Meanf-liq_T=6000_p=323-Fe4O1)


elif [  ${modeTm}   -eq 0    ];
then
   ensembleList=(
   Gmix-liq_T=5500_p=323-0.04 Gmix-liq_T=5500_p=323-0.08
   Gmix-liq_T=5500_p=323-0.12 Gmix-liq_T=5500_p=323-0.16
   Meanf-liq_T=6000_p=323-Fe4O1
   )
fi

mkdir -p pot-ite


timeLog="$workDir/time.log"

# 可选：记录脚本启动时间
echo "=== fit.sh start: $(date '+%F %T')  workDir=$workDir ===" >> "$timeLog"


cd $workDir
i=1
potPath=$workDir/pot-ite/pot_combine*.fs
G_Path=$workDir/G-v${i}

cp -r G-init G-v${i}
cp  $potPath G-v${i}/
cd G-v${i}
./pd.sh $LMP
cd $workDir
G_Path=$workDir/G-v${i}/


echo  "fs potentail of version${i} has generated!"
echo  "******************************************"
cd  $workDir
# cp  $potPath    ./pot-ite
##fit from the CMD-atom-gen with no melting point

for i in {1..5}
do  
    iter_start_epoch=$(date +%s)
    iter_start_ts=$(date '+%F %T')
    ##gen work environment
    mkdir -p ite${i}
    mkdir -p ite${i}/rerun
    if [ ! -d "ite${i}/rerun/init" ]; then
        cp -r init/rerun/init ite${i}/rerun/init
    fi
    cp -Rf init/rerun/atom ite${i}/rerun
    cp -Rf "fit" ite${i}/rerun
    cp -Rf "target" ite${i}/rerun
    cp -r $CMDGenPath ite${i}
    cp -Rf "basis-fun"/pair_json     ite${i}/rerun/init/
    cp -Rf "basis-fun"/pot           ite${i}/rerun/init/
    
    # ## run md to get atomic structure

    cd $workDir
    cd ite${i} 
    cd "atom-gen"/"NVT"
    cp $potPath "pot_combine.fs"
    ./gen-NVT.sh $LMP >md.out

    cd $workDir
    cd ite${i} 
    cd "atom-gen"/"NPT"
    cp $potPath "pot_combine.fs"
    ./gen-NPT.sh $LMP >md.out  


    ## gen atom which can be read by rerun code
    cd $workDir
    cd ite${i}/rerun/atom
    ./gen-atom.sh          ##从xdat生成dump.atom文件

    ## rerun the property 
    cd ../init
    ./gen_property.sh   ${ensembleList[*]} >    rerun.log          ##对基函数生成所有性质
    ./get_value-v2.sh   ${ensembleList[*]} >>   rerun.log          ##gen database

    echo  "database has generated!"

    ## fit and gen  potential
    cd ../fit
    cp $G_Path/Gmix.dat .
    mkdir pair_json
    mkdir tmp
    python fit-v1-G.py  > fit.log    
    potPath=$PWD/"pot_combine.fs"                                    ## gen-pot file

    cd $workDir
    j=$((1+$i))
    cp -r G-init G-v${j}
    cp  $potPath G-v${j}/
    cd G-v${j}
    ./pd.sh $LMP
    cd $workDir
    G_Path=$workDir/G-v${j}/

    echo  "fs potentail of version${i} has generated!"
    echo  "******************************************"
    cd  $workDir
    cp  $potPath    ./pot-ite/"pot_combine-v${i}.fs"
    iter_end_epoch=$(date +%s)
    iter_end_ts=$(date '+%F %T')
    iter_elapsed=$((iter_end_epoch - iter_start_epoch))
    echo "iter=$i start=$iter_start_ts end=$iter_end_ts elapsed=${iter_elapsed}s" >> "$timeLog"
done

echo "=== fit.sh end: $(date '+%F %T') ===" >> "$timeLog"
