# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# we must undefine any fix ave/* fix before using reset_timestep
if "$(is_defined(fix,avall))" then "unfix avall"
reset_timestep 0 
# Choose potential
pair_style     eam/fs
pair_coeff     * * ${pot_file} ${element1} ${element2}

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup output

fix avall all ave/time 1 1 1  v_myenthalpy v_mype c_pressk v_myvol  ave running
variable hTmp equal f_avall[1]
variable peTmp equal f_avall[2]
variable pressTmp equal f_avall[3]
variable volTmp equal f_avall[4]

thermo		10000
thermo_style custom step temp pe press pxx pyy pzz pyz pxz pxy 
thermo_modify norm no





