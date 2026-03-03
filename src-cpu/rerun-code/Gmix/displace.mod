include potential.mod

rerun           ${atomdir}/dump_${dir}.atom  first 0 every 0 &
dump x y z purge yes replace no add yes 

# Obtain new thermo data


# store data
group           a1 type 1
variable        na1 equal count(a1)
group           a2 type 2
variable        na2 equal count(a2)



variable ${dir} equal ${na2}/(${na1}+${na2})
variable h_${dir} equal ${hTmp}/$(atoms)
variable pe_${dir} equal ${peTmp}/$(atoms)
variable press_${dir} equal ${pressTmp}/10000
variable vol_${dir} equal ${volTmp}/$(atoms)





variable dir delete
variable hTmp delete
variable peTmp delete
variable pressTmp delete
variable volTmp delete