include init.mod


variable        myenthalpy equal enthalpy
variable        mypress equal press
variable        mype equal pe
variable        myvol equal vol

compute pressk all pressure NULL virial
## liquid cal
variable dir string x
include displace.mod

## solid cal

variable dir string y
include displace.mod




variable dlambda equal ${pe_x}-${pe_y}*${x}/${y}
variable dH equal ${h_x}-${h_y}*${x}/${y}


print """{
  "pe_x": ${pe_x},
  "pe_y": ${pe_y},
  "press_x": ${press_x},
  "press_y": ${press_y},
  "vol_x": ${vol_x},
  "vol_y": ${vol_y},
  "dlambda":${dlambda},
  "x": ${x},
  "y": ${y},
  "dH": ${dH}
}""" file ${thermo_file} screen no
