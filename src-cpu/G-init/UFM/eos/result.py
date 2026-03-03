#!/usr/bin/env python3
  
import glob
import json

files=sorted(glob.glob('T*/thermo.json'))

print("P T V H  E  natom n_Fe")
for f in files:
    para=json.load(open(f))
    print(para['press']*0.0001,para['temp'],para['vol'],para['enthalpy'],para['energy'],para['natom'],para['na1'])


