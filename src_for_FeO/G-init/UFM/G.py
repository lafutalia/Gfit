#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import json
font = {
        'size'   : 18}
plt.rc('font', **font)


# In[2]:


TI_dat_file=f"TI/TI.dat"
##读取F0所在行的值
with open(TI_dat_file) as f:
    for line in f:
        if 'F0' in line:
            F0=float(line.split()[1])
            break
eos_json=f"eos/thermo.json"
with open(eos_json) as f:
    eos_json=json.load(f)
vol=eos_json['vol'];natom=eos_json['natom'];press=eos_json['press']/10000
temp=eos_json['temp']
vol_per_atom=vol/natom
print(F0,vol_per_atom)
G=F0+press*vol_per_atom*0.006241509074
output = {
    "temp": temp,
    "press": press,
    "G": G
}
with open("result.json", "w") as f:
    json.dump(output, f, indent=4)

