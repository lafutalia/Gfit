#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import json
font = {
        'size'   : 18}
plt.rc('font', **font)


# In[ ]:


x_list=[0,0.04,0.08,0.12,0.16,0.20];F0_list=[];vol_list=[];G_list=[]
for x in x_list:
    print(x)
    if x==0:
        G_dat_file=f"T=5500_p=323_liq_x={x}/result.json"
    else:
        G_dat_file=f"T=5500_p=323_liq_x={x:.2f}/result.json"
    data_json=json.load(open(G_dat_file,'r'))
    G=data_json['G']
    G_list.append(G)
G_mix=np.array(G_list)-((1-np.array(x_list)/0.2)*G_list[0]+np.array(x_list)/0.2*G_list[-1])
plt.plot(x_list,G_mix,'o-')
np.savetxt('Gmix.dat',np.array([x_list,G_mix]).T,header='x Gmix')


# In[ ]:




