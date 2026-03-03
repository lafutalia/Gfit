#!/usr/bin/env python3

import sys
import numpy as np
import os
import datetime
import json
import argparse
import glob
import matplotlib.pylab as plt

parser = argparse.ArgumentParser()
parser.add_argument("-t","--target", help="target gr",action='store')
parser.add_argument("-i","--input",  help="input gr's",action='store') ## read gr from md out
args = parser.parse_args()
with open(args.input,'r') as f:
    gr_raw=np.loadtxt(f,skiprows=4)     #

## pair number
pair_number=int(gr_raw.shape[1]/2-1)  ## pair number    2*npair+1
r=gr_raw[:,1]                    ## r
gr_now=gr_raw[:,2::2]            ##all gr file


## 要对关联的值进行平均


with open(args.target,'r') as f:
    target_raw=np.loadtxt(f,skiprows=1)
target=target_raw[:,1:]
r_target=target_raw[:,0]
plt.figure(figsize=(20,10))
for i in range(gr_now.shape[1]):
    msq = np.sqrt( np.sum( (target[:,i]-gr_now[:,i])**2 )/target.shape[0])
    plt.subplot(1,pair_number,i+1)
    plt.scatter(r,gr_now[:,i], label=f'{i}, Rg=%.5f'%(msq),c='r')
    plt.plot(r_target,target[:,i],label='target')
    print("--------------------------------")
    print(f"Comparing target and input gr")
    plt.rcParams.update({'font.size': 16})
    plt.xlabel("r ($\AA$)")
    plt.ylabel("gr")
    plt.ylim([-0.2,3.3])
    plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig("gr-diff.png")
# plt.show()

