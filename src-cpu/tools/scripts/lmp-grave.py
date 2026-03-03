#!/usr/bin/env python3

import sys
import numpy as np
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-s","--show", help="plot the function", action='store_true')
args = parser.parse_args()

fin=open("gr.all","r")
npair=1

for k in range(3):
   fin.readline()
ndata=int(fin.readline().split()[1])
fin.seek(0)
for k in range(3):
   fin.readline()

gr_data=[ [0 for i in range(ndata)] for j in range(2*npair+1)]
nc=0
for line in fin:
   for i in range(ndata):
      ll=fin.readline().split()
      for j in range(2*npair+1):  # r gr Sum(gr)
         gr_data[j][i]+=float(ll[j+1])
   nc+=1
fin.close()

gdata=np.array(gr_data)/nc  # average
dr = gdata[0][1] - gdata[0][0]

print(f"gr.all read; averaged for {nc} snapshots")

# dln(gr)/dr
lngr=np.log(gdata[1,:])
#print(lngr)
lngr[lngr == -np.inf] = np.nan 
#print(lngr)
MF = np.gradient( lngr , dr)
#print(MF)

fout=open('gr.ave','w+')
print("#r gr lngr dlngr/dr",file=fout)
for i in range(ndata):
   print(gdata[0][i], gdata[1][i], lngr[i],  MF[i], file=fout)
fout.close()

if(args.show):
   import pylab as plt
   plt.figure(figsize=(10,8))
   plt.rcParams.update({'font.size': 12})
   
   plt.subplot(2,2,1)
   plt.plot(gdata[0,:], gdata[1,:], label='g(r)')
   plt.axhline(0,ls='--',c='k')
   plt.xlabel("r ($\AA$)")
   plt.xlim(0,)
   plt.legend()
   
   plt.subplot(2,2,2)
   plt.plot(gdata[0,:], -lngr, label='-ln[g(r)]')
   plt.axhline(0,ls='--',c='k')
   plt.xlabel("r ($\AA$)")
   plt.xlim(0,)
   plt.legend()
   
   plt.subplot(2,2,3)
   plt.plot(gdata[0,:], MF, label='d(ln[g(r)])/dr')
   plt.axhline(0,ls='--',c='k')
   plt.xlabel("r ($\AA$)")
   plt.xlim(0,)
   plt.legend()
   
   plt.tight_layout()
   plt.show()
