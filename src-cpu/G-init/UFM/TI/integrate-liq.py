#!/usr/bin/env python3

import numpy as np
import ufgenerator
import sys

# constants
kB = 8.61733034e-5 # eV/K
hbar = 6.582119569e-16 # eV.s
h = 4.13566766e-15 # eV.s
Na = 6.02214076e23
e = 1.602176634e-19


def F_idealgas (temp, rho, natom, mass, comp):
   beta=1/(kB*temp)
   Omega=np.sqrt(beta*h*h/2/np.pi/mass * 1e3 * e * Na)*1e10   # angstrom
   FE = 0
   for k in range(len(mass)):
      if(comp[k]==0):
         FE += 0
      else:
         FE += comp[k]* ( 3*np.log(Omega[k]) + np.log(rho) -1 + np.log(comp[k]) )
#   FE += 1/(2*natom)*np.log(2*np.pi*natom)
   return 1/beta*FE


def F_UF (temp, rho, p, sigma):
   x = (0.5*(np.pi*sigma*sigma)**1.5)*rho
   press,fe = ufgenerator.get_UF(p,x)
   beta = 1/(kB*temp)
   return fe/beta # eV/atom

# parameters
temp = float(sys.argv[1])    # K
mass = np.array([55.845,15.999]) # g/mol
p = 50
sigma = 1.5

# read lammps output for natom and vol
natom = int(sys.argv[2])
vol = -999  # A^3 
fin=open("md.out","r")
for line in fin:
   if("summary:" in line):
      ll=line.split()
      vol=float(ll[1])  # A^3
fin.close()

comp0 = float(sys.argv[3])  # 0<x<1
comp = np.array([1-comp0, comp0])  # 0<x<1
rho = natom/vol
#print(comp)

data_f=np.loadtxt("forward.dat",skiprows=1)
data_b=np.loadtxt("backward.dat",skiprows=1)

U_eam_f=data_f[:,0]
U_mod_f=data_f[:,1]
lmbda_f=data_f[:,2]

U_eam_b=data_b[:,0]
U_mod_b=data_b[:,1]
lmbda_b=data_b[:,2]

dU_f=U_eam_f - U_mod_f
dU_b=U_eam_b - U_mod_b
W_if = np.trapz(dU_f,lmbda_f)
W_fi = np.trapz(dU_b,lmbda_b)

W = 0.5*(W_if - W_fi)
F_ig = F_idealgas (temp, rho, natom, mass, comp)
F_ufm = F_UF(temp, rho, p, sigma)
F = F_ig + F_ufm  - W  # forward 1->0

import pylab as plt
fig=plt.figure(figsize=(8,6))
plt.rcParams.update({'font.size': 16})

plt.plot(lmbda_f, dU_f, label='Forward')
plt.plot(lmbda_b, dU_b, label='Backward')
plt.xlabel(r"$\lambda$")
plt.ylabel(r"$\Delta U=U_{eam}-U_{mod}$")
plt.legend()
plt.text(0.4, max(dU_f)-0.25*(max(dU_f)-min(dU_f)), 'W:%.6f\nF_ref:%.6f\nF0:%.6f'%(W,F_ig+F_ufm,F))
plt.tight_layout()
plt.savefig("intg.png")
#plt.show()

#print(F)
print("x in Fe1-xOx = ", comp0)
print("rho:", rho)
print("W:", W)
print("F_ig:", F_ig)
print("F_UF:", F_ufm)
print("F_ref:", F_ig+F_ufm)
print("F0:", F)
