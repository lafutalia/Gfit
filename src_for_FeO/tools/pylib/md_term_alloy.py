#!/usr/bin/env python3

import numpy as np
import os
import json
import numpy as np
from copy import deepcopy
from sklearn.linear_model import LinearRegression, Ridge
import datetime
import matplotlib.pylab as plt
import basis_fun
import importlib
importlib.reload(basis_fun)
from basis_fun import gen_pot_fun
from basis_fun import cut_off
import re
import random
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
import numpy as np 
from scipy.optimize import curve_fit
## 定义一个元素类，以拓展到更多元素的fit
import numpy as np
from scipy.optimize import curve_fit

kB = 8.61733326E-5  # Boltzmann 常数

def fit_Gx_single(x_data, a1, a2,T):
    """
      G_mix(x) = G + a0*x + a1*x*(1-x) + a2*x*(1-x)*(2*x-1) + T*kB*( x*ln(x) + (1-x)*ln(1-x) )
    其中:
      G = 0
      0.2*a0 + 0.16*a1 - 0.096*a2 + T*kB*[0.2*ln(0.2) + 0.8*ln(0.8)] = 0

    => a0 = ( -0.16*a1 + 0.096*a2 - c ) / 0.2
    其中 c = T*kB*(0.2*ln(0.2) + 0.8*ln(0.8))
    """
    # 计算常数 c
    c = T * kB * (0.2 * np.log(0.2) + 0.8 * np.log(0.8))
    # 由约束解出 a0
    a0 = (-0.16 * a1 + 0.096 * a2 - c) / 0.2

    # 计算理想混合项
    # 注意 x=0 时 x*log(x) -> 0, 本示例中假设 x_data 都大于 0
    ideal_mix = T * kB * (x_data * np.log(x_data) + (1 - x_data) * np.log(1 - x_data))

    # 计算 G_mix(x)
    G_mix = a0 * x_data \
          + a1 * x_data * (1 - x_data) \
          + a2 * x_data * (1 - x_data) * (2 * x_data - 1) \
          + ideal_mix
    return G_mix

# def Gx_para(x_data,y_data,T_data):
#     combined_data = (x_data, T_data)

#     params, cov = curve_fit(fit_Gx_single, combined_data, y_data)

#     a1_fit, a2_fit = params
#     print("拟合结果: a1 =", a1_fit, ", a2 =", a2_fit)

#     c = T_data * kB * (0.2*np.log(0.2) + 0.8*np.log(0.8))
#     a0_fit = (-0.16*a1_fit + 0.096*a2_fit - c) / 0.2
#     return a0_fit, a1_fit, a2_fit

##根据a0,a1,a2计算Gx-Gx'*x的值
# def G_xGp(para,x):
#     a0,a1,a2,T=para
#     Gx=a0*x+a1*x*(1-x)+a2*x*(1-x)*(2*x-1)+T*kB*(x*np.log(x)+(1-x)*np.log(1-x))
#     Gxpx=a0+a1*(1-2*x)+a2*((1-x)*(2*x-1)+x*(1-2*x)+2*x*(1-x))+T*kB*(np.log(x)-np.log(1-x))
#     return Gx-Gxpx*x


def F_rho_basis(x_rho):
    F_rho=-np.sqrt(x_rho)
    return F_rho
class Elements():
    def __init__(self,tp,mass,num,latt_type,latt_a):
        self.tp=tp
        self.mass=mass
        self.num=num
        self.latt_type=latt_type
        self.latt_a=latt_a
        ## 需要定义term_type


    @classmethod
    def load_from_para(cls,para_dict):
        #从para.json加载元素
        ele_list=[]
        for i,element in enumerate(para_dict['elements']):
            tp=element["type"]
            mass=element["mass"]
            num=element["atom_num"]
            latt_type=element["latt_type"]
            latt_a=element["latt_a"]
            ele_list+=[Elements(tp=tp,mass=mass,num=num,latt_type=latt_type,latt_a=latt_a)]
        Nelements=len(ele_list)
        return Nelements,ele_list


def B2_T(x_r,u):
    plt.figure(figsize=(10,10))
    for T in [5500,6500,7500,8500,9500]:
        kB=8.61733326E-5
        beta=T*kB
        inner=x_r**2*(np.exp(-beta*u)-1)
        B2T=-2*np.pi*np.trapz(x_r,inner)
        rho=np.linspace(5,8,1000)
        plt.plot(1/rho,rho+B2T*rho**2)
        plt.savefig("P-V-T")
def pot_writer(pair='pair.json', para='para.json', out_name='pot.fs',mode=None,if_fit=None,pair_cut=[1.5]):
## 加载元素信息
    with open(para,'r') as f:
        data=json.load(f)
    Nelements,ele_list=Elements.load_from_para(data)
    Nrho = data["Nrho"]
    drho = data["drho"]
    Nr = data["Nr"]
    dr = data["dr"]
    cut = data["cutoff"]
    x_rho=np.array([drho * k for k in range(Nrho)])
    x_r = np.array([dr * k for k in range(Nr)])
    path_now = os.getcwd()
    pot_name = out_name
    pot_dir = os.path.join(path_now)
    path_pot = os.path.join(path_now, pot_name)
    fout = open(path_pot, "w+")
    print("FS potential", file=fout)
    date = datetime.date.today().strftime("%m/%d/%y")
    current_time = datetime.datetime.now().strftime("%H:%M:%S")
    print("Generated on", date, current_time, file=fout)
    print(f"Path:{path_pot}", file=fout)
    # L.4 Element info
    print(Nelements,*[ele.tp for ele in ele_list], file=fout)
    print("%8d %24.14E %8d %24.14E %24.14E" % (Nrho, drho, Nr, dr, cut), file=fout)

    with open(pair,'r') as f:
        pot=json.load(f)

    if if_fit is not None:
        pot=pot['pair']

    pair_info=pot['pair']
    emb_info=pot['emb']
    rho_info=pot['rho']

    ## 写入F(rho)
    for i,elem_i in enumerate(ele_list):
        print("%3d %24.14E %24.14E %6s" % (elem_i.num, elem_i.mass, elem_i.latt_a, elem_i.latt_type),
            file=fout)
        index=1
        emb_dict=emb_info[f'{i+1}']
        F=gen_pot_fun(x_rho,emb_dict)
        cut_off_F=cut_off(x_rho,left_value=None,type=f"emb-{index}")
        # print(f"writting embedding term of atom{i+1}")
        if emb_dict['need_c']==1:
            F=F+cut_off_F
        # print(F)
        for count, number in enumerate(F):
            print("%24.14E" % (number), end=" ", file=fout)
            count += 1
            if (count % 5 == 0):
                print("", file=fout)

        ##写入rho
        for j in range(Nelements):
            rho_dict=rho_info[f'{j+1}{i+1}']
            rho=gen_pot_fun(x_r,rho_dict)
            # print(f"writting rho term of atom pair {j+1}{i+1}")
            for count, number in enumerate(rho):
                print("%24.14E" % (number), end=" ", file=fout)
                count += 1
                if (count % 5 == 0):
                    print("", file=fout)
        ##写入rphi
    count_pair_cut=0
    for i,elem_i in enumerate(ele_list):
        for j in range(0,i+1):
            index=1
            print(f"writting pair term of atom pair {i+1}{j+1}")
            pair_dict=pair_info[f'{i+1}{j+1}']
            left_value_prime,left_value,phi=gen_pot_fun(x_r,pair_dict,need_left=True,need_left_prime=True)
            # print(left_value,left_value_prime)
            # print(f"writting pair term of atom pair {j+1}{i+1}")
            if pair_dict['need_c']==1:
                pair_cut_i=pair_cut[count_pair_cut]
                cutoff_phi=cut_off(x_r,left_value=left_value,left_value_prime=left_value_prime,type=f'pair-{index}',innercut=pair_cut_i)
                phi=phi+cutoff_phi
            phi[0]=1e133
            # B2T=B2_T(x_r,phi)
            for count, number in enumerate(x_r*phi):
                print("%24.14E" % (number), end=" ", file=fout)
                count += 1
                if (count % 5 == 0):
                    print("", file=fout)
            count_pair_cut+=1
    fout.close()  
def pot_writer_alloy(pair='pair.json', para='para.json', out_name='pot.fs',mode=None,if_fit=None,pair_cut=[1.5],if_cd=None):
## 加载元素信息
    with open(para,'r') as f:
        data=json.load(f)
    Nelements,ele_list=Elements.load_from_para(data)
    Nrho = data["Nrho"]
    drho = data["drho"]
    Nr = data["Nr"]
    dr = data["dr"]
    cut = data["cutoff"]
    x_rho=np.array([drho * k for k in range(Nrho)])
    x_r = np.array([dr * k for k in range(Nr)])
    path_now = os.getcwd()
    pot_name = out_name
    pot_dir = os.path.join(path_now)
    path_pot = os.path.join(path_now, pot_name)
    fout = open(path_pot, "w+")
    print("FS potential", file=fout)
    date = datetime.date.today().strftime("%m/%d/%y")
    current_time = datetime.datetime.now().strftime("%H:%M:%S")
    print("Generated on", date, current_time, file=fout)
    print(f"Path:{path_pot}", file=fout)
    # L.4 Element info
    print(Nelements,*[ele.tp for ele in ele_list], file=fout)
    print("%8d %24.14E %8d %24.14E %24.14E" % (Nrho, drho, Nr, dr, cut), file=fout)

    with open(pair,'r') as f:
        pot=json.load(f)

    if if_fit is not None:
        pot=pot['pair']

    pair_info=pot['pair']
    emb_info=pot['emb']
    rho_info=pot['rho']
    
    ## 写入F(rho)
    for i,elem_i in enumerate(ele_list):
        print("%3d %24.14E %24.14E %6s" % (elem_i.num, elem_i.mass, elem_i.latt_a, elem_i.latt_type),
            file=fout)
        index=1
        emb_dict=emb_info[f'{i+1}']
        F=gen_pot_fun(x_rho,emb_dict)
        cut_off_F=cut_off(x_rho,left_value=None,type=f"emb-{index}")
        # print(f"writting embedding term of atom{i+1}")
        if emb_dict['need_c']==1:
            F=F+cut_off_F
        # print(F)
        for count, number in enumerate(F):
            print("%24.14E" % (number), end=" ", file=fout)
            count += 1
            if (count % 5 == 0):
                print("", file=fout)

        ##写入rho;对于eam/alloy类型的势函数。需要调用rho_info的i,i
        
        rho_dict=rho_info[f'{i+1}{i+1}']
        rho=gen_pot_fun(x_r,rho_dict)
        # print(f"writting rho term of atom pair {j+1}{i+1}")
        for count, number in enumerate(rho):
            print("%24.14E" % (number), end=" ", file=fout)
            count += 1
            if (count % 5 == 0):
                print("", file=fout)
        ##写入rphi
    count_pair_cut=0
    for i,elem_i in enumerate(ele_list):
        for j in range(0,i):
            index=1
            pair_dict=pair_info[f'{i+1}{j+1}']
            left_value_prime,left_value,phi=gen_pot_fun(x_r,pair_dict,need_left=True,need_left_prime=True)
            # print(left_value,left_value_prime)
            # print(f"writting pair term of atom pair {j+1}{i+1}")
            if pair_dict['need_c']==1:
                pair_cut_i=pair_cut[count_pair_cut]
                cutoff_phi=cut_off(x_r,left_value=left_value,left_value_prime=left_value_prime,type=f'pair-{index}',innercut=pair_cut_i)
                phi=phi+cutoff_phi
            phi[0]=1e133
            # B2T=B2_T(x_r,phi)
            for count, number in enumerate(x_r*phi):
                print("%24.14E" % (number), end=" ", file=fout)
                count += 1
                if (count % 5 == 0):
                    print("", file=fout)
            count_pair_cut+=1
    ##写入h
    ## 打印注释行
    if if_cd is not None:
        h_info=pot['h']
        N=h_info['N']
        coeffs=h_info['coeffs']
        if len(coeffs) != N+1:
            raise ValueError(f"Length of coeffs should be {N+1}, but got {len(coeffs)}")
        print(f"# h(x) function, coefficients from 0 to {N} th order.", file=fout)
        print("%d" % N, end=" ", file=fout)
        for count, number in enumerate(coeffs):
            print("%24.14E" % (number), end=" ", file=fout)
        fout.close() 
    else:
        fout.close()
def pot_writer_rho(pair='pair.json', para='para.json', out_name='pot.fs',mode=None,if_fit=None,pair_cut=1.5):
## 加载元素信息
    with open(para,'r') as f:
        data=json.load(f)
    Nelements,ele_list=Elements.load_from_para(data)
    Nrho = data["Nrho"]
    drho = data["drho"]
    Nr = data["Nr"]
    dr = data["dr"]
    cut = data["cutoff"]
    x_rho=np.array([drho * k for k in range(Nrho)])
    x_r = np.array([dr * k for k in range(Nr)])
    path_now = os.getcwd()
    pot_name = out_name
    pot_dir = os.path.join(path_now)
    path_pot = os.path.join(path_now, pot_name)
    fout = open(path_pot, "w+")
    print("FS potential", file=fout)
    date = datetime.date.today().strftime("%m/%d/%y")
    current_time = datetime.datetime.now().strftime("%H:%M:%S")
    print("Generated on", date, current_time, file=fout)
    print(f"Path:{path_pot}", file=fout)
    # L.4 Element info
    print(Nelements,*[ele.tp for ele in ele_list], file=fout)
    print("%8d %24.14E %8d %24.14E %24.14E" % (Nrho, drho, Nr, dr, cut), file=fout)

    with open(pair,'r') as f:
        pot=json.load(f)

    if if_fit is not None:
        pot=pot['pair']

    pair_info=pot['pair']
    emb_info=pot['emb']
    rho_info=pot['rho']

    ## 写入F(rho)
    rho_save={}
    for i,elem_i in enumerate(ele_list):
        print("%3d %24.14E %24.14E %6s" % (elem_i.num, elem_i.mass, elem_i.latt_a, elem_i.latt_type),
            file=fout)
        index=1
        emb_dict=emb_info[f'{i+1}']
        F=gen_pot_fun(x_rho,emb_dict)
        cut_off_F=cut_off(x_rho,left_value=None,type=f"emb-{index}")
        # print(f"writting embedding term of atom{i+1}")
        if emb_dict['need_c']==1:
            F=F+cut_off_F
        # print(F)
        for count, number in enumerate(np.zeros_like(F)):
            print("%24.14E" % (number), end=" ", file=fout)
            count += 1
            if (count % 5 == 0):
                print("", file=fout)

        ##写入rho
        for j in range(Nelements):
            rho_dict=rho_info[f'{j+1}{i+1}']
            rho=gen_pot_fun(x_r,rho_dict)
            rho_save[f'{j+1}{i+1}']=rho
            # print(f"writting rho term of atom pair {j+1}{i+1}")
            for count, number in enumerate(rho):
                print("%24.14E" % (number), end=" ", file=fout)
                count += 1
                if (count % 5 == 0):
                    print("", file=fout)
        ##写入rphi
        # print(rho_save.keys())
    for i,elem_i in enumerate(ele_list):
        for j in range(i,Nelements):
            index=1
            pair_dict=pair_info[f'{j+1}{i+1}']
            left_value_prime,left_value,phi=gen_pot_fun(x_r,pair_dict,need_left=True,need_left_prime=True)
            # print(left_value,left_value_prime)
            # print(f"writting pair term of atom pair {j+1}{i+1}")
            if pair_dict['need_c']==1:
                cutoff_phi=cut_off(x_r,left_value=left_value,left_value_prime=left_value_prime,type=f'pair-{index}',innercut=pair_cut)
                phi=phi+cutoff_phi
            phi[0]=1e133
            # B2T=B2_T(x_r,phi)

            for count, number in enumerate(x_r*rho_save[f'{j+1}{i+1}']):
                print("%24.14E" % (number), end=" ", file=fout)
                count += 1
                if (count % 5 == 0):
                    print("", file=fout)
    return rho_save
    fout.close()  
def pot_writer_rho2(pair='pair.json', para='para.json', out_name='pot.fs',mode=None,if_fit=None,pair_cut=1.5):
## 加载元素信息
    with open(para,'r') as f:
        data=json.load(f)
    Nelements,ele_list=Elements.load_from_para(data)
    Nrho = data["Nrho"]
    drho = data["drho"]
    Nr = data["Nr"]
    dr = data["dr"]
    cut = data["cutoff"]
    x_rho=np.array([drho * k for k in range(Nrho)])
    x_r = np.array([dr * k for k in range(Nr)])
    path_now = os.getcwd()
    pot_name = out_name
    pot_dir = os.path.join(path_now)
    path_pot = os.path.join(path_now, pot_name)
    fout = open(path_pot, "w+")
    print("FS potential", file=fout)
    date = datetime.date.today().strftime("%m/%d/%y")
    current_time = datetime.datetime.now().strftime("%H:%M:%S")
    print("Generated on", date, current_time, file=fout)
    print(f"Path:{path_pot}", file=fout)
    # L.4 Element info
    print(Nelements,*[ele.tp for ele in ele_list], file=fout)
    print("%8d %24.14E %8d %24.14E %24.14E" % (Nrho, drho, Nr, dr, cut), file=fout)

    with open(pair,'r') as f:
        pot=json.load(f)

    if if_fit is not None:
        pot=pot['pair']

    pair_info=pot['pair']
    emb_info=pot['emb']
    rho_info=pot['rho']

    ## 写入F(rho)
    for i,elem_i in enumerate(ele_list):
        print("%3d %24.14E %24.14E %6s" % (elem_i.num, elem_i.mass, elem_i.latt_a, elem_i.latt_type),
            file=fout)
        index=1
        emb_dict=emb_info[f'{i+1}']
        F=gen_pot_fun(x_rho,emb_dict)
        cut_off_F=cut_off(x_rho,left_value=None,type=f"emb-{index}")
        # print(f"writting embedding term of atom{i+1}")
        if emb_dict['need_c']==1:
            F=F+cut_off_F
        # print(F)
        for count, number in enumerate(x_rho):
            print("%24.14E" % (number), end=" ", file=fout)
            count += 1
            if (count % 5 == 0):
                print("", file=fout)

        ##写入rho
        rho_save={}
        for j in range(Nelements):
            rho_dict=rho_info[f'{j+1}{i+1}']
            rho=gen_pot_fun(x_r,rho_dict)
            rho_save[f'{j+1}{i+1}']=rho
            # print(f"writting rho term of atom pair {j+1}{i+1}")
            for count, number in enumerate(rho):
                print("%24.14E" % (number), end=" ", file=fout)
                count += 1
                if (count % 5 == 0):
                    print("", file=fout)
        ##写入rphi
    for i,elem_i in enumerate(ele_list):
        for j in range(i,Nelements):
            index=1
            pair_dict=pair_info[f'{j+1}{i+1}']
            left_value_prime,left_value,phi=gen_pot_fun(x_r,pair_dict,need_left=True,need_left_prime=True)
            # print(left_value,left_value_prime)
            # print(f"writting pair term of atom pair {j+1}{i+1}")
            if pair_dict['need_c']==1:
                cutoff_phi=cut_off(x_r,left_value=left_value,left_value_prime=left_value_prime,type=f'pair-{index}',innercut=pair_cut)
                phi=phi+cutoff_phi
            phi[0]=1e133
            # B2T=B2_T(x_r,phi)

            for count, number in enumerate(0*x_r*rho_save[f'{j+1}{i+1}']):
                print("%24.14E" % (number), end=" ", file=fout)
                count += 1
                if (count % 5 == 0):
                    print("", file=fout)
    return rho_save




def plot_fun(pot,pltType='all',select=None,if_alloy=None):
    """
    单纯通过pot文件画出两个元素之间的所有EAM相互作用
    """
    def read_data(ndata,fin):
        ct=0; data=[]
        for line in fin:
            ll=line.split()
            ct+=len(ll)
            for d in ll:
                data.append(float(d))
            if(ct==ndata):
                break
            elif(ct>ndata):
                print("Error: read too many data; check your potential file!")
        return data
    fin=open(pot,'r')


    # L.1-3 header
    for k in range(3):
        print(fin.readline().strip('\n'))
    
    # L.4 elements
    ll=fin.readline().split()
    Nelem=int(ll[0])
    elem=[ll[k] for k in range(1,Nelem+1)]
    # print("Elements:",end=" ")
    fin.close

    # L.5
    ll=fin.readline().split()
    Nrho=int(ll[0]); drho=float(ll[1])
    Nr=int(ll[2]); dr=float(ll[3])
    cutoff=float(ll[4])
    # print("  Nrho     drho     Nr       dr   cutoff")
    # print("%6d %8.5f %6d %8.5f %8.3f"%(Nrho, drho, Nr, dr, cutoff))
    # print("------------------")

    x_rho=[drho*k for k in range(Nrho)]
    x_r=[dr*k for k in range(Nr)]
    leng=[]
    if pltType=='all':
        for i in range(Nelem):
            fin.readline()
            func_F=read_data(Nrho,fin)
            
            # print F(rho)
            plt.figure(1)
            # plt.subplot(1,Nelem,i+1)
            plt.plot(x_rho,func_F,label=f"F_{elem[i]}")
            plt.legend()
            plt.xlabel("rho"); plt.ylabel("F(rho)")
            plt.ylim(-100,   100) 
            plt.xlim(0,200)     
            # print rho(r)
            if if_alloy is not None:
                rho=read_data(Nr,fin)
                plt.figure(2)
                plt.plot(x_r, rho,label=f"rho, {elem[i]}-{elem[i]}")
                plt.xlabel("r"); plt.ylabel("rho(r)")
                leng+=[f"rho, {elem[i]}-{elem[i]}"]
                plt.ylim(0,30)
                plt.xlim(0,6)
                
            plt.legend()
            if if_alloy is not None:
                continue
            for j in range(Nelem):
                rho=read_data(Nr,fin)
                plt.figure(2)
                plt.plot(x_r, rho,label=f"rho, {elem[j]}-{elem[i]}")
                plt.xlabel("r"); plt.ylabel("rho(r)")
                leng+=[f"rho, {elem[j]}-{elem[i]}"]
                plt.ylim(0,30)
                plt.xlim(0,6)
            plt.legend(leng)
            # print phi
        for i in range(Nelem):
            for j in range(max(0,i-1),i+1):
                print
                phi=read_data(Nr,fin)
                for k in range(Nr):
                    if(x_r[k]!=0):
                        phi[k]=phi[k]/x_r[k]  # lammps use r*phi
                    # phi[x_r==0]=1e13
                plt.figure(3)
                plt.ylim(-50,100)
                plt.xlim(0,6)
                plt.plot(x_r,phi,label=f"phi, {elem[i]}-{elem[j]}")
                plt.ylim(-50,100)
                plt.xlabel("r"); plt.ylabel("phi(r)")
                plt.legend()
                print()
                
                
    if pltType=='emb':
        for i in range(Nelem):
            fin.readline()
            func_F=read_data(Nrho,fin)
            
            # print F(rho)
            plt.figure(1)
            # plt.subplot(1,Nelem,i+1)
            plt.plot(x_rho,func_F,label=f"F_{elem[i]}")
            plt.legend()
            plt.xlabel("rho"); plt.ylabel("F(rho)")
            plt.ylim(-100,   100) 
            plt.xlim(0,200)     
            # print rho(r)
            if if_alloy is not None:
                rho=read_data(Nr,fin)
                plt.figure(2)
                plt.plot(x_r, rho,label=f"rho, {elem[i]}-{elem[i]}")
                plt.xlabel("r"); plt.ylabel("rho(r)")
                leng+=[f"rho, {elem[i]}-{elem[i]}"]
                plt.ylim(0,30)
                plt.xlim(0,6)
                plt.legend()
                continue
            for j in range(Nelem):
                rho=read_data(Nr,fin)
                # plt.figure(2)
                # plt.plot(x_r, rho,label=f"rho, {elem[j]}-{elem[i]}")
                # plt.xlabel("r"); plt.ylabel("rho(r)")
                leng+=[f"rho, {elem[j]}-{elem[i]}"]
                # plt.ylim(0,30)
                # plt.xlim(0,6)
            # plt.legend(leng)
            # print phi
        for i in range(Nelem):
            for j in range(max(0,i-1),i+1):
                phi=read_data(Nr,fin)
                for k in range(Nr):
                    if(x_r[k]!=0):
                        phi[k]=phi[k]/x_r[k]  # lammps use r*phi
                    # phi[x_r==0]=1e13
                # plt.figure(3)
                # plt.ylim(-50,100)
                # plt.xlim(0,6)
                # plt.plot(x_r,phi,label=f"phi, {elem[i]}-{elem[j]}")
                # plt.xlabel("r"); plt.ylabel("phi(r)")
                # plt.legend()
                print()
    if pltType=='pair':
        for i in range(Nelem):
            fin.readline()
            func_F=read_data(Nrho,fin)
            
            # print F(rho)
            # plt.figure(1)
            # plt.subplot(1,Nelem,i+1)
            # plt.plot(x_rho,func_F,label=f"F_{elem[i]}")
            # plt.legend()
            # plt.xlabel("rho"); plt.ylabel("F(rho)")
            # plt.ylim(-100,   100) 
            # plt.xlim(0,200)     
            # print rho(r)
            if if_alloy is not None:
                rho=read_data(Nr,fin)
                plt.figure(2)
                plt.plot(x_r, rho,label=f"rho, {elem[i]}-{elem[i]}")
                plt.xlabel("r"); plt.ylabel("rho(r)")
                leng+=[f"rho, {elem[i]}-{elem[i]}"]
                plt.ylim(0,30)
                plt.xlim(0,6)
                plt.legend()
                continue
            for j in range(Nelem):
                rho=read_data(Nr,fin)
            #     plt.figure(2)
            #     plt.plot(x_r, rho,label=f"rho, {elem[j]}-{elem[i]}")
            #     plt.xlabel("r"); plt.ylabel("rho(r)")
            #     leng+=[f"rho, {elem[j]}-{elem[i]}"]
            #     # plt.ylim(0,30)
            #     plt.xlim(0,6)
            # plt.legend(leng)
            # print phi
        for i in range(Nelem):
            for j in range(max(0,i-1),i+1):
                phi=read_data(Nr,fin)
                for k in range(Nr):
                    if(x_r[k]!=0):
                        phi[k]=phi[k]/x_r[k]  # lammps use r*phi
                    # phi[x_r==0]=1e13
                plt.figure(3)
                plt.ylim(-3,10)
                plt.xlim(0,6)
                if select is not None:
                    # print(f"{i+1}{j+1}",select)
                    if f"{i+1}{j+1}" == select:
                        plt.plot(x_r,phi,label=f"phi, {elem[i]}-{elem[j]}")
                        plt.xlabel("r"); plt.ylabel("phi(r)")
                        plt.legend()
                        print()
                else:
                    plt.plot(x_r,phi,label=f"phi, {elem[i]}-{elem[j]}")
                    plt.xlabel("r"); plt.ylabel("phi(r)")
                    plt.legend()
                    print()
    if pltType=='rho':
        for i in range(Nelem):
            fin.readline()
            func_F=read_data(Nrho,fin)
            
            # print F(rho)
            # plt.figure(1)
            # plt.subplot(1,Nelem,i+1)
            # plt.plot(x_rho,func_F,label=f"F_{elem[i]}")
            # plt.legend()
            # plt.xlabel("rho"); plt.ylabel("F(rho)")
            # plt.ylim(-20,100)    
            # print rho(r)
            if if_alloy is not None:
                rho=read_data(Nr,fin)
                plt.figure(2)
                plt.plot(x_r, rho,label=f"rho, {elem[i]}-{elem[i]}")
                plt.xlabel("r"); plt.ylabel("rho(r)")
                leng+=[f"rho, {elem[i]}-{elem[i]}"]
                plt.ylim(0,30)
                plt.xlim(0,6)
                plt.legend()
                continue            
            for j in range(Nelem):
                rho=read_data(Nr,fin)
                plt.figure(2)
                plt.plot(x_r, rho,label=f"rho, {elem[j]}-{elem[i]}")
                plt.xlabel("r"); plt.ylabel("rho(r)")
                leng+=[f"rho, {elem[j]}-{elem[i]}"]
                plt.ylim(0,30)
            plt.legend(leng)
            # print phi
        for i in range(Nelem):
            for j in range(max(0,i-1),i+1):
                phi=read_data(Nr,fin)
                for k in range(Nr):
                    if(x_r[k]!=0):
                        phi[k]=phi[k]/x_r[k]  # lammps use r*phi
                    phi[x_r==0]=1e13
                # plt.figure(3)
                # plt.plot(x_r,phi,label=f"phi, {elem[i]}-{elem[j]}")
                # plt.xlabel("r"); plt.ylabel("phi(r)")
                # plt.legend()
                # plt.ylim(-2,5)    
def divide_random(list,num):
    sampled_list=random.sample(list,num)
    not_sampled_list = [element for element in list if element not in sampled_list]
    return sampled_list,not_sampled_list
def scaled_MD_data(X_martix,max=350,min=0):
    scaler=(X_martix[:,0]/max).reshape(-1,1)
    return X_martix/scaler,scaler.ravel()
def read_json_file(file,target_json):       ##后处理动态数据
    """

    :param file: 读取文件
    :param A: targat_json，根据要拟合的目标来读取MD数据
    :return:
    """
    with open(target_json,'r') as f:
        A=json.load(f)
    with open(file) as f:
        B=json.load(f)
    data,md_array=read_from_template(B,A)
    # if data['pair']['type']!="combination":
    #     r0_n = [data['pair']['values']['r0'], data['pair']['values']['n']]
    # if data['pair']['type']=="combination":
    #     r0_n=[data['pair']]

    ##
    return md_array,{'pair':B['pair'],'data':data}


def read_from_template(B,A):
    """
    :param source_dict:
    :param template_dict:
    :return:
    """
    source_dict=B['data']
    template_dict=A['data']
    result=copy_dict(template_dict,source_dict)
    array=get_array(result,[])
    result['pair']=B['pair']
    return result,array

def copy_dict(template_dict,source_dict):
    result = {}
    for key, value in template_dict.items():
        if isinstance(value, dict):
            result[key] = copy_dict(value, source_dict[key])
        else:
            result[key] = source_dict[key]
    return result
def get_array(input_json,out_put):
    for key,value in input_json.items():
        if isinstance(value, dict):
            # print(f"key '{key}' is a dict")
            get_array(value,out_put)
        else:
            # print(f"key '{key}' is a value")
            out_put+=[value]
    return out_put

def get_dict(my_dict,my_list):
    my_tuples = [(k, v) for k, v in my_dict.items()]

    # 使用列表中的值按顺序替换元组中的键
    for i, tup in enumerate(my_tuples):
        if isinstance(tup[1], dict):
            # 如果值是字典，则递归替换其键
            my_tuples[i] = (tup[0], get_dict(tup[1], my_list))
        else:
            # 否则，只需将其键替换为列表中的下一个值
            my_tuples[i] = (tup[0], my_list.pop(0))
    out_put_dict = dict(my_tuples)

    return out_put_dict


class MD_term():
    data_source="lammps"
    def __init__(self,term_dict,md_dict={},file_name="test.json"):   ## 初始化一个MD单项的类

        self.init_property = md_dict
        for key, value in md_dict.items():
            setattr(self, key, value)
        self.name=MD_term.get_basis_name(file_name)
        self.term=term_dict
        self.interact=['pair','emb']
        self.fix=['rho']
        self.addTermType()
    def addTermType(self):
        interact_list=self.interact
        type_count=[]
        for interact_type in interact_list:
            count=0
            term_interact=self.term[interact_type]
            for element,term_info in term_interact.items():
                count+=len(term_interact[element]['term'])
                count+=term_interact[element]['need_c']
            if count:
                type_count+=[1]
                self.type=interact_type
            elif count==0:
                type_count+=[0]

        if sum(type_count)>=2:
            self.type='combine'
        if sum(type_count)==0:
            self.type=None

        ## 传入一个term，解析其类型
    @staticmethod
    def get_basis_name(filename):
        # 使用正则表达式提取文件名中的特定部分
        pattern = r"pot-(.*?)\.json"
        match = re.search(pattern, filename)
        # print(match)
        if match:
            extracted_part = match.group(1)
            # print(extracted_part)
            return extracted_part
        else:
            # print("未找到匹配的部分")
            return "test"

    @classmethod   ## 类方法，用于构造函数
    def creat_by_data(cls,data_main_d,X,file_name='test.json'):
        """
        :param data_main_d: MD数据
        :param X: MD array
        :return:
        """
        md_dict={"md_data":X}
        term_dict=data_main_d['pair']
        #(term_dict,md_dict={},file_name="test.json",cutoff=None)
        return cls(term_dict,md_dict=md_dict,file_name=file_name)
    @classmethod
    def creat_by_file(cls,input_data,target_json="target.json"):
        c1, c2= read_json_file(input_data, target_json)
        file_name=os.path.basename(input_data)
        ## c1:md_array,c2:原来的字典
        return cls.creat_by_data(c2,np.array(c1),file_name=file_name)
    @classmethod
    def combine(cls,term_list_d):
        """
        :param term_list: 一个已经完全实例化的term类型的列表
        :return: 一个新的合并之后的类
        """
        term_list=deepcopy(term_list_d)
        example_dict=term_list[0].init_property
        example_interact=term_list[0].interact
        example_fix=term_list[0].fix
        new_term_dict=cls.plus_term([term_one.term for term_one in term_list],example_interact,example_fix)
        
        new_md_dict={}

        for key, value in example_dict.items():
            new_md_dict[key]=cls.plus_MD_data([term_one.init_property[key] for term_one in term_list])
            
        return cls(term_dict=new_term_dict,md_dict=new_md_dict,file_name="combine.json")
    @classmethod
    def combine_by_para(cls,term_list,para):
        new_term_list=[]
        if len(para)!=len(term_list):
            print("shape not coisistent")
        else:
            for i,size in enumerate(para):
                new_term_list+=[cls.reshape(term_list[i],size)]
            return cls.combine(new_term_list)
    
    @classmethod
    def reshape(cls,term,size):
        """
        :param term: 传入一个term，生成一个新的term，该term是reshape过的
        :return:
        """
        example_dict = term.init_property
        example_interact=term.interact
        example_fix=term.fix
        new_term_dict=cls.para_reshape(term.term,size,example_interact,example_fix)
        new_md_dict = {}

        for key, value in example_dict.items():
            new_md_dict[key] = cls.md_reshape(term.init_property[key],size)

        return cls(term_dict=new_term_dict,md_dict=new_md_dict,file_name="combine.json")
    def get_a(self):
        
        
        
        
        return self.term['pair']['values']['a']
    @staticmethod
    def para_reshape(term,size,interact_list,fix_list):
        """

        :param list: 传入一个系数列表，将系数reshape一个size
        :return: 一个新的reshape后的list
        """
        new_term_dict={}
        for interact_type in interact_list:
            term_interact=term[interact_type]
            new_term_dict.update({interact_type:{}})
            for element,term_info in term_interact.items():
                dict_new=MD_term.basis_reshape(term_interact[element],size)
                new_term_dict[interact_type].update({element:dict_new})

        for fix_type in fix_list:
            term_interact=term[fix_type]
            new_term_dict.update({fix_type:{}})
            for element,term_info in term_interact.items():
                new_term_dict[fix_type].update({element:term[fix_type][element]}) 
        return new_term_dict   
    
    @staticmethod
    def basis_reshape(term,size):
        new_term=deepcopy(term)
        if term['need_c']==1:
            print("this term can not be reshape because it need constant")
        elif term['need_c']==0:
            for dict1 in new_term['term']:
                dict1['a']*=size
        return new_term




    @staticmethod
    def plus_term(term_list,interact_list,fix_list):    ## 给出一堆由r0,n,a的字典的list  term，匹配一样的r0,n，相加a，如果a=0则剔除a
        new_term_dict={}
        for interact_type in interact_list:
            example_term=term_list[0][interact_type]
            new_term_dict.update({interact_type:{}})
            for element,term_info in example_term.items():
                dict_new=MD_term.basis_combine([term[interact_type][element] for term in term_list])
                new_term_dict[interact_type].update({element:dict_new})

        for fix_type in fix_list:
            example_term=term_list[0][fix_type]
            new_term_dict.update({fix_type:{}})
            for element,term_info in example_term.items():
                # print(example_term)
                dict_new=example_term[element]
                new_term_dict[fix_type].update({element:dict_new})
            
        return new_term_dict
    @staticmethod
    def basis_combine(term_list):
        need_c=0
        dict_map = {}
        for term in term_list:
            need_c+=term['need_c']
            for dict1 in term['term']:
                dict_rm_a=deepcopy(dict1)
                # print(dict_rm_a)
                dict_rm_a.pop('a')
                key=tuple(dict_rm_a.values())
                # print(key)
                if key in dict_map:
                    dict_map[key]['a']+=dict1['a']
                else:
                    dict_map[key]=dict1
        new_list=list(dict_map.values())
        if need_c:
            new_need_c=1
        else:
            new_need_c=0
        # print(dict_map)
        return {"term":new_list,"need_c":new_need_c}

    @staticmethod
    def plus_MD_data(MD_list):
        return np.sum(MD_list,axis=0)


    @staticmethod
    def md_reshape(data,size):
        return data*size

    # @classmethod
    def json_writer(self,target_dict='target.json',output_json='test.json',fitted_valued=True):
        path = os.getcwd()
        date = datetime.date.today().strftime("%m/%d/%y")
        current_time = datetime.datetime.now().strftime("%H:%M:%S")
        with open(target_dict,'r') as f:
            target_json = json.load(f)
        if fitted_valued:
            with open(output_json, "w") as f:
                data1={"pair":self.term,
                       "inforamtion": f'{date} {current_time}',
                       "fitted_value":get_dict(deepcopy(target_json['data']),deepcopy((self.md_data).tolist()))
                       }
                dict_data={}
                # for key,value in self.init_property.items():
                #     if key!="md_data" and key!="mf":
                #         dict_data[key]=value
                # data1['all_data']=dict_data
                json.dump(data1,f,indent=4, separators=(',', ':'))
                f.write('\n')
        else:
            with open(output_json, "w") as f:
                data1={"pair":self.pair_term,"embedding":self.embedding_term,"type":"combination",
                       "need_fc": self.fc,
                       "need_pc": self.pc,
                       "inforamtion": f'{date} {current_time}',
                       }
                dict_data={}
                data1['all_data']=dict_data
                json.dump(data1,f,indent=4, separators=(',', ':'))
                f.write('\n')

    @classmethod ## 传进来很多个term类，进行拟合后释放新的一个term
    def fit_progress(cls,term_list,target,type=None,weight=None,cf_straint=None,lambda_reg=None):
        X_selected=np.vstack([term.get_fitdata(type=type) for term in term_list])
        if cf_straint is None:
            para=cls.linear_regression(X_selected,target,weight=weight,lambda_reg=lambda_reg)
        elif cf_straint is not None:
            X_selected=np.hstack([X_selected,np.eye(X_selected.shape[0])])
            target=np.hstack([target,np.ones(X_selected.shape[0])])
            para = cls.linear_regression(X_selected,target,weight=weight,lambda_reg=lambda_reg)
        return cls.combine_by_para(term_list,para),para

    def get_fitdata(self,type=None):     ## 这个要改写成给定需要的字典来提取
        data_list =[]
        if type==None:
            return self.md_data
        elif type is not None :
            for key in type:
                data_list += [self.init_property[key]]

            return np.hstack(data_list)
    @classmethod ## 传进来很多个term类，进行拟合后释放新的一个term
    def fit_progress_with_bound(cls,term_list,target,bound=None,X0=None,type="MD"):
        X_selected=np.vstack([term.get_fitdata(type=type) for term in term_list])
        # print(X_selected.shape)
        result,para=cls.linear_regression_with_bound(X_selected,target,bound,X0=X0)
        # print(len(term_list),len(para))
        return cls.combine_by_para(term_list,para),para

    @staticmethod
    def linear_regression(X_selected,target,weight=None,lambda_reg=None):
        Y = target
        # 如有权重则先做加权
        if weight is not None:
            X_selected = X_selected * weight
            Y = Y * weight

        X_train = X_selected.T

        # 根据 lambda_reg 决定使用 Ridge 或 LinearRegression
        if lambda_reg is not None and lambda_reg > 0:
            reg = Ridge(alpha=lambda_reg, fit_intercept=False)
        else:
            reg = LinearRegression(fit_intercept=False)

        reg.fit(X_train, Y)
        para_a = reg.coef_
        return para_a
    def get_cost(self,target):
        md_array=deepcopy(self.md_data)
        target_array=deepcopy(target)
        target_array[target==0]+=10
        md_array[target==0]+=10
        self.cost=self.rmse_cost(md_array,target_array)

    def information(self):
        print(f"type:{self.type}\n")
        print(f"pair:{self.pair_term}\n")
        print(f"embedding:{self.embedding_term}\n")
        # print(f"cost with target:{self.cost}")
        print(f"md data:{self.md_data}")

    @classmethod
    def read_from_database(cls,data_dir, target_json='target.json',select_json=None):
        """
        :param data_dir: 数据库的目录
        :return: 返回数据矩阵X和原始的json文件
        """
        term_list = []
        if select_json is not None:
            with open(select_json) as f:
                data = json.load(f)
            string_list = data['pot']
            need_all=data['all']
            if need_all:
                data_list = [filename for filename in os.listdir(data_dir)]
                for file_name in data_list:
                    file_path = os.path.join(data_dir, file_name)
                    if os.path.exists(file_path):
                        term_list += [cls.creat_by_file(file_path, target_json=target_json)]
            else:
                for s in string_list:
                    # 构造文件名
                    file_name = s + ".json"
                    # print(file_name)
                    # 获取文件路径
                    file_path = os.path.join(data_dir, file_name)
                    # 如果文件存在，输出文件路径
                    if os.path.exists(file_path):
                        term_list += [cls.creat_by_file(file_path, target_json=target_json)]

        return term_list


    @classmethod
    def minus_term(cls,term1,term2):
        """
        返回两个单项的差，项的系数相减，重复的项合并，剔除为0的项
        :param term1:
        :param term2:
        :return:
        """
        para=np.array([1,-1])
        term_list=[term1,term2]
        return cls.combine_by_para(term_list,para)
    
    def add_lambda_data(self,lambda_file,dH, name,option=None):  ##输入多组列表？
        lambda_data,alpha_data=MD_term.analysis_lambda(lambda_file)
        if option is None:
            extra_dict = {f"{name}": lambda_data / dH,f"{name}_alpha":alpha_data}
        elif option =="set-zero":
            extra_dict = {f"{name}": 0}
        self.init_property.update(extra_dict)
        for key, value in extra_dict.items():
            setattr(self, key, value)
    def add_G_data(self,lambda_file,name,option=None):  ##输入多组列表？
        lambda_data,alpha_data=MD_term.analysis_lambda(lambda_file)
        if option is None:
            extra_dict = {f"{name}": lambda_data}
        elif option =="set-zero":
            extra_dict = {f"{name}": 0}
        self.init_property.update(extra_dict)
        for key, value in extra_dict.items():
            setattr(self, key, value)
    
    
    def add_cutoff(self,rol_max,rol_mean,rol_int=[10]):
        ## 如果是emb的类型，我们才会加上这个东西，然后还是对emb的进行操作
        interact_list=['emb']
        for interact_type in interact_list:
            example_term=self.term[interact_type] 
            i=0
            for element,term_info in example_term.items():
                cutoff_dict={f"maxCut_{element}":0,f"prime1_{element}":0,f"prime2_{element}":0}
                ##新增键int_{element}，用于存储积分值
                cutoff_dict.update({f"int_{element}":0})
                for basis_term in example_term[element]['term']:
                    a=basis_term['a']
                    r0=basis_term['r0']
                    n=basis_term['n']
                    cutoff_dict[f"int_{element}"]+=a/(n+1)*(rol_int[i]-r0)**(n+1)*((rol_int[i]-r0)>0)   ##仅考虑rho_int大于r0的情况
                    cutoff_dict[f"maxCut_{element}"]+=a*(rol_max[i]-r0)**n
                    cutoff_dict[f"prime1_{element}"]+=a*n*(rol_mean[i]-r0)**(n-1)*((rol_mean[i]-r0)>0)   ##仅考虑rho_mean大于r0的情况
                    cutoff_dict[f"prime2_{element}"]+=a*n*(n-1)*(rol_mean[i]-r0)**(n-2)*((rol_mean[i]-r0)>0)
                self.init_property.update(cutoff_dict)
                for key, value in cutoff_dict.items():
                    setattr(self, key, value)
            i+=1
    def add_F_basis(self,x_rho,element_sel=None):
        ## 如果是emb的类型，我们才会加上这个东西，然后还是对emb的进行操作
        interact_list=['emb']
        for interact_type in interact_list:
            example_term=self.term[interact_type] 
            i=0
            for element,term_info in example_term.items():
                cutoff_dict={f"F_basis_{element}":np.zeros_like(x_rho)}
                for basis_term in example_term[element]['term']:
                    if basis_term["type"]!="polynomial-2":
                        continue
                    a=basis_term['a']
                    r0=basis_term['r0']
                    n=basis_term['n']
                    add_array=a*(x_rho-r0)**n*(((x_rho-r0)>0))   ##仅考虑rho_int大于r0的情况
                    cutoff_dict[f"F_basis_{element}"]=cutoff_dict[f"F_basis_{element}"]+add_array
                self.init_property.update(cutoff_dict)
                for key, value in cutoff_dict.items():
                    setattr(self, key, value)
            i+=1
    def add_phi_basis(self,x_r,element_sel=None):
        ## 如果是emb的类型，我们才会加上这个东西，然后还是对emb的进行操作
        interact_list=['pair']
        for interact_type in interact_list:
            example_term=self.term[interact_type] 
            i=0
            for element,term_info in example_term.items():
                # print(element,element_sel)
                if element_sel is  not None and element not in element_sel:
                    # print("jump")
                    continue
                cutoff_dict={f"phi_basis_{element}":np.zeros_like(x_r)}
                for basis_term in example_term[element]['term']:
                    a=basis_term['a']
                    r0=basis_term['r0']
                    n=basis_term['n']
                    add_array=a*(r0-x_r)**n*(((r0-x_r)>0))   ##仅考虑rho_int大于r0的情况
                    # add_force_array=a*n*(r0-x_r)**(n-1)*(((r0-x_r)>0))   ##仅考虑rho_int大于r0的情况
                    cutoff_dict[f"phi_basis_{element}"]=cutoff_dict[f"phi_basis_{element}"]+add_array
                self.init_property.update(cutoff_dict)
                for key, value in cutoff_dict.items():
                    setattr(self, key, value)
            i+=1    
    def add_force_basis(self,x_r,element_sel=None):
        ## 如果是emb的类型，我们才会加上这个东西，然后还是对emb的进行操作
        interact_list=['pair']
        for interact_type in interact_list:
            example_term=self.term[interact_type] 
            i=0
            for element,term_info in example_term.items():
                # print(element,element_sel)
                if element_sel is  not None and element not in element_sel:
                    # print("jump")
                    continue
                cutoff_dict={f"force_basis_{element}":np.zeros_like(x_r)}
                for basis_term in example_term[element]['term']:
                    a=basis_term['a']
                    r0=basis_term['r0']
                    n=basis_term['n']
                    # add_array=a*(r0-x_r)**n*(((r0-x_r)>0))   ##仅考虑rho_int大于r0的情况
                    add_force_array=-a*n*(r0-x_r)**(n-1)*(((r0-x_r)>0))   ##仅考虑rho_int大于r0的情况
                    cutoff_dict[f"force_basis_{element}"]=cutoff_dict[f"force_basis_{element}"]+add_force_array
                self.init_property.update(cutoff_dict)
                for key, value in cutoff_dict.items():
                    setattr(self, key, value)
            i+=1 
    def add_cutoff_pair(self,innercut_dict,select=None):
        interact_list=['pair']
        for interact_type in interact_list:
            example_term=self.term[interact_type] 
            for element,term_info in example_term.items():
                # print("deeling with element",element)
                if select is not None and element not in select:
                    # print(f"skip {element}")
                    continue
                ##如果innercut_dict的类型不是字典而是一个数值，那么我们就直接用这个数值
                if not isinstance(innercut_dict,dict):
                    innercut=innercut_dict
                elif element in innercut_dict:
                    innercut=innercut_dict[element]
                
                cutoff_dict={f"innerprime1_{element}":0}
                for basis_term in example_term[element]['term']:
                    type=basis_term['type']
                    if type=="polynomial-1":
                        a=basis_term['a']
                        r0=basis_term['r0']
                        n=basis_term['n']
                        # print(a,r0,n,cutoff_dict[f"innerprime1_{element}"],innercut)
                        cutoff_dict[f"innerprime1_{element}"]+=-a*n*(r0-innercut)**(n-1)
                self.init_property.update(cutoff_dict)
                for key, value in cutoff_dict.items():
                    setattr(self, key, value)

        
        
    def add_cutoff_mf(self,mfp):
            ## 如果是emb的类型，我们才会加上这个东西，然后还是对emb的进行操作
            interact_list=['emb']
            for interact_type in interact_list:
                example_term=self.term[interact_type] 
                i=0
                for element,term_info in example_term.items():
                    data=self.init_property[f"{mfp}"][0]
                    # print(data)
                    cutoff_dict={f"mf_{element}":0}
                    for basis_term in example_term[element]['term']:
                        a=basis_term['a']
                        r0=basis_term['r0']
                        n=basis_term['n']
                        cutoff_dict[f"mf_{element}"]+=data
                    self.init_property.update(cutoff_dict)
                    for key, value in cutoff_dict.items():
                        setattr(self, key, value)
                i+=1
            self.init_property[f"{mfp}"]=np.zeros(len(self.init_property[f"{mfp}"]))
    def overfit(self): 
        interact_list=['pair']
        for interact_type in interact_list:
            example_term=self.term[interact_type] 
            i=0
            for element,term_info in example_term.items():
                cutoff_dict={f"coff_{element}":0}
                for basis_term in example_term[element]['term']:
                    a=basis_term['a']
                    r0=basis_term['r0']
                    n=basis_term['n']
                    cutoff_dict[f"coff_{element}"]+=a
                    
                self.init_property.update(cutoff_dict)
                for key, value in cutoff_dict.items():
                    setattr(self, key, value)
        
        
    
        
    @staticmethod
    def get_H_T(solid_json,liquid_json):       ## using this to get the target Tm
        with open(solid_json,'r') as f:
            data=json.load(f)
            temp_solid=data['temp']
            h_solid=data['enthalpy']
            e_solid=data['pe']
        with open(liquid_json,'r') as f:
            data=json.load(f)
            temp_liquid=data['temp']
            h_liquid = data['enthalpy']
            e_liquid=data['pe'] 
            
            
            TmX0=np.mean([temp_liquid,temp_solid])
            dHX0=h_liquid-h_solid
            deX0=e_liquid-e_solid
            a=-dHX0/TmX0

            return TmX0,dHX0,a,deX0
    @staticmethod
    def get_U_T(solid_json,liquid_json):       ## using this to get the target Tm
        with open(solid_json,'r') as f:
            data=json.load(f)
            temp_solid=data['temp']
            e_solid=data['pe']
        with open(liquid_json,'r') as f:
            data=json.load(f)
            temp_liquid=data['temp']
            e_liquid=data['pe'] 
            TmX0=np.mean([temp_liquid,temp_solid])
            deX0=e_liquid-e_solid

            return TmX0,deX0
    @staticmethod
    def get_U_Gmix(x_json,y_json):       ## using this to get the target Tm
        with open(x_json,'r') as f:
            data=json.load(f)
            e_x=data['pe']
            x=data["x"]
        with open(y_json,'r') as f:
            data=json.load(f)
            e_y=data['pe']
            y=data["x"]
            return e_x-x/y*e_y
    @staticmethod
    def get_U_Gxpx(x_json,y_json,x0_json):       ## using this to get the target Tm
        with open(x_json,'r') as f:
            data=json.load(f)
            e_x=data['pe']
            x=data["x"]
        with open(y_json,'r') as f:
            data=json.load(f)
            e_y=data['pe']
            y=data["x"]
        with open(x0_json,'r') as f:
            data=json.load(f)
            e_x0=data['pe']
            x0=data["x"]
        # print(e_x,e_y,e_x0)
        return e_x-x/y*e_y-(1-x/y)*e_x0
    @staticmethod
    def read_H_thermo(thermo): 
        with open(thermo,'r') as f:
            data=json.load(f)
        enthalpy=data['enthalpy']
        temp=data['temp']
        vol=data['vol']
        press=data['press'] ## press in GPa,temp in K ,vol in A^3
        kB=1.380649e-23
        eV=1.60217662e-19
        p_ke=data["presske"]   ## ke contribute to press in GPa
        ke=data["ke"]
        h_ke=5/3*ke
        return p_ke,h_ke,enthalpy
    def read_Mf_thermo(thermo): 
        with open(thermo,'r') as f:
            data=json.load(f)
        # temp=data['temp']
        # vol=data['vol']
        # press=data['press'] ## press in GPa,temp in K ,vol in A^3
        presske=data['presske']
        return presske
    @staticmethod
    def analysis_lambda(lambda_file):
        with open(lambda_file,'r') as f:
            data=json.load(f)
            lambda_data=data['dlambda']
            alpha_data=data['dalpha']
        return lambda_data,alpha_data
    @staticmethod
    def analysis_thermo(thermo_press_file):
        with open(thermo_press_file,'r') as f:
            data = json.load(f)
            press = data["press"]
            pe=data["pe"]
        return press,pe
    
    def add_mf_data(self,meanforce_file,info,option=None):  ##我们要输入的就是
        km_list=info['km']
        name=info['name']
        mf,mfSum=MD_term.analysis_mf(meanforce_file,km_list)
        if option==None:
            extra_dict={f"{name}-{i}":mf[i] for i,km in enumerate(km_list)}
            # mfsum_dict={f"{name}-{i}Sum":mfSum for i,km in enumerate(km_list)}
            # extra_dict.update(mfsum_dict)
        self.init_property.update(extra_dict)
        for key, value in extra_dict.items():
            setattr(self, key, value) 
    @staticmethod
    def Gx_para(x_data,y_data,T_data):
        # print(x_data,y_data,T_data)
        new_fit=lambda x,a1,a2: fit_Gx_single(x,a1,a2,T_data)
        params, cov = curve_fit(new_fit, x_data, y_data)
        ##打印拟合误差
        error=np.sqrt(np.mean((new_fit(x_data,*params)-y_data)**2))/np.mean(y_data)
        
        a1_fit, a2_fit = params
        # print("拟合结果: a1 =", a1_fit, ", a2 =", a2_fit)

        c = T_data * kB * (0.2*np.log(0.2) + 0.8*np.log(0.8))
        a0_fit = (-0.16*a1_fit + 0.096*a2_fit - c) / 0.2
        return a0_fit, a1_fit, a2_fit,error
    @staticmethod
    def G_xGp(para,x):
        a0,a1,a2,T=para
        Gx=a0*x+a1*x*(1-x)+a2*x*(1-x)*(2*x-1)+T*kB*(x*np.log(x)+(1-x)*np.log(1-x))
        Gxpx=a0+a1*(1-2*x)+a2*((1-x)*(2*x-1)+x*(1-2*x)+2*x*(1-x))+T*kB*(np.log(x)-np.log(1-x))
        return Gx-Gxpx*x
    @staticmethod
    def Gx(para,x):
        a0,a1,a2,T=para
        Gx=a0*x+a1*x*(1-x)+a2*x*(1-x)*(2*x-1)+T*kB*(x*np.log(x)+(1-x)*np.log(1-x))
        return Gx
    @staticmethod
    def Gp(para,x):
        ## return the prime of G according to parameter
        a0,a1,a2,T=para
        Gxpx=a0+a1*(1-2*x)+a2*((1-x)*(2*x-1)+x*(1-2*x)+2*x*(1-x))+T*kB*(np.log(x)-np.log(1-x))
        return Gxpx
    def add_mu(self,file,name,x_list,T,xt):
        ##提取self的不同x下的G值
        G=MD_term.analysis_G_xGp(file,name,x_list)
        ##拟合a0,a1,a2的值 
        a0,a1,a2,error=MD_term.Gx_para(x_list,G,T)
        if error>0.01:
            print(f"error={error} for term {self.name} is so large,a0={a0:.2f},a1={a1:.2f},a2={a2:.2f}")
        # print(a0,a1,a2,self.name)
        para=(a0,a1,a2,T)
        ##计算G-Gpx的值
        mu=MD_term.G_xGp(para,xt)
        extra_dict={f"mu-{name}":mu}
        para_dict={f"Gmix-{name}":para}             ##Gmix fitting parameter
        self.init_property.update(extra_dict)
        # self.init_property.update(para_dict)
        for key, value in extra_dict.items():
            setattr(self, key, value) 
        for key, value in para_dict.items():
            setattr(self, key, value)               ##Gmix fitting parameter add to the term
    def add_dmu(self,nameA,nameB,xA,xB):
        ## add property about the difference of chemical potential between A and B
        paraA=getattr(self,f"Gmix-{nameA}")
        paraB=getattr(self,f"Gmix-{nameB}")
        muA=MD_term.G_xGp(paraA,xA)
        muB=MD_term.G_xGp(paraB,xB)
        dmu=muA-muB
        extra_dict={f"dmu-{nameA}-{nameB}":dmu}
        self.init_property.update(extra_dict)
        for key, value in extra_dict.items():
            setattr(self, key, value)
    def add_Gp(self,nameA,nameB,xA,xB):
        ## add property about free energy derivative between A and B
        paraA=getattr(self,f"Gmix-{nameA}")
        paraB=getattr(self,f"Gmix-{nameB}")
        GpxA=MD_term.Gp(paraA,xA)
        GpxB=MD_term.Gp(paraB,xB)
        dGp=GpxA-GpxB
        extra_dict={f"dGp-{nameA}-{nameB}":dGp}
        self.init_property.update(extra_dict)
        for key, value in extra_dict.items():
            setattr(self, key, value)
        
          
    @staticmethod
    def analysis_mf(meanforce_file,km_all):                     
        ## km as a list
        kmax=-12
        data_all=np.loadtxt(meanforce_file,skiprows=1)        ##计算meanforce 后需要将其积分起来
        basis_all=[];sumMf_all=[]
        for i,km in enumerate(km_all):
            data=data_all[:,i+1]
            data[np.isnan(data)]=0                            ## 避开nan值，赋值为0
            basis=data[km:kmax]                                 ## 从km开始关心g(r)的值，比如从0.1开始
            basis_all+=[basis]
        #     r=data[km:,0];mfVal=data[km:,1]
        # dr=np.zeros_like(r);dr[1:]=np.diff(r)
        # sumMf=np.cumsum(dr*mfVal)                         ## 作为求和 我觉得可以直接去改cpp代码，直接使用
        
        return basis_all,sumMf_all
    @staticmethod
    def analysis_G_xGp(database_file,name,x_list):
        with open(database_file,'r') as f:
            data=json.load(f)['data']
        G_list=[]
        for x in x_list:
            G=data[f"{name}-{x}"]['dUmix']
            G_list+=[G]
        
        return np.array(G_list)
    @staticmethod
    def get_effective_gr(temp,tar_file,grcut):
        kB=8.61733326E-5
        gr_tar = np.loadtxt(tar_file, skiprows=1)
        num_gr = gr_tar.shape[0] ## 读取gr的行数据
        gmin = grcut
        rall=gr_tar[:,0]
        rfit_all=[];mf_tar_all=[];kmin_all=[];gr_tar_all=[]
        kmax=-12
        for gtype in range(1,gr_tar.shape[1]):
            count = 0;kmin = 9999;rm = 0
            gr_one=gr_tar[:,gtype]
            for k in range(num_gr):
                if (gr_one[k] > gmin):  ## 如果g大于最小值则计入要拟合的值
                    count += 1  ## 则拟合
                    if (k < kmin):  ##设定要拟合的阈值
                        kmin = k
                        rm = rall[kmin]
            print(f"Total gr data: {num_gr}; fitting data: {count}; from {rm} (k={kmin})")
            rfit = gr_tar[kmin:kmax, 0]
            a_gr=np.polyfit(rfit,np.log(gr_one[kmin:kmax]),30)
            plt.figure()
            # plt.plot(rfit,np.log(gr_one[kmin:kmax]))
            # plt.plot(rfit,np.polyval(a_gr,rfit))
            a_mf=np.polyder(a_gr)
            # plt.plot(rfit,np.polyval(a_mf,rfit))
            mf_tar=np.polyval(a_mf,rfit)*kB*temp
            rfit_all+=[rfit]
            kmin_all+=[kmin]
            mf_tar_all+=[mf_tar]
            gr_tar_all+=[gr_one[kmin:kmax]]
            npair=gr_tar.shape[1]-1
            
        return rfit_all,kmin_all,mf_tar_all,gr_tar_all,npair
    def plt_term_fun(self,out_name=False,target_json='target.json',filename=None,type='all',pair_cut=1.5):
        if filename==None:
            if out_name:
                term_json_file=f"{out_name}.json"
                self.json_writer(output_json=term_json_file,target_dict=target_json)
                pot_writer(pair=term_json_file,out_name=f"{out_name}",if_fit=1,pair_cut=pair_cut)
                plot_fun(pot=f"{out_name}",pltType=type)
            elif out_name ==False:
                term_json_file=f"combine.json"
                self.json_writer(output_json=term_json_file,target_dict=target_json)
                pot_writer(pair=term_json_file,out_name=f"pot_combine.fs",if_fit=1,pair_cut=pair_cut)
                plot_fun(pot=f"pot_combine.fs",pltType=type)
        elif filename is not None:
            term_json_file = filename
            out_name="pot.fs"
            pot_writer(pair=term_json_file, out_name=f"{out_name}",if_fit=1,pair_cut=pair_cut)

    def addConstant(self):
        ##为所需要的东西增加cutoff
        interact_list=['pair']
        for interact_type in interact_list:
            example_term=self.term[interact_type]
            for element,term_info in example_term.items():
                example_term[element]['need_c']=1

    def addGr(self,gr_cut,rfit,name,temp):
        mf_fit=getattr(self,f"{name}")        ##获得已经拟合的gr数据
        print(len(mf_fit))
        lngr_int=[]
        kB=8.61733326E-5
        for i,r in enumerate(rfit):
            lngr_int+=[np.trapz(mf_fit[:i+1],rfit[:i+1])]
        lngr_int=np.array(lngr_int)
        gr_val=gr_cut*np.exp(lngr_int/(kB*temp))
        setattr(self, f"gr-{name}", gr_val)
    @classmethod
    def pot_writer_from_json(cls,json_file,pot_name='pot.fs'):
        pot_writer(pair=json_file,out_name=f"{pot_name}")

    @classmethod
    def slect_by_name(cls,term_list,name):
        i=[]
        for term in term_list:
            if term.name==name:
                i+=[]
                return term
        if len(i)==0:
            print(f"no term name {name}")
            # return None
    @classmethod
    def rm_by_name(cls,term_list,name_list):
        select_term=[cls.slect_by_name(term_list,name) for name in name_list]
        list1=term_list
        list2=select_term
        return list(set(list1) - set(list2))


