#!/usr/bin/env python
# coding: utf-8

# In[4]:


#!/usr/bin/env python3
import matplotlib.pyplot as plt
import md_term_alloy
import importlib
importlib.reload(md_term_alloy)
from md_term_alloy import MD_term
from md_term_alloy import pot_writer
from md_term_alloy import get_array
from md_term_alloy import plot_fun
import numpy as np
import json
import os
from copy import deepcopy
import argparse
import matplotlib
font = {
        'size'   : 18}
plt.rc('font', **font)
# matplotlib.use('Agg')
parser = argparse.ArgumentParser()
parser.add_argument("-o","--output", help="output potential",action='store')

current_directory = os.getcwd()
print(current_directory)
parent_directory = os.path.dirname(current_directory)
data_base=os.path.join(parent_directory,"init","database")
pair_json=os.path.join(current_directory,"pair_json")

## 按照文件名匹配
dir_target=os.path.join(parent_directory,"target")

dir_atom=os.path.join(parent_directory,"atom")
dir_init=os.path.join(parent_directory,"init")
target_solid=os.path.join(dir_target,"target_solid.json")
with open(target_solid,'r') as f:
    target_json=json.load(f)


mf_list=["Meanf-liq_T=6000_p=323-Fe4O1"]
mfp_list=["mf_liq_T=6000_p=323-Fe4O1"]

xl_list=[0.04,0.08,0.12,0.16]
Gmix_list=[f"Gmix-liq_T=5500_p=323-{x}" for x in xl_list]
Gmixp_list=[f"Gmix_liq_T_5500_p_323_{x}" for x in xl_list]
dU_X0_G=[]
for i,G in enumerate(Gmix_list):
    G_x=os.path.join(dir_atom,G,"thermo_x.json")
    G_y=os.path.join(dir_atom,G,"thermo_y.json")
    dU_X0_G.append(MD_term.get_U_Gmix(G_x,G_y))

##定义perturbation
b=0
args = parser.parse_args(['-o',f'pot_combine.fs'])
def G_perturbation(b,x):
    return b*x*(1-x/0.2)


Gmix_now_l=np.array([-0.033716, -0.039292, -0.03376 , -0.019931])
##load Gmix from Gmix.dat。寻找和xl_list对应的Gmix,Gmix.dat中包含其他x的Gmix，所以需要挑选

data_Gmix=np.loadtxt("Gmix.dat")
Gmix_now_l=data_Gmix[1:5,1]
Gmix_now=np.concatenate([Gmix_now_l])
Gmix_target=np.zeros_like(Gmix_now)
for i,G in enumerate(Gmix_list):
    Gmix_target[i]=target_json['data'][G]['dlambda']+G_perturbation(b,xl_list[i])
    target_json['data'][G]['dlambda']+=(dU_X0_G[i]-Gmix_now[i]+G_perturbation(b,xl_list[i]))



temp_list=[6000]
info_list=[]
target_mf_list=[]
rfit_list=[]
gr_target=[]
km_list=[]
npair_list=[]
weight_mf_list=[]

x_rho=np.arange(25,100,1)
p_ke_all=[]
for i,mf in enumerate(mf_list):
    thermo_file=os.path.join(dir_atom,mf,"thermo.json")
    target_liquid=os.path.join(dir_atom,mf,f"gr-target.dat")
    rfit,km,target_mf_one,gr_temp,npair=MD_term.get_effective_gr(temp_list[i],target_liquid,grcut=0.1)
    p_ke=MD_term.read_Mf_thermo(thermo_file)
    p_ke_all+=[p_ke]
    target_mf_list+=[target_mf_one]
    rfit_list+=[rfit]
    index_long=np.where(rfit[1]>4)
    gr_target+=[gr_temp]
    info_list+=[{"km":km,"name":mfp_list[i]}]
    km_list+=[km]
    npair_list+=[npair]
    weight_mf_list+=[[np.ones_like(target_mf) for target_mf in target_mf_one]]
    target_json['data'][mf]['press']-=p_ke
weight_mf_list[0][1][index_long]*=2
weight_mf_list[0][0]*=1
weight_mf_list[0][2]*=1
target_MD=np.array(get_array(target_json['data'],[]))

good_pot_json=os.path.join(dir_target,"good_pot.json")
tt=MD_term.read_from_database(data_base,select_json=good_pot_json,target_json=target_solid)
print("reading terms")
for j,term in enumerate(tt):
    out_put_name=os.path.join(pair_json,f"pair_{term.name}.json")
    # term.add_cutoff([200],[100])
    for i,mfp in enumerate(mfp_list):
        dir_meanforce=os.path.join(parent_directory,"init",mf_list[i])
        meanforce_file=os.path.join(dir_meanforce,f"meanforce_{term.name}.dat")
        term.add_mf_data(meanforce_file,info_list[i])
    term.add_cutoff([100],[39])
    term.add_cutoff_pair(1.35)
    term.add_F_basis(x_rho)
    tt[j]=term
with open("target_new.json",'w') as f:
    json.dump(target_json,f,indent=4, separators=(',', ':'))



# In[ ]:


print(Gmix_target)
print(target_MD)


# In[ ]:


weight_MD=1/(np.abs(target_MD)+1)*1
weight_MD[0]*=10


weight_MD[1:5]=np.ones_like(weight_MD[1:5])*100
# weight_MD[11]=np.ones_like(weight_MD[11])*80
mf_type=[f"{mfp}-{pair}"   for i,mfp in enumerate(mfp_list) for pair in range(npair_list[i])]
print(target_MD*weight_MD)

weight_mf_all=np.hstack([weight_mf_list[i][pair]   for i,mfp in enumerate(mfp_list) for pair in range(npair_list[i])])
target_mf=np.hstack([target_mf_list[i][j]  for i in range(len(target_mf_list)) for j in range(len(target_mf_list[i]))])
weight_mf=np.ones_like(target_mf)/len(target_mf)*40*weight_mf_all
weight_mf[target_mf<2]*=2
weight_rho=np.ones_like(x_rho)*0.5
target_rho=np.zeros_like(x_rho)

# target_all=np.hstack([target_MD,target_mf,8000,3,-60])
# type_list=['md_data']+mf_type+["maxCut_1"]+[f"prime1_1"]+[f"innerprime1_21"]
# weight_all=np.hstack([weight_MD,weight_mf,1/8000,1,0.1])
target_all=np.hstack([target_MD,target_mf,-60,-60,100,10,target_rho])

type_list=['md_data']+mf_type+[f"innerprime1_21"]+[f"innerprime1_22"]+["maxCut_2"]+[f"prime2_2"]+[f"F_basis_2"]
weight_all=np.hstack([weight_MD,weight_mf,0.001,0.01,0.001,0,weight_rho])



# In[ ]:


weight_MD


# In[ ]:


##读取data_base文件夹下的所有文件名，使用正则表达式匹配文件名，得到所有的文件名
print(data_base)
# data_base=os.path.join(parent_directory,"init","database")
file_list=os.listdir(data_base)
##匹配文件名中含有pair的
file_list_pair=[file for file in file_list if "pair" in file]
##去除文件名的pot-前缀和.json后缀
file_list_pair=[file[4:-5] for file in file_list_pair]
##匹配文件名中含有emb的
file_list_emb=[file for file in file_list if "emb" in file]
##去除文件名的pot-前缀和.json后缀
file_list_emb=[file[4:-5] for file in file_list_emb]
print(file_list_pair)
print(file_list_emb)


# 

# In[ ]:


ite_all=file_list_pair+file_list_emb
term_single=MD_term.slect_by_name(tt,"0-type-single-elements-single") 
term_list_ite=[MD_term.slect_by_name(tt,term) for term in ite_all]

target_0=target_all-term_single.get_fitdata(type=type_list)
fit_pot_1,para_1=MD_term.fit_progress(term_list=term_list_ite,target=target_0,type=type_list,weight=weight_all)
fit_pot=MD_term.combine([term_single,fit_pot_1])

fit_pot.addConstant()
fit_pot.json_writer(target_dict=target_solid,output_json="pair_fitted_0.json")
plt.figure(2222,figsize=(10,10))
fit_pot.plt_term_fun(target_json=target_solid,out_name=args.output,pair_cut=[1.35,1.35,1.35])


plt.ylim(-2,8)
plt.savefig('pair.png')


plt.figure(25,figsize=(10,10))
count=1
test_term=test_term=tt[24]

# for i,mfpName in enumerate(mf_list):
for j in range(npair_list[0]):
    plt.subplot(3,3,count)
    print(count)
    plt.plot(rfit_list[0][j],getattr(fit_pot,f"{mfp_list[0]}-{j}"))
    plt.plot(rfit_list[0][j],target_mf_list[0][j])
    count+=1
    plt.ylim(-2,16)
plt.savefig("mf-fit")

plt.figure(26,figsize=(20,10))
count=1
for i,mfpName in enumerate(mf_list):
    for j in range(npair_list[i]):
        plt.subplot(1,6,count)
        # plt.plot(rfit_list[i][j],getattr(fit_pot,f"{mfp_list[i]}-{j}"))
        # plt.plot(rfit_list[i][j],target_mf_list[i][j])
        count+=1
        plt.plot(rfit_list[i][j],getattr(test_term,f"{mfp_list[i]}-{j}"),label=f"{test_term.name}")
        # plt.legend()
# plt.legend(['fit','target','MIM'])
plt.savefig("mf-test")
fit_pot.init_property["prime2_2"]


# In[ ]:


fitted_data=json.load(open(f"pot_combine.fs.json",'r'))
fit_data_G=[fitted_data["fitted_value"][G]["dlambda"] for G in Gmix_list]

fit_data_G_l=np.array(fit_data_G[0:4])
Gmix_list_l=np.array(Gmix_list[0:4])

G_fit_list=[]
for i,G in enumerate(Gmix_list_l):
    G_fit=fit_data_G[i]-dU_X0_G[i]+Gmix_now[i]
    print(f"{G} {G_fit}");G_fit_list+=[G_fit]
plt.plot(xl_list,Gmix_now[0:len(xl_list)],label="Gmix_now")
plt.plot(xl_list,Gmix_target[0:len(xl_list)],label="Gmix_target")
plt.plot(xl_list,G_fit_list,label="Gmix_fit")
plt.legend()
plt.savefig("Gmix.png")



