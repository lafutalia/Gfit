import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
##去除divide by zero的警告
np.seterr(divide='ignore')
## def own basis fun for certain destination
def H0(x):
    return (np.array(x > 0)).astype(int)

def H(x):
    return (np.array(x >= 0)).astype(int)

def ufm_pair(x,term):
    p=term["p"]
    beta=term["beta"]
    sigma=term["sigma"]
    a=term["a"]
    fun=-p/beta*np.log(1-np.exp(-(x/sigma)**2))
    ##将x=0处的值设为1e10
    fun[x==0]=1e10
    return 0,0,fun
def ufm_pair_prime_sigma(x,term):
    p=term["p"]
    beta=term["beta"]
    sigma=term["sigma"]
    a=term["a"]
    fun=2*x**2*p/beta/sigma**3*np.exp(-(x/sigma)**2)/(1-np.exp(-(x/sigma)**2))
    ##将x=0处的值设为1e10
    fun[x==0]=1e10
    return 0,0,fun    
    
    
def zero_rho(x,term):
    return 0,0, np.zeros_like(x)
def zero_emb(x,term):
    return 0,0, np.zeros_like(x)
def TBM_rho(x,term):
    N=term["N"]
    epsilon=term["epsilon"]
    out_cut=term["out_cut"]
    need_shift=term["need_shift"]
    fun=(N*x**3*np.exp(-x*epsilon))**2
    fun_outcut=(N*out_cut**3*np.exp(-out_cut*epsilon))**2
    if need_shift:
        fun=fun-fun_outcut*x/out_cut
    fun[x==0]=1e10
    # print(fun.shape)
    return 0,0,fun


def Alfe_pair(x,term):
    A=term["A"]
    a=term["a"]
    n=term["n"]
    epsilon=term["epsilon"]
    innercut=term["innercut"]
    need_shift=term["need_shift"]
    out_cut=term["out_cut"]
    fun=epsilon*(A/x)**n*H(x-innercut)
    fun_outcut=epsilon*(A/out_cut)**n
    # print(fun_outcut,need_shift)
    if need_shift:
        # print(1)
        fun=(fun-fun_outcut)*H(out_cut-x)
    else:
        fun=(fun)*H(out_cut-x)
    left_value=epsilon*(A/innercut)**n
    left_value_prime=-n*epsilon*A**n/innercut**(n+1)
    if a!=1:
        print("WARNING: try to change the value of the coef in Alfe potential")
    return left_value_prime,left_value,fun

def Alfe_rho(x,term):
    A=term["A"]
    m=term["m"]
    out_cut=term["out_cut"]
    need_shift=term["need_shift"]
    fun=(A/x)**m
    fun_outcut=(A/out_cut)**m
    if need_shift:
        fun=(fun-fun_outcut)*H(out_cut-x)
    else:
        fun=(fun)*H(out_cut-x)
    fun[x==0]=1e10

    return 0,0,fun

def Alfe_emb(x,term):
    a=term["a"]
    epsilon=term["epsilon"]
    C=term["C"]
    fun=-epsilon*C*x**(1/2)
    if a!=1:
        print("WARNING: try to change the value of the coef in Alfe potential")
    return 0,0,fun

def polynomial_1(x,term):
    a=term['a']
    n=term['n']
    x0=term['r0']
    innercut=term['innercut']
    left_value=a*(x0-innercut)**n
    left_value_prime=-n*a*(x0-innercut)**(n-1)
    return left_value_prime,left_value,a*(x0-x)**n*H(x0-x)*H(x-innercut)

def polynomial_2(x,term):
    a=term['a']
    n=term['n']
    x0=term['r0']
    return 0,a*(x-x0)**n*H(x-x0)

def exp_polynomial_1(x,term):
    k=term['k']
    n=term['n']
    a=term['a']
    x0=term['r0']
    return 0,a*np.exp(k*x)*(np.abs(x-x0))**n


def exp_polynomial_1_mixing(x,term):
    k=term['k']
    n=term['n']
    a=term['a']
    x0=term['r0']
    # print(k,n,a,x0)
    value=np.zeros_like(x)
    for i,kt in enumerate(k):
        # print(k)
        tmp=a[i]*np.exp(k[i]*x)*(np.abs(x-x0[i]))**n[i]
        value=value+tmp
    return 0,value/2

def gen_pot_fun(x,term_dict,need_left=False,need_left_prime=False):
    term_extra_dict={"Alfe_pair":Alfe_pair,"Alfe_rho":Alfe_rho,"Alfe_emb":Alfe_emb,
                     "zero_rho":zero_rho,"zero_emb":zero_emb,
                     "ufm_pair":ufm_pair,"ufm_pair_prime_sigma":ufm_pair_prime_sigma,
                     "TBM_rho":TBM_rho}
    term_list=term_dict['term']
    need_c=term_dict['need_c']
    term_value=[]
    left_value=[];left_value_prime=[]
    if len(term_list)==0:
        left_all=0
        left_all_prime=0
        term_all=np.zeros_like(x)
    elif len(term_list)!=0:
        for i,term in enumerate(term_list):
            if term['type'] in term_extra_dict.keys():
                left_one_value_prime,left_one_value,one_value=term_extra_dict[term['type']](x,term)
                left_value+=[left_one_value]
                term_value+=[one_value]
                left_value_prime+=[left_one_value_prime]
            if term['type']=='polynomial-1':
                left_one_value_prime,left_one_value,one_value=polynomial_1(x,term)
                left_value+=[left_one_value]
                term_value+=[one_value]
                left_value_prime+=[left_one_value_prime]
            if term['type']=='polynomial-2':
                left_one_value,one_value=polynomial_2(x,term)
                left_value+=[left_one_value]
                term_value+=[one_value]
                left_value_prime+=[0]
            if term['type']=='exp-polynomial-1':
                left_one_value,one_value=exp_polynomial_1(x,term)
                left_value+=[left_one_value]
                term_value+=[one_value]
                left_value_prime+=[0]
            if term['type']=='exp-polynomial-1-mixing':
                left_one_value,one_value=exp_polynomial_1_mixing(x,term)
                left_value+=[left_one_value]
                term_value+=[one_value]
                left_value_prime+=[0]
           
        left_all=np.sum(np.array(left_value))
        left_all_prime=np.sum(np.array(left_value_prime))
        term_all=np.sum(np.array(term_value),axis=0)
    if need_left==True and need_left_prime==False:
        return left_all,term_all
    if need_left==False and need_left_prime==True:
        return left_all_prime,term_all
    if need_left==True and need_left_prime==True:
        return left_all_prime,left_all,term_all
    if need_left==False and need_left_prime==False:
        return term_all

def constant_pair_1(r,value,value_prime,innercut=1.5):
    np.seterr(divide='ignore')
    # phi_r = (9.7342365892908E+03 / r) * (0.1818 * np.exp(-2.8616724320005E+01 * r) + 0.5099 * np.exp(
    #     -8.4267310396064E+00 * r) + 0.2802 * np.exp(-3.6030244464156E+00 * r) + 0.02817 * np.exp(
    #     -1.8028536321603E+00 * r)) * H0(1.0000 - r) + np.exp(
    #     1.0969029696852E-01 + 1.7897391251247E+01 * r - 1.7775009177332E+01 * r ** 2 + 4.5602908011289E+00 * r ** 3) * H(
    #     r - 1.0000) * H0(innercut - r)+value*H0(innercut - r)
    # phi_3 = np.exp(
    #     1.0969029696852E-01 + 1.7897391251247E+01 * innercut - 1.7775009177332E+01 * innercut ** 2 + 4.5602908011289E+00 * innercut ** 3) * H(
    #     innercut - 1.0000) * H(innercut - r)
    
    # t=phi_r - phi_3
    # return t
    
    def cut_off_fun(r,var,para):
        x,f_value,f_value_prime,k1,k2=para
        A1,A2=var
        return 1/r*(A2*np.exp(-k2*r)+A1*np.exp(-k1*r))*H(x-r)+10*(x-r)*1/r**2*H(x-r)
    def cut_off_f1(var,para):
        x,f_value,f_value_prime,k1,k2=para
        A1,A2=var
        fun_value=1/x*(A2*np.exp(-k2*x)+A1*np.exp(-k1*x))
        fun_value_prime=-1/x**2*(A2*np.exp(-k2*x)+A1*np.exp(-k1*x))-(1/x)*(k2*A2*np.exp(-k2*x)+A1*k1*np.exp(-k1*x))
        return fun_value-f_value,fun_value_prime-f_value_prime
    k1,k2=(-0.1,-0.5)
    A1,A2=fsolve(lambda var : cut_off_f1(var,(innercut,value,value_prime,k1,k2)),[np.random.randn(),np.random.randn()])
    # print(A1,A2)
    var=(A1,A2)
    para=(innercut,value,value_prime,k1,k2)
    # print(cut_off_fun(innercut,var,para))
    return cut_off_fun(r,var,para)

def constant_pair_2(r,value,value_prime,innercut=1.5):
    np.seterr(divide='ignore')
    
    n=0.5
    def cut_off_fun(r,var,para):
        x,f_value,f_value_prime,k1,k2=para
        A1,A2=var
        # print(innercut)
        return 1/(r)**(n)*(A2*np.exp(-k2*r)+A1*np.exp(-k1*r))*H(x-r)+0*(x-r)*1/r**2*H(x-r)
    def cut_off_f1(var,para):
        x,f_value,f_value_prime,k1,k2=para
        A1,A2=var
        fun_value=(1/(x)**n)*(A2*np.exp(-k2*x)+A1*np.exp(-k1*x))
        fun_value_prime=-(n)/x**(n+1)*(A2*np.exp(-k2*x)+A1*np.exp(-k1*x))-(1/(x)**n)*(k2*A2*np.exp(-k2*x)+A1*k1*np.exp(-k1*x))
        return fun_value-f_value,fun_value_prime-f_value_prime
    k1,k2=(-0.1,-0.2)
    # print(value_prime)
    A1,A2=fsolve(lambda var : cut_off_f1(var,(innercut,value,value_prime,k1,k2)),[np.random.randn(),np.random.randn()])
    # print(A1,A2)
    var=(A1,A2)
    para=(innercut,value,value_prime,k1,k2)
    # print(cut_off_fun(innercut,var,para))
    return cut_off_fun(r,var,para)    



def constant_emb_1(rol):
    return -np.sqrt(rol)

def cut_off(x,innercut=1.5,left_value=None,left_value_prime=None,type=None):
    if type=='pair-1':
        if left_value==None or left_value_prime==None:
           ## 提出警告
              print('Warning: The left value and left value prime are not given')
              return 
        else:
            return constant_pair_2(x,left_value,left_value_prime,innercut=innercut)
    if type=='emb-1':
        return constant_emb_1(x)
    if type==None:
        return 0