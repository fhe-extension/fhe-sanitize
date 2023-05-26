from sage.all import *
import math,ssl
from estimator import *
from estimator.nd import NoiseDistribution, stddevf
from estimator.lwe_parameters import LWEParameters


def bitSizeOf(x):
    try: return math.log2(x)
    except: return 0

def norm2 (dim,variance):
 return math.sqrt(2*math.pi*variance*dim)

def eLen(dim):
    v=math.log(2 * dim * (1+1.0/(2**(-80))))
    return v/(2*math.pi*math.pi)


n=538      # LWE dimension

N = 1024   #RLWE polynomial dimension
d=1         # dimension of TLWE
B = 2**5 #Gadget basis decomposition
ell= 7    #precision of the decomposition
Bks = 2**1     #keyswitch basis decomposition
t=11
dprime=1

q= B**ell   #ciphertext modulus
sds = 2**(-33.8)        # noise standard deviation
sde = 2**(-33.8)
sd_ks = 2**(-13.6)

var_bk= sde*sde
var_sk= sds*sds
var_ks = sd_ks*sd_ks

var_e=var_bk

norm_e= norm2(N*dprime,var_e)
norm_s= norm2(N*d,q*q*var_sk)

bound_var_r= (1/(q**2))*((1.0+ q*norm_e)**2)*eLen(dprime*N)
bound_var_eprime= (1/(q**2))*((1.0+norm_s)**2)*eLen(d*N)

var_r=2**(int(10*bitSizeOf(bound_var_r))/10)
var_eprime=var_r

bound_var_e2= var_eprime * norm_s**2 + bound_var_eprime

var_e2=2**(int(10*bitSizeOf(bound_var_e2))/10)

norm_e_bk= norm2(N*(d+1)*ell,var_e)
bound_var_cdot= (1/(q**2))*(1.0+B**2)*(1.0+q*norm_e_bk)**2*eLen((d+1)*ell*N)

var_cdot=2**(int(10*bitSizeOf(bound_var_cdot))/10)

vpk= var_e2 + norm_s**2*var_eprime + (q*norm_e)**2*var_r
vbr= n*norm_e_bk**2*var_cdot*q**2
vks_round= d*N*Bks**(-2*t)/4
vks_err= d*N*t*Bks*Bks*var_ks/4
vks= vks_round + vks_err

final_stdev= math.sqrt(var_cdot+vpk+vbr+vks)
proba_err= bitSizeOf(math.erfc(1/(4*math.sqrt(2)*final_stdev)))

dists = NoiseDistribution.DiscreteGaussian(stddev=sds*q)
diste = NoiseDistribution.DiscreteGaussian(stddev=sde*q)
params = LWEParameters(n=N,q=q,Xs=dists,Xe=diste)
LWE.estimate.rough(params)
LWE.estimate(params)




# correctness for message space Z_{t}
def correctness_lwe(t,stdev,mean): 
    #bound = q/(2.0 * t) - float(error_params['mean'])
    bound = 1/(2.0 * t) - float(mean)
    error_out = math.erfc(bound/(math.sqrt(2) * stdev)) 
    print("[stddev: " + str(stdev) + "]", "[mean: " + str(mean) + "]", "[proba error of output LWE: " + str(error_out) + "]")




#lemma 11 - bound_on_r=boound for min(r,r') 
def bound_Ginv (eps_lem11):
    dim=(d+1)*ell*N
    norm_of_e=norm2(dim,var_bk)
    print( "norm of BK error :" + str(norm_of_e))
    left_=(math.sqrt(1+B*B))*(q*norm_of_e+1)
    e_= 1.0/(eps_lem11)+1.0
    right_=math.sqrt((math.log(2*dim*e_)/math.pi))
    return left_*right_



def bound_for_vpk1():
    norm_e= norm2(N*d,var_e)
    bound1= (1/(q**2))*((1.0+ q*norm_e)**2)*eLen(dprime*N)
    print("[min(var_r and v_e2) " + str(min(var_r,var_e2)), "[bound is: " + str(bound1) + "]")

def bound_for_vpk2():
    norm_sk= norm2(N*dprime*q,var_sk)
    bound2= (1/(q**2))*((1.0+norm_sk)**2)*eLen(d*N)
    var_u= var_e2- norm_sk*norm_sk*var_eprime
    print("[min(var_eprime and var_u) " + str(min(var_eprime,var_u)), "[bound is: " + str(bound2) + "]")


def circuit_priv_sd():
    eeps_lem11 = 2**(-80)
    print("bound on (r,r')  " + str(bound_Ginv(eeps_lem11)))
    

bound_for_vpk1()
bound_for_vpk2()
circuit_priv_sd()


# variance of keyswitching
def var_keyswitch(vks):
    return (N/4.0)*d*Bks**(-2*t) + N*d*Bks*Bks*vks/(4.0) 

def vpk():
    norm_pow2_sk= norm2(N*dprime*q,var_sk)**2
    vpk=(var_e2+ 2*math.pi*N*dprime*(q**2) *var_r* var_e+ norm_pow2_sk*var_eprime)
    print("vpk is "+ str(vpk))
    return vpk

# VARIANCE AFTER SANITIZATION ALGORITHM #
def variance_cprivacy(n_dim,vbk,vks,vpk,r,t):
    varks=var_keyswitch(vks)
    var_output= r*r/(q*q*2*math.pi) +n_dim*(d+1)*ell*N*r*r*vbk + vpk +varks
    return var_output

def stdev_cprivacy(n_dim,vbk,vks,vpk,r,t):
    var_output=variance_cprivacy(n_dim,vbk,vks,vpk,r,t)
    return math.sqrt(var_output)


def variance_priv():
    #_m=350      #number of samples >= 310 or 350 for eps=2**(-130) - bound on number of samples
    a_ks=2**(-20) # for n=612, security ok
    v_ks=a_ks*a_ks
    #v_pk=v_bk
    v_=0
    _r=  919012
    print("Error variance after circuit private bootstrapping: "+ str(variance_cprivacy(n,var_bk,v_ks,v_,_r,t)))
    print("Error stdev after circuit bootstrapping: "+ str(stdev_cprivacy(n,var_bk,v_ks,v_,_r,t)))
   

 
variance_priv()

correctness_lwe(2,0.026,0) 
