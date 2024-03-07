# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 23:07:54 2023

@author: taomihog
"""

import numpy as np

#constants
Kelvin = 298.15
Joul = 4184
Avogadro = 6.022E+23
Boltzmann=0.0138065
kT=Kelvin*Boltzmann

#pseudo constatns
LPDS = 51.97#dsDNA persistence length
KDS = 1318#dsDNA elastic modulus
L0DS = 0.338#dsDNA contour length per bp
LPSS = 0.765#ssDNA persistence length
KSS = 470#ssDNA elastic modulus  
L0SS = 0.554#ssDNA contour length per nt

#experiment conditions
jds = 2200
salt_mM =100
kpillar = 0.07406 #z stiffness of AOT (Here, power = 200 mW)

def getEDNA(seq, salt_mM, kT):
    LUT_mapping = {
        'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3,
        'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
        'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11,
        'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15
    }
    LUTdH = [-7.28, -5.8, -5.21, -4.63, -8.96, -8.57, -9.66, -5.21, -8.16, -10.1, -8.57, -5.8, -8.31, -8.16, -8.96, -7.28]
    LUTdS = [-20.28, -14.46, -12.89, -11.62, -24.48, -22.3, -24.43, -12.89, -22.46, -25.96, -22.3, -14.46, -25.06, -22.46, -24.48, -20.28]
    LUTm = [0.145, 0.099, 0.07, 0.117, 0.091, 0.063, 0.132, 0.07, 0.155, 0.079, 0.063, 0.099, 0.091, 0.155, 0.091, 0.145]
    
    energy = np.zeros(len(seq), dtype=np.float)
    eff_salt = np.log(salt_mM * 0.001) / 298
    
    for i in range(len(seq) - 1):
        j = LUT_mapping.get(seq[i:i+2], -1)
        if j != -1:
            energy[i + 1] = energy[i] + LUTdH[j] - (LUTm[j] * eff_salt + LUTdS[j] * 0.001) * Kelvin
    
    energy = -kT * energy * 1e21 * Joul / Avogadro / Boltzmann / Kelvin
    print("energy:", energy[-1])
    return energy

with open('..\\data\\NEB_H5alpha_Accessory_colonization_factor_AcfD.txt', 'r') as file:
    sequence = file.read().replace('\n', '')
EDNA = getEDNA(sequence,salt_mM,kT)
print()

#%%utility function
#basic physics models, import from my c++ code
def alpha2phi_Odijk95(alpha, k0_eff): 
    if alpha < 0.25:
        return 0
    res = 1.0 - 0.5 / np.sqrt(alpha) + alpha / k0_eff
    if np.isnan(res):
        print ("alpha2phi_Odijk95", res)
    return res

def integ_alphadphi_Odijk95(alpha, k0_eff):
    res = alpha * alpha2phi_Odijk95(alpha, k0_eff) - \
        alpha - np.sqrt(alpha) + 0.5 * alpha * alpha / k0_eff
    if np.isnan(res):
        print ("integ_alphadphi_Odijk95", res)
    return res

def alpha2phi_Smith95(alpha, k0_eff):
    res = (np.cosh(2.0 * alpha) / np.sinh(2.0 * alpha) - 0.5 / alpha) + alpha / k0_eff
    if np.isnan(res):
        print ("alpha2phi_Smith95", alpha, res)
    return res

def integ_alphadphi_Smith95(alpha, k0_eff):
    #integ actually starts from 1, but it's OK since it is for partition function calculation
    res = alpha * alpha2phi_Smith95(alpha, k0_eff) -\
        0.5 * np.log(0.5 * np.sinh(2.0 * alpha) / alpha) + 0.5 * alpha * alpha / k0_eff
    if np.isnan(res):
        print ("integ_alphadphi_Smith95", res)
    return res

#energy and force from DNA model
def lz_ss (force, jss): #ssDNA's length per base
    return 2 * jss * L0SS * alpha2phi_Smith95(force * LPSS / kT, KSS * LPSS / kT)

def lz_ds (force): #dsDNA's length per base
    return jds * L0DS * alpha2phi_Odijk95(force * LPDS / kT, KDS * LPDS / kT)

def le_ss (force, jss):#function version of ssDNA's energy per bp:
    return 2.0 * jss * kT * L0SS * integ_alphadphi_Smith95(force * LPSS / kT, KSS * LPSS / kT) / LPSS

def le_ds (force):#function version of dsDNA's energy per bp:
    return jds * kT * L0DS * integ_alphadphi_Odijk95(force * LPDS / kT, KDS * LPDS / kT) / LPDS

def delta_ext(F,z0,jss):
    return z0-(F/kpillar + lz_ss(F, jss) + lz_ds(F));
#%%

tor_binary_search = 1.0e-7

fmin = 0.1
fmax = 500

def FindForce(jss,z0):
    f1 = fmin
    f2 = fmax

    y1 = delta_ext(f1,z0,jss)
    y10 = y1
    y2 = delta_ext(f2,z0,jss)
    y20 = y2

    if (y1 * y2 >= 0):
        if ( y1 > 0) :
            return fmin-1 #f is too small
        else:
            return fmax+1

    cnt = 0 #in case the root is not found
    while (cnt < 10000):
        cnt += 1
        
        fm = (f1 + f2) * 0.5
        ym = delta_ext(fm,z0,jss)
        
        if (abs(ym) <= tor_binary_search):
            return fm

        if (y1 < y2 and ym > 0 or y1 > y2 and ym < 0):
            f2 = fm
            y2 = ym
        elif (y1 < y2 and ym < 0 or y1 > y2 and ym > 0):
            f1 = fm
            y1 = ym
        
        else:
            print("error1: ", y1, y2 ,ym, jss,z0, y10, y20)
            return np.inf #means some weird error
    
    print("error2: ", jss,z0)
    return np.inf

vectorized_FindForce = np.vectorize(FindForce)


def getEtot (jss,force):
    if force <= fmin or force >= fmax:
        return np.inf
    
    temp = EDNA[jss] + 0.5 * force * force / kpillar + le_ss (force, jss) + le_ds (force)
    if (np.isnan(temp)):
        print("Nan jss and force: ", jss, force)
    return temp

vectorized_getEtot = np.vectorize(getEtot)
    
    
z0max = 5100
z0min = 600
z0step = 500
z0_res = np.arange(z0min, z0max, z0step)
f_res = np.zeros((z0_res.size))

jssmax = EDNA.size
jss_arr = np.arange(EDNA.size)
e_arr = np.zeros((jssmax))
f_arr = np.zeros((jssmax))

for i,z0 in enumerate(z0_res):
    f_arr = vectorized_FindForce(jss_arr, z0)
    e_arr = vectorized_getEtot(jss_arr, f_arr)
    
    emin = np.min(e_arr)
    e_arr -= emin
    
    p_arr = np.exp(-e_arr / kT)
    fp_arr = f_arr * p_arr
    #jp_arr = jss_arr * p_arr
    #...
    psum = np.sum(p_arr)
    fpsum = np.sum(fp_arr)
    f_res[i] = fpsum / psum
    print(z0, f_res[i])

z0_res = z0_res - f_res / kpillar

import matplotlib.pyplot as plt

fig,ax = plt.subplots(figsize=(16, 4))
ax.plot(z0_res, f_res, '-b')
import pandas as pd
pyResult = pd.read_csv("standard.txt",delimiter='\t')
ax.plot(pyResult['# zMean'],pyResult['FMean'], "k", linewidth = 2, label = "python code")
plt.xlabel('extension, nm')
plt.ylabel('force, pN')
plt.show()
