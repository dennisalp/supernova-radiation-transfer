'''
2019-06-10, Dennis Alp, dalp@kth.se

Show energy fractions of direct line emission, down-scattered continuum, and absorption. Requested by Anders Jerkstrand.

mod=b15_lmc
python src/energy_budget.py ${mod}_lum 0 1000 ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
mod=n20_lmc
python src/energy_budget.py ${mod}_lum 0 1000 ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
mod=l15
python src/energy_budget.py ${mod}_lum 0 1000 ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
mod=w15
python src/energy_budget.py ${mod}_lum 0 1000 ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
mod=iib
python src/energy_budget.py ${mod}_lum 0 1000 ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
mod=m157b-2
python src/energy_budget.py ${mod}_lum 0 1000 ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
mod=m167b-2
python src/energy_budget.py ${mod}_lum 0 1000 ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)



################################################################
# Parameters
Msun = 1.989e33 # g
uu = 1.660539040e-24 # g
kev2erg = 1.60218e-9 # erg keV-1
tau_ni56 = 6.075/np.log(2) # milne04, junde99, junde11 nuclear data sheets
tau_co56 = 77.236/np.log(2) # milne04, junde99, junde11 nuclear data sheets
tau_ti44 = 58.9*365.25/np.log(2) # ahmad06
tau_ni57 = 35.60/24./np.log(2) # bhat98
tau_co57 = 271.74/np.log(2) # bhat98
ratio_57_56 = 2*0.02309436100878436 # twice (kurfess92, fransson93) solar (lodders03), table 6
nlc = 101
ph_per_decay = np.array([4.80106, 3.213, 2.90798, 1.056681])



################################################################
ene_ti44 = np.array([67.875, 78.337, 146.212, 1157.031, 1499.43,2144.2, 2656.41, 3301.2])
p_ti44 = np.array([94.4e-2, 96e-2, 0.089e-2, 99.9e-2, 0.912e-2, 0.0069e-2, 0.115e-2, 0.0031e-2])

ene_ni56 = np.array([158.38, 269.50, 480.44, 749.95,811.85,1561.80])
p_ni56 = np.array([98.8e-2, 36.5e-2, 36.5e-2, 49.5e-2, 86.0e-2, 14.0e-2])

ene_co56 = np.array([263.410, 411.380, 486.540, 655.000, 674.700, 733.511, 787.742, 846.771, 852.780, 896.531, 977.373, 997.330, 1037.840, 1089.030, 1140.404, 1159.920, 1175.102, 1198.780, 1238.282, 1272.200, 1335.589, 1360.215, 1442.750, 1462.340, 1640.404, 1771.351, 1810.772, 1963.714, 2015.181, 2034.755, 2113.123,  2212.933, 2276.360, 2373.700, 2598.459, 2657.450, 3009.596, 3201.962, 3253.416, 3272.990, 3451.152, 3547.930, 3600.700, 3611.800, 510.9989461])
p_co56 = np.array([0.022e-2, 0.025e-2, 0.055e-2, 0.038e-2, 0.038e-2, 0.200e-2, 0.310e-2, 100.000e-2, 0.050e-2, 0.070e-2, 1.430e-2, 0.112e-2, 13.990e-2, 0.050e-2, 0.150e-2, 0.100e-2, 2.270e-2, 0.050e-2, 67.600e-2, 0.020e-2, 0.120e-2, 4.330e-2, 0.200e-2, 0.077e-2, 0.060e-2, 15.600e-2, 0.640e-2, 0.720e-2, 3.080e-2, 7.880e-2, 0.380e-2, 0.350e-2, 0.110e-2, 0.080e-2, 17.280e-2, 0.025e-2, 1.049e-2, 3.240e-2, 7.930e-2, 1.889e-2, 0.953e-2, 0.198e-2, 0.018e-2, 0.010e-2, 38.00e-2])

ene_co57 = np.array([14.41300, 122.0614, 136.4743, 230.29, 339.54,352.36, 366.75, 569.92, 692.03, 706.40])
p_co57 = np.array([9.16e-2, 85.60e-2, 10.68e-2, 0.0004e-2, 0.0139e-2,0.0132e-2, 0.0013e-2, 0.017e-2, 0.157e-2, 0.0253e-2])



################################################################
# Parse input
print('Expects arguments: lab t_min t_max 44ti/sc 56ni 56co 57co\nPass "none" to skip isotope.')
if not len(sys.argv) == 8:
    print('ERROR. Wrong number of input arguments')
    sys.exit()
    
lab = sys.argv[1] # label for the output
t_min = float(sys.argv[2])
t_max = float(sys.argv[3])

fid = []
for ii in range(4, 8):
    fid.append(sys.argv[ii]) # ID of dataset
    if not sys.argv[ii] == 'none':
        mod = '_'.join(sys.argv[ii].split('_')[:-1])
print('MODEL: ' + mod)    



################################################################
# Help functions
def d2s(d):
    return d*24.*60.*60.
def s2d(s):
    return s/(24.*60.*60.)
 
def print_lc(out, meta, dat, wei, lum):
    loc_wei = dat['ee']/1e3*kev2erg*wei
    idx = (dat['ns'] == 0)

    out[:, int(meta['src_abu'])+1] = lum # (erg~s$^{-1}$)

    flx, _ = np.histogram(dat['tt'][idx], tt, weights=loc_wei[idx])
    out[:, int(meta['src_abu'])+1+4] = flx/np.diff(tt) # (erg~s$^{-1}$)

    flx, _ = np.histogram(dat['tt'][~idx], tt, weights=loc_wei[~idx])
    out[:, int(meta['src_abu'])+1+2*4] = flx/np.diff(tt) # (erg~s$^{-1}$)
    
def bateman_equation(N0, t1, t2, tt):
    t1, t2 = d2s(t1), d2s(t2)
    return N0/(t1-t2)*(np.exp(-tt/t1)-np.exp(-tt/t2))

def get_meta(ff):
    meta = {}
    with open('out/' + ff + '.txt') as ff:
        for line in ff:
            (key, val) = line.split('\t')
            meta[key] = val.strip()
    return meta

def get_ini_num(meta):
    if meta['src_abu']=='0':
        return 1.5e-4*(51.2/50)**2*Msun/(43.9596901*uu)
    elif meta['src_abu']=='1' or meta['src_abu']=='2':
        return 0.07*Msun/(55.942132*uu)
    elif meta['src_abu']=='3':
        return 0.07*Msun/(55.942132*uu)*ratio_57_56

def get_wei(meta, dat, ini_num, tt, dt):
    if meta['src_abu'] == '0':
        wei = dt*ini_num*np.exp(-np.array(dat['ti'])/d2s(tau_ti44))/d2s(tau_ti44)
        avg_ene = np.average(ene_ti44, weights=p_ti44)*kev2erg
        lum = avg_ene*ph_per_decay[0]*ini_num*np.exp(-(tt[:-1]+tt[1:])/2/d2s(tau_ti44))/d2s(tau_ti44)
    elif meta['src_abu'] == '1':
        wei = dt*ini_num*np.exp(-np.array(dat['ti'])/d2s(tau_ni56))/d2s(tau_ni56)
        avg_ene = np.average(ene_ni56, weights=p_ni56)*kev2erg
        lum = avg_ene*ph_per_decay[1]*ini_num*np.exp(-(tt[:-1]+tt[1:])/2/d2s(tau_ni56))/d2s(tau_ni56)
    elif meta['src_abu'] == '2':
        wei = dt*bateman_equation(ini_num, tau_ni56, tau_co56, np.array(dat['ti']))
        avg_ene = np.average(ene_co56, weights=p_co56)*kev2erg
        lum = avg_ene*ph_per_decay[2]*bateman_equation(ini_num, tau_ni56, tau_co56, (tt[:-1]+tt[1:])/2)
    elif meta['src_abu'] == '3':
        wei = dt*bateman_equation(ini_num, tau_ni57, tau_co57, np.array(dat['ti']))
        avg_ene = np.average(ene_co57, weights=p_co57)*kev2erg
        lum = avg_ene*ph_per_decay[3]*bateman_equation(ini_num, tau_ni57, tau_co57, (tt[:-1]+tt[1:])/2)
    return wei, lum


################################################################
tt = d2s(np.linspace(t_min, t_max, nlc))
out = np.zeros((nlc-1, 1+3*4))
out[:, 0] = tt[:-1]

for ii, ff in enumerate(fid):
    if ff == 'none':
        continue

    print(ff)
    dat = pd.read_csv('out/' + ff + '.esc', delim_whitespace=True, names=['ei', 'ee', 'ox', 'oy', 'oz', 'ns', 'tt', 'ti'])
    meta = get_meta(ff)
    print(meta)
    if not int(meta['src_abu']) == ii:
        print('ERROR. File order does not match src_abu.')
        sys.exit(1)
    ini_num = get_ini_num(meta)

    # Compute photons packet^-1 (weights)
    n_dec = float(meta['n_pho'])/ph_per_decay[int(meta['src_abu'])]
    m_tmin = float(meta['tmin'])
    m_tmax = float(meta['tmax'])
    dt = d2s((m_tmax-m_tmin)/n_dec)
    wei, lum = get_wei(meta, dat, ini_num, tt, dt)

    header = 'From Alp et al. (2019). First column: Time at the early side of the bin in seconds. Then follows 3 sets of 4 luminosities in units of erg s^-1. The 3 sets are total luminosity, direct line luminosity, and down-scattered continuum luminosity. The 4 luminosities are for different isotopes; 44Ti+44Sc, 56Ni, 56Co, 57Co'
    print_lc(out, meta, dat, wei, lum)



################################################################
tt = s2d(out[:,0])
np.savetxt('dat/' + lab + '.txt', out, header = header)
fig = plt.figure(figsize=(5, 3.75))

tot = np.sum(out[:,1:5], axis=1)
lin = np.sum(out[:,5:9], axis=1)
con = np.sum(out[:,9:], axis=1)
pha = tot-lin-con

plt.semilogy(tt, tot)
plt.semilogy(tt, lin)
plt.semilogy(tt, con)
plt.semilogy(tt, pha)

plt.xlabel('Time (d)')
plt.ylabel('Luminosity (erg~s$^{-1}$)')
plt.legend(['Total', 'Direct line emision', 'Down-scattered continuum', 'Deposited'])

fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/' + lab + '.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
#plt.show()
#db()
