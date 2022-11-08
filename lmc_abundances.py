'''
0000-00-00, Dennis Alp, dalp@kth.se

Description of this file.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

# Constants, cgs
cc = 2.99792458e10 # cm s-1
GG = 6.67259e-8 # cm3 g-1 s-2
hh = 6.6260755e-27 # erg s
DD = 51.2 # kpc
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
Lsun = 3.828e33 # erg s-1
Rsun = 6.957e10 # cm
Tsun = 5772 # K
uu = 1.660539040e-24 # g
SBc = 5.670367e-5 # erg cm-2 K-4 s-1
kB = 1.38064852e-16 # erg K-1

ato_mas = {
    'ar36': 35.96754510600*uu,
    'c12' : 12.00000000000*uu,
    'ca40': 39.96259098000*uu,
    'ca44': 43.95548180000*uu,
    'co56': 55.93983930000*uu,
    'cr48': 47.95403200000*uu,
    'fe52': 51.94811400000*uu,
    'fe56': 55.93493630000*uu,
    'he4' : 4.002603254150*uu,
    'mg24': 23.98504170000*uu,
    'n'   : 1.008664915880*uu,
    'ne20': 19.99244017540*uu,
    'ni56': 55.94213200000*uu,
    'o16' : 15.99491461956*uu,
    'p'   : 1.007276466879*uu,
    's32' : 31.97207100000*uu,
    'sc44': 43.95940280000*uu,
    'si28': 27.97692653250*uu,
    'ti44': 43.95969010000*uu,
    'x56' : 55.93493630000*uu}

# Note that number abundances are normalized to hydrogen whereas mass fractions are normalized!
all_masses = np.array([ato_mas['p'], ato_mas['he4'], ato_mas['c12'], ato_mas['o16'], ato_mas['ne20'], ato_mas['mg24'], ato_mas['si28'], ato_mas['s32'], ato_mas['ar36'], ato_mas['ca40'], ato_mas['ti44'], ato_mas['cr48'], ato_mas['fe52']])

ele = 'H He C O Ne Mg Si S Ar Ca Ti Cr Fe'.split()



################################################################
# From lodders03
lodders = [12, 10.899, 8.39, 8.69, 7.87, 7.55, 7.54, 7.19, 6.55, 6.34, 4.92, 5.65, 7.47]
#          H   He      C     O     Ne    Mg    Si    S     Ar    Ca    Ti    Cr    Fe
solar_num = 10**np.array(lodders)

#russell90, kurt98
He = np.average(np.array([10.96, 10.91]), weights=1/np.array([0.06, 0.05]))
#hunter07b, russell89, russell90, kurt98
C = np.average(np.array([7.75, 8.04, 7.66, 7.81]))
#kurt98, 2*russell90, 2*trundle07, hunter07b
O = np.average(np.array([8.37, 8.37, 8.25, 8.33, 8.40, 8.35]))
#kurt98, 2*russell90
Ne = np.average(np.array([7.55, 7.68, 7.50]))
#2*trundle07, hunter07b
Mg = np.average(np.array([7.06, 7.08, 7.05]))
#2*trundle07, hunter07b
Si = np.average(np.array([7.19, 7.21, 7.20]))
#2*russell90
S = np.average(np.array([6.87, 6.62]), weights=1/np.array([0.14, 0.23]))
#2*russell90
Ar = np.average(np.array([6.07, 6.59]), weights=1/np.array([0.25, 0.07]))
#russell89, russell90
Ca = np.average(np.array([5.89, 6.05]), weights=1/np.array([0.16, 0.03]))
#russell92
Ti = 4.81
#russell89, russell90
Cr = np.average(np.array([5.47, 5.35]), weights=1/np.array([0.14, 0.22]))
#russell89, russell90, 2*trundle07
Fe = np.average(np.array([7.23, 7.22, 7.23, 7.24]), weights=1/np.array([0.17, 0.09, 0.10, 0.12]))
lmc_num = 10**np.array([12, He, C, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe])

print('He {:4.2f}\nC  {:4.2f}\nO  {:4.2f}\nNe {:4.2f}\nMg {:4.2f}\nSi {:4.2f}\nS  {:4.2f}\nAr {:4.2f}\nCa {:4.2f}\nTi {:4.2f}\nCr {:4.2f}\nFe {:4.2f}'.format(lodders[1], lodders[2], lodders[3], lodders[4], lodders[5], lodders[6], lodders[7], lodders[8], lodders[9], lodders[10], lodders[11], lodders[12]))

print('He {:4.2f}\nC  {:4.2f}\nO  {:4.2f}\nNe {:4.2f}\nMg {:4.2f}\nSi {:4.2f}\nS  {:4.2f}\nAr {:4.2f}\nCa {:4.2f}\nTi {:4.2f}\nCr {:4.2f}\nFe {:4.2f}'.format(He, C, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe))

print('He {:4.2f}\nC  {:4.2f}\nO  {:4.2f}\nNe {:4.2f}\nMg {:4.2f}\nSi {:4.2f}\nS  {:4.2f}\nAr {:4.2f}\nCa {:4.2f}\nTi {:4.2f}\nCr {:4.2f}\nFe {:4.2f}'.format(He-lodders[1], C-lodders[2], O-lodders[3], Ne-lodders[4], Mg-lodders[5], Si-lodders[6], S-lodders[7], Ar-lodders[8], Ca-lodders[9], Ti-lodders[10], Cr-lodders[11], Fe-lodders[12]))

################################################################
# Cross sections
ver = np.loadtxt('sig_verner96.txt')
ene = np.logspace(0, np.log10(4e6), 501)
cs_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15]
sig_sol = []
norm_sol = (solar_num*all_masses).sum()
e0 = np.argmax(ene>3e4)-1
for ii, sn in enumerate(solar_num):
    plt.loglog(ene/1e3, sn*ver[cs_idx[ii], :]/norm_sol)
    sig_sol.append(sn*ver[cs_idx[ii], e0]/norm_sol)
plt.legend(ele)

plt.figure()
sig_lmc = []
norm = (lmc_num*all_masses).sum()
for ii, ln in enumerate(lmc_num):
    plt.loglog(ene/1e3, ln*ver[cs_idx[ii], :]/norm)
    sig_lmc.append(ln*ver[cs_idx[ii], e0]/norm)

plt.legend(ele)
sig_sol = np.array(sig_sol)
sig_lmc = np.array(sig_lmc)
print('Metallicity based on photoabsorption opacity at 30 keV = {0:7.5f} Z'.format(sig_lmc.sum()/sig_sol.sum()))

plt.figure()
plt.plot(lmc_num/solar_num)
plt.show()
db()
