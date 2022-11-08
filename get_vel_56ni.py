'''
2019-07-09, Dennis Alp, dalp@kth.se

Get the mass-weighted average radial velocity of the fastest 1% of
  56Ni. Response to referee report.
'''

import os
import sys
import pdb
from glob import glob

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

def get_dat(dat, lab):
    tmp = np.array(dat[lab])
    tmp[:,:,-1] = tmp[:, :, -2]
    return tmp

# Gravitational to baryonic mass using Eq. (36) of lattimer01
def mg2mb(Mg, RR, Mb=0.):
    beta = GG*Mg*Msun/(RR*1e5*cc**2)
    Eb = Mg*0.6*beta/(1-0.5*beta)
    return Mg+Eb-Mb

def print_masses(mas, leg):
    for ii, ll in enumerate(leg):
        print(ll, mas[ii])

# Constants, cgs
cc = 2.99792458e10 # cm s-1
hh = 6.6260755e-27 # erg s
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
GG = 6.67259e-8 # cm3 g-1 s-2
Lsun = 3.828e33 # erg s-1
uu = 1.660539040e-24  # g


################################################################
# Parameter
dm = 0.2 # Msun
dv = 100 # km s-1
ff = sys.argv[2]
path = sys.argv[1] + ff + '.h5'
print ('Current file:', ff)

tt = {'B15': 13489751.4508,
      'B15_lmc': 13489751.4508,
      'B15_87a': 13489751.4508,
      'B15_lmc_min': 13489751.4508,
      'B15_lmc_max': 13489751.4508,
      'M157b-2_min': 86401.5,
      'M157b-2_max': 86401.5,
      'N20': 12488928.8787,
      'N20_lmc': 12488928.8787,
      'L15': 12630166.8468,
      'W15': 12755626.5895,
      'IIb': 1565943.04402,
      'BHH': 13489751.4508,
      'B1D': 13489751.4508,
      'BEO': 13489751.4508,
      'BEN': 13489751.4508,
      'HMM': 14832200.0,
      'B15_1d': 85858.3109787,
      'M15': 3785.8960366549186}

dat = h5py.File(path, 'r')
rad = np.array(dat['radius']).astype('double')
tht = np.array(dat['theta']).astype('double')
phi = np.array(dat['phi']).astype('double')
vex = get_dat(dat, 'vex').astype('double')
den = get_dat(dat, 'den').astype('double')


################################################################
# Check volume and density
ato_idx = {'ar36':  0,
           'c12' :  1,
           'ca40':  2,
           'ca44':  3,
           'co56':  4,
           'cr48':  5,
           'fe52':  6,
           'fe56':  7,
           'he4' :  8,
           'mg24':  9,
           'n'   : 10,
           'ne20': 11,
           'ni56': 12,
           'o16' : 13,
           'p'   : 14,
           's32' : 15,
           'sc44': 16,
           'si28': 17,
           'ti44': 18,
           'x56' : 19}

ato_mas = uu*np.array([35.96754510600, 12.00000000000, 39.96259098000, 43.95548180000, 55.93983930000, 47.95403200000, 51.94811400000, 55.93493630000, 4.002603254150, 23.98504170000, 1.008664915880, 19.99244017540, 55.94213200000, 15.99491461956, 1.007276466879, 31.97207100000, 43.95940280000, 27.97692653250, 43.95969010000, 55.93493630000])

mid_rad = (rad[:-1]+rad[1:])/2.
mid_phi = (phi[:-1]+phi[1:])/2.
mid_tht = (tht[:-1]+tht[1:])/2.

mid_phi = mid_phi[:, np.newaxis, np.newaxis]
mid_tht = mid_tht[np.newaxis, :, np.newaxis]
mid_rad = mid_rad[np.newaxis, np.newaxis, :]

dif_phi = np.diff(phi)[:, np.newaxis, np.newaxis]
dif_tht = np.diff(tht)[np.newaxis, :, np.newaxis]
dif_rad = np.diff(rad)[np.newaxis, np.newaxis, :]

# Print composition
vol = mid_rad**2*np.sin(mid_tht)*dif_rad*dif_tht*dif_phi
mas = np.sum(vol*den)/Msun
print('Check volume and (ejecta) mass:', np.sum(vol)/(4*np.pi*(rad[-1]**3-rad[0]**3)/3.), mas)

# Compute some useful quantities
if ff in tt.keys():
    rad2vel = 1/(tt.get(ff)*1e5)
else:
    rad2vel = 1/(dat['time'][0]*1e5)

Ek = vol*den*(mid_rad*rad2vel*1e5)**2/2.
print('Energy (homologous):', np.sum(Ek))
Ek = vol*den*vex**2/2.
print('Energy (vex):', np.sum(Ek))
#pdb.set_trace()
    
den_ele = den*(get_dat(dat, 'fe56')+get_dat(dat, 'co56')+get_dat(dat, 'ni56')+0.5*get_dat(dat, 'x56' ))
mas_tot = np.sum(vol*den_ele, axis=(0, 1))/Msun
mas_cum = np.cumsum(mas_tot)
mas_coo = np.cumsum(np.sum(vol*den, axis=(0, 1))/Msun)

################################################################
# Get the mass-weighted average radial velocity of the fastest 1% of 56Ni
ii = np.where(mas_cum > 0.99*mas_cum[-1])
vel = mid_rad[0,0]*rad2vel
vel_coo = np.average(vel[ii], weights=mas_tot[ii])
mas_coo = mas_coo[np.argmax(vel > vel_coo)]
print(ff + ' <v>_1%(56Ni):{0:9.3f} km s-1,{1:7.3f} Msun'.format(vel_coo, mas_coo))
#pdb.set_trace()
