'''
0000-00-00, Dennis Alp, dalp@kth.se

Overplot the sensitivity of some insturments ontop of the model predictions.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as units
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

# Constants, cgs
cc = 2.99792458e10 # cm s-1
hh = 6.6260755e-27 # erg s
DD = 51.2 # kpc
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
Lsun = 3.828e33 # erg s-1
uu = 1.660539040e-24 # g

def nu2lam(nu):
    return cc/nu

def nu2kev(nu):
    return hh*nu/kev2erg
    
def lam2kev(lam):
    return hh*cc/lam/kev2erg

def lam2nu(lam):
    return cc/lam

def kev2lam(ev):
    return hh*cc/ev/kev2erg

def kev2nu(ev):
    return ev/hh*kev2erg

def spec_smooth(ee, spec):
    print('WARNING Spectra are smoothed!')
    smo = gaussian_filter(spec, 5)
    lines = gaussian_filter(np.abs(smo/spec-1),3)
    smo = np.where((lines>0.06) & (ee>50) & (ee<3100), spec, smo)
    return  gaussian_filter(smo, 1)
    
def lc_smooth(lc):
    print('WARNING Light curves are smoothed! This MUST be tuned by hand!')
    return 10**gaussian_filter(np.log10(np.where(lc<10**-5.3, 10**-5.3, lc)), 3)

# 847 keV line sensitivity
lim_spi = 4e-5 # ph s-1 cm-2, from figure 18 of roques03, 3sig, 1 Ms
lim_eas = 5e-6 # ph s-1 cm-2, from table 1 of de_angelis17, 3sig, 1 Ms
NGC_6946 = 7000 # kpc, Fireworks Galaxy
non_stripped = 3000 # kpc
stripped = non_stripped # kpc
non_stripped_line = 200 # kpc
stripped_line = non_stripped_line # kpc
zz = z_at_value(cosmo.luminosity_distance, non_stripped*units.kpc)
print('Non-stripped at distance {0:.1f} Mpc (z = {1:.5f}).'.format(non_stripped/1000, zz))
zz = z_at_value(cosmo.luminosity_distance, stripped*units.kpc)
print('Stripped at distance {0:.1f} Mpc (z = {1:.5f}).'.format(stripped/1000, zz))
zz = z_at_value(cosmo.luminosity_distance, non_stripped_line*units.kpc)
print('Non-stripped (line) at distance {0:.1f} Mpc (z = {1:.5f}).'.format(non_stripped_line/1000, zz))
zz = z_at_value(cosmo.luminosity_distance, stripped_line*units.kpc)
print('Stripped (line) at distance {0:.1f} Mpc (z = {1:.5f}).'.format(stripped_line/1000, zz))
plt2b15 = (51.2/non_stripped)**2
plt2iib = (51.2/stripped)**2
plt2b15_line = (51.2/non_stripped_line)**2
plt2iib_line = (51.2/stripped_line)**2

#Integral: Figure 17 and 18 from \citet{roques03}, delta E/E=0.5, 1 Ms, 3 sigma
#NuSTAR: Figure 1 koglin05, delta E/E=0.5, 1 Ms, 3 sigma
#Chandra: Figure 6 takahashi10, delta E/E=0.5, 1 Ms, ? sigma
#e-ASTROGAM: Figure 1 de_angelis17, delta E/E=?, 50 ks, 3 sigma
cxo = np.loadtxt('/Users/silver/box/sci/lib/t/takahashi10/chandra_sensitivity.txt')
nus = np.loadtxt('/Users/silver/box/sci/lib/k/koglin05/nustar_sensitivity.txt')
spi = np.loadtxt('/Users/silver/box/sci/lib/r/roques03/integral_spi_sensitivity.txt')
eas = np.loadtxt('/Users/silver/box/sci/lib/d/de_angelis17/fig_5.txt')
eas[:,1] = eas[:,1]/eas[:,0]
eas[:,0] = nu2kev(eas[:,0])
eas[:,1] = eas[:,1]/(eas[:,0]*kev2erg)*1/(hh/kev2erg)
ea2 = np.loadtxt('/Users/silver/box/sci/lib/d/de_angelis17/fig_4.txt')
ea2[:,1] = ea2[:,1]/ea2[:,0]
ea2[:,0] = ea2[:,0]*1.e3
ea2[:,1] = ea2[:,1]/(ea2[:,0]**2*kev2erg)*((8*24*3600)/1e6)**0.5
ea3 = np.loadtxt('/Users/silver/box/sci/lib/d/de_angelis17/fig_7.txt')
ea3[:,1] = ea3[:,1]/(ea3[:,0]/1e6)
ea3[:,0] = ea3[:,0]/1.e3
ea3[:,1] = ea3[:,1]/(ea3[:,0]**2*kev2erg)*((365*24*3600)/1e6)**0.5

#iib_spec = [np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/iib_spec_50d.txt'),
#            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/iib_spec_125d.txt'),
#            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/iib_spec_250d.txt')]
iib_spec = [np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/iib_spec_100d.txt')]
iib_llc = np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/iib_llc847.txt')
b15_spec = [np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/b15_spec_300d.txt'),
            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/l15_spec_300d.txt'),
            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/w15_1z_spec_300d.txt')
#            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/ic_spec_50d.txt'),
#            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/w7_spec_20d.txt')
#            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/w7_spec_50d.txt'),
#            np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/w7_spec_70d.txt')
            ]
b15_llc = np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/b15_llc847.txt')
l15_llc = np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/l15_llc847.txt')
w15_llc = np.loadtxt('/Users/silver/box/phd/pro/sne/grt/bin/dat/w15_1z_llc847.txt')

################################################################
fig = plt.figure(figsize=(5, 3.75))
#colors = ['#1f77b4', '#ff7f0e', '#d62728', 'r', 'g', 'b', '#2ca02c']
colors = ['#1f77b4', '#ff7f0e', '#d62728', '#2ca02c']
for ii, b15 in enumerate(b15_spec[:]):
#    plt.loglog(b15[:,0]/1e3, (b15[:,3])*plt2b15, color='gray')
    plt.loglog(b15[:,0]/1e3, spec_smooth(b15[:,0]/1e3, np.sum(b15[:,1:5], 1))*plt2b15, color=colors[ii])
for iib in iib_spec:
#    plt.loglog(iib[:,0]/1e3, (iib[:,3])*plt2iib, color='gray')
    plt.loglog(iib[:,0]/1e3, spec_smooth(iib[:,0]/1e3, np.sum(iib[:,1:5], 1))*plt2iib, color=colors[ii+1])

#b15 = b15_spec[-1]
#plt.loglog(b15[:,0]/1e3, b15[:,3]*plt2b15, 'gray', zorder=-99)
plt.legend(['B15 ($Z_\mathrm{eff}=0.03$ Z$_{\mathrm{eff,}\odot}$)', 'L15 ($Z_\mathrm{eff}=0.30$ Z$_{\mathrm{eff,}\odot}$)', 'W15 ($Z_\mathrm{eff}=1.12$ Z$_{\mathrm{eff,}\odot}$)', 'IIb ($Z_\mathrm{eff}=0.36$ Z$_{\mathrm{eff,}\odot}$)'])

plt.loglog(cxo[:,0], cxo[:,1], 'k:')
plt.loglog(nus[:,0], nus[:,1], 'k--')
plt.loglog(spi[:,0], spi[:,1], 'k-')
#plt.loglog(eas[:,0], eas[:,1], 'k-.')
#plt.loglog(ea2[:,0], ea2[:,1], 'k-.')
plt.loglog(ea3[:,0], ea3[:,1], 'k-.')

plt.xlabel('Energy (keV)')
plt.ylabel('Flux density at ' + str(int(stripped/1e3)) + ' Mpc ($\gamma$~cm$^{-2}$~s$^{-1}$~keV$^{-1}$)')
plt.xlim([4, 4e3])
plt.ylim([1e-9, 1e-5])
fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/sensitivity_continuum.pdf', bbox_inches='tight', pad_inches=0.1)

################################################################
fig = plt.figure(figsize=(5, 3.75))
plt.semilogy(b15_llc[:,0]/(3600*24), lc_smooth(b15_llc[:,1])*plt2b15_line, color=colors[0])
plt.semilogy(l15_llc[:,0]/(3600*24), lc_smooth(l15_llc[:,1])*plt2b15_line, color=colors[1])
plt.semilogy(w15_llc[:,0]/(3600*24), lc_smooth(w15_llc[:,1])*plt2b15_line, color=colors[2])
#plt.semilogy(b15_llc[:,0]/(3600*24), (b15_llc[:,3])*plt2b15_line, color=colors[0])
#plt.semilogy(l15_llc[:,0]/(3600*24), (l15_llc[:,3])*plt2b15_line, color=colors[1])
#plt.semilogy(w15_llc[:,0]/(3600*24), (w15_llc[:,3])*plt2b15_line, color=colors[2])
plt.semilogy(iib_llc[:,0]/(3600*24), iib_llc[:,1]*plt2iib_line, color=colors[3])
plt.axhline(y=lim_spi, color='k')
plt.axhline(y=lim_eas, color='k', ls='-.')

#plt.legend(['B15 at ' + str(non_stripped_line) + ' kpc', 'IIb at ' + str(int(stripped_line/1e3)) + ' Mpc'])
plt.legend(['B15 ($Z_\mathrm{eff}=0.03$ Z$_{\mathrm{eff,}\odot}$)', 'L15 ($Z_\mathrm{eff}=0.30$ Z$_{\mathrm{eff,}\odot}$)', 'W15 ($Z_\mathrm{eff}=1.12$ Z$_{\mathrm{eff,}\odot}$)', 'IIb ($Z_\mathrm{eff}=0.36$ Z$_{\mathrm{eff,}\odot}$)'])

plt.xlim([0, 1000])
plt.ylim([2e-6, 4e-3])
plt.xlabel('Time (day)')
plt.ylabel('Flux at ' + str(int(stripped_line)) + ' kpc ($\gamma$~cm$^{-2}$~s$^{-1}$)')
fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/sensitivity_line.pdf', bbox_inches='tight', pad_inches=0.1)

plt.show()
db()
