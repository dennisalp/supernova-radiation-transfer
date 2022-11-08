'''
2020-10-16, Dennis Alp, dalp@kth.se

Takes data from get_dat and creates spectra.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter
import matplotlib.font_manager

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)



################################################################
# Parameters
# Constants, cgs
cc = 2.99792458e10 # cm s-1
hh = 6.6260755e-27
DD = 51.2 # kpc
pc = 3.086e18 # cm
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
kev2erg = 1.60218e-9 # erg keV-1
Msun = 1.989e33 # g
Lsun = 3.828e33 # erg s-1
uu = 1.660539040e-24 # g



################################################################
# Help functions
def d2s(d):
    return d*24.*60.*60.
def s2d(s):
    return s/(24.*60.*60.)

def fix_lim():
    if mode == 'lc':
        plt.ylim(bottom=0)
        plt.gca().grid(True, which='both', ls=":")
    elif sys.argv[2][:11] == 'z_spec_300d':
        plt.xlim([10, 100])
        plt.ylim([1e-5, 3e-4])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
def get_label(labs):
    ret = []
    for lab in labs:
        if lab == 'b15_01z_spec_300d' and sys.argv[2][:11] == 'z_spec_300d':
            ret.append('B15 ($Z_\mathrm{eff}=0.12$ Z$_{\mathrm{eff,}\odot}$)')
        else:
            ret.append(lab.replace('_', ' '))

    return ret

# def spec_smooth(ee, spec, sig_smo=0):
#     emin_smo = 50
#     post_sig = 1
#     if sys.argv[2] == 'oth_spec':
#         tmp = np.where(spec<1e-6, 1e-10, spec)
#         spec = np.where(ee<30, 10**gaussian_filter(np.log10(tmp),10), spec)
#         emin_smo = 50
#         lines_lim = 0.1
#         sig_smo = 5
#     else:
#         return spec

#     print('WARNING Spectra are smoothed!')
#     smo = gaussian_filter(spec, sig_smo)
#     lines = gaussian_filter(np.abs(smo/spec-1),3)
#     smo = np.where((lines>lines_lim) & (ee>emin_smo) & (ee<3100), spec, smo)
#     return  gaussian_filter(smo, post_sig)
    
# def lc_smooth(lc, tm2=0):
#     if sys.argv[2] == '87a_lc_45-105_kev':
#         print('WARNING Light curves are smoothed! This MUST be tuned by hand!')
#         return gaussian_filter(lc, 2)
#     else:
#         return lc
    
def spec_smooth(ee, spec):
    # ii = np.where((ee > 0.8*511) & (ee < 1.2*511))[0]
    nn = 10
    out = np.convolve(spec, np.ones(nn)/nn, 'same')
    # out[ii] = spec[ii]
    return out
    spec = np.where(spec<1e-10, 1e-10, spec)
    return 10**gaussian_filter(np.log10(spec),5)

def lc_smooth(lc, tm2=0):
    if sys.argv[2] == '87a_lc_45-105_kev':
        print('WARNING Light curves are smoothed! This MUST be tuned by hand!')
        return gaussian_filter(lc, 2)
    else:
        return lc

def plt_lc(ff, norm):
    dat = np.loadtxt('../dat/' + ff + '.txt')
    tt = s2d(dat[:,0])
    lc = dat[:,1]

    if 'dod_dir' in sys.argv[2]:
        if ii == 3:
            plt.plot(tt, lc, lw=2, color='k', zorder=100, drawstyle='steps-post')
        else:
            # order colors, hard coded
            # col_idx = (np.trapz(smo, tt)-0.2*743.42492230869625)/(1.*2159.757475649814-0.2*743.42492230869625)
            # col_idx = cm.get_cmap('viridis', 256)(col_idx)
            # plt.plot(tt, lc_smooth(smo,tmp), lw=1., color=col_idx)
            plt.plot(tt, lc, lw=2, color=viridis[ii-4], drawstyle='steps-post')
    else:
        plt.plot(tt, lc, lw=1, drawstyle='steps-post')
    
    plt.xlabel('Time (d)')
    plt.ylabel('Luminosity (erg~s$^{-1}$)')
    plt.xlim([tt[0], tt[-1]])    

def mk_log_lc():
    if mode == 'lc':
        out = '/Users/silver/box/phd/pro/mag/nus/fig/' + sys.argv[2] + '_log.pdf'
        plt.yscale('log')
        plt.xlim(left=3)
        plt.xscale('log')
        ylim = plt.gca().get_ylim()
        plt.ylim(bottom=1e-5*ylim[1], top=2*ylim[1])
        fig.savefig(out, bbox_inches='tight', pad_inches=0.1, dpi=300)    

def plt_spec(ff):
    dat = np.loadtxt('../dat/' + ff + '.txt')
    xx = dat[:,0]/1.e3
    yy = dat[:,1]
    yy = spec_smooth(xx, yy)
    
    if 'dod_dir' in sys.argv[2]:
        if ii == 3:
            plt.loglog(xx, yy, drawstyle='steps-post', lw=2, color='k', zorder=100)
        else:
            plt.loglog(xx, yy, drawstyle='steps-post', lw=2, color=viridis[ii-4])
    else:
        plt.loglog(xx, yy, drawstyle='steps-post', lw=2)

    plt.xlabel('Energy (keV)')
    plt.ylabel('Number luminosity density ($\gamma$~s$^{-1}$~keV$^{-1}$)')
    plt.gca().grid(True, which='both', ls=":")
    # plt.ylim([1e-7, 1e-3])




################################################################
# Main
nn = len(sys.argv)-3
print('Number of data sets:', nn)
mode = sys.argv[1]
viridis = cm.get_cmap('viridis', 256)
viridis = viridis(np.linspace(0, 1, nn-1))

fig = plt.figure(figsize=(5, 3.75))
for ii in range(3, len(sys.argv)):
    print(sys.argv[ii])
    if mode == 'spec':
        plt_spec(sys.argv[ii])
    elif mode == 'lc':
        plt_lc(sys.argv[ii], False)

if nn > 1 and not 'dod' in sys.argv[2]:
    plt.legend(get_label(sys.argv[3:]), loc='best')

fix_lim()
out = '/Users/silver/box/phd/pro/mag/nus/fig/' + sys.argv[2] + '.pdf'
fig.savefig(out, bbox_inches='tight', pad_inches=0.1, dpi=300)
mk_log_lc()

# plt.show()
