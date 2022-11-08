'''
2018-10-03, Dennis Alp, dalp@kth.se

Takes data from get_dat and creates spectra.
'''

from __future__ import division, print_function
import os
import pdb
import sys
from glob import glob

import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter

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
v_lmc = 287 # km s-1

tau_ni56 = 6.075/np.log(2) # milne04, junde99, junde11 nuclear data sheets
tau_co56 = 77.236/np.log(2) # milne04, junde99, junde11 nuclear data sheets
tau_ti44 = 58.9*365.25/np.log(2) # ahmad06
tau_ni57 = 35.60/24./np.log(2) # bhat98
tau_co57 = 271.74/np.log(2) # bhat98
ratio_57_56 = 2*0.02309436100878436 # twice (kurfess92, fransson93) solar (lodders03), table 6

################################################################
# Help functions
def d2s(d):
    return d*24.*60.*60.
def s2d(s):
    return s/(24.*60.*60.)

def e2v(ee, e0):
    return -(ee/e0-1.)*cc/1e5

def gauss(ee, mag, mu, sig):
    return mag/(sig*np.sqrt(2*np.pi))*np.exp(-0.5*((ee-mu)/sig)**2)

def delmc(xx):
    return xx/(1-v_lmc*1e5/cc)

def fwhm2sig(fwhm):
    return fwhm/(2*np.sqrt(2*np.log(2)))

def bateman_equation(N0, t1, t2, tt):
    t1, t2 = d2s(t1), d2s(t2)
    return N0/(t1-t2)*(np.exp(-tt/t1)-np.exp(-tt/t2))

def get_ini_num(nickel_mass=0.07):
    return nickel_mass*Msun/(55.94213200000*uu)

def fix_lim():
    if sys.argv[2][:11] == 'z_spec_300d':
        plt.xlim([10, 100])
        plt.ylim([1e-5, 3e-4])
    elif sys.argv[2] == 'oth_llc847':
        plt.xlim([0, 1000])
        plt.ylim([1e-4, 1e-1])
    elif sys.argv[2] == 'oth_spec':
        plt.xlim([10, 4000])
        plt.ylim([1e-6, 4e-3])
    elif sys.argv[2] == '87a_llc847+1238' or sys.argv[2] == 'b15_m157b-2_llc847+1238':
        plt.gca().set_yscale("linear", nonposy='clip')
        plt.ylim([-0.1, 5])
        plt.xlim([0, 1000])
    elif sys.argv[2] == 'b15_llc847+1238':
#        plt.gca().set_yscale("log", nonposy='clip')
#        plt.ylim([1e-4, 1e-2])
        plt.gca().set_yscale("linear", nonposy='clip')
        plt.ylim([-1e-4, 5e-3])
        plt.xlim([0, 1000])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    elif sys.argv[2] == '87a_spec_300d':
        plt.xlim([10, 500])
        plt.ylim([1e-6, 3e-4])
    elif sys.argv[2] == 'b15_lmc_spec_300d' or sys.argv[2] == 'b15_m157b-2_spec_300d':
        plt.xlim([10, 4000])
        plt.ylim([1e-6, 3e-4])
    elif sys.argv[2] == 'sample_spec':
        plt.xlim([10, 4000])
        plt.ylim([1e-6, 3e-4])
    elif sys.argv[2] == 'b15_lmc_lc_45-105_kev':
        plt.xlim([0, 1000])
        plt.ylim([-2e-6, 1e-4])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    elif sys.argv[2] == '87a_lc_45-105_kev' or sys.argv[2] == 'b15_m157b-2_lc_45-105_kev':
        plt.xlim([0, 1000])
        plt.ylim([-0.3, 8.3])        
    elif sys.argv[2] == 'lmc_87a_lc_45-105_kev':
        plt.xlim([0, 1000])
        plt.ylim([0, 9])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    elif '_lc_45-105_kev_dod_dir' in sys.argv[2]:
        plt.xlim([0, 1000])
        plt.ylim([0, 8])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    elif sys.argv[2] == 'm157b_all_lc_45-105_kev':
        plt.xlim([0, 1000])
        plt.ylim([0, 8])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    elif sys.argv[2] == 'm157b_all_llc847+1238':
        plt.ylim([0, 1.7])
        plt.xlim([0, 1000])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        
def get_label(labs):
    ret = []
    for lab in labs:
        if lab == 'b15_01z_spec_300d' and sys.argv[2][:11] == 'z_spec_300d':
            ret.append('B15 ($Z_\mathrm{eff}=0.12$ Z$_{\mathrm{eff,}\odot}$)')
        elif lab == 'b15_025z_spec_300d' and sys.argv[2][:11] == 'z_spec_300d':
            ret.append('B15 ($Z_\mathrm{eff}=0.27$ Z$_{\mathrm{eff,}\odot}$)')
        elif lab == 'b15_87a_spec_300d' and sys.argv[2][:11] == 'z_spec_300d':
            ret.append('B15 ($Z_\mathrm{eff}=0.28$ Z$_{\mathrm{eff,}\odot}$)')
        elif lab == 'b15_lmc_spec_300d' and sys.argv[2][:11] == 'z_spec_300d':
            ret.append('B15 ($Z_\mathrm{eff}=0.55$ Z$_{\mathrm{eff,}\odot}$)')
        elif lab == 'b15_1z_spec_300d' and sys.argv[2][:11] == 'z_spec_300d':
            ret.append('B15 ($Z_\mathrm{eff}=1.00$ Z$_{\mathrm{eff,}\odot}$)')
        elif lab == 'l15_llc847':
            ret.append('L15')
        elif lab == 'w15_llc847':
            ret.append('W15')
        elif lab == 'iib_llc847':
            ret.append('IIb')
        elif lab == 'l15_spec_300d':
            ret.append('L15 300 d')
        elif lab == 'w15_spec_300d':
            ret.append('W15 300 d')
        elif lab == 'iib_spec_125d':
            ret.append('IIb 125 d')
        elif lab == 'iib_spec_100d':
            ret.append('IIb 100 d')
        elif lab == 'b15_lmc_llc847+1238':
            ret.append('B15')
        elif lab == 'n20_lmc_llc847+1238':
            ret.append('N20')
        elif lab == 'hmm_llc847+1238':
            ret.append('10HMM')
        elif lab == 'b15_lmc_spec_300d' and sys.argv[2] == 'oth_spec':
            ret.append('B15 300 d')
        elif lab[:17] == 'b15_lmc_spec_300d':
            ret.append('B15')
        elif lab == 'n20_lmc_spec_300d':
            ret.append('N20')
        elif lab == 'b15_lmc_llc847':
            ret.append('B15')
        elif lab == 'hmm_spec_300d':
            ret.append('10HMM')
        elif lab == 'b15_lmc_lc_45-105_kev' and sys.argv[2] == 'lmc_87a_lc_45-105_kev':
            ret.append('B15 LMC')
        elif lab == 'b15_87a_lc_45-105_kev' and sys.argv[2] == 'lmc_87a_lc_45-105_kev':
            ret.append('B15 SN 1987A')
        elif lab == 'b15_lmc_lc_45-105_kev':
            ret.append('B15')
        elif lab == 'n20_lmc_lc_45-105_kev':
            ret.append('N20')
        elif lab == 'hmm_lc_45-105_kev':
            ret.append('10HMM')
        elif lab[:15] == 'm157b-2_spec_300' or lab[:16] == 'm157b-2_spec_300':
            ret.append('M15-7b')
        elif lab[:20] == 'm157b-2_lc_45-105_kev' or lab[:21] == 'm157b-2_lc_45-105_kev':
            ret.append('M15-7b')
        elif lab[:18] == 'm157b-2_llc847+1238' or lab[:19] == 'm157b-2_llc847+1238':
            ret.append('M15-7b')
        elif lab[:15] == 'm167b+_spec_300' or lab[:16] == 'm167b-2_spec_300':
            ret.append('M16-7b')
        elif lab[:20] == 'm167b+_lc_45-105_kev' or lab[:21] == 'm167b-2_lc_45-105_kev':
            ret.append('M16-7b')
        elif lab[:18] == 'm167b+_llc847+1238' or lab[:19] == 'm167b-2_llc847+1238':
            ret.append('M16-7b')
        else:
            ret.append(lab.replace('_', ' '))

    if sys.argv[2] == 'oth_llc847':
        return ret[:4]
    elif sys.argv[2] == 'b15_m157b-2_llc847+1238':
        return ret[:2]
    elif sys.argv[2] == 'm157b_all_lc_45-105_kev' or sys.argv[2] == 'm157b_all_llc847+1238':
        return ['$1.40\\times{}10^{51}$~erg', '$1.43\\times{}10^{51}$~erg', '$1.44\\times{}10^{51}$~erg', '$1.76\\times{}10^{51}$~erg']
    elif 'm157b_all' in sys.argv[2]:
        return ret[:4]
    elif sys.argv[2] == 'lmc_87a_lc_45-105_kev':
        return ret[:2]
    
    return ret

def spec_smooth(ee, spec, sig_smo=0):
    emin_smo = 50
    post_sig = 1
    if sys.argv[2] == 'oth_spec':
        tmp = np.where(spec<1e-6, 1e-10, spec)
        spec = np.where(ee<30, 10**gaussian_filter(np.log10(tmp),10), spec)
        emin_smo = 50
        lines_lim = 0.1
        sig_smo = 5
    elif sys.argv[2][:11] == 'z_spec_300d':
        sig_smo = 20
        lines_lim = 0.06
    elif sys.argv[2] == '87a_spec_300d':
        sig_smo = 20 if sig_smo == 0 else sig_smo
        lines_lim = 99999
    elif sys.argv[2] == 'b15_lmc_spec_300d' or sys.argv[2] == 'b15_m157b-2_spec_300d' or sys.argv[2] == 'sample_spec':
        sig_smo = 10 if sig_smo == 0 else sig_smo
        lines_lim = 0.06 if sig_smo == 0 else 0.32
        post_sig = post_sig if sig_smo == 0 else 2.5
    elif 'dod_dir' in sys.argv[2]:
        sig_smo = 10 if sig_smo == 0 else sig_smo
        lines_lim = 0.06 if sig_smo == 0 else 0.32
        post_sig = post_sig if sig_smo == 0 else 2.5        
    else:
        return spec
        sig_smo = 5
        lines_lim = 0.06

    print('WARNING Spectra are smoothed!')
    smo = gaussian_filter(spec, sig_smo)
    lines = gaussian_filter(np.abs(smo/spec-1),3)
    smo = np.where((lines>lines_lim) & (ee>emin_smo) & (ee<3100), spec, smo)
    return  gaussian_filter(smo, post_sig)
    
def lc_smooth(lc, tm2=0):
    if sys.argv[2] == '87a_lc_45-105_kev':
        print('WARNING Light curves are smoothed! This MUST be tuned by hand!')
        return gaussian_filter(lc, 2)
    elif sys.argv[2] == 'b15_lmc_lc_45-105_kev':
        print('WARNING Light curves are smoothed! This MUST be tuned by hand!')
        return gaussian_filter(lc, 2)
    elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev':
        print('WARNING Light curves are smoothed! This MUST be tuned by hand!')
        return gaussian_filter(lc, 2)
    elif sys.argv[2] == 'oth_llc847':
        tmp = -4.6
        tm2 = 1 if tm2 == 0 else tm2
    elif sys.argv[2] == '87a_llc847+1238' or sys.argv[2] == 'b15_m157b-2_llc847+1238':
        tmp = -9.9 if tm2 == 0 else -9.9
        tm2 = 1 if tm2 == 0 else tm2
        lc[:10] = gaussian_filter(lc[:10], 2)
    elif sys.argv[2] == 'INSERT MISSING HERE':
        tmp = -5
        tm2 = 3
    elif 'dod_dir' in sys.argv[2]:
        return gaussian_filter(lc, tm2)
    elif 'm157b_all_lc_45-105_kev' in sys.argv[2]:
        return gaussian_filter(lc, tm2)
    elif 'm157b_all_llc847+1238' in sys.argv[2]:
        return gaussian_filter(lc, tm2)
    elif 'lmc_87a_lc_45-105_kev' in sys.argv[2]:
        return gaussian_filter(lc, tm2)
    else:
        print('asddsaasd')
        return lc

    print('WARNING Light curves are smoothed! This MUST be tuned by hand!')
    return 10**gaussian_filter(np.log10(np.where(lc<10**tmp, 10**tmp, lc)), tm2)

#def get_ini_num2():
#    if 'b15_1d' in sys.argv[ii]:
#        path = 'kicks/B15_1d_ini_num.csv'
#    elif 'b15_lmc' in sys.argv[ii]:
#        path = 'kicks/B15_LMC_ini_num.csv'
#    elif 'b15_01z' in sys.argv[ii]:
#        path = 'kicks/B15_01Z_ini_num.csv'
#    elif 'b15_025z' in sys.argv[ii]:
#        path = 'kicks/B15_025Z_ini_num.csv'
#    elif 'b15_1z' in sys.argv[ii]:
#        path = 'kicks/B15_1Z_ini_num.csv'
#    elif 'b15_km' in sys.argv[ii]:
#        path = 'kicks/B15_km_ini_num.csv'
#    elif 'b15' in sys.argv[ii]:
#        path = 'kicks/B15_ini_num.csv'
#    elif 'n20_lmc' in sys.argv[ii]:
#        path = 'kicks/N20_LMC_ini_num.csv'
#    elif 'n20' in sys.argv[ii]:
#        path = 'kicks/N20_ini_num.csv'
#    elif 'l15' in sys.argv[ii]:
#        path = 'kicks/L15_ini_num.csv'
#    elif 'w15_1z' in sys.argv[ii]:
#        path = 'kicks/W15_1Z_ini_num.csv'
#    elif 'w15' in sys.argv[ii]:
#        path = 'kicks/W15_ini_num.csv'
#    elif 'iib' in sys.argv[ii]:
#        path = 'kicks/IIb_ini_num.csv'
#    elif 'bhh' in sys.argv[ii]:
#        path = 'kicks/BHH_ini_num.csv'
#    elif 'b1d' in sys.argv[ii]:
#        path = 'kicks/B1D_ini_num.csv'
#    elif 'beo' in sys.argv[ii]:
#        path = 'kicks/BEO_ini_num.csv'
#    elif 'ben' in sys.argv[ii]:
#        path = 'kicks/BEN_ini_num.csv'
#    elif 'hmm' in sys.argv[ii]:
#        path = 'kicks/HMM_ini_num.csv'
#    elif 'm15' in sys.argv[ii]:
#        path = 'kicks/M15_ini_num.csv'
#
#    ini_nums = {}
#    with open(path) as ff:
#        for line in ff:
#            line = line.split()
#            ini_nums[line[0]] = float(line[1])
#
#    return ini_nums['feconix']

def ovrplt_spec(tt):
    if tt=='none':
        return
    else:
        tt = int(tt)

    n_avg = 3
    gin_low = np.loadtxt('/Users/silver/box/sci/lib/i/inoue91/ginga_6_16.txt')
    gin_low[:,1:] = gin_low[:,1:]/(11*kev2erg)
    gin_hig = np.loadtxt('/Users/silver/box/sci/lib/i/inoue91/ginga_16_28.txt')
    gin_hig[:,1:] = gin_hig[:,1:]/(24*kev2erg)
    gin_low_idx = np.argsort(np.abs(tt-gin_low[:,0]))[:n_avg]
    gin_hig_idx = np.argsort(np.abs(tt-gin_hig[:,0]))[:n_avg]
    gin_dat = np.zeros((2,2))
    gin_dat[0,0] = np.mean(gin_low[gin_low_idx, 1])
    gin_dat[0,1] = np.sqrt(np.sum(gin_low[gin_low_idx,2]**2))/3
    gin_dat[1,0] = np.mean(gin_hig[gin_hig_idx, 1])
    gin_dat[1,1] = np.sqrt(np.sum(gin_hig[gin_hig_idx,2]**2))/3

    if tt==185:
        dat = np.loadtxt('/Users/silver/box/sci/lib/s/syunyaev90/hexe_spec_fig_1a_185d.txt')
    elif tt==250:
        dat = np.loadtxt('/Users/silver/box/sci/lib/s/syunyaev90/hexe_spec_fig_1a_250d.txt')
    elif tt==320 and not (sys.argv[2][:11] == 'z_spec_300d' or sys.argv[2] == 'b15_m157b-2_spec_300d'):
        dat = np.loadtxt('/Users/silver/box/sci/lib/s/syunyaev90/hexe_spec_fig_1a_320d.txt')
        px1 = np.loadtxt('/Users/silver/box/sci/lib/s/syunyaev90/px1_spec_fig_1a_320d.txt')
        msfc = np.loadtxt('/Users/silver/box/sci/lib/p/pendleton95/msfc_spec_fig_7_248d.txt')
        grip = np.loadtxt('/Users/silver/box/sci/lib/p/palmer93/caltech_spec_tab_3_d268.txt')
        grip[:,1] = 1e-5*grip[:,1]
        grip[:,4:6] = 1e-5*grip[:,4:6]
        grip[:,4] = np.where(-grip[:,4] > grip[:,1], -grip[:,1]+1e-10, grip[:,4])
        plt.errorbar(msfc[:,0], msfc[:,1], xerr=np.abs(msfc[:,2:4]).T, yerr=np.abs(msfc[:,4:6]).T, fmt='none', ecolor='#8c564b', zorder=100)
        plt.errorbar(px1[:,0], px1[:,1], xerr=np.abs(px1[:,2:4]).T, yerr=np.abs(px1[:,4:6]).T, fmt='none', ecolor='#e377c2', zorder=100)
        plt.errorbar(grip[:,0], grip[:,1], xerr=np.abs(grip[:,2:4]).T, yerr=np.abs(grip[:,4:6]).T, fmt='none', ecolor='#bcbd22', zorder=100)
    elif tt==320:
        dat = np.loadtxt('/Users/silver/box/sci/lib/s/syunyaev90/hexe_spec_fig_1a_320d.txt')
    elif tt==425:
        dat = np.loadtxt('/Users/silver/box/sci/lib/s/syunyaev90/hexe_spec_fig_1a_425d.txt')
    else:
        print('ERROR. No spectrum available at requested time.')
        sys.exit(1)
    plt.errorbar(dat[:,0], dat[:,1], xerr=np.abs(dat[:,2:4]).T, yerr=np.abs(dat[:,4:6]).T, fmt='none', ecolor=ovr_col, zorder=100)
#    lim = gin_dat[0,0]+sts.norm.ppf(1-sts.norm.sf(3)*2)*gin_dat[0,1]
    yerr = np.array([np.min((gin_dat[0,1], gin_dat[0,0]-1e-16)), gin_dat[0,1]])[:, np.newaxis]
    plt.errorbar(11, gin_dat[0,0], xerr=5, yerr=yerr, ecolor='#7f7f7f', ms=2, zorder=99)
    plt.errorbar(22, gin_dat[1,0], xerr=6, yerr=gin_dat[1,1], ecolor='#7f7f7f', ms=2, zorder=99)
    
def ovrplt_lp(ff, ee):
    dat = np.loadtxt('/Users/silver/box/sci/lib/t/tueller91b/line_profiles_tab_1.txt')
    
    if '847' in ff:
        eline = 846.771
        idx = (dat[:,4] > 827) & (dat[:,4] < 867)
    elif '1238' in ff:
        eline = 1238.282
        idx = (dat[:,4] > 1218) & (dat[:,4] < 1258)

    for ii in range(0, np.sum(idx)):
        # Plot Gaussian profiles
        mag = dat[idx][ii,1]
        mu = delmc(dat[idx][ii,4])
        sig = fwhm2sig(dat[idx][ii,7])
        profile = gauss(ee, mag, mu, sig)
        plt.plot(e2v(ee, eline), profile)

        # Show the uncertainty
        emax = ee[np.argmax(profile)]
        pmax = np.amax(profile)
        xerr = e2v(np.array([-dat[idx][ii,5]+eline, dat[idx][ii,6]+eline]), eline)
        yerr = np.array([(-dat[idx][ii,2]/dat[idx][ii,1])*pmax, (dat[idx][ii,3]/dat[idx][ii,1])*pmax])
        plt.errorbar(e2v(emax, eline), pmax, xerr=xerr[:,np.newaxis], yerr=yerr[:,np.newaxis], fmt='none', ecolor=(0.7, 0.7, 0.7), zorder=-10)
        plt.text(e2v(mu, eline)+300, pmax, 'Day ' + str(int(dat[idx][ii,0])), color='k')

def ovrplt_lc(ee):
    if sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' or sys.argv[2] == '87a_lc_45-105_kev':
        scale_factor = 1.e5
    else:
        scale_factor = 1.

    if 'hexe' in ee:
        dat = np.loadtxt('/Users/silver/box/sci/lib/s/syunyaev90/hexe_lc_table_1.txt')
        dat[:, 2:] = dat[:, 2:]*1e-6*scale_factor
        tmid = (dat[:,0]+dat[:,1])/2
        xerr = np.c_[tmid-dat[:,0], dat[:,1]-tmid].T
        
    if ee=='hexe_low':
        plt.errorbar(tmid, dat[:, 2], xerr=xerr, yerr=dat[:,3], fmt='none', ecolor=ovr_col)
    elif ee=='hexe_mid':
        plt.errorbar(tmid, dat[:, 4], xerr=xerr, yerr=dat[:,5], fmt='none', ecolor=ovr_col)
    elif ee=='hexe_high':
        plt.errorbar(tmid, dat[:, 6], xerr=xerr, yerr=dat[:,7], fmt='none', ecolor=ovr_col)
    elif ee=='hexe_all':
        yerr = np.sqrt((dat[:,3]*30)**2+(60*dat[:,5])**2+(95*dat[:,7])**2)/185
        plt.errorbar(tmid, (dat[:, 2]*30+dat[:, 4]*60+dat[:, 6]*95)/185, xerr=xerr, yerr=yerr, fmt='none', ecolor=ovr_col)
    elif ee=='ginga_soft':
        dat = np.loadtxt('/Users/silver/box/sci/lib/i/inoue91/ginga_6_16.txt')
        plt.errorbar(dat[:,0], scale_factor*dat[:,1], yerr=scale_factor*dat[:,2], fmt='ko', ms=2, ecolor=ovr_col, zorder=100)
    elif ee=='ginga_hard':
        dat = np.loadtxt('/Users/silver/box/sci/lib/i/inoue91/ginga_16_28.txt')
        plt.errorbar(dat[:,0], scale_factor*dat[:,1], yerr=scale_factor*dat[:,2], fmt='ko', ms=2, ecolor=ovr_col, zorder=100)
        
def ovrplt_llc(mode):
    if sys.argv[2] == '87a_llc847+1238' or sys.argv[2] == 'b15_m157b-2_llc847+1238':
        scale_factor = 1.e3
    else:
        scale_factor = 1.
    
    # SMM
    dat = np.loadtxt('/Users/silver/box/sci/lib/l/leising90/smm_lines_tab_1.txt')
    dat[:, 2:] = dat[:, 2:]*1e-4
    tmid = (dat[:,0]+dat[:,1])/2
    ttt = np.linspace(np.amin(tmid), np.amax(2*tmid), 1000)
    xerr = np.c_[tmid-dat[:,0], dat[:,1]-tmid].T
    if '847+1238' in mode:
        plt.errorbar(tmid, scale_factor*(dat[:, 2]+dat[:, 4]), xerr=xerr, yerr=scale_factor*np.sqrt(dat[:,3]**2+dat[:,5]**2), fmt='none', ecolor=ovr_col, zorder=-9)
        lin_yie = 1.676
    elif '847' in mode:
        plt.errorbar(tmid, scale_factor*dat[:, 2], xerr=xerr, yerr=scale_factor*dat[:,3], fmt='none', ecolor=ovr_col, zorder=-9)
        lin_yie = 1.000
    elif '1238' in mode:
        plt.errorbar(tmid, scale_factor*dat[:, 4], xerr=xerr, yerr=scale_factor*dat[:,5], fmt='none', ecolor=ovr_col, zorder=-9)
        lin_yie = 0.676
    
    # Balloons
    markers = {'287': '*', #mahoney88
               '249': '^', #sandie88
               '320': 'x', #rester89
               '433': 'o', #tueller91b
               '613': 'o'} #tueller91b
    balloon = np.loadtxt('/Users/silver/box/sci/lib/t/tueller91b/line_fluxes_tab_1.txt')
    for ii, ball in enumerate(balloon):
        if ii < balloon.shape[0]-1 and balloon[ii,0] == balloon[ii+1,0]:
            continue
        elif ii > 0 and balloon[ii-1, 0] == balloon[ii,0]:
            bal2 = balloon[ii-1]
            yer1 = np.array([-ball[2], ball[3]])[:, np.newaxis].T
            yer2 = np.array([-bal2[2], bal2[3]])[:, np.newaxis].T
            fmt = markers[str(int(ball[0]))]
            plt.errorbar(ball[0], scale_factor*(ball[1]+bal2[1]), yerr=scale_factor*np.sqrt(yer1**2+yer2**2), fmt=fmt, ecolor='#7f7f7f', color='#7f7f7f', zorder=-10)
        else:
            if ball[4] < 1000:
                scale = 1.676
            else:
                scale = 1.676/0.676
            fmt = markers[str(int(ball[0]))]
            plt.errorbar(ball[0], scale_factor*scale*ball[1], yerr=scale_factor*scale*np.array([-ball[2], ball[3]])[:, np.newaxis].T, fmt=fmt, ecolor='#7f7f7f', color='#7f7f7f', zorder=-10)

    # UCR
    ucr = np.loadtxt('/Users/silver/box/sci/lib/a/ait_ouamer92/1238_kev.txt')
    scale = 1.676/0.676
    plt.errorbar(ucr[0], scale_factor*scale*ucr[1], yerr=scale_factor*scale*np.array([-ucr[2], ucr[3]])[:, np.newaxis].T, fmt='d', ecolor='#7f7f7f', color='#7f7f7f', zorder=-10)

    # GRIP
    grip1 = np.loadtxt('/Users/silver/box/sci/lib/p/palmer93/847_kev.txt')
    grip2 = np.loadtxt('/Users/silver/box/sci/lib/p/palmer93/1238_kev.txt')
    for ii in range(grip1.shape[0]):
        g1 = grip1[ii]
        g2 = grip2[ii]
        yer1 = np.array([-g1[2], g1[3]])[:, np.newaxis].T
        yer2 = np.array([-g2[2], g2[3]])[:, np.newaxis].T
        plt.errorbar(g1[0], scale_factor*(g1[1]+g2[1]), yerr=scale_factor*np.sqrt(yer1**2+yer2**2), fmt='s', ecolor='#7f7f7f', color='#7f7f7f', zorder=-10)
        
    # Stuff
    flx = lin_yie*bateman_equation(get_ini_num(), tau_ni56, tau_co56, d2s(ttt))
    flx = flx/(4*np.pi*(DD*kpc)**2)
    lims = plt.gca().get_ylim()
    plt.plot(ttt, scale_factor*flx, color='#7f7f7f', zorder=-10)
    plt.ylim((0.7*lims[0], 0.7*lims[1]))
    
def plt_lc(ff, norm):
    dat = np.loadtxt('dat/' + ff + '.txt')
    ne = (dat.shape[1]-1)
    dat[:,0] = s2d(dat[:,0])
    tt = (dat[1:,0]+dat[:-1,0])/2
        
    if nn == 1:
        for jj in range(1,ne+1):
            plt.plot(tt, dat[:-1,jj], alpha=0.3)
        smo = gaussian_filter(np.sum(dat[:-1,1:ne+1], axis=1), 0.001)
        smo = np.sum(dat[:-1,1:ne+1], axis=1)
        plt.plot(tt, smo, 'k')
        plt.legend(['$^{44}$Ti', '$^{56}$Ni', '$^{56}$Co', '$^{57}$Co', 'Sum'])

    else:
        smo = gaussian_filter(np.sum(dat[:-1,1:ne+1], axis=1), 0.0001)
        smo = np.sum(dat[:-1,1:ne+1], axis=1)
        # Normalize to maximum flux
        if norm:
            smo = smo/np.amax(smo)
#        plt.plot(tt, smo, color='gray')
        if sys.argv[2] == 'b15_lmc_lc_45-105_kev' and ii == 4:
            plt.plot(tt, lc_smooth(smo), color='#1f77b4', ls='--', alpha=0.5)
        elif sys.argv[2] == 'b15_lmc_lc_45-105_kev' and ii == 5:
            plt.plot(tt, lc_smooth(smo), color='#1f77b4', ls='-.', alpha=0.5)
        elif sys.argv[2] == 'b15_lmc_lc_45-105_kev' and ii == 3:
            plt.plot(tt, lc_smooth(smo))
        elif sys.argv[2] == '87a_lc_45-105_kev':
            plt.plot(tt, 1e5*lc_smooth(smo))
        elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' and ii == 3:
            plt.plot(tt, 1e5*lc_smooth(smo), color='#1f77b4')
        elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' and ii == 4:
            plt.plot(tt, 1e5*lc_smooth(smo), color='#ff7f0e')
        elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' and ii == 5:
            plt.plot(tt, 1e5*lc_smooth(smo), color='#1f77b4', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' and ii == 6:
            plt.plot(tt, 1e5*lc_smooth(smo), color='#1f77b4', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' and ii == 7:
#            tmp[1000:] = gaussian_filter(tmp[1000:],30)
            plt.plot(tt, 1e5*lc_smooth(smo), color='#ff7f0e', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' and ii == 8:
            plt.plot(tt, 1e5*lc_smooth(smo), color='#ff7f0e', ls='-.', alpha=0.5, zorder=-3)
        elif 'dod_dir' in sys.argv[2]:
            smo = 1.e5*smo
            tmp = 2.
            if ii == 3:
                plt.plot(tt, lc_smooth(smo,tmp), lw=1.5, color='k', zorder=100)
            elif ii == 4:
                plt.plot(tt, lc_smooth(smo,tmp), lw=1.5, color='k', zorder=100, ls='--')
            elif ii == 5:
                plt.plot(tt, lc_smooth(smo,tmp), lw=1.5, color='k', zorder=100, ls='-.')
            elif ii >= 26:
                plt.plot(tt, lc_smooth(smo,tmp), lw=1.5, color='#ff7f0e', zorder=102)
            else:
                if sys.argv[2] == 'm157b-2_lc_45-105_kev_dod_dir':
                    col_idx = (np.trapz(smo, tt)-0.2*743.42492230869625)/(1.*2159.757475649814-0.2*743.42492230869625)
                    col_idx = cm.get_cmap('viridis', 256)(col_idx)
                    plt.plot(tt, lc_smooth(smo,tmp), lw=1., color=col_idx)
                elif sys.argv[2] == 'b15_lmc_lc_45-105_kev_dod_dir':
                    col_idx = (np.trapz(smo, tt)-0.7*2087.07018421475)/(1.*3012.3011256424803-0.7*2087.07018421475)
                    col_idx = cm.get_cmap('viridis', 256)(col_idx)
                    plt.plot(tt, lc_smooth(smo,tmp), lw=1., color=col_idx)
                else:
                    plt.plot(tt, lc_smooth(smo,tmp), lw=2, color=viridis[ii])
        elif 'lmc_87a_' == sys.argv[2][:8]:
            smo = 1.e5*smo
            tmp = 2
            if ii == 3:
                plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4')
            elif ii == 5:
                plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', ls='--', alpha=0.5)
            elif ii == 6:
                plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', ls='-.', alpha=0.5)
            elif ii == 4:
                plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e')
            elif ii == 7:
                plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', ls='--', alpha=0.5)
            elif ii == 8:
                plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', ls='-.', alpha=0.5)
                
        elif 'm157b_all_lc_45-105_kev' == sys.argv[2]:
            smo = 1.e5*smo
            tmp = 2
            smo[:4] = 0.
            if ii in [3]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', zorder=1200)
            if ii in [4]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', zorder=1200)
            if ii in [5]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#2ca02c', zorder=1200)
            if ii in [6]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#d62728', zorder=1200)
            if ii in [7]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', zorder=100, alpha=0.5, ls='--')
            if ii in [9]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', zorder=100, alpha=0.5, ls='--')
            if ii in [11]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#2ca02c', zorder=100, alpha=0.5, ls='--')
            if ii in [13]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#d62728', zorder=100, alpha=0.5, ls='--')
            if ii in [8]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', zorder=100, alpha=0.5, ls='-.')
            if ii in [10]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', zorder=100, alpha=0.5, ls='-.')
            if ii in [12]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#2ca02c', zorder=100, alpha=0.5, ls='-.')
            if ii in [14]:
                plt.plot(tt, lc_smooth(smo,tmp), color='#d62728', zorder=100, alpha=0.5, ls='-.')
                
        else:
            plt.plot(tt, smo)
    
    plt.xlabel('Time (d)')
    if norm:
        plt.ylabel('Normalized flux density')
    elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev' or sys.argv[2] == '87a_lc_45-105_kev' or 'dod_dir' in sys.argv[2] or 'm157b_all_lc_45-105_kev' == sys.argv[2] or sys.argv[2] == 'lmc_87a_lc_45-105_kev':
        plt.ylabel('Flux density at 51.2 kpc ($10^{-5}$~$\gamma$~s$^{-1}$~keV$^{-1}$~cm$^{-2}$)')        
    else:
        plt.ylabel('Flux density at 51.2 kpc ($\gamma$~s$^{-1}$~keV$^{-1}$~cm$^{-2}$)')
    if 'iib' in ff.lower() or 'ib' in ff.lower() or 'ic' in ff.lower():
        plt.xlim([0, 300])
    else:
        plt.xlim([0, 1000])
    
def plt_llc(ff, mode):
    dat = np.loadtxt('dat/' + ff + '.txt')
    ne = (dat.shape[1]-1)
    dat[:,0] = s2d(dat[:,0])
    tt = (dat[1:,0]+dat[:-1,0])/2

    if '847' in mode:
        eline = 846.771
    elif '1238' in mode:
        eline = 1238.282
        
    smo = gaussian_filter(np.sum(dat[:-1,1:ne+1], axis=1), 0.001)
    if sys.argv[2][:4] == 'oth_':
#        plt.semilogy(tt, smo, color='gray')
        if sys.argv[2] == 'oth_llc847' and ii == 7:
            plt.semilogy(tt, lc_smooth(smo, 3), color='#1f77b4', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'oth_llc847' and ii == 8:
            plt.semilogy(tt, lc_smooth(smo, 3), color='#1f77b4', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'oth_llc847' and ii == 9:
            plt.semilogy(tt, lc_smooth(smo, 3), color='#ff7f0e', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'oth_llc847' and ii == 10:
            plt.semilogy(tt, lc_smooth(smo, 3), color='#ff7f0e', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'oth_llc847' and ii == 11:
            plt.semilogy(tt, lc_smooth(smo), color='#2ca02c', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'oth_llc847' and ii == 12:
            plt.semilogy(tt, lc_smooth(smo), color='#2ca02c', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'oth_llc847' and ii == 6:
            plt.semilogy(tt, lc_smooth(smo, 3), color='#7f7f7f', ls='-', zorder=-10)
        else:
            plt.semilogy(tt, lc_smooth(smo, 2))
    elif 'dod_dir' in sys.argv[2]:
        if ii == 3:
            plt.plot(tt, lc_smooth(smo,1), lw=2, color='k', zorder=100)
        elif ii == 4:
            plt.plot(tt, lc_smooth(smo,1), lw=2, color='k', zorder=100)
        elif ii == 5:
            plt.plot(tt, lc_smooth(smo,1), lw=2, color='k', zorder=100)
        else:
            plt.plot(tt, lc_smooth(smo,1), lw=2, color=viridis[ii])
    elif 'm157b_all_llc847+1238' == sys.argv[2]:
        smo = 1.e3*smo
        tmp = 2.
        if ii in [3]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', zorder=1200)
        if ii in [4]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', zorder=1200)
        if ii in [5]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#2ca02c', zorder=1200)
        if ii in [6]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#d62728', zorder=1200)
        if ii in [7]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', zorder=100, alpha=0.5, ls='--')
        if ii in [9]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', zorder=100, alpha=0.5, ls='--')
        if ii in [11]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#2ca02c', zorder=100, alpha=0.5, ls='--')
        if ii in [13]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#d62728', zorder=100, alpha=0.5, ls='--')
        if ii in [8]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#1f77b4', zorder=100, alpha=0.5, ls='-.')
        if ii in [10]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#ff7f0e', zorder=100, alpha=0.5, ls='-.')
        if ii in [12]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#2ca02c', zorder=100, alpha=0.5, ls='-.')
        if ii in [14]:
            plt.plot(tt, lc_smooth(smo,tmp), color='#d62728', zorder=100, alpha=0.5, ls='-.')
    elif 'lmc_87a_' == sys.argv[2][:8]:
        if ii == 3:
            plt.plot(tt, smo, color='#1f77b4')
        elif ii == 5:
            plt.plot(tt, smo, color='#1f77b4', ls='--', alpha=0.5)
        elif ii == 6:
            plt.plot(tt, smo, color='#1f77b4', ls='-.', alpha=0.5)
        elif ii == 4:
            plt.plot(tt, smo, color='#ff7f0e')
        elif ii == 7:
            plt.plot(tt, smo, color='#ff7f0e', ls='--', alpha=0.5)
        elif ii == 8:
            plt.plot(tt, smo, color='#ff7f0e', ls='-.', alpha=0.5)
    else:
#        plt.plot(tt, smo, color='gray')
        if sys.argv[2] == 'b15_m157b-2_llc847+1238' and ii == 5:
            plt.semilogy(tt, 1.e3*lc_smooth(smo, 2), color='#1f77b4', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_llc847+1238' and ii == 6:
            plt.semilogy(tt, 1.e3*lc_smooth(smo, 2), color='#1f77b4', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_llc847+1238' and ii == 7:
            plt.semilogy(tt, 1.e3*lc_smooth(smo, 2), color='#ff7f0e', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_llc847+1238' and ii == 8:
            plt.semilogy(tt, 1.e3*lc_smooth(smo, 2), color='#ff7f0e', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_llc847+1238':
            plt.semilogy(tt, 1.e3*lc_smooth(smo))
        elif sys.argv[2] == '87a_llc847+1238':
            plt.semilogy(tt, 1.e3*lc_smooth(smo))
        else:
            plt.plot(tt, smo)
        
    plt.xlabel('Time (d)')
    if sys.argv[2] == 'm157b_all_llc847+1238' or sys.argv[2] == '87a_llc847+1238' or sys.argv[2] == 'b15_m157b-2_llc847+1238':
        plt.ylabel('Flux at 51.2 kpc ($10^{-3}$~$\gamma$~s$^{-1}$~cm$^{-2}$)')
    else:
        plt.ylabel('Flux at 51.2 kpc ($\gamma$~s$^{-1}$~cm$^{-2}$)')
    plt.xlim([0, 1000])

    if ii == len(sys.argv)-2 and not (sys.argv[2][:4] == 'oth_' or 'm157b_all_llc847+1238' == sys.argv[2]):
        ovrplt_llc(mode)
    
def plt_spec(ff):
    dat = np.loadtxt('dat/' + ff + '.txt')
    ne = (dat.shape[1]-1)
    mid_ene = np.zeros(2*dat.shape[0])
    mid_flx = np.zeros((2*dat.shape[0], ne))
    mid_ene[::2] = dat[:,0]/1e3
    mid_ene[1::2] = dat[:,0]/1e3

    for jj in range(1,ne+1):
        mid_flx[0, jj-1] = 1e-99
        mid_flx[2:-1:2, jj-1] = dat[:-1, jj]
        mid_flx[1:-2:2, jj-1] = dat[:-1, jj]
        mid_flx[-1, jj-1] = 1e-99
    
    if nn == 1:
        for jj in range(0, ne):
            plt.loglog(mid_ene, spec_smooth(mid_ene, mid_flx[:, jj]), alpha=0.3)
        plt.loglog(mid_ene, spec_smooth(mid_ene, np.sum(mid_flx, axis=1)), 'k')
        plt.legend(['$^{44}$Ti', '$^{56}$Ni', '$^{56}$Co', '$^{57}$Co', 'Sum'])
        
#    elif nn == 1:
#        for jj in range(1,ne+1):
#            plt.loglog(mid_ene, mid_flx[:, jj-1], alpha=0.3)
#        plt.loglog(mid_ene, np.sum(mid_flx, axis=1), 'k')
#        plt.legend(['$^{44}$Ti', '$^{56}$Ni', '$^{56}$Co', '$^{57}$Co', 'Sum'])

    else:
        tmp = np.sum(mid_flx, axis=1)
        if sys.argv[2] == 'oth_spec' and ii == 6:
            plt.semilogy(mid_ene, spec_smooth(mid_ene, tmp), color='#7f7f7f', ls='-', zorder=-43)
        elif sys.argv[2] == 'b15_lmc_spec_300d' and ii == 4:
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), color='#1f77b4', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_lmc_spec_300d' and ii == 5:
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), color='#1f77b4', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_spec_300d' and ii == 5:
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), color='#1f77b4', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_spec_300d' and ii == 6:
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), color='#1f77b4', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_spec_300d' and ii == 7:
            tmp[1000:] = gaussian_filter(tmp[1000:],30)
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), color='#ff7f0e', ls='--', alpha=0.5, zorder=-3)
        elif sys.argv[2] == 'b15_m157b-2_spec_300d' and ii == 8:
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 40), color='#ff7f0e', ls='-.', alpha=0.5, zorder=-3)
        elif sys.argv[2] == '87a_spec_300d' and ii == 6:
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 70))
        elif 'dod_dir' in sys.argv[2]:
            if ii == 3:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp), lw=2, color='k', zorder=100, alpha=0.5)
            elif ii == 4:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), lw=2, color='k', zorder=100, alpha=0.5)
            elif ii == 5:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), lw=2, color='k', zorder=100, alpha=0.5)
            else:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), lw=2, color=viridis[ii])
        elif 'm157b_all_spec_300d' == sys.argv[2]:
            if ii in [3]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp), lw=1, color='#1f77b4', zorder=1200)
            if ii in [4]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp), lw=1, color='#ff7f0e', zorder=1200)
            if ii in [5]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp), lw=1, color='#2ca02c', zorder=1200)
            if ii in [6]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp), lw=1, color='#d62728', zorder=1200)
            if ii in [7,8]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), lw=1, color='#1f77b4', zorder=100, alpha=0.5)
            if ii in [9,10]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), lw=1, color='#ff7f0e', zorder=100, alpha=0.5)
            if ii in [11,12]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), lw=1, color='#2ca02c', zorder=100, alpha=0.5)
            if ii in [13,14]:
                plt.loglog(mid_ene, spec_smooth(mid_ene, tmp, 20), lw=1, color='#d62728', zorder=100, alpha=0.5)
        else:
            plt.loglog(mid_ene, spec_smooth(mid_ene, tmp))

    plt.xlabel('Energy (keV)')
    plt.ylabel('Flux density at 51.2 kpc ($\gamma$~s$^{-1}$~keV$^{-1}$~cm$^{-2}$)')
    plt.ylim([1e-7, 1e-3])

def plt_lp(ff):
    dat = np.loadtxt('dat/' + ff + '.txt')
    ne = (dat.shape[1]-1)
    mid_ene = np.zeros(2*dat.shape[0])
    mid_flx = np.zeros((2*dat.shape[0], ne))
    mid_ene[::2] = dat[:,0]/1e3
    mid_ene[1::2] = dat[:,0]/1e3

    for jj in range(1,ne+1):
        mid_flx[::2, jj-1] = dat[:, jj]
        mid_flx[1::2, jj-1] = dat[:, jj]

    if '847' in ff:
        eline = 846.771
    elif '1238' in ff:
        eline = 1238.282

    mid_v = e2v(mid_ene, eline)
    if 'dod_dir' in sys.argv[2]:
        if ii == 3:
            plt.plot(mid_v, np.sum(mid_flx, axis=1), lw=2, color='k', zorder=100, alpha=0.5)
        elif ii == 4:
            plt.plot(mid_v, np.sum(mid_flx, axis=1), lw=2, color='k', zorder=100, alpha=0.5)
        elif ii == 5:
            plt.plot(mid_v, np.sum(mid_flx, axis=1), lw=2, color='k', zorder=100, alpha=0.5)
        else:
            plt.plot(mid_v, np.sum(mid_flx, axis=1), lw=2, color=viridis[ii])
    else:
        plt.plot(mid_v, np.sum(mid_flx, axis=1))

    plt.xlabel('Velocity (km~s$^{-1}$)')
    plt.ylabel('Flux density at 51.2 kpc ($\gamma$~s$^{-1}$~keV$^{-1}$~cm$^{-2}$)')
    plt.xlim([np.amin(mid_v), np.amax(mid_v)])
#    plt.ylim([1e-7, 1e-3])

    if ii == len(sys.argv)-2 and not sys.argv[-1]=='disable':
        ovrplt_lp(ff, mid_ene)
        
################################################################
# Main
nn = len(sys.argv)-4
print('Number of data sets:', nn)
mode = sys.argv[1]
ovr_col = 'gray' if nn == 1 else 'k'
viridis = cm.get_cmap('viridis', 256)
viridis = viridis(np.linspace(0, 1, len(sys.argv)-1))

fig = plt.figure(figsize=(5, 3.75))
for ii in range(3, len(sys.argv)-1):
    print(sys.argv[ii])
    if mode == 'spec':
        plt_spec(sys.argv[ii])
    elif 'lp' in mode:
        plt_lp(sys.argv[ii])
    elif 'llc' in mode:
        plt_llc(sys.argv[ii], mode)
    elif mode == 'norm':
        plt_lc(sys.argv[ii], True)
    elif mode == 'lc':
        plt_lc(sys.argv[ii], False)

if nn > 1:
    if sys.argv[2] == '87a_llc847+1238' or sys.argv[2] == 'b15_m157b-2_llc847+1238':
        plt.legend(get_label(sys.argv[3:-1]), loc='upper right')
    elif sys.argv[2] == 'b15_lmc_spec_300d':
        pass
    elif sys.argv[2] == 'b15_m157b-2_spec_300d':
        plt.legend(get_label(sys.argv[3:5]), loc='best')
    elif sys.argv[2] == 'b15_m157b-2_lc_45-105_kev':
        plt.legend(get_label(sys.argv[3:5]), loc='best')
    elif sys.argv[2] == 'b15_lmc_lc_45-105_kev' or 'dod_dir' in sys.argv[2]:
        pass
    else:
        plt.legend(get_label(sys.argv[3:-1]), loc='best')

fix_lim()
        
if mode == 'spec' and not sys.argv[2] == 'sample_spec':
    ovrplt_spec(sys.argv[-1])
elif mode == 'lc':
    ovrplt_lc(sys.argv[-1])
#elif 'lp' in mode:
#    ovrplt_lp(int(sys.argv[-1]), mode)

fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/' + sys.argv[2] + '.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
#plt.show()
