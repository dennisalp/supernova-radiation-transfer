'''
2020-10-16, Dennis Alp, dalp@kth.se

Visualizes the output from grt.f, i.e. takes the event lists from the gamma-ray transfer and visualizes them.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap
import pandas as pd
import h5py

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)





################################################################
# Constants, cgs
cc = 2.99792458e10 # cm s-1
GG = 6.67259e-8 # cm3 g-1 s-2
hh = 6.6260755e-27 # erg s
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
mp = 1.67262192369e-24 # g





################################################################
# Parameters
nn = 4001
nlc = 101
ndir = 20
nphi = 72
ntht = 36
#nphi = 180
#ntht = 90
sig = np.deg2rad(5) # Sigma of the convolution kernel on the sphere
opening = np.deg2rad(30) # Half-opening angle of cone for directional quantities






################################################################
# Help functions
def d2s(d):
    return d*24.*60.*60.
def s2d(s):
    return s/(24.*60.*60.)

def sph2car(phi, tht):
    xx = np.cos(phi)*np.cos(tht)
    yy = np.sin(phi)*np.cos(tht)
    zz = np.sin(tht)
    return xx, yy, zz

def car2sph(xx, yy, zz):
    phi = np.arctan2(yy, xx)
    tht = np.arctan(zz/np.sqrt(xx**2+yy**2))
    return phi, tht

def get_sol_ang(tht):
    return 2*np.pi*(1-np.cos(tht))

def kernel(sky, alpha, area):
    return np.sum(area*sky*1/(2*np.pi*sig**2)*np.exp(-0.5*(alpha/sig)**2))

# Vertices of a dodecahedron
def get_dod_dir(jj):
    gr = (1+np.sqrt(5))/2
    if jj == 0:
        phi, tht = car2sph(1,1,1)
    elif jj == 1:
        phi, tht = car2sph(1,1,-1)
    elif jj == 2:
        phi, tht = car2sph(1,-1,1)
    elif jj == 3:
        phi, tht = car2sph(1,-1,-1)
    elif jj == 4:
        phi, tht = car2sph(-1,1,1)
    elif jj == 5:
        phi, tht = car2sph(-1,1,-1)
    elif jj == 6:
        phi, tht = car2sph(-1,-1,1)
    elif jj == 7:
        phi, tht = car2sph(-1,-1,-1)
    elif jj == 8:
        phi, tht = car2sph(0,gr,1/gr)
    elif jj == 9:
        phi, tht = car2sph(0,gr,-1/gr)
    elif jj == 10:
        phi, tht = car2sph(0,-gr,1/gr)
    elif jj == 11:
        phi, tht = car2sph(0,-gr,-1/gr)
    elif jj == 12:
        phi, tht = car2sph(1/gr,0,gr)
    elif jj == 13:
        phi, tht = car2sph(1/gr,0,-gr)
    elif jj == 14:
        phi, tht = car2sph(-1/gr,0,gr)
    elif jj == 15:
        phi, tht = car2sph(-1/gr,0,-gr)
    elif jj == 16:
        phi, tht = car2sph(gr,1/gr,0)
    elif jj == 17:
        phi, tht = car2sph(gr,-1/gr,0)
    elif jj == 18:
        phi, tht = car2sph(-gr,1/gr,0)
    elif jj == 19:
        phi, tht = car2sph(-gr,-1/gr,0)
    return np.array([phi, tht])
        
def smooth(sky):
    # Allocate arrays
    nphi, ntht = sky.shape
    smo_map = np.zeros((nphi, ntht))
    bin_phi = np.linspace(-np.pi, np.pi, nphi+1)
    mid_phi = (bin_phi[:-1]+bin_phi[1:])/2
    bin_tht = np.linspace(-np.pi/2, np.pi/2, ntht+1)
    mid_tht = (bin_tht[:-1]+bin_tht[1:])/2
    area = (bin_phi[1]-bin_phi[0])*(bin_tht[1]-bin_tht[0])*np.cos(mid_tht)
    # Mass weighted thetea coordinate (int xcos(x) dx / int cos(x) dx)
    mid_tht = np.diff(np.cos(bin_tht)+bin_tht*np.sin(bin_tht))/np.diff(np.sin(bin_tht))
    bin_phi, bin_tht = np.meshgrid(bin_phi, bin_tht, indexing='ij')
    mid_phi, mid_tht = np.meshgrid(mid_phi, mid_tht, indexing='ij')

    # Convolution
    xx, yy, zz = sph2car(mid_phi, mid_tht)
    xyz = np.dstack((xx, yy, zz))
    for ii in range(nphi):
        for jj in range(ntht):
            tmp = np.array([xx[ii,jj], yy[ii, jj], zz[ii,jj]])
            cos_angle = np.sum(tmp*xyz, axis=2)
            cos_angle = np.where(cos_angle > 1, 1., cos_angle)
            cos_angle = np.where(cos_angle < -1, -1., cos_angle)
            smo_map[ii, jj] = kernel(sky, np.arccos(cos_angle), area)
    return smo_map

def print_lc(out, meta, dat, wei, t_min, t_max, e_min, e_max, nlc, erg_or_cts, obs_dir):
    # Unit selection
    if erg_or_cts == 'erg':
        loc_wei = dat['ee']/1e3*kev2erg*wei
    elif erg_or_cts == 'cts':
        loc_wei = wei
    else:
        print('ERROR. erg_or_cts must be "erg" or "cts".')
        sys.exit(1)

    # For direction selection
    if obs_dir is None:
        sol_ang = 4*np.pi
        loc_ope = np.pi
        loc_dir = np.array([0., 0., 0.])
    else:
        sol_ang = get_sol_ang(opening)
        loc_ope = opening
        loc_dir = sph2car(obs_dir[0], obs_dir[1])
    cos_ang = np.sum(loc_dir*np.c_[dat['ox'],dat['oy'],dat['oz']], axis=1)
        
    idx = (dat['ee'] > e_min) & (dat['ee'] < e_max) & (cos_ang > np.cos(loc_ope))
    tt = np.linspace(t_min, t_max, nlc)
    flx, tmp = np.histogram(dat['tt'][idx], tt, weights=loc_wei[idx]*4*np.pi/sol_ang)
    flx = flx/np.diff(tt)
    
    flx = np.concatenate((flx, [flx[-1]+(flx[-1]-flx[-2])])) # Extrapolate one bin for cosmetics
    out[:, 0] = tt
    out[:, 1] = flx

def print_spec(out, meta, dat, wei, t_min, t_max, e_min, e_max, nn, obs_dir):
    # For direction selection
    if obs_dir is None:
        sol_ang = 4*np.pi
        loc_ope = np.pi
        loc_dir = np.array([0., 0., 0.])
    else:
        sol_ang = get_sol_ang(opening)
        loc_ope = opening
        loc_dir = sph2car(obs_dir[0], obs_dir[1])
    cos_ang = np.sum(loc_dir*np.c_[dat['ox'],dat['oy'],dat['oz']], axis=1)
        
    # Get the complete spectrum
    idx = (dat['tt'] > t_min) & (dat['tt'] < t_max) & (cos_ang > np.cos(loc_ope))
    ee = np.logspace(np.log10(e_min), np.log10(e_max), nn)
    flx, tmp = np.histogram(dat['ee'][idx], ee, weights=wei[idx]*4*np.pi/sol_ang)
    flx = flx/(np.diff(ee)/1e3*(t_max-t_min))

    flx = np.concatenate((flx, [flx[-1]+(flx[-1]-flx[-2])])) # Extrapolate one bin for cosmetics
    out[:, 0] = ee
    out[:, 1] = flx

def print_map(out, meta, dat, wei, t_min, t_max, e_min, e_max, nn, erg_or_cts):
    idx = (dat['ee'] > e_min) & (dat['ee'] < e_max) & (dat['tt'] > t_min) & (dat['tt'] < t_max)
    phi, tht = car2sph(dat['ox'][idx], dat['oy'][idx], dat['oz'][idx])
    solid_angle = 360/nphi*180/ntht # square degrees
    # Unit selection
    if erg_or_cts == 'erg':
        loc_wei = dat['ee']/1e3*kev2erg*wei
    elif erg_or_cts == 'cts':
        loc_wei = wei

    sky, _, _ = np.histogram2d(np.array(phi), np.array(tht), bins=(bin_phi, bin_tht), weights=loc_wei[idx]/np.cos(tht))
    sky = sky/(solid_angle*(t_max-t_min)) # unit conversion
    out[:,:] = out+sky

def get_dat(ff):
    npy = '../out/' + ff + '.npy'
    if os.path.isfile(npy):
        print('Reading binary')
        dat = pd.DataFrame(np.load(npy), columns=['ei', 'ee', 'ox', 'oy', 'oz', 'ns', 'tt', 'ti'])
        print('Loaded')
    else:
        print('Reading CSV')
        dat = pd.read_csv('../out/' + ff + '.esc', delim_whitespace=True, names=['ei', 'ee', 'ox', 'oy', 'oz', 'ns', 'tt', 'ti'])
        print('Loaded')
        np.save(npy, dat)
    return dat

def get_meta(ff):
    meta = {}
    with open('../out/' + ff + '.txt') as ff:
        for line in ff:
            (key, val) = line.split('\t')
            meta[key] = val.strip()
    return meta

def get_wei(dat, meta):
    def get_lum(tt):
        mm = 1.4*Msun
        rr = 1.2e6
        tmp = 'Assuming M={0:.1f} Msun, R={1:.0f} km, P={2:.1f} ms, B={3:.2f}*10^14 G'
        tmp = tmp.format(mm/Msun, rr/1e5, p0*1.e3, bb/1e14)
        print(tmp)

        tp = 3/40*mm*cc**3/(np.pi**2*rr**4)*p0**2/bb**2
        ll = 2**5/3*np.pi**4/cc**3*rr**6/p0**4*bb**2*1/(1+tt/tp)**2
        return eta*ll

    def get_num(ll):
        E = np.logspace(np.log10(e0), np.log10(e1), 100001)

        # Compute spectrum normalization
        if np.abs(gam-2) < 1.e-6:
            n0 = ll/(np.log(e1)-np.log(e0))
        else:
            g2 = -gam+2
            n0 = ll*g2/(e1**g2-e0**g2)
            # int n0*E**g1 dE = n0/xg2*(e_max**g2-e_min**g2) = ll

        check = np.trapz(n0[0]*E**(-gam+1), E)
        check = np.abs(ll[0]-check)/ll[0]
        if check > 1.e-4:
            print('Integration of spectrum normalization failed.')
            sys.exit(1)


        # Compute photon luminosity
        if np.abs(gam-1) < 1.e-6:
            nn = n0*(np.log(e1)-np.log(e0))
        else:
            g1 = -gam+1
            nn = n0/g1*(e1**g1-e0**g1)
        
        check = np.trapz(n0[0]*E**(-gam+1), E)/np.trapz(n0[0]*E**-gam, E)
        check = np.abs(ll[0]-check*nn[0])/ll[0]
        if check > 1.e-4:
            print('Integration of spectrum normalization failed.')
            sys.exit(1)
            
        return nn

    def get_time_wei(nn):
        pps = n_pho/(m_tmax-m_tmin)
        return nn/pps

    def get_spec_wei(wei):
        sim_gam = float(meta['sim_gam'])

        s1 = -sim_gam+1
        g1 = -gam+1
        if np.abs(s1) < 1.e-6:
            n_sim = np.log(e1)-np.log(e0)
        else:
            n_sim = (e1**s1-e0**s1)/s1
        if np.abs(g1) < 1.e-6:
            n_gam = np.log(e1)-np.log(e0)
        else:
            n_gam = (e1**g1-e0**g1)/g1

        ei = np.array(dat['ei'])*kev2erg/1.e3
        wei = wei*(ei**-gam/n_gam)/(ei**-sim_gam/n_sim)
        return wei
    

    e0 = float(meta['e_min'])*kev2erg/1.e3
    e1 = float(meta['e_max'])*kev2erg/1.e3
    n_pho = float(meta['n_pho'])
    m_tmin = d2s(float(meta['tmin']))
    m_tmax = d2s(float(meta['tmax']))

    ll = get_lum(np.array(dat['ti']))
    nn = get_num(ll)
    wei = get_time_wei(nn)
    wei = get_spec_wei(wei)
    return wei


def get_hea():
    if mode == 'lc':
        hea = 'See Alp et al. (2018; 2019). First column: Time at the early side of the bin in seconds. Luminosity in units of erg s^-1'
    elif mode == 'spec':
        hea = 'See Alp et al. (2018; 2019). First column: Energy at the lower end of the bin in eV. Number luminosity density in units of photons s^-1 keV^-1.'
    elif mode == 'map':
        hea = 'See Alp et al. (2018; 2019). Spherical projection with a resolution of 72 phi bins times 36 theta bins. This files contains a 2D array of the flux in units of photons s^-1 deg^-2. Theta is defined from 0 (north pole) to pi (south pole).'        
    return hea








################################################################
# Parse input
print('Expects arguments: lab e_min e_max t_min t_max mode data_set P/ms B/1e14G Gamma efficiency')
if not len(sys.argv) == 13:
    print('ERROR. Wrong number of input arguments')
    sys.exit()
    
lab = sys.argv[1] # label for the output
e_min = float(sys.argv[2])*1e3
e_max = float(sys.argv[3])*1e3
t_min = float(sys.argv[4])
t_max = float(sys.argv[5])
mode = sys.argv[6]

ff = sys.argv[7] # ID of dataset
mod = '_'.join(ff.split('_')[:-1])
print('MODEL: ' + mod)    

p0 = float(sys.argv[8])/1000
bb = float(sys.argv[9])*1e14
gam = float(sys.argv[10])
eta = float(sys.argv[11])
erg_or_cts = sys.argv[12]




################################################################
# Allocate
bin_phi = np.linspace(-np.pi, np.pi, nphi+1)
mid_phi = (bin_phi[:-1]+bin_phi[1:])/2
bin_tht = np.linspace(-np.pi/2, np.pi/2, ntht+1)
mid_tht = (bin_tht[:-1]+bin_tht[1:])/2

if mode == 'lc':
    out = np.zeros((nlc, 2))
    out_dod = list(np.zeros((ndir, nlc, 2)))
elif mode == 'spec':
    out = np.zeros((nn, 2))
    out_dod = list(np.zeros((ndir, nn, 2)))
elif mode == 'map':
    out = np.zeros((nphi, ntht))
else:
    print('ERROR. mode must be "lc", "spec", or "map".')
    sys.exit(1)


# Load packets and meta data
print(ff)
dat = get_dat(ff)
meta = get_meta(ff)
print(meta)

# Compute photons packet^-1 (weights)
wei = get_wei(dat, meta)





################################################################
# Do stuff
if mode == 'lc' or mode == 'spec':
    if mode == 'lc':
        print_lc(out, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nlc, erg_or_cts, None)
        for jj in range(ndir):
            dod_dir = get_dod_dir(jj)
            print_lc(out_dod[jj], meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nlc, erg_or_cts, dod_dir)

    elif mode == 'spec':
        print_spec(out, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, None)
        for jj in range(ndir):
            dod_dir = get_dod_dir(jj)
            print_spec(out_dod[jj], meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, dod_dir)
            
    for jj in range(ndir):
        np.savetxt('../dat/' + lab + '_dir_' + "{:03d}".format(jj) + '.txt', out_dod[jj], header=get_hea())

elif mode == 'map':
    print_map(out, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, erg_or_cts)
    smo = smooth(out)
    np.savetxt('../dat/' + lab + '_smo.txt', smo, header=get_hea())
    
np.savetxt('../dat/' + lab + '.txt', out, header=get_hea())
