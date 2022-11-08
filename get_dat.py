'''
2018-10-03, Dennis Alp, dalp@kth.se

Visualizes the output from dalp_continuous, i.e. takes the event lists from the gamma-ray transfer and visualizes them.
'''

from __future__ import division, print_function
import os
import pdb
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
# Parameters
Msun = 1.989e33 # g
uu = 1.660539040e-24 # g
kev2erg = 1.60218e-9 # erg keV-1
kpc = 3.086e21 # cm
mpc = 3.086e24 # cm
DD = 51.2*kpc # cm
tau_ni56 = 6.075/np.log(2) # milne04, junde99, junde11 nuclear data sheets
tau_co56 = 77.236/np.log(2) # milne04, junde99, junde11 nuclear data sheets
tau_ti44 = 58.9*365.25/np.log(2) # ahmad06
tau_ni57 = 35.60/24./np.log(2) # bhat98
tau_co57 = 271.74/np.log(2) # bhat98
ratio_57_56 = 2*0.02309436100878436 # twice (kurfess92, fransson93) solar (lodders03), table 6
nn = 4001
nlc = 101
nllc = 101
nlp = 101
ndir = 20
nphi = 72
ntht = 36
#nphi = 180
#ntht = 90
sig = np.deg2rad(5) # Sigma of the convolution kernel on the sphere
opening = np.deg2rad(30) # Opening angle of cone for directional quantities
ph_per_decay = np.array([4.80106, 3.213, 2.90798, 1.056681])



################################################################
# Parse input
print('Expects arguments: lab e_min e_max t_min t_max erg_or_cts mode 44ti/sc 56ni 56co 57co\nPass "none" to skip isotope.')
if not len(sys.argv) == 12:
    print('ERROR. Wrong number of input arguments')
    sys.exit()
    
lab = sys.argv[1] # label for the output
e_min = float(sys.argv[2])*1e3
e_max = float(sys.argv[3])*1e3
t_min = float(sys.argv[4])
t_max = float(sys.argv[5])
erg_or_cts = sys.argv[6]
mode = sys.argv[7]

fid = []
for ii in range(8, 12):
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

def get_min_max_dir():
    tmp_map = np.zeros((nphi, ntht))
    long_runs = ['020', '021', '015', '016']
    for ii, ff in enumerate(fid):
        ff = mod + '_' + long_runs[ii]
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
        wei = get_wei(meta, dat, ini_num, dt)
        print_map(tmp_map, meta, dat, wei, d2s(0), d2s(1000), 0, 4e6, nn, 'cts')

    # Projection-convolution
    smo_map = smooth(tmp_map)
    phi_min, tht_min = np.unravel_index(smo_map.argmin(), smo_map.shape)
    phi_max, tht_max = np.unravel_index(smo_map.argmax(), smo_map.shape)

    phi_min, tht_min = mid_phi[phi_min], mid_tht[tht_min]
    phi_max, tht_max = mid_phi[phi_max], mid_tht[tht_max]

    min_dir = np.array([phi_min, tht_min])
    max_dir = np.array([phi_max, tht_max])
    print('min_dir:', min_dir)
    print('max_dir:', max_dir)
    return min_dir, max_dir

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
    flx = flx/(np.diff(tt)*4*np.pi*DD**2*(e_max-e_min)/1e3) # flux at 51.2 kpc (s$^{-1}$~cm$^{-2}$)
    out[:, 0] = tt[:-1]
    out[:, int(meta['src_abu'])+1] = flx

#    idx = (dat['ee'] > e_min) & (dat['ee'] < e_max) & (dat['ns'] == 0) & (cos_ang > np.cos(loc_ope))
#    flx, tmp = np.histogram(dat['tt'][idx], tt, weights=loc_wei[idx]*4*np.pi/sol_ang)
#    flx = flx/(np.diff(tt)*4*np.pi*DD**2*(e_max-e_min)/1e3) # flux at 51.2 kpc (s$^{-1}$~cm$^{-2}$)
#    out[:, int(meta['src_abu'])+5] = flx

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
    flx = flx/(np.diff(ee)/1e3*(t_max-t_min)*4*np.pi*DD**2)
    out[:, 0] = ee[:-1]
    out[:, int(meta['src_abu'])+1] = flx

#    # Get spectrum of direct escapes
#    idx = (dat['tt'] > t_min) & (dat['tt'] < t_max) & (dat['ns'] == 0) & (cos_ang > np.cos(loc_ope))
#    flx, tmp = np.histogram(dat['ee'][idx], ee, weights=wei[idx]*4*np.pi/sol_ang)
#    flx = flx/(np.diff(ee)/1e3*(t_max-t_min)*4*np.pi*DD**2)
#    out[:, int(meta['src_abu'])+5] = flx

def print_llc(out, meta, dat, wei, t_min, t_max, mode, nllc, obs_dir):
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
    if '847' in mode:
        eline = 846.771*1e3
    elif '1238' in mode:
        eline = 1238.282*1e3
    else:
        print('ERROR. Invalid line energy.')
        sys.exit(1)

    ew = 1.02
#    e_min = eline/ew
#    e_max = eline*ew
    if '847+1238' in mode:
        eline1 = 846.771*1e3
        eline2 = 1238.282*1e3
        idx = ((np.abs(dat['ei']-eline1) < 1) | ((np.abs(dat['ei']-eline2) < 1)))
        idx = (cos_ang > np.cos(loc_ope)) & idx & (dat['ns'] == 0)
    else:
        idx = (cos_ang > np.cos(loc_ope)) & (np.abs(dat['ei']-eline) < 1) & (dat['ns'] == 0)
    tt = np.linspace(t_min, t_max, nllc)
    flx, tmp = np.histogram(dat['tt'][idx], tt, weights=wei[idx]*4*np.pi/sol_ang)
    flx = flx/(np.diff(tt)*4*np.pi*DD**2) # flux at 51.2 kpc (s$^{-1}$~cm$^{-2}$)
    out[:, 0] = tt[:-1]
    out[:, 1] = flx

def print_lp(out, meta, dat, wei, t_min, t_max, mode, nlp, obs_dir):
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
    if '847' in mode:
        eline = 846.771*1e3
    elif '1238' in mode:
        eline = 1238.282*1e3
    else:
        print('ERROR. Invalid line energy.')
        sys.exit(1)

    ew = 1.02
    e_min = eline/ew
    e_max = eline*ew
    dir_esc = (np.abs(dat['ei']-eline) < 1) & (dat['ns'] == 0)
    idx = (dat['tt'] > t_min) & (dat['tt'] < t_max) & (cos_ang > np.cos(loc_ope)) & dir_esc
    ee = np.linspace(e_min, e_max, nlp)
    flx, tmp = np.histogram(dat['ee'][idx], ee, weights=wei[idx]*4*np.pi/sol_ang)
    flx = flx/(np.diff(ee)/1e3*(t_max-t_min)*4*np.pi*DD**2)
    out[:, 0] = ee[:-1]
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
#    path = '/Users/silver/dat/sne/' + meta['model'] + '.h5'
#    dat = h5py.File(path, 'r')
#
#    rad = np.array(dat['radius'])
#    tht = np.array(dat['theta'])
#    phi = np.array(dat['phi'])
#    den = np.array(dat['den'])
#    den[:, :, -1] = den[:, :, -2]
#
#    ################################################################
#    # Check volume and density
#    ato_idx = {'ar36': 0, 'c12': 1, 'ca40': 2, 'ca44': 3, 'co56': 4, 'cr48': 5, 'fe52': 6, 'fe56': 7, 'he4' : 8, 'mg24': 9, 'n': 10, 'ne20': 11, 'ni56': 12, 'o16': 13, 'p': 14, 's32': 15, 'sc44': 16, 'si28': 17, 'ti44': 18, 'x56': 19}
#
#    ato_mas = uu*np.array([35.96754510600, 12.00000000000, 39.96259098000, 43.95548180000, 55.93983930000, 47.95403200000, 51.94811400000, 55.93493630000, 4.002603254150, 23.98504170000, 1.008664915880, 19.99244017540, 55.94213200000, 15.99491461956, 1.007276466879, 31.97207100000, 43.95940280000, 27.97692653250, 43.95969010000, 55.93493630000])
#
#    mid_rad = (rad[:-1]+rad[1:])/2.
#    mid_phi = (phi[:-1]+phi[1:])/2.
#    mid_tht = (tht[:-1]+tht[1:])/2.    
#    mid_phi = mid_phi[:, np.newaxis, np.newaxis]
#    mid_tht = mid_tht[np.newaxis, :, np.newaxis]
#    mid_rad = mid_rad[np.newaxis, np.newaxis, :]
#    dif_phi = np.diff(phi)[:, np.newaxis, np.newaxis]
#    dif_tht = np.diff(tht)[np.newaxis, :, np.newaxis]
#    dif_rad = np.diff(rad)[np.newaxis, np.newaxis, :]
#
#    vol = mid_rad**2*np.sin(mid_tht)*dif_rad*dif_tht*dif_phi
#    if meta['src_abu']=='0':
#        tmp = dat['ca44']/ato_mas[ato_idx['ca44']]+dat['sc44']/ato_mas[ato_idx['sc44']]+dat['ti44']/ato_mas[ato_idx['ti44']]
#        tmp[:,:,-1] = tmp[:,:,-2]
#        ini_num =  np.sum(vol*den*tmp)
#        print_num = ini_num*ato_mas[ato_idx['ti44']]/Msun
#    elif meta['src_abu']=='1' or meta['src_abu']=='2':
#        tmp = dat['fe56']/ato_mas[ato_idx['fe56']]+dat['co56']/ato_mas[ato_idx['co56']]+dat['ni56']/ato_mas[ato_idx['ni56']]+0.5*(dat['x56']/ato_mas[ato_idx['x56']])
#        tmp[:,:,-1] = tmp[:,:,-2]
#        ini_num =  np.sum(vol*den*tmp)
#        print_num = ini_num*ato_mas[ato_idx['x56']]/Msun
#    elif meta['src_abu']=='3':
#        tmp = dat['fe56']/ato_mas[ato_idx['fe56']]+dat['co56']/ato_mas[ato_idx['co56']]+dat['ni56']/ato_mas[ato_idx['ni56']]+0.5*(dat['x56']/ato_mas[ato_idx['x56']])
#        tmp[:,:,-1] = tmp[:,:,-2]
#        ini_num = ratio_57_56*np.sum(vol*den*tmp)
#        print_num = ini_num*ato_mas[ato_idx['x56']]/Msun
#
#    print('Mass of ' + ['44Ti', '56Ni', '56Co', '57Co'][int(meta['src_abu'])] + ':', print_num)
#    print('Should match:', ini_num)

# pre-computed #     # Newer faster pre-computed version
# pre-computed #     path = 'kicks/' + meta['model'] + '_ini_num.csv'
# pre-computed #     ini_nums = {}
# pre-computed #     with open(path) as ff:
# pre-computed #         for line in ff:
# pre-computed #             line = line.split()
# pre-computed #             ini_nums[line[0]] = float(line[1])
# pre-computed # 
# pre-computed #     if meta['src_abu']=='0':
# pre-computed #         ini_num = ini_nums['cascti']
# pre-computed #     elif meta['src_abu']=='1' or meta['src_abu']=='2':
# pre-computed #         ini_num = ini_nums['feconix']
# pre-computed #     elif meta['src_abu']=='3':
# pre-computed #         ini_num = ini_nums['feconix']*ratio_57_56
# pre-computed # #    print('Should match:', ini_num)
# pre-computed #     return ini_num
    if meta['src_abu']=='0':
        return 1.5e-4*(51.2/50)**2*Msun/(43.9596901*uu)
    elif meta['src_abu']=='1' or meta['src_abu']=='2':
        return 0.07*Msun/(55.942132*uu)
    elif meta['src_abu']=='3':
        return 0.07*Msun/(55.942132*uu)*ratio_57_56

def get_wei(meta, dat, ini_num, dt):
    if meta['src_abu'] == '0':
        wei = dt*ini_num*np.exp(-np.array(dat['ti'])/d2s(tau_ti44))/d2s(tau_ti44)
    elif meta['src_abu'] == '1':
        wei = dt*ini_num*np.exp(-np.array(dat['ti'])/d2s(tau_ni56))/d2s(tau_ni56)
    elif meta['src_abu'] == '2':
        wei = dt*bateman_equation(ini_num, tau_ni56, tau_co56, np.array(dat['ti']))
    elif meta['src_abu'] == '3':
        wei = dt*bateman_equation(ini_num, tau_ni57, tau_co57, np.array(dat['ti']))
    return wei



################################################################
# Prepare
bin_phi = np.linspace(-np.pi, np.pi, nphi+1)
mid_phi = (bin_phi[:-1]+bin_phi[1:])/2
bin_tht = np.linspace(-np.pi/2, np.pi/2, ntht+1)
mid_tht = (bin_tht[:-1]+bin_tht[1:])/2

if mode == 'lc':
    out = np.zeros((nlc-1, 5))
    out_min = np.zeros((nlc-1, 5))
    out_max = np.zeros((nlc-1, 5))
    out_dod = list(np.zeros((ndir, nlc-1, 5)))
elif 'llc' in mode:
    out = np.zeros((nlc-1, 2))
    out_min = np.zeros((nlc-1, 2))
    out_max = np.zeros((nlc-1, 2))
    out_dod = list(np.zeros((ndir, nlc-1, 5)))
elif mode == 'spec':
    out = np.zeros((nn-1, 5))
    out_min = np.zeros((nn-1, 5))
    out_max = np.zeros((nn-1, 5))
    out_dod = list(np.zeros((ndir, nn-1, 5)))
elif 'lp' in mode:
    out = np.zeros((nlp-1, 2))
    out_min = np.zeros((nlp-1, 2))
    out_max = np.zeros((nlp-1, 2))
    out_dod = list(np.zeros((ndir, nlp-1, 5)))
elif mode == 'map':
    out = np.zeros((nphi, ntht))
else:
    print('ERROR. mode must be "lc", "spec", or "map".')
    sys.exit(1)

# Load everything once to get min and max dir
if os.path.isfile('min_max_dir/'+mod+'.txt'):
    tmp = np.loadtxt('min_max_dir/'+mod+'.txt')
    min_dir = tmp[0]
    max_dir = tmp[1]
else:
    min_dir, max_dir = get_min_max_dir()
    np.savetxt('min_max_dir/'+mod+'.txt', np.vstack((min_dir, max_dir)))

# Load packets
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
    wei = get_wei(meta, dat, ini_num, dt)

    if mode == 'lc' or mode == 'spec' or 'lp' in mode or 'llc' in mode:
        if mode == 'lc':
            header = 'From Alp et al. (2019). First column: Time at the early side of the bin in seconds. Then follows four fluxes in units of photons s^-1 keV^-1 cm^-2 (multiply by the energy width in keV for total flux) scaled to a distance of 51.2 kpc (SN 1987A): 44Ti+44Sc, 56Ni, 56Co, 57Co'
            print_lc(out, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nlc, erg_or_cts, None)
            print_lc(out_min, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nlc, erg_or_cts, min_dir)
            print_lc(out_max, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nlc, erg_or_cts, max_dir)
            for jj in range(ndir):
                dod_dir = get_dod_dir(jj)
                print_lc(out_dod[jj], meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nlc, erg_or_cts, dod_dir)

        elif mode == 'spec':
            header = 'From Alp et al. (2019). First column: Energy at the lower end of the bin in eV. Then follows four fluxes in units of photons s^-1 keV^-1 cm^-2 scaled to a distance of 51.2 kpc (SN 1987A): 44Ti+44Sc, 56Ni, 56Co, 57Co'
            print_spec(out, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, None)
            print_spec(out_min, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, min_dir)
            print_spec(out_max, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, max_dir)
            for jj in range(ndir):
                dod_dir = get_dod_dir(jj)
                print_spec(out_dod[jj], meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, dod_dir)

        elif 'lp' in mode:
            header = 'From Alp et al. (2019). First column: Energy at the lower end of the bin in eV. Second column is the line flux in units of photons s^-1 keV^-1 cm^-2 scaled to a distance of 51.2 kpc (SN 1987A).'
            print_lp(out, meta, dat, wei, d2s(t_min), d2s(t_max), mode, nlp, None)
            print_lp(out_min, meta, dat, wei, d2s(t_min), d2s(t_max), mode, nlp, min_dir)
            print_lp(out_max, meta, dat, wei, d2s(t_min), d2s(t_max), mode, nlp, max_dir)
            for jj in range(ndir):
                dod_dir = get_dod_dir(jj)
                print_lp(out_dod[jj], meta, dat, wei, d2s(t_min), d2s(t_max), mode, nlp, dod_dir)

        elif 'llc' in mode:
            header = 'From Alp et al. (2019). First column: Time at the early side of the bin in seconds. Second column is the line flux in units of photons s^-1 cm^-2 scaled to a distance of 51.2 kpc (SN 1987A).'
            print_llc(out, meta, dat, wei, d2s(t_min), d2s(t_max), mode, nllc, None)
            print_llc(out_min, meta, dat, wei, d2s(t_min), d2s(t_max), mode, nllc, min_dir)
            print_llc(out_max, meta, dat, wei, d2s(t_min), d2s(t_max), mode, nllc, max_dir)
            for jj in range(ndir):
                dod_dir = get_dod_dir(jj)
                print_llc(out_dod[jj], meta, dat, wei, d2s(t_min), d2s(t_max), mode, nllc, dod_dir)

                
        np.savetxt('dat/' + lab + '_min_dir.txt', out_min, header = header)
        np.savetxt('dat/' + lab + '_max_dir.txt', out_max, header = header)
        for jj in range(ndir):
            np.savetxt('dat/' + lab + '_dir_' + "{:03d}".format(jj) + '.txt', out_dod[jj], header = header)
    elif mode == 'map':
        print_map(out, meta, dat, wei, d2s(t_min), d2s(t_max), e_min, e_max, nn, erg_or_cts)
        smo = smooth(out)
        header = 'From Alp et al. (2019). Spherical projection with a resolution of 72 phi bins times 36 theta bins. This files contains a 2D array of the flux in units of photons s^-1 deg^-2. Theta is defined from 0 (north pole) to pi (south pole).'
        np.savetxt('dat/' + lab + '_smo.txt', smo, header = header)
    
np.savetxt('dat/' + lab + '.txt', out, header = header)
