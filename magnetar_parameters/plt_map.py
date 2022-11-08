'''
2018-10-03, Dennis Alp, dalp@kth.se

Takes data from get_dat and creates spectra.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
import re
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap
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
sig = np.deg2rad(5) # Sigma of the convolution kernel on the sphere
opening = np.deg2rad(30) # Opening angle of cone for directional quantities



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

def fix_ann(lab):
    if '_lmc' in lab:
        return lab.replace('_lmc', '')
    else:
        return lab

def smooth(sky):
    # Allocate arrays
    smo_map = np.zeros((nphi, ntht))
    bin_phi = np.linspace(-np.pi, np.pi, nphi+1)
    mid_phi = (bin_phi[:-1]+bin_phi[1:])/2
    bin_tht = np.linspace(-np.pi/2, np.pi/2, ntht+1)
    mid_tht = (bin_tht[:-1]+bin_tht[1:])/2
    area = (bin_phi[1]-bin_phi[0])*(bin_tht[1]-bin_tht[0])*np.cos(mid_tht)
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

#def get_min_max_dir():
#    sky = np.loadtxt('dat/' + sys.argv[1].replace('_smo','') + '.txt')
#    # Projection-convolution
#    smo_map = smooth(sky)
#
#    phi_min, tht_min = np.unravel_index(smo_map.argmin(), smo_map.shape)
#    phi_max, tht_max = np.unravel_index(smo_map.argmax(), smo_map.shape)
#
#    phi_min, tht_min = mid_phi[phi_min], mid_tht[tht_min]
#    phi_max, tht_max = mid_phi[phi_max], mid_tht[tht_max]
#
#    min_dir = np.array([phi_min, tht_min])
#    max_dir = np.array([phi_max, tht_max])
#    print('min_dir:', min_dir)
#    print('max_dir:', max_dir)
#    return min_dir, max_dir



################################################################
# Main
sky = np.loadtxt('../dat/' + sys.argv[1] + '.txt')
mod = sys.argv[2]
print('MODEL: ' + mod)    
nphi, ntht = sky.shape
bin_phi = np.linspace(-np.pi, np.pi, nphi+1)
mid_phi = (bin_phi[:-1]+bin_phi[1:])/2
bin_tht = np.linspace(-np.pi/2, np.pi/2, ntht+1)
mid_tht = (bin_tht[:-1]+bin_tht[1:])/2

fig = plt.figure()
plt.gca().axis('off')
mm = Basemap(projection='hammer', lon_0 = 0, llcrnrlon=True, llcrnrlat=True)
mm.drawmapboundary(fill_color='w', linewidth=0)
fig.set_size_inches(5, 3.75)

bin_phi, bin_tht = np.rad2deg(bin_phi), np.rad2deg(bin_tht)
bin_phi, bin_tht = np.meshgrid(bin_phi, bin_tht, indexing='ij')

plt_dat = np.zeros((sky.shape[0]+1, sky.shape[1]+1))
plt_dat[:-1,:-1] = sky
vmin = np.percentile(sky, 1)
vmax = np.percentile(sky,99)
mesh = mm.pcolormesh(bin_phi, bin_tht, plt_dat, latlon=True, edgecolors='face', lw=0, cmap='viridis', vmin=vmin, vmax=vmax)
cbar = mm.colorbar(mesh, location='bottom', pad="5%")
cbar.set_label('Angular Number Flux (s$^{-1}$~deg$^{-2}$)')


annotation = sys.argv[1].split('_')[0].upper() + ' ' + re.findall('([0-9]{0,9})d', sys.argv[1])[0] + ' d'
plt.text(0, 17e6, fix_ann(annotation), color='k', fontsize=10)

fig.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)
fig.savefig('/Users/silver/box/phd/pro/mag/nus/fig/' + sys.argv[1] + '.pdf', pad_inches=0.1, dpi=300)
# plt.show()
