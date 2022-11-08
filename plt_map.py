'''
2018-10-03, Dennis Alp, dalp@kth.se

Takes data from get_dat and creates spectra.
'''

from __future__ import division, print_function
import os
import pdb
import sys
import re
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap
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
    elif 'M157B-2' in lab:
        return lab.replace('M157B-2', 'M15-7b')
    elif lab[:7] == 'feconix':
        return lab.replace('feconix', '$^{56}$Ni')
    elif lab[:3] == 'min':
        return lab.replace('min', 'Min')
    elif lab[:3] == 'max':
        return lab.replace('max', 'Max')
    elif lab[:2] == 'ns':
        return lab.replace('ns', 'NS')
    else:
        return lab

def shift_lon_lat(mm, lon, lat, lab):
    b15 = sys.argv[1][:12] == 'b15_lmc_map_'
    n20 = sys.argv[1] == 'n20_lmc_map_300d_smo'
    l15 = sys.argv[1] == 'l15_map_300d_smo'
    w15 = sys.argv[1] == 'w15_map_300d_smo'
    m157b = sys.argv[1] == 'm157b-2_map_300d_smo'
    iib = sys.argv[1] == 'iib_map_100d_smo'
    print(sys.argv[1], m157b, lab)
    
    if lon < 0 and lat < 0: # Bottom left
        lon, lat = mm(lon, lat)
        if w15 and lab[:3] == 'min':
            return lon+4.2e5, lat-1.6e6
        elif l15 and lab[:2] == 'ns':
            return lon-3e6, lat+3e5
        elif b15 and lab[:3] == 'min':
            return lon+3.7e5, lat+6e5
        elif iib and lab[:3] == 'min':
            return lon+4.5e5, lat-1.4e6
        return lon+3.7e5, lat+6e5
    elif lon > 0 and lat < 0: # Bottom right
        lon, lat = mm(lon, lat)
        if l15 and lab[:3] == 'min':
            return lon-2.7e6, lat+3e5
        return lon-3e6, lat+3e5
    elif lon < 0 and lat > 0: # Top left
        lon, lat = mm(lon, lat)
        if l15 and lab[:3] == 'max':
            return lon+4e5, lat+4e5
        elif l15 and lab[:7] == 'feconix':
            return lon+2e5, lat+2e5
        elif m157b and lab[:7] == 'feconix':
            return lon+5.e5, lat-1.4e6
        elif m157b and lab[:3] == 'max':
            return lon+5.e5, lat+4e5
        return lon+3.7e5, lat-6e5
    else:                   # Top right
        lon, lat = mm(lon, lat)
        if w15 and lab[:7] == 'feconix':
            return lon-5e6, lat-1.5e6
        elif w15 and lab[:3] == 'max':
            return lon-3.2e6, lat
        elif n20 and lab[:7] == 'feconix':
            return lon-5.2e6, lat-1.1e6
        elif n20 and lab[:3] == 'max':
            return lon+5.e5, lat-1.2e6
        elif b15 and lab[:7] == 'feconix':
            return lon+5.3e5, lat-1.2e6
        elif b15 and lab[:3] == 'max':
            return lon-3e6, lat+4e5
        elif iib and lab[:7] == 'feconix':
            return lon-5.6e6, lat+5.2e5
        elif iib and lab[:3] == 'max':
            return lon-3e6, lat-1e6
        return lon-4.5e6, lat-2e5

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

def ovrplt(mm):
    # Parse the input
    kicks = {}
    if 'iib' in sys.argv[1]:
        tmp = 'kicks/IIb_mas_cen.csv'
    elif 'b15_lmc_' in sys.argv[1]:
        tmp = 'kicks/B15_LMC_mas_cen.csv'
    else:
        tmp = 'kicks/' + sys.argv[1].split('_')[0].upper() + '_mas_cen.csv'
    with open(tmp) as ff:
        for line in ff:
            line = line.split()
            buf = np.array([float(line[1]), float(line[2]), float(line[3])])
            kicks[line[0]] = buf

    # Plots the points
    tmp = np.loadtxt('min_max_dir/'+mod+'.txt')
    min_dir = tmp[0]
    max_dir = tmp[1]
    for ii in range(3, len(sys.argv)):
        if sys.argv[ii] == 'min' or sys.argv[ii] == 'max':
            if sys.argv[ii] == 'min':
                lon, lat = mm(np.rad2deg(min_dir[0]), np.rad2deg(min_dir[1]))
                lo2, la2 = shift_lon_lat(mm, np.rad2deg(min_dir[0]), np.rad2deg(min_dir[1]),sys.argv[ii])
                col = 'w'
            elif sys.argv[ii] == 'max':
                lon, lat = mm(np.rad2deg(max_dir[0]), np.rad2deg(max_dir[1]))
                lo2, la2 = shift_lon_lat(mm, np.rad2deg(max_dir[0]), np.rad2deg(max_dir[1]),sys.argv[ii])
                col = 'k'

            mm.scatter(lon, lat)
            annotation = sys.argv[ii]
            plt.text(lo2, la2, fix_ann(annotation), color=col, fontsize=10)
            continue
        
        kick = kicks[sys.argv[ii]]
        phi = np.arctan2(kick[1], kick[0])
        tht = np.arctan(kick[2]/np.sqrt(kick[0]**2+kick[1]**2))
        lon, lat = mm(np.rad2deg(phi), np.rad2deg(tht))
        print(sys.argv[ii], phi, tht)
        mm.scatter(lon, lat)
        annotation = sys.argv[ii] + ' ' + str(int(np.round(np.sqrt(np.sum(kick**2)))))
        lon, lat = shift_lon_lat(mm, np.rad2deg(phi), np.rad2deg(tht), annotation)
        if annotation[:2] == 'ns':
            col = 'w'
        else:
            col = 'k'
        plt.text(lon, lat, fix_ann(annotation), color=col, fontsize=10)

        

################################################################
# Main
sky = np.loadtxt('dat/' + sys.argv[1] + '.txt')
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
if 'iib' in sys.argv[1]:
    annotation =  'IIb ' + re.findall('([0-9]{0,9})d', sys.argv[1])[0] + ' d'
else:
    annotation = sys.argv[1].split('_')[0].upper() + ' ' + re.findall('([0-9]{0,9})d', sys.argv[1])[0] + ' d'
plt.text(0, 17e6, fix_ann(annotation), color='k', fontsize=10)
ovrplt(mm)
fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/' + sys.argv[1] + '.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
#plt.show()
#pdb.set_trace()
