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

ns_mas = {'B15': 1.24,
          'B15_lmc': 1.24,
          'B15_87a': 1.24,
          'B15_lmc_min': 1.24,
          'B15_lmc_max': 1.24,
          'N20': 1.45,
          'L15': 1.58,
          'W15': 1.37,
          'IIb': 1.37,
          'BHH': 1.24,
          'B1D': 1.24,
          'BEO': 1.24,
          'BEN': 1.24,
          'HMM': 1.4,
          'B15_1d': 1.24,
          'M157b': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b+': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b++': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b-1': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b-2': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b-2_min': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b-2_max': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b-3': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M157b-4': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M158b': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M158b-': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M164a': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M164a-': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M167b': fsolve(mg2mb, 1.4, args=(12, 1.55))[0], # ONLY THIS ONE IN MENON19
          'M167b+': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M167b-1': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M167b-2': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M177a': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M177a+': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M178a': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M178a+': fsolve(mg2mb, 1.4, args=(12, 1.55))[0],
          'M15': 1.4,
          'W7': 1.e-5}

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
    
## All
#ne = 11
#den_ele = np.zeros(den.shape + (ne,))
#den_ele[:,:,:,  0] = den*(get_dat(dat, 'p'   ))
#den_ele[:,:,:,  1] = den*(get_dat(dat, 'he4' ))
#den_ele[:,:,:,  2] = den*(get_dat(dat, 'c12' ))
#den_ele[:,:,:,  3] = den*(get_dat(dat, 'o16' ))
#den_ele[:,:,:,  4] = den*(get_dat(dat, 'ne20'))
#den_ele[:,:,:,  5] = den*(get_dat(dat, 'mg24'))
#den_ele[:,:,:,  6] = den*(get_dat(dat, 'si28'))
#den_ele[:,:,:,  7] = den*(get_dat(dat, 's32' ))
#den_ele[:,:,:,  8] = den*(get_dat(dat, 'ar36'))
#den_ele[:,:,:,  9] = den*(get_dat(dat, 'ca40')+get_dat(dat, 'ca44')+get_dat(dat, 'cr48')+get_dat(dat, 'ti44')+get_dat(dat, 'sc44'))
#den_ele[:,:,:, 10] = den*(get_dat(dat, 'fe56')+get_dat(dat, 'fe52')+get_dat(dat, 'co56')+get_dat(dat, 'ni56')+get_dat(dat, 'x56' ))
# wongwathanarat15
ne = 7
den_ele = np.zeros(den.shape + (ne,))
den_ele[:,:,:, 0] = den*(get_dat(dat, 'p'   ))
den_ele[:,:,:, 1] = den*(get_dat(dat, 'he4' ))
den_ele[:,:,:, 2] = den*(get_dat(dat, 'c12' ))
den_ele[:,:,:, 3] = den*(get_dat(dat, 'o16' ))
den_ele[:,:,:, 4] = den*(get_dat(dat, 'si28'))
den_ele[:,:,:, 5] = den*(get_dat(dat, 'ca44')+get_dat(dat, 'ti44')+get_dat(dat, 'sc44'))
#den_ele[:,:,:, 5] = den*(get_dat(dat, 'ti44'))
den_ele[:,:,:, 6] = den*(get_dat(dat, 'fe56')+get_dat(dat, 'co56')+get_dat(dat, 'ni56')+get_dat(dat, 'x56' ))
plt_leg = ['H', 'He', 'C', 'O', 'Si', 'Ti', 'Ni+X']
mas_rad = np.sum(vol[:,:,:,np.newaxis]*den_ele, axis=(0, 1))/Msun
mas_tot = np.cumsum(np.sum(vol*den, axis=(0, 1)))/Msun+ns_mas[ff]
mas_nor = np.sum(mas_rad, axis=0)
print_masses(mas_nor, plt_leg)



################################################################
# Rebin uniformly in mass
nm = np.int(np.ceil(mas_tot[-1]/dm))
step = 0
reb = np.zeros((den.shape[2], 4))
for ii in range(0, den.shape[2]):
    while mas_tot[ii] > step*dm:
        step += 1
    step -= 1
    
    if ii == 0:
        reb[ii, :2] = (step, 1.)
    elif mas_tot[ii-1] >= (step)*dm:
        reb[ii, :2] = (step, 1.)
    else:
        lin_int=(step*dm-mas_tot[ii-1])/(mas_tot[ii]-mas_tot[ii-1])
        reb[ii, :2] = (step-1, lin_int)
        reb[ii, 2:] = (step, 1-lin_int)


rad_dis = np.zeros((ne, nm))
for ii in range(0, den.shape[2]):
    idx = reb[ii,0].astype('int')
    frac = reb[ii,1]
    rad_dis[:, idx] = rad_dis[:, idx] + frac*mas_rad[ii, :]/(dm*mas_nor)
    idx = reb[ii,2].astype('int')
    frac = reb[ii,3]
    rad_dis[:, idx] = rad_dis[:, idx] + frac*mas_rad[ii, :]/(dm*mas_nor)

xx = np.repeat(np.linspace(0, nm*dm, nm+1), 2)
hlp_dis = np.zeros((ne, 2*nm+2))
hlp_dis[:, 1:-2:2] = rad_dis
hlp_dis[:, 2:-1:2] = rad_dis
fig = plt.figure(figsize=(5, 3.75))
for ii in range(0, ne):
    plt.semilogy(xx, hlp_dis[ii])

plt.legend(plt_leg, ncol=2)
plt.ylim([5e-4, 3e0])
plt.ylabel('$(\Delta M_I/M_I)/\Delta M(r)$ ([0.2~M$_\odot$]$^{-1}$)')
plt.xlabel('Mass (M$_\odot$)')
np.savetxt('dat/' + ff.lower() + '_rad_dis_mas.txt', np.r_[np.linspace(0, nm*dm, nm+1)[np.newaxis, 1:], rad_dis].T, header='outer_mass, ' + ', '.join(plt_leg))
fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/' + ff.lower() + '_rad_dis_mas.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)



################################################################
# Rebin uniformly in velocity
nv = np.int(np.ceil(rad2vel*rad[-1]/dv))
step = 0
reb = np.zeros((den.shape[2], 4))
for ii in range(0, den.shape[2]):
    while rad2vel*rad[ii] > step*dv:
        step += 1
    step -= 1
    
    if ii == 0:
        reb[ii, :2] = (step, 1.)
    elif rad2vel*rad[ii-1] >= (step)*dv:
        reb[ii, :2] = (step, 1.)
    else:
        lin_int=(step*dv-rad2vel*rad[ii-1])/(rad2vel*rad[ii]-rad2vel*rad[ii-1])
        reb[ii, :2] = (step-1, lin_int)
        reb[ii, 2:] = (step, 1-lin_int)


rad_dis = np.zeros((ne, nv))
for ii in range(0, den.shape[2]):
    idx = reb[ii,0].astype('int')
    frac = reb[ii,1]
    rad_dis[:, idx] = rad_dis[:, idx] + frac*mas_rad[ii, :]/(dv*mas_nor)
    idx = reb[ii,2].astype('int')
    frac = reb[ii,3]
    rad_dis[:, idx] = rad_dis[:, idx] + frac*mas_rad[ii, :]/(dv*mas_nor)

xx = np.repeat(np.linspace(0, nv*dv, nv+1), 2)
hlp_dis = np.zeros((ne, 2*nv+2))
hlp_dis[:, 1:-2:2] = rad_dis
hlp_dis[:, 2:-1:2] = rad_dis
fig = plt.figure(figsize=(5, 3.75))
for ii in range(0, ne):
    plt.semilogy(xx/1e3, hlp_dis[ii])

plt.legend(plt_leg, ncol=2)
plt.xlim([-0.2, 7.5])
plt.ylim([3e-9, 3e-3])
plt.ylabel('$(\Delta M_I/M_I)/\Delta v(r)$ ([100~km~s$^{-1}$]$^{-1}$)')
plt.xlabel('Velocity (1000~km~s$^{-1}$)')
np.savetxt('dat/' + ff.lower() + '_rad_dis_vel.txt', np.r_[np.linspace(0, nv*dv, nv+1)[np.newaxis, 1:], rad_dis].T, header='outer_velocity, ' + ', '.join(plt_leg))
fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/' + ff.lower() + '_rad_dis_vel.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)



################################################################
# Centers of mass
mas_ele = vol[:,:,:,np.newaxis]*den_ele
vv0 = np.zeros((mid_phi.size, mid_tht.size, mid_rad.size , 1, 3))
vv1 = np.zeros((mid_phi.size, mid_tht.size, mid_rad.size , 1, 3))
vv0[:,:,:,0,0] = mid_rad*rad2vel*np.cos(mid_phi)*np.sin(mid_tht)
vv0[:,:,:,0,1] = mid_rad*rad2vel*np.sin(mid_phi)*np.sin(mid_tht)
vv0[:,:,:,0,2] = mid_rad*rad2vel*np.cos(mid_tht)*mid_phi/mid_phi
vv1[:,:,:,0,0] = vex/1e5*np.cos(mid_phi)*np.sin(mid_tht)
vv1[:,:,:,0,1] = vex/1e5*np.sin(mid_phi)*np.sin(mid_tht)
vv1[:,:,:,0,2] = vex/1e5*np.cos(mid_tht)*mid_phi/mid_phi

ns0 = -np.sum(den[:,:,:,np.newaxis]*vol[:,:,:,np.newaxis]*vv0[:,:,:,0,:], axis=(0,1,2))/(Msun*ns_mas[ff])
ns1 = -np.sum(den[:,:,:,np.newaxis]*vol[:,:,:,np.newaxis]*vv1[:,:,:,0,:], axis=(0,1,2))/(Msun*ns_mas[ff])
pp0_ele = np.sum(mas_ele[:,:,:,:, np.newaxis]*vv0, axis=(0,1,2))/(Msun*mas_nor)[:,np.newaxis]
pp1_ele = np.sum(mas_ele[:,:,:,:, np.newaxis]*vv1, axis=(0,1,2))/(Msun*mas_nor)[:,np.newaxis]
pp0_ele = np.r_[pp0_ele, ns0[np.newaxis,:]]
pp1_ele = np.r_[pp1_ele, ns1[np.newaxis,:]]
pp0_nor = np.append(np.sqrt(np.sum(pp0_ele**2, axis=1)), np.sqrt(np.sum(ns0**2)))
pp1_nor = np.append(np.sqrt(np.sum(pp1_ele**2, axis=1)), np.sqrt(np.sum(ns1**2)))

# Plot
plt_leg.append('NS')

fig = plt.figure()
#plt.gca().axis('off')
mm = Basemap(projection='hammer', lon_0 = 0, llcrnrlon=True, llcrnrlat=True)
mm.drawmapboundary(fill_color='w', linewidth=1)
mm.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0]) # draw parallels
mm.drawmeridians(np.arange(-180.,180.,60.)) # draw meridians
fig.set_size_inches(5, 3.75)

for ii, kick in enumerate(pp0_ele):
    phi = np.arctan2(kick[1], kick[0])
    tht = np.arctan(kick[2]/np.sqrt(kick[0]**2+kick[1]**2))
    lon, lat = mm(np.rad2deg(phi), np.rad2deg(tht))
    mm.scatter(lon, lat)

    ann = plt_leg[ii] + ' ' + str(int(np.round(pp0_nor[ii])))
    plt.text(lon, lat, ann, color='k', fontsize=10)

np.savetxt('dat/' + ff.lower() + '_mas_cen.txt', np.c_[np.array(plt_leg), pp0_ele], header='element,           x (km s-1),           y (km s-1),           z (km s-1)', fmt='%9.7s, %20.18s, %20.18s, %20.18s')
fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/' + ff.lower() + '_mas_cen.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

fig = plt.figure()
#plt.gca().axis('off')
mm = Basemap(projection='hammer', lon_0 = 0, llcrnrlon=True, llcrnrlat=True)
mm.drawmapboundary(fill_color='w', linewidth=1)
mm.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0]) # draw parallels
mm.drawmeridians(np.arange(-180.,180.,60.)) # draw meridians
fig.set_size_inches(5, 3.75)

for ii, kick in enumerate(pp1_ele):
    phi = np.arctan2(kick[1], kick[0])
    tht = np.arctan(kick[2]/np.sqrt(kick[0]**2+kick[1]**2))
    lon, lat = mm(np.rad2deg(phi), np.rad2deg(tht))
    mm.scatter(lon, lat)

    ann = plt_leg[ii] + ' ' + str(int(np.round(pp1_nor[ii])))
    plt.text(lon, lat, ann, color='k', fontsize=10)

np.savetxt('dat/' + ff.lower() + '_mas_ce2.txt', np.c_[np.array(plt_leg), pp1_ele], header='element,           x (km s-1),           y (km s-1),           z (km s-1)', fmt='%9.7s, %20.18s, %20.18s, %20.18s')
fig.savefig('/Users/silver/box/phd/pro/sne/grt/fig/' + ff.lower() + '_mas_ce2.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

#pdb.set_trace()
#plt.show()
