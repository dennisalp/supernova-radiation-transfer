#time python -i src/toy.py /Users/silver/dat/sne/ B15
import os
import sys
import pdb
from glob import glob

import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D

def get_dat(dat, lab):
    tmp = np.array(dat[lab])
    tmp[:,:,-1] = tmp[:, :, -2]
    return tmp

def sph2car(phi, tht):
    xx = np.cos(phi)*np.sin(tht)
    yy = np.sin(phi)*np.sin(tht)
    zz = np.cos(tht)
    return xx, yy, zz

def reduce_hydrogen(xx):
    global den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    frac = 1-(1-xx)*p
    den  = den*frac
    ar36 = ar36/frac
    c12  = c12 /frac
    ca40 = ca40/frac
    ca44 = ca44/frac
    co56 = co56/frac
    cr48 = cr48/frac
    fe52 = fe52/frac
    fe56 = fe56/frac
    he4  = he4 /frac
    mg24 = mg24/frac
    n    = n   /frac
    ne20 = ne20/frac
    ni56 = ni56/frac
    o16  = o16 /frac
    s32  = s32 /frac
    sc44 = sc44/frac
    si28 = si28/frac
    ti44 = ti44/frac
    x56  = x56 /frac
    p    = p   /frac*xx
def reduce_helium(xx):
    global den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    frac = 1-(1-xx)*he4
    den  = den*frac
    ar36 = ar36/frac
    c12  = c12 /frac
    ca40 = ca40/frac
    ca44 = ca44/frac
    co56 = co56/frac
    cr48 = cr48/frac
    fe52 = fe52/frac
    fe56 = fe56/frac
    p    =    p/frac
    mg24 = mg24/frac
    n    = n   /frac
    ne20 = ne20/frac
    ni56 = ni56/frac
    o16  = o16 /frac
    s32  = s32 /frac
    sc44 = sc44/frac
    si28 = si28/frac
    ti44 = ti44/frac
    x56  = x56 /frac
    he4  = he4 /frac*xx

def get_dir():
    def hlp(xx, ii, jj, wei):
        return np.tile(np.average(xx, axis=(0,1), weights=wei), (nphi, ntht, 1))

    half_opening_angle = np.deg2rad(30)
    
    area = (phi[1]-phi[0])*(theta[1]-theta[0])*np.sin(mid_tht)
    bin_phi, bin_tht = np.meshgrid(phi, theta, indexing='ij')
    mid_phi_loc, mid_tht_loc = np.meshgrid(mid_phi, mid_tht, indexing='ij')
    xx, yy, zz = sph2car(mid_phi_loc, mid_tht_loc)
    xyz = np.dstack((xx, yy, zz))

    tmp = np.loadtxt('min_max_dir/'+ ff +'.txt')
    if 'min' in out:
        phi_dir, tht_dir = tmp[0]
    elif 'max' in out:
        phi_dir, tht_dir = tmp[1]
    else:
        print('ERROR Must be min or max direction')
        sys.exit(1)
#    phi_dir = float(sys.argv[4])
#    phi_dir = int(np.floor((phi_dir/(2*np.pi)+0.5)*(phi.size-1)))
#    tht_dir = float(sys.argv[5])
#    tht_dir = int(np.floor((tht_dir/np.pi+0.5)*(theta.size-1)))

    tht_dir = np.pi/2.-tht_dir
#    phi_dir = int(np.floor((phi_dir/(2*np.pi)+0.5)*(phi.size-1)))
#    tht_dir = int(np.floor((tht_dir/np.pi)*(theta.size-1)))

    xyz_dir = np.array([sph2car(phi_dir, tht_dir)])
    cos_angle = np.sum(xyz_dir*xyz, axis=2)
    cos_angle = np.where(cos_angle > np.cos(half_opening_angle), 1., 0.)
    wei = np.tile(cos_angle[:,:,np.newaxis], (1, 1, nrad))
    
    global den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    den  = hlp(den , phi_dir, tht_dir, wei)
    ar36 = hlp(ar36, phi_dir, tht_dir, wei)
    c12  = hlp(c12 , phi_dir, tht_dir, wei)
    ca40 = hlp(ca40, phi_dir, tht_dir, wei)
    ca44 = hlp(ca44, phi_dir, tht_dir, wei)
    co56 = hlp(co56, phi_dir, tht_dir, wei)
    cr48 = hlp(cr48, phi_dir, tht_dir, wei)
    fe52 = hlp(fe52, phi_dir, tht_dir, wei)
    fe56 = hlp(fe56, phi_dir, tht_dir, wei)
    he4  = hlp(he4 , phi_dir, tht_dir, wei)
    mg24 = hlp(mg24, phi_dir, tht_dir, wei)
    n    = hlp(n   , phi_dir, tht_dir, wei)
    ne20 = hlp(ne20, phi_dir, tht_dir, wei)
    ni56 = hlp(ni56, phi_dir, tht_dir, wei)
    o16  = hlp(o16 , phi_dir, tht_dir, wei)
    s32  = hlp(s32 , phi_dir, tht_dir, wei)
    sc44 = hlp(sc44, phi_dir, tht_dir, wei)
    si28 = hlp(si28, phi_dir, tht_dir, wei)
    ti44 = hlp(ti44, phi_dir, tht_dir, wei)
    x56  = hlp(x56 , phi_dir, tht_dir, wei)
    p    = hlp(p   , phi_dir, tht_dir, wei)
    
def make_1d():
    def hlp(xx, wei):
        return np.tile(np.average(xx, axis=(0,1), weights=wei), (nphi, ntht, 1))
    
    global den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    wei = np.tile(np.sin(mid_tht), (nphi, 1, nrad))
    we2  = den * wei
    den  = hlp(den , wei)
    ar36 = hlp(ar36, we2)
    c12  = hlp(c12 , we2)
    ca40 = hlp(ca40, we2)
    ca44 = hlp(ca44, we2)
    co56 = hlp(co56, we2)
    cr48 = hlp(cr48, we2)
    fe52 = hlp(fe52, we2)
    fe56 = hlp(fe56, we2)
    he4  = hlp(he4 , we2)
    mg24 = hlp(mg24, we2)
    n    = hlp(n   , we2)
    ne20 = hlp(ne20, we2)
    ni56 = hlp(ni56, we2)
    o16  = hlp(o16 , we2)
    s32  = hlp(s32 , we2)
    sc44 = hlp(sc44, we2)
    si28 = hlp(si28, we2)
    ti44 = hlp(ti44, we2)
    x56  = hlp(x56 , we2)
    p    = hlp(p   , we2)
    
def add_oxygen(xx):
    global den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    xtra = xx*p*den
    den  = den+xtra
    frac = (den-xtra)/den
    ar36 = ar36*frac
    c12  = c12 *frac
    ca40 = ca40*frac
    ca44 = ca44*frac
    co56 = co56*frac
    cr48 = cr48*frac
    fe52 = fe52*frac
    fe56 = fe56*frac
    he4  = he4 *frac
    mg24 = mg24*frac
    n    = n   *frac
    ne20 = ne20*frac
    ni56 = ni56*frac
    s32  = s32 *frac
    sc44 = sc44*frac
    si28 = si28*frac
    ti44 = ti44*frac
    x56  = x56 *frac
    p    = p   *frac
    o16  = o16 *frac+xtra/den

def add_nickel(xx):
    global den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    xtra = xx*p*den
    den  = den+xtra
    frac = (den-xtra)/den
    ar36 = ar36*frac
    c12  = c12 *frac
    ca40 = ca40*frac
    ca44 = ca44*frac
    co56 = co56*frac
    cr48 = cr48*frac
    fe52 = fe52*frac
    fe56 = fe56*frac
    he4  = he4 *frac
    mg24 = mg24*frac
    n    = n   *frac
    ne20 = ne20*frac
    ni56 = ni56*frac+xtra/den
    o16  = o16 *frac
    s32  = s32 *frac
    sc44 = sc44*frac
    si28 = si28*frac
    ti44 = ti44*frac
    x56  = x56 *frac
    p    = p   *frac
 
def mix_n_move(mix, move):
    global den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    tmp_radius = mix*radius
    tmp_mid_rad = (tmp_radius[:-1]+tmp_radius[1:])/2.
    tmp_xx = tmp_mid_rad*np.sin(mid_tht)*np.cos(mid_phi)
    tmp_yy = tmp_mid_rad*np.sin(mid_tht)*np.sin(mid_phi)
    tmp_zz = tmp_mid_rad*np.cos(mid_tht)*mid_phi/mid_phi # get shape right

    max_dir = np.array([0.7776338150935392, 1.07111704269264+np.pi/2]) # b15_lmc_map_300d_smo, feconix
    shift_move = move*1e5*tt[sys.argv[2]]
    shift_xx = shift_move*np.sin(max_dir[1])*np.cos(max_dir[0])
    shift_yy = shift_move*np.sin(max_dir[1])*np.sin(max_dir[0])
    shift_zz = shift_move*-np.cos(max_dir[1]) # get shape right
    tmp_xx += shift_xx
    tmp_yy += shift_yy
    tmp_zz += shift_zz

    tmp_radius = np.sqrt(tmp_xx**2+tmp_yy**2+tmp_zz**2)
    tmp_phi = np.arctan2(tmp_yy,tmp_xx)
    tmp_tht = np.arccos(tmp_zz/tmp_radius)
    tmp_fe56 = np.zeros(fe56.shape)
    tmp_co56 = np.zeros(co56.shape)
    tmp_ni56 = np.zeros(ni56.shape)
    tmp_x56 = np.zeros(x56.shape)

    rad_ind = np.searchsorted(radius, tmp_radius)-1
    rad_ind = np.where(rad_ind >= 1200, 1199, rad_ind)
    rad_ind = np.where(rad_ind < 0, 0, rad_ind)
    tht_ind = np.searchsorted(theta, tmp_tht)-1
    tht_ind = np.where(tht_ind >= 90, 89, tht_ind)
    tht_ind = np.where(tht_ind < 0, 0, tht_ind)
    phi_ind = np.searchsorted(phi, tmp_phi)-1
    phi_ind = np.where(phi_ind >= 180, 179, phi_ind)
    phi_ind = np.where(phi_ind < 0, 0, phi_ind)

    feconix = fe56+co56+ni56+x56
    new_den = den*(1-feconix)
    frac = ar36+c12+ca40+ca44+cr48+fe52+he4+mg24+n+ne20+o16+s32+sc44+si28+ti44+p
    ar36 = ar36/frac
    c12  = c12 /frac
    ca40 = ca40/frac
    ca44 = ca44/frac
    cr48 = cr48/frac
    fe52 = fe52/frac
    he4  = he4 /frac
    mg24 = mg24/frac
    n    = n   /frac
    ne20 = ne20/frac
    o16  = o16 /frac
    s32  = s32 /frac
    sc44 = sc44/frac
    si28 = si28/frac
    ti44 = ti44/frac
    p    = p   /frac
    for ii in tqdm(range(rad_ind.shape[0])):
        for jj in range(rad_ind.shape[1]):
            for kk in range(rad_ind.shape[2]):
                ind = (phi_ind[ii,jj,kk], tht_ind[ii,jj,kk], rad_ind[ii,jj,kk])

                den2add = vol[ii, jj, kk]*den[ii, jj, kk]*feconix[ii, jj, kk]/vol[ind]
                new_den[ind] = new_den[ind]+den2add

                frac = den2add/new_den[ind]
                ar36[ind] = ar36[ind]*(1-frac)
                c12 [ind] = c12 [ind]*(1-frac)
                ca40[ind] = ca40[ind]*(1-frac)
                ca44[ind] = ca44[ind]*(1-frac)
                cr48[ind] = cr48[ind]*(1-frac)
                fe52[ind] = fe52[ind]*(1-frac)
                he4 [ind] = he4 [ind]*(1-frac)
                mg24[ind] = mg24[ind]*(1-frac)
                n   [ind] = n   [ind]*(1-frac)
                ne20[ind] = ne20[ind]*(1-frac)
                o16 [ind] = o16 [ind]*(1-frac)
                s32 [ind] = s32 [ind]*(1-frac)
                sc44[ind] = sc44[ind]*(1-frac)
                si28[ind] = si28[ind]*(1-frac)
                ti44[ind] = ti44[ind]*(1-frac)
                p   [ind] = p   [ind]*(1-frac)
                tmp_fe56[ind] = tmp_fe56[ind]*(1-frac)
                tmp_co56[ind] = tmp_co56[ind]*(1-frac)
                tmp_ni56[ind] = tmp_ni56[ind]*(1-frac)
                tmp_x56 [ind] = tmp_x56 [ind]*(1-frac)
                tmp_fe56[ind] = tmp_fe56[ind]+frac*fe56[ii,jj,kk]/feconix[ii,jj,kk]
                tmp_co56[ind] = tmp_co56[ind]+frac*co56[ii,jj,kk]/feconix[ii,jj,kk]
                tmp_ni56[ind] = tmp_ni56[ind]+frac*ni56[ii,jj,kk]/feconix[ii,jj,kk]
                tmp_x56 [ind] = tmp_x56 [ind]+frac*x56[ii,jj,kk]/feconix[ii,jj,kk]

    return new_den, tmp_fe56, tmp_co56, tmp_ni56, tmp_x56
    
################################################################
tt = {'B15': 13489751.4508,
      'B15_LMC': 13489751.4508,
      'N20': 12488928.8787,
      'L15': 12630166.8468,
      'W15': 12755626.5895,
      'IIb': 1565943.04402}

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

out = sys.argv[3]
ff = sys.argv[2]
path = sys.argv[1] + ff + '.h5'
print('Current file:', ff)

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



################################################################
# Load data
dat = h5py.File(path, 'r')
radius = np.array(dat['radius']).astype('double')
theta  = np.array(dat['theta' ]).astype('double')
phi    = np.array(dat['phi'   ]).astype('double')
vex    = np.array(dat['vex'   ]).astype('double')
den    = np.array(dat['den'   ]).astype('double')
ar36   = np.array(dat['ar36'  ]).astype('double')
c12    = np.array(dat['c12'   ]).astype('double')
ca40   = np.array(dat['ca40'  ]).astype('double')
ca44   = np.array(dat['ca44'  ]).astype('double')
co56   = np.array(dat['co56'  ]).astype('double')
cr48   = np.array(dat['cr48'  ]).astype('double')
fe52   = np.array(dat['fe52'  ]).astype('double')
fe56   = np.array(dat['fe56'  ]).astype('double')
he4    = np.array(dat['he4'   ]).astype('double')
mg24   = np.array(dat['mg24'  ]).astype('double')
n      = np.array(dat['n'     ]).astype('double')
ne20   = np.array(dat['ne20'  ]).astype('double')
ni56   = np.array(dat['ni56'  ]).astype('double')
o16    = np.array(dat['o16'   ]).astype('double')
p      = np.array(dat['p'     ]).astype('double')
s32    = np.array(dat['s32'   ]).astype('double')
sc44   = np.array(dat['sc44'  ]).astype('double')
si28   = np.array(dat['si28'  ]).astype('double')
ti44   = np.array(dat['ti44'  ]).astype('double')
x56    = np.array(dat['x56'   ]).astype('double')

mid_rad = (radius[:-1]+radius[1:])/2.
mid_phi = (phi[:-1]+phi[1:])/2.
mid_tht = (theta[:-1]+theta[1:])/2.
mid_phi = mid_phi[:, np.newaxis, np.newaxis]
mid_tht = mid_tht[np.newaxis, :, np.newaxis]
mid_rad = mid_rad[np.newaxis, np.newaxis, :]
dif_phi = np.diff(phi)[:, np.newaxis, np.newaxis]
dif_tht = np.diff(theta)[np.newaxis, :, np.newaxis]
dif_rad = np.diff(radius)[np.newaxis, np.newaxis, :]
xx = mid_rad*np.sin(mid_tht)*np.cos(mid_phi)
yy = mid_rad*np.sin(mid_tht)*np.sin(mid_phi)
zz = mid_rad*np.cos(mid_tht)
vol = mid_rad**2*np.sin(mid_tht)*dif_rad*dif_tht*dif_phi
nphi = mid_phi.shape[0]
ntht = mid_tht.shape[1]
nrad = mid_rad.shape[2]


################################################################
# Do stuff
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1])/Msun
print('Total mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*p[:,:,:-1])/Msun
print('Hydrogen mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*he4[:,:,:-1])/Msun
print('Helium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*c12[:,:,:-1])/Msun
print('Carbon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*o16[:,:,:-1])/Msun
print('Oxygen mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ne20[:,:,:-1])/Msun
print('Neon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*mg24[:,:,:-1])/Msun
print('Magnesium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*si28[:,:,:-1])/Msun
print('Silicon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*s32[:,:,:-1])/Msun
print('Sulphur mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ar36[:,:,:-1])/Msun
print('Argon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ca40[:,:,:-1])/Msun
print('Calcium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*(ti44[:,:,:-1]+sc44[:,:,:-1]+ca44[:,:,:-1]))/Msun
print('Titanium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*(ni56[:,:,:-1]+co56[:,:,:-1]+fe56[:,:,:-1]))/Msun
print('Nickel mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*x56[:,:,:-1])/Msun
print('X mass:', mas)

#reduce_hydrogen(0.)
#reduce_helium(0.)
get_dir()
#make_1d()
#add_oxygen(0.005)
#add_nickel(0.001)
#den, fe56, co56, ni56, x56 = mix_n_move(1., 500)

# Check stuff
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1])/Msun
print('Total mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*p[:,:,:-1])/Msun
print('Hydrogen mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*he4[:,:,:-1])/Msun
print('Helium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*c12[:,:,:-1])/Msun
print('Carbon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*o16[:,:,:-1])/Msun
print('Oxygen mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ne20[:,:,:-1])/Msun
print('Neon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*mg24[:,:,:-1])/Msun
print('Magnesium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*si28[:,:,:-1])/Msun
print('Silicon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*s32[:,:,:-1])/Msun
print('Sulphur mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ar36[:,:,:-1])/Msun
print('Argon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ca40[:,:,:-1])/Msun
print('Calcium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*(ti44[:,:,:-1]+sc44[:,:,:-1]+ca44[:,:,:-1]))/Msun
print('Titanium mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*(ni56[:,:,:-1]+co56[:,:,:-1]+fe56[:,:,:-1]))/Msun
print('Nickel mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*x56[:,:,:-1])/Msun
print('X mass:', mas)
sum_all = ar36+c12+ca40+ca44+co56+cr48+fe52+fe56+he4+mg24+n+ne20+ni56+o16+s32+sc44+si28+ti44+x56+p
print(np.amin(sum_all[:,:,:-1]), np.amax(sum_all[:,:,:-1]))



################################################################
# Print to file
out = h5py.File(path.replace(ff, out), 'w')
out.create_dataset('radius', data=radius)
out.create_dataset('theta' , data=theta )
out.create_dataset('phi'   , data=phi   )
out.create_dataset('vex'   , data=vex   )
out.create_dataset('den'   , data=den   )
out.create_dataset('ar36'  , data=ar36  )
out.create_dataset('c12'   , data=c12   )
out.create_dataset('ca40'  , data=ca40  )
out.create_dataset('ca44'  , data=ca44  )
out.create_dataset('co56'  , data=co56  )
out.create_dataset('cr48'  , data=cr48  )
out.create_dataset('fe52'  , data=fe52  )
out.create_dataset('fe56'  , data=fe56  )
out.create_dataset('he4'   , data=he4   )
out.create_dataset('mg24'  , data=mg24  )
out.create_dataset('n'     , data=n     )
out.create_dataset('ne20'  , data=ne20  )
out.create_dataset('ni56'  , data=ni56  )
out.create_dataset('o16'   , data=o16   )
out.create_dataset('p'     , data=p     )
out.create_dataset('s32'   , data=s32   )
out.create_dataset('sc44'  , data=sc44  )
out.create_dataset('si28'  , data=si28  )
out.create_dataset('ti44'  , data=ti44  )
out.create_dataset('x56'   , data=x56   )
out.close()

#pdb.set_trace()
