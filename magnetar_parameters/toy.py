'''
2020-10-14, Dennis Alp, dalp@kth.se

Make toy SN models.

time python toy.py /Users/silver/dat/sne/ IIb SLSN-I_000 4 3 2
1D version of 001

time python toy.py /Users/silver/dat/sne/ IIb SLSN-I_001 4 3 2
3D version of 000

'''

import os
import sys
from pdb import set_trace as db
from glob import glob

import h5py
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D


################################################################
# Help functions
def get_dat(dat, lab):
    tmp = np.array(dat[lab]).astype('double')
    if len(tmp.shape) == 3:
        tmp[:,:,-1] = tmp[:, :, -2]
    return tmp

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
    
def mk_bubble():
    global den, mas
    for ii in range(100):
        den[:,:,ii] = (ii+1)/100.*den[:,:,ii]
    mas = vol*den
    
def rm_hhe():
    global mas, m_old, ek, den, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p

    den = den*(1-(p+he4))
    p = p*0.
    he4 = he4*0.
    frac = ar36+c12+ca40+ca44+co56+cr48+fe52+fe56+mg24+n+ne20+ni56+o16+s32+sc44+si28+ti44+x56
    
    ar36 = ar36/frac
    c12  = c12 /frac
    ca40 = ca40/frac
    ca44 = ca44/frac
    co56 = co56/frac
    cr48 = cr48/frac
    fe52 = fe52/frac
    fe56 = fe56/frac
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

    mas = vol*den
    ek = (lor-1)*mas*cc**2
    m_old = np.sum(mas)/Msun

def truncate(ii):
    global radius, dif_rad, mid_rad, mas, den, vol, vel, vex, lor, ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p, nrad

    nrad = ii
    print('New time: {0:.5f} s, truncated to {1:d}({2:d}) radial elements\n'.format(tt, ii, ii+1))
    
    mas = mas[:,:,:ii]
    radius = radius[:ii+1]
    mid_rad = mid_rad[:,:,:ii]
    dif_rad = dif_rad[:,:,:ii]
    
    den  = den [:,:,:ii]
    vol  = vol [:,:,:ii]
    vel  = vel [:,:,:ii]
    vex  = vex [:,:,:ii]
    lor  = lor [:,:,:ii]
    ar36 = ar36[:,:,:ii]
    c12  = c12 [:,:,:ii]
    ca40 = ca40[:,:,:ii]
    ca44 = ca44[:,:,:ii]
    co56 = co56[:,:,:ii]
    cr48 = cr48[:,:,:ii]
    fe52 = fe52[:,:,:ii]
    fe56 = fe56[:,:,:ii]
    he4  = he4 [:,:,:ii]
    mg24 = mg24[:,:,:ii]
    n    = n   [:,:,:ii]
    ne20 = ne20[:,:,:ii]
    ni56 = ni56[:,:,:ii]
    o16  = o16 [:,:,:ii]
    s32  = s32 [:,:,:ii]
    sc44 = sc44[:,:,:ii]
    si28 = si28[:,:,:ii]
    ti44 = ti44[:,:,:ii]
    x56  = x56 [:,:,:ii]
    p    = p   [:,:,:ii]
    
def mk_sl(m_new, e_new):
    global tt, den, vex, vel, lor, mas, ek
    den = den*m_new/m_old
    ek = ek*m_new/m_old

    tt = tt/np.sqrt(e_new*1.e51/ek.sum())

    vel = mid_rad/tt
    vex = np.ones(vex.shape)*vel
    truncate(np.argmax(vel > 0.9*cc))
    
    lor = 1/np.sqrt(1-(vel/cc)**2)
    mas = vol*den
    ek = (lor-1)*mas*cc**2

def contrast():
    global den
    old = np.average(den, axis=(0,1))
    den = den**boost
    new = np.average(den, axis=(0,1))
    den = den*old/new

def mk_1d():
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

boost = float(sys.argv[6])
e_new = float(sys.argv[5])
m_new = float(sys.argv[4])
out = sys.argv[3]
ff = sys.argv[2]
path = sys.argv[1] + ff + '.h5'

tmp = '\nCurrent file: {0:s}\nNew mass: {1:.1f} Msun, New energy: {2:.1f} Bethe\n'
tmp = tmp.format(ff, m_new, e_new)
print(tmp)

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
tt = tt[ff]
dat = h5py.File(path, 'r')
radius = get_dat(dat, 'radius')
theta  = get_dat(dat, 'theta' )
phi    = get_dat(dat, 'phi'   )
vex    = get_dat(dat, 'vex'   )
den    = get_dat(dat, 'den'   )
ar36   = get_dat(dat, 'ar36'  )
c12    = get_dat(dat, 'c12'   )
ca40   = get_dat(dat, 'ca40'  )
ca44   = get_dat(dat, 'ca44'  )
co56   = get_dat(dat, 'co56'  )
cr48   = get_dat(dat, 'cr48'  )
fe52   = get_dat(dat, 'fe52'  )
fe56   = get_dat(dat, 'fe56'  )
he4    = get_dat(dat, 'he4'   )
mg24   = get_dat(dat, 'mg24'  )
n      = get_dat(dat, 'n'     )
ne20   = get_dat(dat, 'ne20'  )
ni56   = get_dat(dat, 'ni56'  )
o16    = get_dat(dat, 'o16'   )
p      = get_dat(dat, 'p'     )
s32    = get_dat(dat, 's32'   )
sc44   = get_dat(dat, 'sc44'  )
si28   = get_dat(dat, 'si28'  )
ti44   = get_dat(dat, 'ti44'  )
x56    = get_dat(dat, 'x56'   )


################################################################
# Help variables
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

nphi = mid_phi.shape[0]
ntht = mid_tht.shape[1]
nrad = mid_rad.shape[2]

vol = mid_rad**2*np.sin(mid_tht)*dif_rad*dif_tht*dif_phi
vel = mid_rad/tt
lor = 1/np.sqrt(1-(vel/cc)**2)
mas = vol*den
ek = (lor-1)*mas*cc**2


################################################################
# Do stuff
m_old = np.sum(mas)/Msun
print('Total mass:', m_old)
print('Energy:', ek.sum()/1.e51)
tmp = np.sum(mas*p)/Msun
print('Hydrogen mass:', tmp)
tmp = np.sum(mas*he4)/Msun
print('Helium mass:', tmp)
tmp = np.sum(mas*c12)/Msun
print('Carbon mass:', tmp)
tmp = np.sum(mas*o16)/Msun
print('Oxygen mass:', tmp)
tmp = np.sum(mas*ne20)/Msun
print('Neon mass:', tmp)
tmp = np.sum(mas*mg24)/Msun
print('Magnesium mass:', tmp)
tmp = np.sum(mas*si28)/Msun
print('Silicon mass:', tmp)
tmp = np.sum(mas*s32)/Msun
print('Sulphur mass:', tmp)
tmp = np.sum(mas*ar36)/Msun
print('Argon mass:', tmp)
tmp = np.sum(mas*ca40)/Msun
print('Calcium mass:', tmp)
tmp = np.sum(mas*(ti44+sc44+ca44))/Msun
print('Titanium mass:', tmp)
tmp = np.sum(mas*(ni56+co56+fe56))/Msun
print('Nickel mass:', tmp)
tmp = np.sum(mas*x56)/Msun
print('X mass:', tmp, '\n')

mk_bubble()
rm_hhe()
mk_sl(m_new, e_new)
contrast()
mk_1d()

# Check stuff
tmp = np.sum(mas)/Msun
print('Total mass:', tmp)
print('Energy:', ek.sum()/1.e51)
tmp = np.sum(mas*p)/Msun
print('Hydrogen mass:', tmp)
tmp = np.sum(mas*he4)/Msun
print('Helium mass:', tmp)
tmp = np.sum(mas*c12)/Msun
print('Carbon mass:', tmp)
tmp = np.sum(mas*o16)/Msun
print('Oxygen mass:', tmp)
tmp = np.sum(mas*ne20)/Msun
print('Neon mass:', tmp)
tmp = np.sum(mas*mg24)/Msun
print('Magnesium mass:', tmp)
tmp = np.sum(mas*si28)/Msun
print('Silicon mass:', tmp)
tmp = np.sum(mas*s32)/Msun
print('Sulphur mass:', tmp)
tmp = np.sum(mas*ar36)/Msun
print('Argon mass:', tmp)
tmp = np.sum(mas*ca40)/Msun
print('Calcium mass:', tmp)
tmp = np.sum(mas*(ti44+sc44+ca44))/Msun
print('Titanium mass:', tmp)
tmp = np.sum(mas*(ni56+co56+fe56))/Msun
print('Nickel mass:', tmp)
tmp = np.sum(mas*x56)/Msun
print('X mass:', tmp)
sum_all = ar36+c12+ca40+ca44+co56+cr48+fe52+fe56+he4+mg24+n+ne20+ni56+o16+s32+sc44+si28+ti44+x56+p
print(np.amin(sum_all), np.amax(sum_all))


################################################################
# Plots
# plt.hist(np.log10(den[:,:,200]).ravel(),30, alpha=0.3)
# plt.hist(np.log10(den[:,:,400]).ravel(),30, alpha=0.3)
# plt.hist(np.log10(den[:,:,600]).ravel(),30, alpha=0.3)
# plt.yscale('log')
# plt.show()


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
