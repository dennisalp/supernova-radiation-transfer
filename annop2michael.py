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

def split_iron_dirty(ca40, ti44, fe52, ni56, env_idx):
    ii = env_idx[ff]

    if not ff in fixed:
        lim = ni56[ii[0], ii[1], ii[2]]
        print('Envelope iron mass fraction is', lim)
        fe52 = np.where(ni56 < 1.01*lim, fe52+ni56, fe52+lim)
        ni56 = np.where(ni56 < 1.01*lim, 0, ni56-lim)

    lim = ti44[ii[0], ii[1], ii[2]]
    print('Envelope titanium mass fraction is', lim)
    ca40 = np.where(ti44 < 1.01*lim, ca40+ti44, ca40+lim)
    ti44 = np.where(ti44 < 1.01*lim, 0, ti44-lim)
    return ca40, ti44, fe52, ni56

################################################################

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
print('Current file:', ff, '(' + out + ')')

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
radius = np.array(dat['rad'])
theta  = np.array(dat['tht' ])
phi    = np.array(dat['phi'   ])
velx    = np.array(dat['velx' ])
vely    = np.array(dat['vely' ])
velz    = np.array(dat['velz' ])
den    = np.array(dat['densty'])
ar36   = np.array(dat['ar36'  ])
c12    = np.array(dat['c12'   ])
ca40   = np.array(dat['ca40'  ])
cr48   = np.array(dat['cr48'  ])
fe52   = np.array(dat['fe52'  ])
he4    = np.array(dat['he4'   ])
mg24   = np.array(dat['mg24'  ])
ne20   = np.array(dat['ne20'  ])
ni56   = np.array(dat['ni56'  ])
co56   = np.array(dat['ni56'  ])*1e-9
fe56   = np.array(dat['ni56'  ])*1e-9
o16    = np.array(dat['o16'   ])
p      = np.array(dat['p1'    ])+np.array(dat['h1'    ])
s32    = np.array(dat['s32'   ])
si28   = np.array(dat['si28'  ])
ti44   = np.array(dat['ti44'  ])
sc44   = np.array(dat['ti44'  ])*1e-9
ca44   = np.array(dat['ti44'  ])*1e-9
x56    = np.array(dat['x56'   ])
n      = np.array(dat['x56'   ])*1e-9
tt = dat['time']

env_idx = {
    'M157b-1.281': [50,50,1600],
    'M157b-1.364': [50,50,1600], # Fixed promoridal iron
    'M157b-2.240': [50,50,1600],
    'M157b-2.327': [50,50,1600], # Fixed promoridal iron
    'M157b-3.270': [50,50,1600],
    'M157b-3.345': [50,50,1600], # Fixed promoridal iron
    'M157b-4.318': [50,50,1600], # Fixed promoridal iron
    'M158b-1.248': [50,50,1600],
    'M158b-2.218': [50,50,1600],
    'M164a-1.207': [50,50,1600],
    'M164a-2.216': [50,50,1600],
    'M167b-2.214': [50,50,1600],
    'M167b-2.271': [50,50,1600], # Fixed promoridal iron
    'M167b-1.199': [50,50,1600],
    'M167b-1.313': [50,50,1600], # Fixed promoridal iron
    'M177a-1.179': [50,50,1600],
    'M177a-2.176': [50,50,1600],
    'M178a-3.241': [50,50,1600],
    'M178a-4.214': [50,50,1600]
    }

#pdb.set_trace()
fixed = [
    'M157b-1.364',
    'M157b-2.327',
    'M157b-3.345',
    'M157b-4.318',
    'M167b-2.271',
    'M167b-1.313']
ca40, ti44, fe52, ni56 = split_iron_dirty(ca40, ti44, fe52, ni56, env_idx)

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
vol = mid_rad.astype('double')**2*np.sin(mid_tht)*dif_rad*dif_tht*dif_phi
nphi = mid_phi.shape[0]
ntht = mid_tht.shape[1]
nrad = mid_rad.shape[2]

vex = velx.copy()

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
#make_1d()
#add_oxygen(0.005)
#add_nickel(0.001)
#den, fe56, co56, ni56, x56 = mix_n_move(1., 500)

# Check stuff
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1])/Msun
print('Total mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*p[:,:,:-1])/Msun
print('Hydrogen mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*o16[:,:,:-1])/Msun
print('Oxygen mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*(ni56[:,:,:-1]+co56[:,:,:-1]+fe56[:,:,:-1]))/Msun
print('Nickel mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*x56[:,:,:-1])/Msun
print('X mass:', mas)
sum_all = ar36+c12+ca40+ca44+co56+cr48+fe52+fe56+he4+mg24+ne20+ni56+o16+s32+sc44+si28+ti44+x56+p+n
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
out.create_dataset('time'   , data=tt   )
out.close()

pdb.set_trace()
