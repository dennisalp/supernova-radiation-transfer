#time python -i src/change_metallicity.py /Users/silver/dat/sne/ B15_025Z
import os
import sys
import pdb
from glob import glob

import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_dat(dat, lab):
    tmp = np.array(dat[lab])
    tmp[:,:,-1] = tmp[:, :, -2]
    return tmp

def get_mas_comp():
    tmp = np.array([p, # 0 H
                  he4, # 1 He
                  c12, # 2 C
                  o16, # 3 O
                 ne20, # 4 Ne
                 mg24, # 5 Mg
                 si28, # 6 Si
                  s32, # 7 S
                 ar36, # 8 Ar
                 ca40+ # 9 Ca
                 ca44+
                 sc44+
                 ti44, 
                 cr48, # 10 Ti
                 fe52, # 11 Cr
                 fe56+ # 12 Fe
                 co56+
                 ni56+
                  x56])
    return tmp

def get_comp():
    tmp = np.array([p/ato_mas['p'],    # 0 H
                  he4/ato_mas[ 'he4'], # 1 He
                  c12/ato_mas[ 'c12'], # 2 C
                  o16/ato_mas[ 'o16'], # 3 O
                 ne20/ato_mas['ne20'], # 4 Ne
                 mg24/ato_mas['mg24'], # 5 Mg
                 si28/ato_mas['si28'], # 6 Si
                  s32/ato_mas[ 's32'], # 7 S
                 ar36/ato_mas['ar36'], # 8 Ar
                 ca40/ato_mas['ca40']+ # 9 Ca
                 ca44/ato_mas['ca44']+
                 sc44/ato_mas['sc44']+
                 ti44/ato_mas['ti44'], 
                 cr48/ato_mas['cr48'], # 10 Ti
                 fe52/ato_mas['fe52'], # 11 Cr
                 fe56/ato_mas['fe56']+ # 12 Fe
                 co56/ato_mas['co56']+
                 ni56/ato_mas['ni56']+
                  x56/ato_mas[ 'x56']])
    return tmp/tmp[0,:,:,:]



################################################################
# Constants, cgs
newz = float(sys.argv[3])
ff = sys.argv[2]
path = sys.argv[1] + ff + '.h5'
print ('Current file:', ff)

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

env_idx = {
    'B15': [40,40,670],
    'N20': [20,20,700],
    'L15': [70,70,999],
    'W15': [80,80,900],
    'M15': [50,50,1600],
    'M16': [50,50,1600],
    'HMM': [ 0, 0,600],
    'IIb': [10,10,650]}

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
# Compute abundance patterns
#russell90, kurt98
He = np.average(np.array([10.96, 10.91]), weights=1/np.array([0.06, 0.05]))
#hunter07b, russell89, russell90, kurt98
C = np.average(np.array([7.75, 8.04, 7.66, 7.81]))
#kurt98, 2*russell90, 2*trundle07, hunter07b
O = np.average(np.array([8.37, 8.37, 8.25, 8.33, 8.40, 8.35]))
#kurt98, 2*russell90
Ne = np.average(np.array([7.55, 7.68, 7.50]))
#2*trundle07, hunter07b
Mg = np.average(np.array([7.06, 7.08, 7.05]))
#2*trundle07, hunter07b
Si = np.average(np.array([7.19, 7.21, 7.20]))
#2*russell90
S = np.average(np.array([6.87, 6.62]), weights=1/np.array([0.14, 0.23]))
#2*russell90
Ar = np.average(np.array([6.07, 6.59]), weights=1/np.array([0.25, 0.07]))
#russell89, russell90
Ca = np.average(np.array([5.89, 6.05]), weights=1/np.array([0.16, 0.03]))
#russell92
Ti = 4.81
#russell89, russell90
Cr = np.average(np.array([5.47, 5.35]), weights=1/np.array([0.14, 0.22]))
#russell89, russell90, 2*trundle07
Fe = np.average(np.array([7.23, 7.22, 7.23, 7.24]), weights=1/np.array([0.17, 0.09, 0.10, 0.12]))
lmc_num = 10**np.array([12, He, C, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe])

# From lodders03
solar_num = 10**np.array([12, 10.899, 8.39, 8.69, 7.87, 7.55, 7.54, 7.19, 6.55, 6.34, 4.92, 5.65, 7.47])
#                         H   He      C     O     Ne    Mg    Si    S     Ar    Ca    Ti    Cr    Fe
# Ca = ca40+ca44+sc44+ti44
# Ti = cr48
# Cr = fe52
# Fe = fe56+co56+ni56+x56

# Note that number abundances are normalized to hydrogen whereas mass fractions are normalized!
all_masses = np.array([ato_mas['p'], ato_mas['he4'], ato_mas['c12'], ato_mas['o16'], ato_mas['ne20'], ato_mas['mg24'], ato_mas['si28'], ato_mas['s32'], ato_mas['ar36'], ato_mas['ca40'], ato_mas['ti44'], ato_mas['cr48'], ato_mas['fe52']])
lmc_mas = lmc_num*all_masses
lmc_mas = lmc_mas/np.sum(lmc_mas)
lmc_num = lmc_mas/all_masses
lmc_num = lmc_num/lmc_num[0]

solar_mas = solar_num*all_masses
solar_mas = solar_mas/np.sum(solar_mas)
solar_mas[2:] = newz*solar_mas[2:]
solar_mas[:2] = solar_mas[:2]*(1-solar_mas.sum()+solar_mas[:2].sum())/solar_mas[:2].sum()
solar_num = solar_mas/all_masses
solar_num = solar_num/solar_num[0]

# Cross sections
ver = np.loadtxt('sig_verner96.txt')

################################################################
# Load data
dat = h5py.File(path, 'r')
radius = np.array(dat['radius']).astype('double')
theta  = np.array(dat['theta' ]).astype('double')
phi    = np.array(dat['phi'   ]).astype('double')
vex    = get_dat(dat, 'vex'   ).astype('double')
den    = get_dat(dat, 'den'   ).astype('double')
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

comp = get_mas_comp()
tmp = env_idx[ff[:3]]
ele = 'H He C O Ne Mg Si S Ar Ca Ti Cr Fe'.split()
# These are primordial
# Ca = ca40+ca44+sc44+ti44
# Ti = cr48
# Cr = fe52
# Fe = fe56+co56+ni56+x56
for ii in range(comp.shape[0]):
    el = comp[ii, tmp[0], tmp[1], tmp[2]]
    if 'lmc' in ff.lower() or ff == 'M167b' or ff == 'M157b+':
        if ii == 11:
            print('{0:5s} {1:13.10f}'.format(ele[ii]+'+'+ele[12], (el+comp[12, tmp[0], tmp[1], tmp[2]])/(lmc_mas[ii]+lmc_mas[12])))
        elif ii == 12:
            continue
        else:
            print('{0:5s} {1:13.10f}'.format(ele[ii], el/lmc_mas[ii]))
    else:
        if ii == 11:
            print('{0:5s} {1:13.10f}'.format(ele[ii]+'+'+ele[12], (el+comp[12, tmp[0], tmp[1], tmp[2]])/(solar_mas[ii]+solar_mas[12])))
        elif ii == 12:
            continue
        else:
            print('{0:5s} {1:13.10f}'.format(ele[ii], el/solar_mas[ii]))


print(comp[11:, tmp[0], tmp[1], tmp[2]], solar_mas[11:]/newz)
print('Metallicity based on iron peak = {0:7.5f} Z'.format(newz*comp[11:, tmp[0], tmp[1], tmp[2]].sum()/solar_mas[11:].sum()))
if np.amin(fe52)<0:
    print('ERROR This should never be negative, it could be of how fe52 is assigned by change_metallicity:', np.amin(fe52))

ene = np.logspace(0, np.log10(4e6), 501)
cs_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15]
sig_sol = []
norm_sol = (solar_num[:2]*all_masses[:2]).sum()+(solar_num[2:]*all_masses[2:]/newz).sum()
e0 = np.argmax(ene>3e4)-1
for ii, sn in enumerate(solar_num):
    if ii > 1:
        sn /= newz
    plt.loglog(ene/1e3, sn*ver[cs_idx[ii], :]/norm_sol)
    sig_sol.append(sn*ver[cs_idx[ii], e0]/norm_sol)
plt.legend(ele)

plt.figure()
num_comp = get_comp()
sig_tot = []
norm = (num_comp[:, tmp[0], tmp[1], tmp[2]]*all_masses).sum()
for ii, sn in enumerate(num_comp[:, tmp[0], tmp[1], tmp[2]]):
    plt.loglog(ene/1e3, sn*ver[cs_idx[ii], :]/norm)
    sig_tot.append(sn*ver[cs_idx[ii], e0]/norm)

plt.legend(ele)
sig_sol = np.array(sig_sol)
sig_sol[11] = sig_sol[11]+sig_sol[12]
sig_sol = sig_sol[:-1]
sig_tot = np.array(sig_tot)
sig_tot[11] = sig_tot[11]+sig_tot[12]
sig_tot = sig_tot[:-1]
print('Metallicity based on photoabsorption opacity at 30 keV = {0:7.5f} Z'.format(sig_tot.sum()/sig_sol.sum()))

#plt.show()
#pdb.set_trace()

#python src/check_metallicity.py /Users/silver/dat/sne/ B15 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ B15_LMC 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ B15_01Z 0.10;
#python src/check_metallicity.py /Users/silver/dat/sne/ B15_025Z 0.25;
#python src/check_metallicity.py /Users/silver/dat/sne/ B15_1Z 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ IIb_025Z 0.25;
#python src/check_metallicity.py /Users/silver/dat/sne/ IIb 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ M15_025Z 0.25;
#python src/check_metallicity.py /Users/silver/dat/sne/ M15 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ W15_025Z 0.25;
#python src/check_metallicity.py /Users/silver/dat/sne/ W15_1Z 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ W15 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ L15_025Z 0.25;
#python src/check_metallicity.py /Users/silver/dat/sne/ L15 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ N20_025Z 0.25;
#python src/check_metallicity.py /Users/silver/dat/sne/ N20 1.00;
#python src/check_metallicity.py /Users/silver/dat/sne/ N20_LMC 1.00

#Current file: B15
#H      1.0203884537
#He     0.9724006231
#C      1.5269381405
#O      0.0090817603
#Ne     0.0000261437
#Mg     0.0000001576
#Si     0.0000001382
#S      0.0000002708
#Ar     0.0000010507
#Ca     0.0007008927
#Ti     0.0000366726
#Cr+Fe  0.0000019488
#[1.00004936e-10 2.15932484e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.00000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.02582 Z
#Current file: B15_LMC
#H      0.9962584140
#He     1.0000000000
#C      5.8085708961
#O      1.0000000000
#Ne     1.0000000000
#Mg     1.0000000000
#Si     1.0000000000
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0014628884
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[6.58237560e-04 2.15932484e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.56776 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.54664 Z
#Current file: B15_01Z
#H      0.9958728972
#He     1.0000000000
#C     15.2693814050
#O      1.0000000000
#Ne     1.0000000000
#Mg     1.0000000000
#Si     1.0000000000
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0069935903
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[1.15934433e-04 2.15932484e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.10000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.12461 Z
#Current file: B15_025Z
#H      0.9962999562
#He     1.0000000000
#C      6.1077525620
#O      1.0000000000
#Ne     1.0000000000
#Mg     1.0000000000
#Si     1.0000000000
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0027974361
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[2.89839322e-04 2.15932484e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.25000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.27250 Z
#Current file: B15_1Z
#H      0.9984589762
#He     1.0000000000
#C      1.5269381405
#O      1.0000000000
#Ne     1.0000000000
#Mg     1.0000000000
#Si     1.0000000000
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0006993590
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[1.15936377e-03 2.15932484e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 1.00000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 1.00407 Z
#Current file: IIb_025Z
#H      0.8626245144
#He     1.3733658349
#C     10.1546708995
#O      5.6478874227
#Ne     6.3531743566
#Mg     4.3770430962
#Si     7.5114720954
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0004776658
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[2.89840875e-04 6.06552173e-10] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.25000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.55365 Z
#Current file: IIb
#H      0.8711745525
#He     1.3860818191
#C      2.5386677249
#O      1.4119718557
#Ne     1.5882935892
#Mg     1.0942607740
#Si     1.8778680239
#S      0.0000002708
#Ar     0.0000010507
#Ca     0.0001209521
#Ti     0.0000366708
#Cr+Fe  0.0000006094
#[1.00000000e-10 6.06552173e-10] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.00000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.35703 Z
#Current file: M15_025Z
#H      0.7880294382
#He     1.6476498010
#C      6.2820468622
#O      2.5824199602
#Ne     1.9757889807
#Mg     1.9814609974
#Si     1.5270860382
#S      1.4158471247
#Ar     1.2069201562
#Ca     1.7087908434
#Ti    16.5719787301
#Cr+Fe  1.9619791240
#[3.99383361e-06 5.64669103e-04] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.49049 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.48971 Z
#Current file: M15
#H      0.7953311140
#He     1.6629053784
#C      1.5705117156
#O      0.6456049901
#Ne     0.4939472452
#Mg     0.4953652494
#Si     0.3817715095
#S      0.3539617812
#Ar     0.3017300391
#Ca     0.4271977109
#Ti     4.1429946825
#Cr+Fe  0.4870499381
#[9.72444224e-14 5.64669103e-04] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.48705 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.48297 Z
#Current file: W15_025Z
#H      0.8626739583
#He     1.3732084195
#C     10.1531685346
#O      5.6484871692
#Ne     6.3531746431
#Mg     4.3770432367
#Si     7.5114722682
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0033292674
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[2.89839406e-04 2.07618031e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.25000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.55367 Z
#Current file: W15_1Z
#H      0.8689724333
#He     1.3859229462
#C      2.5382921337
#O      1.4121217923
#Ne     1.5882936608
#Mg     1.0942608092
#Si     1.8778680671
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0008323168
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[1.15936385e-03 2.07618031e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 1.00000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 1.12377 Z
#Current file: W15
#H      0.8712244523
#He     1.3859229462
#C      2.5382921337
#O      1.4121217923
#Ne     1.5882936608
#Mg     1.0942608092
#Si     1.8778680671
#S      0.0000002708
#Ar     0.0000010507
#Ca     0.0008338526
#Ti     0.0000366708
#Cr+Fe  0.0000018770
#[1.00000000e-10 2.07618031e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.00000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.35704 Z
#Current file: L15_025Z
#H      0.8831621440
#He     1.3238078411
#C      3.7761099986
#O      6.0673002040
#Ne     6.1257303111
#Mg     3.3936208827
#Si     3.7764071480
#S      1.0000000000
#Ar     1.0000000000
#Ca     3.8472838053
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[2.89839416e-04 2.06610237e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.25000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.48752 Z
#Current file: L15
#H      0.8918806388
#He     1.3360649682
#C      0.9440274997
#O      1.5168250510
#Ne     1.5314325778
#Mg     0.8484052207
#Si     0.9441017870
#S      0.0000002708
#Ar     0.0000010507
#Ca     0.9618209513
#Ti     0.0000366708
#Cr+Fe  0.0000018684
#[1.00000000e-10 2.06610237e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.00000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.29577 Z
#Current file: N20_025Z
#H      0.6806090456
#He     2.0142118480
#C      1.0000000000
#O      1.0000000000
#Ne     1.0000000000
#Mg     1.0000000000
#Si     1.0000000000
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0025928452
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[2.89839454e-04 2.02801029e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.25000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.26784 Z
#Current file: N20
#H      0.6909690596
#He     2.0328614207
#C      0.0000000455
#O      0.0000000171
#Ne     0.0000000905
#Mg     0.0000001576
#Si     0.0000001382
#S      0.0000002708
#Ar     0.0000010507
#Ca     0.0006497449
#Ti     0.0000366708
#Cr+Fe  0.0000018355
#[1.00000000e-10 2.02801029e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.00000 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.01847 Z
#Current file: N20_LMC
#H      0.6923844007
#He     1.9038308174
#C      1.0000000000
#O      1.0000000000
#Ne     1.0000000000
#Mg     1.0000000000
#Si     1.0000000000
#S      1.0000000000
#Ar     1.0000000000
#Ca     1.0013558998
#Ti     1.0000000000
#Cr+Fe  1.0000000000
#[6.58237691e-04 2.02801029e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.56776 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.54173 Z
#Current file: M167b
#H      0.8620157068
#He     1.4040877457
#C      3.4283854261
#O      0.6362997911
#Ne     0.9500382450
#Mg     1.5374653760
#Si     0.8453921126
#S      0.9306911103
#Ar     0.3619288454
#Ca     0.8936005193
#Ti     5.4021799629
#Cr+Fe  0.8578603524
#[5.6467776e-04 1.1293555e-12] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.48706 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.47192 Z
#Current file: M157b+
#H      0.8755205657
#He     1.3643897097
#C      3.3487850982
#O      0.6556501233
#Ne     0.9500369882
#Mg     1.5374660877
#Si     0.8453894402
#S      0.9306844798
#Ar     0.3619290288
#Ca     0.8936006360
#Ti     5.4021825722
#Cr+Fe  0.8578608830
#[5.6467811e-04 1.1293562e-12] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.48706 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.47208 Z
#Current file: B15_87A
#H      0.7250106148
#He     1.9001608365
#C      1.5269381405
#O      0.1891430390
#Ne     0.4707500369
#Mg     0.2377517357
#Si     0.3332655041
#S      0.2665844897
#Ar     0.3739300458
#Ca     0.5288895878
#Ti     0.5659660270
#Cr+Fe  0.2402221060
#[2.78503165e-04 2.15932484e-09] [1.59753344e-05 1.14339059e-03]
#Metallicity based on iron peak = 0.24022 Z
#Metallicity based on photoabsorption opacity at 30 keV = 0.28476 Z
