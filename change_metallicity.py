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

def change_z():
    global ar36, c12, ca40, ca44, co56, cr48, fe52, fe56, he4, mg24, n, ne20, ni56, o16, s32, sc44, si28, ti44, x56, p
    comp = get_mas_comp()
    for ii in range(nphi):
        for jj in range(ntht):
            for kk in range(nrad):
#solar_num = 10**np.array([12, 10.899, 8.39, 8.69, 7.87, 7.55, 7.54, 7.19, 6.55, 6.34, 4.92, 5.65, 7.47])
#                          H   He      C     O     Ne    Mg    Si    S     Ar    Ca    Ti    Cr    Fe
# Ca = ca40+ca44+sc44+ti44
# Ti = cr48
# Cr = fe52
# Fe = fe56+co56+ni56+x56
                if p[ii,jj,kk] > 0.1:
                    if comp[1,ii,jj,kk] < solar_mas[1]: # compare by number
                        he4[ii,jj,kk] = solar_mas[1] # assign by mass
                    if comp[2,ii,jj,kk] < solar_mas[2]:
                        c12[ii,jj,kk] = solar_mas[2]
                    if comp[3,ii,jj,kk] < solar_mas[3]:
                        o16[ii,jj,kk] = solar_mas[3]
                    if comp[4,ii,jj,kk] < solar_mas[4]:
                        ne20[ii,jj,kk] = solar_mas[4]
                    if comp[5,ii,jj,kk] < solar_mas[5]:
                        mg24[ii,jj,kk] = solar_mas[5]
                    if comp[6,ii,jj,kk] < solar_mas[6]:
                        si28[ii,jj,kk] = solar_mas[6]
                    if comp[7,ii,jj,kk] < solar_mas[7]:
                        s32[ii,jj,kk] = solar_mas[7]
                    if comp[8,ii,jj,kk] < solar_mas[8]:
                        ar36[ii,jj,kk] = solar_mas[8]
                    if comp[9,ii,jj,kk] < solar_mas[9]:
                        ca40[ii,jj,kk] = solar_mas[9]
                    if comp[10,ii,jj,kk] < solar_mas[10]:
                        cr48[ii,jj,kk] = solar_mas[10]
                    if comp[11,ii,jj,kk] < solar_mas[11]:
                        fe52[ii,jj,kk] = solar_mas[11]
                    if comp[12,ii,jj,kk] < solar_mas[12]:
# This one's special because both Cr and Fe are assigned to fe52. This is to avoid the Fortran MC code taking *56 as fresh nickel.
                        fe52[ii,jj,kk] = fe52[ii,jj,kk]+solar_mas[12]-(fe56[ii,jj,kk]+co56[ii,jj,kk]+ni56[ii,jj,kk]+x56[ii,jj,kk])
    p = 1-n-he4-c12-o16-ne20-mg24-si28-s32-ar36-ca40-ca44-sc44-ti44-cr48-fe52-fe56-co56-ni56-x56

    


################################################################
# Constants, cgs
newz = float(sys.argv[4])
out = sys.argv[3]
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
tt = {'B15': 13489751.4508,
      'N20': 12488928.8787,
      'L15': 12630166.8468,
      'W15': 12755626.5895,
      'IIb': 1565943.04402}
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
if 'lmc' in out.lower():
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
    solar_num = 10**np.array([12, He, C, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe])
elif '87a' in out.lower():
    # lundqvist96, mattila10
    He = np.average(np.array([11.40, 11.23]))
    # hunter07b, lundqvist96
    C = np.average(np.array([7.75, 7.51]))
    # hunter07b, dewey08, zhekov09, mattila10
    O = np.average(np.array([8.35, 7.85, 7.85, 8.20, 8.27]))
    # dewey08, zhekov09, mattila10
    Ne = np.average(np.array([7.56, 7.55, 7.93]))
    # 2*trundle07, hunter07b
    Mg = np.average(np.array([7.06, 7.08, 7.05]))
    # 2*trundle07, hunter07b
    Si = np.average(np.array([7.19, 7.21, 7.20]))
    # dewey08, zhekov09, mattila10
    S = np.average(np.array([6.69, 6.68, 7.12]), weights=1/np.array([0.11, 0.09, (0.19+0.34)/2]))
    # zhekov09, mattila10
    Ar = np.average(np.array([6.29, 6.23]))
    # zhekov09, mattila10
    Ca = np.average(np.array([5.89, 6.51]))
    # russell92
    Ti = 4.81
    # russell89, russell90
    Cr = np.average(np.array([5.47, 5.35]), weights=1/np.array([0.14, 0.22]))
    # dewey08, zhekov09, mattila10, dewey12 = 6.98293488
    Fe = np.average(np.array([6.97, 6.97, 6.98, 7.022182518111363]), weights=1/np.array([0.04, 0.02, (0.19+0.34)/2, 0.04]))
    solar_num = 10**np.array([12, He, C, O, Ne, Mg, Si, S, Ar, Ca, Ti, Cr, Fe])    
else:
    # From lodders03
    solar_num = 10**np.array([12, 10.899, 8.39, 8.69, 7.87, 7.55, 7.54, 7.19, 6.55, 6.34, 4.92, 5.65, 7.47])
    #                         H   He      C     O     Ne    Mg    Si    S     Ar    Ca    Ti    Cr    Fe
    # Ca = ca40+ca44+sc44+ti44
    # Ti = cr48
    # Cr = fe52
    # Fe = fe56+co56+ni56+x56

# Note that number abundances are normalized to hydrogen whereas mass fractions are normalized!
solar_mas = solar_num*np.array([ato_mas['p'], ato_mas['he4'], ato_mas['c12'], ato_mas['o16'], ato_mas['ne20'], ato_mas['mg24'], ato_mas['si28'], ato_mas['s32'], ato_mas['ar36'], ato_mas['ca40'], ato_mas['ti44'], ato_mas['cr48'], ato_mas['fe52']])
solar_mas = solar_mas/np.sum(solar_mas)
solar_mas[2:] = newz*solar_mas[2:]
solar_mas[:2] = solar_mas[:2]*(1-solar_mas.sum()+solar_mas[:2].sum())/solar_mas[:2].sum()
solar_num = solar_mas/np.array([ato_mas['p'], ato_mas['he4'], ato_mas['c12'], ato_mas['o16'], ato_mas['ne20'], ato_mas['mg24'], ato_mas['si28'], ato_mas['s32'], ato_mas['ar36'], ato_mas['ca40'], ato_mas['ti44'], ato_mas['cr48'], ato_mas['fe52']])
solar_num = solar_num/solar_num[0]

################################################################
# Load data
dat = h5py.File(path, 'r')
radius = np.array(dat['radius'])
theta  = np.array(dat['theta' ])
phi    = np.array(dat['phi'   ])
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
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ar36[:,:,:-1])/Msun
print('Argon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*fe52[:,:,:-1])/Msun
print('52Fe mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*(ni56[:,:,:-1]+co56[:,:,:-1]+fe56[:,:,:-1]))/Msun
print('Nickel mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*x56[:,:,:-1])/Msun
print('X mass:', mas)
#pdb.set_trace()
# comp = get_comp()
# plt.semilogy(comp[:,0,0,60])
# plt.semilogy(solar_num)
# plt.show()
change_z()

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
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*ar36[:,:,:-1])/Msun
print('Argon mass:', mas)
mas = np.sum(vol[:,:,:-1]*den[:,:,:-1]*fe52[:,:,:-1])/Msun
print('52Fe mass:', mas)
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

#time python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_01Z 0.10;
#time python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_025Z 0.25;
#time python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_1Z 1.00;
#time python src/change_metallicity.py /Users/silver/dat/sne/ IIb IIb_025Z 0.25;
#time python src/change_metallicity.py /Users/silver/dat/sne/ M15 M15_025Z 0.25;
#time python src/change_metallicity.py /Users/silver/dat/sne/ W15 W15_025Z 0.25;
#time python src/change_metallicity.py /Users/silver/dat/sne/ W15 W15_1Z 1.00;
#time python src/change_metallicity.py /Users/silver/dat/sne/ L15 L15_025Z 0.25;
#time python src/change_metallicity.py /Users/silver/dat/sne/ N20 N20_025Z 0.25;
#time python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_LMC 1.00;
#time python src/change_metallicity.py /Users/silver/dat/sne/ N20 N20_LMC 1.00;

#Current file: B15
#Total mass: 14.212189132818908
#Hydrogen mass: 8.20128056129579
#Helium mass: 5.4397401640371
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.17712335936398704
#Argon mass: 0.004088708558056851
#52Fe mass: 0.0014306813292138457
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#Total mass: 14.212189132818908
#Hydrogen mass: 8.112653243844226
#Helium mass: 5.518869898265972
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.18246632507275942
#Argon mass: 0.0041838677873310345
#52Fe mass: 0.002581476940611863
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_01Z 0.10  136.48s user 11.33s system 113% cpu 2:10.77 total
#Current file: B15
#Total mass: 14.212189132818908
#Hydrogen mass: 8.20128056129579
#Helium mass: 5.4397401640371
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.17712335936398704
#Argon mass: 0.004088708558056851
#52Fe mass: 0.0014306813292138457
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#Total mass: 14.212189132818908
#Hydrogen mass: 8.101135355436696
#Helium mass: 5.515048460209737
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.19145696093090797
#Argon mass: 0.004329198274951195
#52Fe mass: 0.004337704994796238
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_025Z 0.25  134.17s user 10.28s system 112% cpu 2:08.19 total
#Current file: B15
#Total mass: 14.212189132818908
#Hydrogen mass: 8.20128056129579
#Helium mass: 5.4397401640371
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.17712335936398704
#Argon mass: 0.004088708558056851
#52Fe mass: 0.0014306813292138457
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#Total mass: 14.212189132818908
#Hydrogen mass: 8.041865209821093
#Helium mass: 5.495995429796493
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.23741420416659492
#Argon mass: 0.005069360828189608
#52Fe mass: 0.013265953188956302
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_1Z 1.00  134.64s user 10.20s system 112% cpu 2:08.36 total
#Current file: IIb
#Total mass: 3.7253966162860634
#Hydrogen mass: 0.514213735728629
#Helium mass: 2.114386414799729
#Carbon mass: 0.1727925896266994
#Oxygen mass: 0.5501511419423326
#Argon mass: 3.7303426945708057e-10
#52Fe mass: 3.7303426945708057e-10
#Nickel mass: 0.0483925166142115
#X mass: 0.07786205832329468
#Total mass: 3.7253966162860634
#Hydrogen mass: 0.5138165259394974
#Helium mass: 2.114386414799729
#Carbon mass: 0.1727925896266994
#Oxygen mass: 0.550196212656793
#Argon mass: 1.9898945717784122e-05
#52Fe mass: 0.0002409221199306109
#Nickel mass: 0.0483925166142115
#X mass: 0.07786205832329468
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ IIb IIb_025Z 0.25  115.04s user 12.53s system 116% cpu 1:49.85 total
#Current file: M15
#Total mass: 19.536772898058363
#Hydrogen mass: 11.312035625245516
#Helium mass: 6.850956779746657
#Carbon mass: 0.09712753396692723
#Oxygen mass: 0.8200203408168988
#Argon mass: 0.005995491118371232
#52Fe mass: 0.002558577916235859
#Nickel mass: 0.05718706729689512
#X mass: 0.12510936046172738
#Total mass: 19.536772898058363
#Hydrogen mass: 11.242498019049219
#Helium mass: 6.92042749084181
#Carbon mass: 0.09712753396692723
#Oxygen mass: 0.8200203408168988
#Argon mass: 0.005995491118371232
#52Fe mass: 0.002625463249024381
#Nickel mass: 0.05718706729689512
#X mass: 0.12510936046172738
#0.9999999999999993 1.0000000000000004
#python src/change_metallicity.py /Users/silver/dat/sne/ M15 M15_025Z 0.25  106.67s user 8.50s system 111% cpu 1:43.03 total
#Current file: W15
#Total mass: 14.007558915146777
#Hydrogen mass: 7.17153171222636
#Helium mass: 5.4238706786554385
#Carbon mass: 0.22701943409785488
#Oxygen mass: 0.7246974480121224
#Argon mass: 1.4061626352615656e-09
#52Fe mass: 1.4061626352615656e-09
#Nickel mass: 0.05347456720050229
#X mass: 0.08382319000775516
#Total mass: 14.007558915146777
#Hydrogen mass: 7.160003002092208
#Helium mass: 5.430720746772509
#Carbon mass: 0.22701943409785488
#Oxygen mass: 0.7246974480121224
#Argon mass: 0.0003304773250897823
#52Fe mass: 0.0029006312998783325
#Nickel mass: 0.05347456720050229
#X mass: 0.08382319000775516
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ W15 W15_025Z 0.25  179.98s user 14.79s system 110% cpu 2:56.58 total
#Current file: W15
#Total mass: 14.007558915146777
#Hydrogen mass: 7.17153171222636
#Helium mass: 5.4238706786554385
#Carbon mass: 0.22701943409785488
#Oxygen mass: 0.7246974480121224
#Argon mass: 1.4061626352615656e-09
#52Fe mass: 1.4061626352615656e-09
#Nickel mass: 0.05347456720050229
#X mass: 0.08382319000775516
#Total mass: 14.007558915146777
#Hydrogen mass: 7.144761071023422
#Helium mass: 5.4303398595564065
#Carbon mass: 0.22701943409785488
#Oxygen mass: 0.7258254086746742
#Argon mass: 0.0013219092643004205
#52Fe mass: 0.012042813479808174
#Nickel mass: 0.05347456720050229
#X mass: 0.08382319000775516
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ W15 W15_1Z 1.00  213.95s user 32.73s system 106% cpu 3:50.93 total
#Current file: L15
#Total mass: 13.683830958370002
#Hydrogen mass: 7.29234855735665
#Helium mass: 5.257746242635396
#Carbon mass: 0.17762957112160016
#Oxygen mass: 0.5765619607146504
#Argon mass: 1.3736937994090285e-09
#52Fe mass: 1.3736937994090285e-09
#Nickel mass: 0.0316283653550527
#X mass: 0.11684156752250163
#Total mass: 13.683830958370002
#Hydrogen mass: 7.287462972539514
#Helium mass: 5.258338951560417
#Carbon mass: 0.17762957112160016
#Oxygen mass: 0.5765619607146504
#Argon mass: 0.00032379950824999146
#52Fe mass: 0.0027033809948652852
#Nickel mass: 0.0316283653550527
#X mass: 0.11684156752250163
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ L15 L15_025Z 0.25  170.56s user 14.54s system 111% cpu 2:45.43 total
#Current file: N20
#Total mass: 14.280358526880576
#Hydrogen mass: 5.667003291724694
#Helium mass: 6.520113533109648
#Carbon mass: 0.09933509897204688
#Oxygen mass: 1.3295919981839628
#Argon mass: 0.0046841515204573666
#52Fe mass: 0.001894863718842457
#Nickel mass: 0.04226743820930174
#X mass: 0.08171661813033346
#Total mass: 14.280358526880576
#Hydrogen mass: 5.637229423974678
#Helium mass: 6.520185015661118
#Carbon mass: 0.10459681936862726
#Oxygen mass: 1.343692224888993
#Argon mass: 0.0049291897643942125
#52Fe mass: 0.004878934832404471
#Nickel mass: 0.04226743820930174
#X mass: 0.08171661813033346
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ N20 N20_025Z 0.25  145.12s user 13.72s system 113% cpu 2:19.61 total
#Current file: B15
#Total mass: 14.212189132818908
#Hydrogen mass: 8.20128056129579
#Helium mass: 5.4397401640371
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.17712335936398704
#Argon mass: 0.004088708558056851
#52Fe mass: 0.0014306813292138457
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#Total mass: 14.212189132818908
#Hydrogen mass: 7.957332517308803
#Helium mass: 5.637094463931908
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.2033663702307689
#Argon mass: 0.00490392536598756
#52Fe mass: 0.008099059728199796
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#0.9999999999999993 1.0000000000000009
#python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_LMC 1.00  147.35s user 19.25s system 112% cpu 2:28.19 total
#Current file: N20
#Total mass: 14.280358526880576
#Hydrogen mass: 5.667003291724694
#Helium mass: 6.520113533109648
#Carbon mass: 0.09933509897204688
#Oxygen mass: 1.3295919981839628
#Argon mass: 0.0046841515204573666
#52Fe mass: 0.001894863718842457
#Nickel mass: 0.04226743820930174
#X mass: 0.08171661813033346
#Total mass: 14.280358526880576
#Hydrogen mass: 5.61605879311802
#Helium mass: 6.520213351801876
#Carbon mass: 0.10486965557618345
#Oxygen mass: 1.3548955470712014
#Argon mass: 0.005506132170290785
#52Fe mass: 0.008695870781190876
#Nickel mass: 0.04226743820930174
#X mass: 0.08171661813033346
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ N20 N20_LMC 1.00  140.16s user 12.61s system 109% cpu 2:19.09 total
#Current file: B15
#Total mass: 14.212189132818908
#Hydrogen mass: 8.20128056129579
#Helium mass: 5.4397401640371
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.17712335936398704
#Argon mass: 0.004088708558056851
#52Fe mass: 0.0014306813292138457
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#Total mass: 14.212189132818908
#Hydrogen mass: 6.0712073651350895
#Helium mass: 7.54533189183694
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.18779504420990406
#Argon mass: 0.004450247729826842
#52Fe mass: 0.004222792246961864
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_87A 1.00  131.78s user 13.38s system 97% cpu 2:29.22 total
#Current file: B15
#Total mass: 14.212189132818908
#Hydrogen mass: 8.20128056129579
#Helium mass: 5.4397401640371
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.17712335936398704
#Argon mass: 0.004088708558056851
#52Fe mass: 0.0014306813292138457
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#Total mass: 14.212189132818908
#Hydrogen mass: 7.967072137492051
#Helium mass: 5.640393514160189
#Carbon mass: 0.12143025755178678
#Oxygen mass: 0.19185974241364437
#Argon mass: 0.004581348017537986
#52Fe mass: 0.005233465125880631
#Nickel mass: 0.03450609215785225
#X mass: 0.06890804520244967
#0.9999999999999993 1.0000000000000007
#python src/change_metallicity.py /Users/silver/dat/sne/ B15 B15_87A_TMP 1.00  125.95s user 12.34s system 97% cpu 2:22.27 total
#
