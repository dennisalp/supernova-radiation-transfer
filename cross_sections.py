from __future__ import division, print_function
import os
import pdb
import sys
from glob import glob

import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt

print('hoeflich92 has switched order of H and He in Table 3')
print('The high-energy asymptote of Mg is a bit odd in Verner96')

ver = np.loadtxt('sig_verner96.txt')
hof = np.loadtxt('sig_hoeflich92.txt')
ene=np.logspace(-3,np.log10(4e3),501)
ele = 'H He C O Ne Mg Si S Ar 40Ca 44Ca Sc Ti Cr 52Fe 56Fe Co Ni X'.split()
for ii in range(ver.shape[0]):
    print(ii)
    plt.loglog(ene,ver[ii,:])
    plt.loglog(ene,hof[ii,:])
    plt.legend(['verner96', 'heoflich92'])
    plt.title(ele[ii])
    plt.show()

pdb.set_trace()
