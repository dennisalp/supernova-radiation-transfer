'''
2020-10-16, Dennis Alp, dalp@kth.se

Plots Compton and pair production cross sections.
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob

import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter
import matplotlib.font_manager

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

xx = np.logspace(-1, 2, 100)
plt.loglog(xx,(1./3.+0.141*xx-0.12*xx**2+(1.+0.5*xx)*(1.+xx)**(-2))*4.989346146e-25)
plt.loglog(xx,((np.log(1.+xx)+0.06)/xx)*4.989346146e-25)
plt.loglog(xx,((np.log(1.+xx)+0.5-1./(2.+0.076*xx))/xx)*4.989346146e-25)
plt.loglog(xx,(0.10063*(xx-1.022))*1.e-27)
plt.loglog(xx,(0.0481+0.301*(xx-1.5))*1.e-27)
plt.show()
