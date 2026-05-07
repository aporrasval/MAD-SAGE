# Adam Stevens, 2022
# Post-analysis of output of PSO with Dark Sage to check fidelity etc

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0, '/home/ajp79/project/fs1/DarkSage/2024DarkSage/Dirty-AstroPy')
from galprops import galread as gr
from galprops import galcalc as gc
from galprops import galplot as gp

import warnings
warnings.filterwarnings("ignore")

fsize = 20
matplotlib.rcParams.update({'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10, 'xtick.major.width': 1, 'ytick.major.width': 1, 'ytick.minor.size': 5, 'xtick.minor.size': 5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 1, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman', 'legend.numpoints': 1, 'legend.columnspacing': 1, 'legend.fontsize': fsize-4, 'lines.markeredgewidth': 1.0, 'errorbar.capsize': 4.0, 'xtick.top': True, 'ytick.right': True})

base_dir = '/home/ajp79/project/fs1/DarkSage/priya_proj/runGUREFT90/DarkSage-astevens-disc/autocalibration/01/'
track_dir = base_dir + 'tracks/'
dump_file = base_dir + 'full_dump.txt'

# lists that will contain the x and y data for the observations and Dark Sage runs that were used to compute the goodness of fit
SMF_x_obs, SMF_y_obs, SMF_y_mod = [], [], []
QLF_x_obs, QLF_y_obs, QLF_y_mod = [], [], []

# HIMF_x_obs, HIMF_y_obs, HIMF_y_mod = [], [], []
# CSFH_x_obs, CSFH_y_obs, CSFH_y_mod = [], [], []

f = open(dump_file, 'r')
lines = f.readlines()
on_a_data_line = False # tracks if the current line in the dump is data of interest
last_line = False
obs = ''
line_start = ''

for l, line in enumerate(lines):
#    print('doing line', l, on_a_data_line, last_line, obs, line_start)

    # a line that specifies a data file tells us which constraint the next data belong to
    if 'weibel24_SMFz6.csv' in line:
        obs = 'SMF'
    elif 'JWST_QLFz6.csv' in line:
        obs = 'QLF'
    # elif 'CSFH_DSILVA+23_Ver_Final.csv' in line:
    #     obs = 'CSFH'
    
    # catches an instance where an array of data we want starts
    if line[:5]=='obs x' or line[:5]=='obs y' or line[:5]=='mod y': 
        line_start = line[:5]
        on_a_data_line = True
        line = line[8:]
        data = np.array([], dtype=np.float32) # start appending data into an array
        
    if on_a_data_line:
        if line[-2:]==']\n':
            last_line = True
            line = line[:-2]
        else:
            last_line = False
            line = line[:-1]

        line_data = np.array(line.split())
        data = np.append(data, line_data.astype(np.float32))
            
        if last_line:
            
            if obs == 'SMF' and line_start == 'obs x':
                SMF_x_obs += [data]
            if obs == 'SMF' and line_start == 'obs y':
                SMF_y_obs += [data]
            if obs == 'SMF' and line_start == 'mod y':
                SMF_y_mod += [data]

            if obs == 'QLF' and line_start == 'obs x':
                QLF_x_obs += [data]
            if obs == 'QLF' and line_start == 'obs y':
                QLF_y_obs += [data]
            if obs == 'QLF' and line_start == 'mod y':
                QLF_y_mod += [data]

            # if obs == 'CSFH' and line_start == 'obs x':
            #     CSFH_x_obs += [data]
            # if obs == 'CSFH' and line_start == 'obs y':
            #     CSFH_y_obs += [data]
            # if obs == 'CSFH' and line_start == 'mod y':
            #     CSFH_y_mod += [data]
            
            on_a_data_line = False
    else:
        last_line = False


print(len(SMF_x_obs), len(SMF_y_obs), len(SMF_y_mod))
print(len(QLF_x_obs), len(QLF_y_obs), len(QLF_y_mod))

# print(len(HIMF_x_obs), len(HIMF_y_obs), len(HIMF_y_mod))
# print(len(CSFH_x_obs), len(CSFH_y_obs), len(CSFH_y_mod))

print("QLF Lbol: ", QLF_x_obs)
print("QLF Phi: ", QLF_y_mod)

fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)
ax1.plot(SMF_x_obs[0], SMF_y_obs[0], 'ko')
ax2.plot(QLF_x_obs[0], QLF_y_obs[0], 'ko')

# ax2.plot(HIMF_x_obs[0], HIMF_y_obs[0], 'ko')
# ax3.plot(CSFH_x_obs[0], CSFH_y_obs[0], 'ko')

for i in range(np.min([len(SMF_x_obs), len(QLF_x_obs)])):
    ax1.plot(SMF_x_obs[i], SMF_y_mod[i], '-')
    ax2.plot(QLF_x_obs[i], QLF_y_mod[i], '-')

    #print("QLF Lbol: ", QLF_x_obs[i])
    #print("QLF Phi: ", QLF_y_mod[i])

    # ax2.plot(HIMF_x_obs[i], HIMF_y_mod[i], '-')
    # ax3.plot(CSFH_x_obs[i], CSFH_y_mod[i], '-')

ax1.axis([7,11.1,-6,0.0])
ax2.axis([43,48,-10,-2.0])
# ax3.axis([0,14,-2.7,-0.5])

ax1.set_xlabel(r'$\log_{10}(m_*~[{\rm M}_\odot])$')
ax1.set_ylabel(r'$\log_{10}(\Phi_*~[{\rm cMpc}^{-3}\,{\rm dex}^{-1}])$')

ax2.set_xlabel(r'$\log_{10}(L_{\rm bol}~[erg \ s^{-1}])$')
ax2.set_ylabel(r'$\log_{10}(\Phi_{\rm HI}~[{\rm cMpc}^{-3}\,{\rm dex}^{-1}])$')

# ax3.set_xlabel(r'Look-back time [Gyr]')
# ax3.set_ylabel(r'$\log_{10}(\bar{\rho}_{\rm SFR}~[{\rm M}_\odot\,{\rm yr}^{-1}\,{\rm cMpc}^{-3}])$')

gp.savepng(base_dir+'reconstructed_constraints', xpixplot=2000, ypixplot=500)
