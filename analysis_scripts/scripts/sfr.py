import gallus as gallus
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument('snapshot')
parser.add_argument('do_fit')
args = parser.parse_args()

snapshot = str(args.snapshot).zfill(3)
do_fit = int(args.do_fit)

folder = '/home/rennehan/Documents/Test/Protocluster'
config = 'protocluster.yml'
snap = '%s/snapshot_%s.hdf5' % (folder, snapshot)

nbins = 300

def over_x(x, a):
    return a / x

def hyp_func(x, a, c):
    return a / np.sinh(c * x)

def exp_func(x, a, c, d):
    return a * np.exp(-c * x) + d

bbox_lim = 1e4
bbox_vec = [-bbox_lim, bbox_lim]
bbox = [bbox_vec, bbox_vec, bbox_vec]

sa = gallus.SnapshotAnalysis(config, snap, bounding_box = bbox)

time_range = [0, float(sa.ds.current_time.in_units('Myr'))]
x = np.linspace(time_range[0], time_range[1], nbins)

initmass = sa.ad[('PartType4', 'StellarInitialMass')]
initmass = sa.ds.arr(initmass, 'code_mass')
initmass = initmass.in_units('Msun')

times = sa.ad[('PartType4', 'StellarFormationTime')]
times = sa.ds.arr(times, 'code_time')
times = times.in_units('Myr')

sfr = sa.get_sfr(initmass, times, time_range, nbins)

# Convert to yr^-1 
sfr *= 1.0e-6

sfr[np.isnan(sfr)] = 0

if do_fit:
    popt, pcov = curve_fit(exp_func, np.delete(x, 0), np.delete(sfr, 0))

    print('Exp. fit')
    print(popt)

plt.figure()
plt.ylabel('SFR [M$_\odot$ yr$^{-1}$]', fontsize = 15)
plt.xlabel('Time [Myr]', fontsize = 15)
plt.ylim([0, 1e2])
plt.xlim([np.amin(x), np.amax(x)])
plt.plot(x, sfr, label = 'SFR', lw = 2, alpha = 0.7, c = 'k')

if do_fit:
    plt.plot(x, exp_func(x, popt[0], popt[1], popt[2]), 
                label ='(%.2E)*exp(- t / (%.2E)) + (%.2E)' % (popt[0], 1.0 / popt[1], popt[2]), 
                lw = 2.5, ls = '--', alpha = 0.9)

plt.legend(loc = 'upper right', fontsize = 11)
plt.savefig('%s/sfr_%s.png' % (folder, snapshot), dpi = 300)
plt.show()
plt.close()

