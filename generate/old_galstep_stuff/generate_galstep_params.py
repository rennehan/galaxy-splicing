#!/usr/bin/env python2
import numpy as np
import argparse as ap
import itertools

parser = ap.ArgumentParser()
parser.add_argument('gas_mass_resol')
args = parser.parse_args()

gas_mass_resol = float(args.gas_mass_resol)
dm_mass_resol = 10.0 * gas_mass_resol
stellar_mass_resol = dm_mass_resol

fac = 1.0e10

total_Ngas = 0
total_Ndm = 0
total_Nstars = 0

labels = itertools.cycle(('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
                            'K', 'L', 'M', 'N'))

halo_params = {'M_halo': 0,
                'N_halo': 0,
                'a_halo': 47,
                'halo_core': False}

disk_params = {'M_disk': 0,
                'N_disk': 0,
                'Rd': 3.5,
                'z0': 0.7,
                'factor': 0.8}

bulge_params = {'M_bulge': 0,
                'N_bulge': 0,
                'a_bulge': 1.5,
                'bulge_core': False}

gas_params = {'M_gas': 1,
                'N_gas': 0,
                'z0_gas': 0.035,
                'Z': 0.001}

global_params = {'N_rho': 256,
                    'rho_max': 300,
                    'Nz': 256,
                    'z_max': 3000}

# 1st col: Mgas, 2nd col: Mdyn
data = np.loadtxt('data.txt', skiprows = 1, usecols = (1, 2))

# row[0]: Mgas, row[1]: Mdyn
for row in data:
    gas_mass = row[0] * fac
    dm_mass = row[1] * fac - 2.0 * gas_mass
    stellar_mass = gas_mass

    Ngas = int(gas_mass / gas_mass_resol)
    Ndm = int(dm_mass / dm_mass_resol)
    Nstars = int(stellar_mass / stellar_mass_resol)
    
    temp_halo_params = halo_params
    temp_disk_params = disk_params
    temp_bulge_params = bulge_params
    temp_gas_params = gas_params
    
    temp_halo_params['M_halo'] = dm_mass / fac
    temp_halo_params['N_halo'] = Ndm

    temp_disk_params['M_disk'] = stellar_mass / (2.0 * fac)
    temp_disk_params['N_disk'] = Nstars / 2

    temp_bulge_params['M_bulge'] = stellar_mass / (2.0 * fac)
    temp_bulge_params['N_bulge'] = Nstars / 2

    temp_gas_params['M_gas'] = gas_mass / fac
    temp_gas_params['N_gas'] = Ngas

    label = labels.next()
    file_name = '%s/params_galaxy.ini' % (label)

    with open(file_name, 'wb') as f:
        f.write('[halo]\n')

        for key in halo_params:
            f.write(key + ' = ' + str(temp_halo_params[key]) + '\n')

        f.write('\n[disk]\n')

        for key in disk_params:
            f.write(key + ' = ' + str(temp_disk_params[key]) + '\n')

        f.write('\n[bulge]\n')

        for key in bulge_params:
            f.write(key + ' = ' + str(temp_bulge_params[key]) + '\n')

        f.write('\n[gas]\n')

        for key in gas_params:
            f.write(key + ' = ' + str(temp_gas_params[key]) + '\n')

        f.write('\n[global]\n')

        for key in global_params:
            f.write(key + ' = ' + str(global_params[key]) + '\n')

    total_Ngas += Ngas
    total_Ndm += Ndm
    total_Nstars += Nstars

total_part = total_Ngas + total_Ndm + total_Nstars

print('TotNgas: %g' % total_Ngas)
print('TotNdm: %g' % total_Ndm)
print('TotNstars: %g' % total_Nstars)
print('TotPart: %g' % total_part)
