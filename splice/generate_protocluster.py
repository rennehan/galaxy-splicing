#!/usr/bin/env python3
import numpy as np
import h5py
import math
import itertools
from scipy.stats import norm
import pickle
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument('resol', help = 'Current resolution of the galaxies (gas resolution). Ex: 2e5, 5e5, 1e6, etc.')
parser.add_argument('--load', help = 'Load a previously stored configuration. Must specify resolution.')
args = parser.parse_args()

# Main base path where everything takes place
folder = '/home/rennehan/ownCloud/Collaboration/protocluster/splice'
mass_resol = str(args.resol) # Make sure this is a string! Old value: '2e5'
do_bulge = False
do_stars = False

mass_table = True       # Is the MassTable present in the header?
density_present = False # Is the gas field Density present?
num_metals = 1          # How many metal fields are there?

fname_out = '%s/spliced_%s.hdf5' % (folder, mass_resol)

# Radius of sphere to place objects in (kpc)
sphere_radius = 65.0

common_fields = ['Coordinates', 'Masses', 'Velocities']
gas_fields = ['Density', 'InternalEnergy', 'Metallicity']

labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']

def rand_rotation_matrix(deflection = 1.0, randnums = None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, 
    competely random rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, 
    they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/
    #gemsiii/rand_rotation.c
    
    if randnums is None:
        randnums = np.random.uniform(size = (3,))
        
    theta, phi, z = randnums
    
    theta = theta * 2.0 * deflection * np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0 * np.pi  # For direction of pole deflection.
    z = z * 2.0 * deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M
    
def sample_velocity():
    v = np.asarray([137.0, 107.0, 830.0, 196.0, 312.0, 623.0, -74.0, -492.0, 
                    537.0, -251.0, 862.0, -147.0, 261.0, 319.0])
    
    mean, std = norm.fit(v)
    v -= mean
    mean = 0
    
    mag = np.abs(np.random.normal(mean, std, 1))
    max_v = np.amax(np.abs(v))

    # Didn't make it within bounds, loop until we do
    while (mag > max_v):
        mag = np.abs(np.random.normal(mean, std, 1))
        
    return mag
    
def generate_random_position(direction = False):
    phi = np.random.uniform(0, 2.0 * np.pi)
    costheta = np.random.uniform(-1, 1)
    u = np.random.uniform(0, 1)

    theta = np.arccos(costheta)
    r = sphere_radius * u

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    if direction == False:
        return x, y, z
    else:
        return x / r, y / r, z / r


def set_new_coords(coords_vec, rotation, rand_vec):
    for i in range(len(coords_vec)):
        temp_coords_vec = np.array([coords_vec[i, 0],
                                    coords_vec[i, 1],
                                    coords_vec[i, 2]])
        new_coords_vec = np.dot(rotation, temp_coords_vec)

        for j in range(0, 3):
            coords_vec[i, j] = new_coords_vec[j] + rand_vec[j]

    return coords_vec


def set_new_vels(vels_vec, rand_vel_vec):
    for i in range(0, 3):
        vels_vec[:, i] += rand_vel_vec[i]

    return vels_vec


def grab_property(f, part_type, field):
    try:
        prop = np.asarray(f['/PartType%d/%s' % (part_type, field)])
    except KeyError:
        prop = np.asarray([])
        #print('KeyError: PartType%d/%s' % (part_type, field))

    return prop

indices = itertools.cycle(np.arange(0, len(labels) + 1, 1))

Ngas = np.zeros(len(labels), dtype = 'uint32')
Ndm = np.copy(Ngas)
Ndisk = np.copy(Ngas)
Nbulge = np.copy(Ngas)
Nstar = np.copy(Ngas)
Nbh = np.copy(Ngas)

gas = {}
dm = {}
disk = {}
bulge = {}
star = {}
bh = {}

for field in common_fields:    
    gas.update({field: np.asarray([])})
    dm.update({field: np.asarray([])})
    disk.update({field: np.asarray([])})
    bulge.update({field: np.asarray([])})
    star.update({field: np.asarray([])})
    bh.update({field: np.asarray([])})

for field in gas_fields:
    gas.update({field: np.asarray([])})
    
    if field == 'Metallicity':
        disk.update({field: np.asarray([])})
        bulge.update({field: np.asarray([])})
        star.update({field: np.asarray([])})
        

for label in labels:
    print('Operating on galaxy %s' % label)
    data_file = '%s/%s/%s.hdf5' % (folder, mass_resol, label)
    log_file = '%s/%s.pkl' % (folder, label)
    old_file = '%s/%s/%s.pkl' % (folder, args.load, label)

    f = h5py.File(data_file, 'r')
    
    label_idx = next(indices)

    gas_coords = grab_property(f, 0, 'Coordinates')
    dm_coords = grab_property(f, 1, 'Coordinates')
    disk_coords = grab_property(f, 2, 'Coordinates')
    bulge_coords = grab_property(f, 3, 'Coordinates')
    star_coords = grab_property(f, 4, 'Coordinates')
    bh_coords = grab_property(f, 5, 'Coordinates')

    Ngas[label_idx] = len(gas_coords)
    Ndm[label_idx] = len(dm_coords)
    Ndisk[label_idx] = len(disk_coords)
    Nbulge[label_idx] = len(bulge_coords)
    Nstar[label_idx] = len(star_coords)
    Nbh[label_idx] = len(bh_coords)

    # The /PartTypeX/Masses dataset does not exist if mass_table is true.
    # In this case we have to set them manually using the MassTable
    # header option.
    if mass_table:
        mass_table_values = f['/Header'].attrs['MassTable']
        
        gas_masses = np.zeros(Ngas[label_idx])
        dm_masses = np.zeros(Ndm[label_idx])
        disk_masses = np.zeros(Ndisk[label_idx])
        bulge_masses = np.zeros(Nbulge[label_idx])
        star_masses = np.zeros(Nstar[label_idx])
        bh_masses = np.zeros(Nbh[label_idx])

        gas_masses[:] = mass_table_values[0]
        dm_masses[:] = mass_table_values[1]
        disk_masses[:] = mass_table_values[2]
        bulge_masses[:] = mass_table_values[3]
        star_masses[:] = mass_table_values[4]
        bh_masses[:] = mass_table_values[5]
    else:
        gas_masses = grab_property(f, 0, 'Masses')
        dm_masses = grab_property(f, 1, 'Masses')
        disk_masses = grab_property(f, 2, 'Masses')
        bulge_masses = grab_property(f, 3, 'Masses')
        star_masses = grab_property(f, 4, 'Masses')
        bh_masses = grab_property(f, 5, 'Masses')

    gas_vels = grab_property(f, 0, 'Velocities')
    dm_vels = grab_property(f, 1, 'Velocities')
    disk_vels = grab_property(f, 2, 'Velocities')
    if do_bulge:
        bulge_vels = grab_property(f, 3, 'Velocities')
    if do_stars:
        star_vels = grab_property(f, 4, 'Velocities')
    bh_vels = grab_property(f, 5, 'Velocities')

    if not args.load:
        rand_x, rand_y, rand_z = generate_random_position()
    
        rand_vec = [rand_x, rand_y, rand_z]

        # Get random velocity magnitude
        rand_mag_vel = sample_velocity()
        rand_rot_vel = rand_rotation_matrix()
    
        # Pick a unit vector
        unit_vec = [0, 0, 1]
        rand_vel = rand_mag_vel * np.dot(rand_rot_vel, unit_vec)
    
        rand_vel_x = rand_vel[0]
        rand_vel_y = rand_vel[1]
        rand_vel_z = rand_vel[2]
   
        rand_vel_vec = [rand_vel_x, rand_vel_y, rand_vel_z]

        rotation = rand_rotation_matrix()
    
        # Save the information for later
        log_data = {'RandomVelocityX':      rand_vel_x,
                    'RandomVelocityY':      rand_vel_y,
                    'RandomVelocityZ':      rand_vel_z,
                    'RandomPositionX':      rand_x,
                    'RandomPositionY':      rand_y,
                    'RandomPositionZ':      rand_z,
                    'RandomRotationMatrix': rotation}
    
        with open(log_file, 'wb') as lf:
            pickle.dump(log_data, lf, protocol = pickle.HIGHEST_PROTOCOL)
    else:
        rand_x = rand_y = rand_z = 0
        rand_vel_x = rand_vel_y = rand_vel_z = 0
        rotation = []

        with open(old_file, 'rb') as lf:
            old_data = pickle.load(lf)

        rand_x = old_data['RandomPositionX']
        rand_y = old_data['RandomPositionY']
        rand_z = old_data['RandomPositionZ']

        rand_vec = [rand_x, rand_y, rand_z]

        rand_vel_x = old_data['RandomVelocityX']
        rand_vel_y = old_data['RandomVelocityY']
        rand_vel_z = old_data['RandomVelocityZ']

        rand_vel_vec = [rand_vel_x, rand_vel_y, rand_vel_z]

        rotation = old_data['RandomRotationMatrix']

    gas_coords = set_new_coords(gas_coords, rotation, rand_vec)
    dm_coords = set_new_coords(dm_coords, rotation, rand_vec)
    disk_coords = set_new_coords(disk_coords, rotation, rand_vec)
    if do_bulge:
        bulge_coords = set_new_coords(bulge_coords, rotation, rand_vec)
    if do_stars:
        star_coords = set_new_coords(star_coords, rotation, rand_vec)
    bh_coords = set_new_coords(bh_coords, rotation, rand_vec)

    gas_vels = set_new_vels(gas_vels, rand_vel_vec)
    dm_vels = set_new_vels(dm_vels, rand_vel_vec)
    disk_vels = set_new_vels(disk_vels, rand_vel_vec)
    if do_bulge:
        bulge_vels = set_new_vels(bulge_vels, rand_vel_vec)
    if do_stars:
        star_vels = set_new_vels(star_vels, rand_vel_vec)
    bh_vels = set_new_vels(bh_vels, rand_vel_vec)

    if label == 'A': 
        gas['Coordinates'] = gas_coords
        gas['Masses'] = gas_masses
        gas['Velocities'] = gas_vels
        
        if density_present:
            gas['Density'] = grab_property(f, 0, 'Density')
            
        gas['InternalEnergy'] = grab_property(f, 0, 'InternalEnergy')
        gas['Metallicity'] = grab_property(f, 0, 'Metallicity')
        
        dm['Coordinates'] = dm_coords
        dm['Masses'] = dm_masses
        dm['Velocities'] = dm_vels
        
        disk['Coordinates'] = disk_coords
        disk['Masses'] = disk_masses
        disk['Velocities'] = disk_vels
       
        if do_bulge:
            bulge['Coordinates'] = bulge_coords
            bulge['Masses'] = bulge_masses
            bulge['Velocities'] = bulge_vels

        if do_stars:
            star['Coordinates'] = star_coords
            star['Masses'] = star_masses
            star['Velocities'] = star_vels

        bh['Coordinates'] = bh_coords
        bh['Masses'] = bh_masses
        bh['Velocities'] = bh_vels
    else:
        gas['Coordinates'] = np.vstack((gas['Coordinates'], gas_coords))
        gas['Masses'] = np.concatenate((gas['Masses'], gas_masses))
        gas['Velocities'] = np.vstack((gas['Velocities'], gas_vels))
        if density_present:
            gas['Density'] = np.concatenate((gas['Density'], grab_property(f, 0, 'Density')))
            
        gas['InternalEnergy'] = np.concatenate((gas['InternalEnergy'], grab_property(f, 0, 'InternalEnergy')))
        
        if num_metals > 1:
            gas['Metallicity'] = np.vstack((gas['Metallicity'], grab_property(f, 0, 'Metallicity')))
        else:
            gas['Metallicity'] = np.concatenate((gas['Metallicity'], grab_property(f, 0, 'Metallicity')))
        
        dm['Coordinates'] = np.vstack((dm['Coordinates'], dm_coords))
        dm['Masses'] = np.concatenate((dm['Masses'], dm_masses))
        dm['Velocities'] = np.vstack((dm['Velocities'], dm_vels))
        
        disk['Coordinates'] = np.vstack((disk['Coordinates'], disk_coords))
        disk['Masses'] = np.concatenate((disk['Masses'], disk_masses))
        disk['Velocities'] = np.vstack((disk['Velocities'], disk_vels))
        
        if do_bulge:
            bulge['Coordinates'] = np.vstack((bulge['Coordinates'], bulge_coords))
            bulge['Masses'] = np.concatenate((bulge['Masses'], bulge_masses))
            bulge['Velocities'] = np.vstack((bulge['Velocities'], bulge_vels))
    
        if do_stars:
            star['Coordinates'] = np.vstack((star['Coordinates'], star_coords))
            star['Masses'] = np.concatenate((star['Masses'], star_masses))
            star['Velocities'] = np.vstack((star['Velocities'], star_vels))

        bh['Coordinates'] = np.vstack((bh['Coordinates'], bh_coords))
        bh['Masses'] = np.concatenate((bh['Masses'], bh_masses))
        bh['Velocities'] = np.vstack((bh['Velocities'], bh_vels))

    f.close()

tot_gas_mass = np.sum(gas['Masses'])
tot_dm_mass = np.sum(dm['Masses'])
tot_disk_mass = np.sum(disk['Masses'])
tot_bulge_mass = np.sum(bulge['Masses'])
tot_star_mass = np.sum(star['Masses'])
tot_bh_mass = np.sum(bh['Masses'])

tot_mass = tot_gas_mass + tot_dm_mass + tot_disk_mass + tot_bulge_mass + tot_star_mass + tot_bh_mass

# We need to remove the bulk velocity of the system. This is important so that
# the group of galaxies doesn't fly away rapidly.
#
# We start by mass weighting the velocities of all of the components, and then
# subtracting away the result of the sum.
print('Correcting for bulk velocities.')
for i in range(0, 3):
    mw_gv = np.sum(gas['Masses'] * gas['Velocities'][:, i])
    mw_dmv = np.sum(dm['Masses'] * dm['Velocities'][:, i])
    mw_dv = np.sum(disk['Masses'] * disk['Velocities'][:, i])
    if do_bulge:
        mw_bv = np.sum(bulge['Masses'] * bulge['Velocities'][:, i])
    else:
        mw_bv = 0
    if do_stars:
        mw_sv = np.sum(star['Masses'] * star['Velocities'][:, i])
    else:
        mw_sv = 0
    mw_bhv = np.sum(bh['Masses'] * bh['Velocities'][:, i])

    bv = (mw_gv + mw_dmv + mw_dv + mw_bv + mw_sv + mw_bhv) / tot_mass
   
    print('Correcting for bulk velocity (comp. %d): %g' % (i, bv))

    gas['Velocities'][:, i] -= bv
    dm['Velocities'][:, i] -= bv
    disk['Velocities'][:, i] -= bv
    if do_bulge:
        bulge['Velocities'][:, i] -= bv
    if do_stars:
        star['Velocities'][:, i] -= bv
    bh['Velocities'][:, i] -= bv

tot_Ngas = int(np.sum(Ngas))
tot_Ndm = int(np.sum(Ndm))
tot_Ndisk = int(np.sum(Ndisk))
tot_Nbulge = int(np.sum(Nbulge))
tot_Nstar = int(np.sum(Nstar))
tot_Nbh = int(np.sum(Nbh))

tot_part = int(tot_Ngas + tot_Ndm + tot_Ndisk + tot_Nbulge + tot_Nstar + tot_Nbh)

Npart = np.asarray([tot_Ngas, tot_Ndm, tot_Ndisk, tot_Nbulge, tot_Nstar, tot_Nbh],
                   dtype = 'uint32')

new_ids = np.arange(1, tot_part + 1, 1, dtype = 'uint32')
  
print('Generating HDF5 file.')

# Prepare new IC file
fp = h5py.File(fname_out, 'w')

header = fp.create_group('Header')

# Set your headers here
header.attrs['MassTable'] = np.zeros(6);
header.attrs['Time'] = 0.0;  # initial time
header.attrs['Redshift'] = 0.0; # initial redshift
header.attrs['BoxSize'] = 1.0e4; # box size
header.attrs['NumFilesPerSnapshot'] = 1; # number of files 
header.attrs['Omega0'] = 0.; # z=0 Omega_matter
header.attrs['OmegaLambda'] = 0.; # z=0 Omega_Lambda
header.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
header.attrs['Flag_Sfr'] = 1; # flag indicating whether star formation is on or off
header.attrs['Flag_Cooling'] = 1; # flag indicating whether cooling is on or off
header.attrs['Flag_StellarAge'] = 1; # flag indicating whether stellar ages are to be saved
header.attrs['Flag_Metals'] = 11; # flag indicating whether metallicity are to be saved
header.attrs['Flag_Feedback'] = 1; # flag indicating whether some parts of springel-hernquist model are active
header.attrs['Flag_DoublePrecision'] = 1; # flag indicating whether ICs are in single/double precision
header.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs

# Running index for the new ParticleIDs
running_idx = 0
next_ids = new_ids[running_idx:running_idx + tot_Ngas]

# Gas
print('Adding gas particles.')
p = fp.create_group('PartType0')

p.create_dataset('ParticleIDs', data = next_ids, dtype = 'uint32')
p.create_dataset('Coordinates', data = gas['Coordinates'], dtype = 'float64')
p.create_dataset('Velocities', data = gas['Velocities'], dtype = 'float64')
p.create_dataset('Masses', data = gas['Masses'], dtype = 'float64')
if density_present:
    p.create_dataset('Density', data = gas['Density'], dtype = 'float64')
p.create_dataset('Metallicity', data = gas['Metallicity'], dtype = 'float64')
p.create_dataset('InternalEnergy', data = gas['InternalEnergy'], 
                 dtype = 'float64')
running_idx += tot_Ngas

next_ids = new_ids[running_idx:running_idx + tot_Ndm]

print('Adding dark matter particles.')
p = fp.create_group('PartType1')
p.create_dataset('ParticleIDs', data = next_ids, dtype = 'uint32')
p.create_dataset('Coordinates', data = dm['Coordinates'], dtype = 'float64')
p.create_dataset('Velocities', data = dm['Velocities'], dtype = 'float64')
p.create_dataset('Masses', data = dm['Masses'], dtype = 'float64')
running_idx += tot_Ndm

next_ids = new_ids[running_idx:running_idx + tot_Ndisk]

print('Adding disk particles.')
p = fp.create_group('PartType2')
p.create_dataset('ParticleIDs', data = next_ids, dtype = 'uint32')
p.create_dataset('Coordinates', data = disk['Coordinates'], dtype = 'float64')
p.create_dataset('Velocities', data = disk['Velocities'], dtype = 'float64')
p.create_dataset('Masses', data = disk['Masses'], dtype = 'float64')
running_idx += tot_Ndisk

if do_bulge:
    print('Adding bulge particles.')
    next_ids = new_ids[running_idx:running_idx + tot_Nbulge]

    p = fp.create_group('PartType3')
    p.create_dataset('ParticleIDs', data = next_ids, dtype = 'uint32')
    p.create_dataset('Coordinates', data = bulge['Coordinates'], dtype = 'float64')
    p.create_dataset('Velocities', data = bulge['Velocities'], dtype = 'float64')
    p.create_dataset('Masses', data = bulge['Masses'], dtype = 'float64')
    running_idx += tot_Nbulge

if do_stars:
    print('Adding star particles.')
    next_ids = new_ids[running_idx:running_idx + tot_Nstars]

    p = fp.create_group('PartType4')
    p.create_dataset('ParticleIDs', data = next_ids, dtype = 'uint32')
    p.create_dataset('Coordinates', data = stars['Coordinates'], dtype = 'float64')
    p.create_dataset('Velocities', data = stars['Velocities'], dtype = 'float64')
    p.create_dataset('Masses', data = stars['Masses'], dtype = 'float64')
    running_idx += tot_Nstars

next_ids = new_ids[running_idx:running_idx + tot_Nbh]

print('Adding black hole particles.')
p = fp.create_group('PartType5')
p.create_dataset('ParticleIDs', data = next_ids, dtype = 'uint32')
p.create_dataset('Coordinates', data = bh['Coordinates'], dtype = 'float64')
p.create_dataset('Velocities', data = bh['Velocities'], dtype = 'float64')
p.create_dataset('Masses', data = bh['Masses'], dtype = 'float64')
running_idx += tot_Nbh

print('running_idx: %d' % running_idx)
print('tot_part: %d' % tot_part)
assert running_idx == tot_part


header.attrs['NumPart_ThisFile'] = Npart;
header.attrs['NumPart_Total'] = Npart;
header.attrs['NumPart_Total_HighWord'] = 0 * Npart;

fp.close()
