#!/usr/bin/env python2
import numpy as np
import h5py
import argparse as ap
import os.path

parser = ap.ArgumentParser()
parser.add_argument('folder')
parser.add_argument('start')
parser.add_argument('end')
args = parser.parse_args()

start = int(args.start)
end = int(args.end)

# Asplund 2009
Zsun = 0.02
Zdef_bulge = 1.0e-3 * Zsun
Zdef_disk = 1.0e-2 * Zsun

part_types = {'Gas': 0,
              'Halo': 1,
              'Disk': 2,
              'Bulge': 3,
              'Stars': 4}

fields = {'Gas': 
                {'Coordinates':         'float64',
                 'Density':             'float64',
                 'InternalEnergy':      'float64',
                 'Masses':              'float64',
                 'Metallicity':         'float64',
                 'Velocities':          'float64',
                 'ParticleIDs':         'uint32',
                 'ElectronAbundance':   'float64',
                 'SmoothingLength':     'float64',
                 'StarFormationRate':   'float64'},
          'Halo': 
                {'Coordinates':     'float64',
                 'Masses':          'float64',
                 'ParticleIDs':     'uint32',
                 'Velocities':      'float64'},
          'Disk': 
                {'Coordinates':             'float64',
                 'Masses':                  'float64',
                 'ParticleIDs':             'uint32',
                 'StellarFormationTime':    'float64',
                 'Metallicity':             'float64',
                 'Velocities':              'float64'},
          'Bulge':
                 {'Coordinates':            'float64',
                 'Masses':                  'float64',
                 'ParticleIDs':             'uint32',
                 'StellarFormationTime':    'float64',
                 'Metallicity':             'float64',
                 'Velocities':              'float64'},
          'Stars':
                  {'Coordinates':           'float64',
                 'Masses':                  'float64',
                 'ParticleIDs':             'uint32',
                 'StellarFormationTime':    'float64',
                 'Metallicity':             'float64',
                 'Velocities':              'float64'}}
    

for snap in range(start, end + 1):
    snapstr = str(snap).zfill(3)
    
    infile = '%s/snapshot_%s.hdf5' % (args.folder, snapstr)
    outfile = '%s/restructured/snapshot_%s.hdf5' % (args.folder, snapstr)
    
    #if os.path.isfile(outfile):
    #    continue
    
    f_in = h5py.File(infile, 'r+')

    old_header = f_in['/Header']

    values = {}
    
    # Construct the values array similarly to the fields array
    for part_type in part_types:
        values.update({part_type: {}})
        
        for field in fields[part_type]:
            field_path = '/PartType%d/%s' % (part_types[part_type], field)
            
            try:
                field_value = np.asarray(f_in[field_path])
            except KeyError:
                field_value = np.zeros(1)
            
            values[part_type].update({field: field_value})

    Npart = np.zeros(len(part_types), dtype = 'uint32')
    
    # We need the total number of particles to assign new particle IDs
    for part_type in part_types:
        # Coordinates is common to every particle, can use to count particles
        Npart[part_types[part_type]] = len(values[part_type]['Coordinates'])
    
    Npart_total = np.sum(Npart)
    new_ids = np.arange(1, Npart_total + 1, 1, dtype = 'uint32')
    
    # Prepare new IC file
    f_out = h5py.File(outfile, 'w')
    
    header = f_out.create_group('Header')

    # Copy over the old header
    for key in old_header.attrs:
        header.attrs[key] = old_header.attrs[key]

        # Only allow 1 metal field for now
        if key == 'Flag_Metals':
            header.attrs[key] = 1

    # Keep track of the index
    running_idx = 0
    for part_type in part_types:
        id_stretch = Npart[part_types[part_type]]
        next_ids = new_ids[running_idx:running_idx + id_stretch]
        
        p = f_out.create_group('PartType' + str(part_types[part_type]))

        for field in fields[part_type]:
            if len(values[part_type][field]) > 1:
                if field == 'Metallicity':
                    values[part_type][field] = values[part_type][field][:, 0]

                if field == 'ParticleIDs':
                    values[part_type][field] = next_ids
            else:
                # Pre-existing stars need a default metallicity and age
                if field == 'Metallicity':
                    values[part_type][field] = np.zeros(len(next_ids))
                    
                    if part_type == 'Bulge':
                        values[part_type][field][:] = Zdef_bulge
                    else:
                        values[part_type][field][:] = Zdef_disk
                if field == 'StellarFormationTime':
                    values[part_type][field] = np.zeros(len(next_ids))
            
            p.create_dataset(field, data = values[part_type][field], 
                                 dtype = fields[part_type][field])
            
        running_idx += id_stretch
        
    f_out.close()
    f_in.close()
