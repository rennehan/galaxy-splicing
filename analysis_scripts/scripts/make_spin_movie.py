from visualization import movie_maker_merger as mm
from visualization import image_maker as imaker
import sys
import numpy as np

smaster = '/home/rennehan/scratch/EddyDiffusion/MFM/TurbOff/TurbOff_Protocluster_HighRes/data/'
#smaster = '/home/rennehan/scratch/EddyDiffusion/MFM/TurbOff/TurbOff_Protocluster_Makegalaxy/data/'
snapshot_dir = smaster
output_dir = '/home/rennehan/analysis/protocluster/images/spin/'
cosmo = 0.
size = 150.
type_to_make = 'star'
snum = 31
frames_per_gyr = 800

frame_start = 0
frame_end = 240

# Start the camera opening angle here
fov_0 = 45.0
fov_final = 60.0
dfov = float((fov_final - fov_0) / frame_end)

# Angular quantities
phi_0 = 240.0
theta_0 = 90.0
rot_angle = 360.0

nframe = frame_end - frame_start + 1
frame_grid = np.arange(0, nframe)

# Angular quantities
dphi = float(rot_angle / frame_end)

phi_grid = dphi * frame_grid
fov_grid = dfov * frame_grid

frame_list = np.arange(frame_start, frame_end + 1, 1)
    
for frame in frame_list:
    phi = phi_0 + phi_grid[frame]
    fov = fov_final #fov_0 + fov_grid[frame] 

    mm.movie_maker(xmax_of_box = size,
                    snapdir = snapshot_dir,
                    outputdir_master = output_dir,
                    show_gasstarxray = type_to_make,
                    cosmological = cosmo,
                    frames_per_gyr = frames_per_gyr,
                    i_snap_min = snum, 
                    i_snap_max = snum,
                    skip_bh = 1,
                    set_fixed_center = [0.,0.,0.],
                    theta_0 = theta_0, 
                    phi_0 = phi, 
                    delta_phi = 0.,
                    add_extension = '_%s' % str(frame).zfill(5),
                    camera_opening_angle = fov,
                    pixels = 512,
                    z_to_add = 0.01,
                    use_old_extinction_routine = 0,
                    frame_min = snum * 2,
                    frame_max = snum * 2)

