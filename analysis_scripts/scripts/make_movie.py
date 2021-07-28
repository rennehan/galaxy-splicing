#from visualization import movie_maker as mm
from visualization import movie_maker_merger as mm
from visualization import image_maker as imaker
import sys
import numpy as np

smaster = '/home/rennehan/scratch/EddyDiffusion/MFM/TurbOff/TurbOff_Protocluster_Makegalaxy/data/'
snapshot_dir = smaster
#sdir='blah'
sdir=''
#output_dir='/mnt/home/chayward/ceph/SPT2349/images/v1/'
output_dir = '/home/rennehan/analysis/protocluster/images/'
cosmo = 0.
size = 150. # old size
#size = 100
type_to_make = 'star'

mm.movie_maker(xmax_of_box=size,snapdir=snapshot_dir,outputdir_master=output_dir,\
      show_gasstarxray=type_to_make,cosmological=cosmo,frames_per_gyr=800.,
      i_snap_min = 0, i_snap_max = 81,
      skip_bh=1,
      set_fixed_center=[0.,0.,0.],\
      theta_0=90.,phi_0=240.,delta_phi=0.,
      add_extension='_v2',\
      camera_opening_angle=60.,\
      pixels=512,\
      z_to_add=0.01,\
      use_old_extinction_routine=0)#,\
      #h_max=0.01)

#      frame_min=np.int(sys.argv[1]),frame_max=np.int(sys.argv[2]))
