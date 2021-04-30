import yt

bbox_lim = 1e5
bbox_vec = [-bbox_lim, bbox_lim]
bbox = [bbox_vec, bbox_vec, bbox_vec]

ds = yt.load('spliced_1e6.hdf5', bounding_box = bbox)

pp = yt.ProjectionPlot(ds, 'z', ('gas', 'density'))
pp.set_width((135, 'kpc'))
pp.set_unit('density', 'Msun/kpc**2')
pp.save('1e6_proj.png')

