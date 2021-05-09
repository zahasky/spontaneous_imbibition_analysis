# -*- coding: utf-8 -*-
"""
SI_model_generation_2D_2phase.py
Date: 4/13/2021
@author: Christopher Zahasky 

Horizontal 2D model with temperatures and pressures specified at
each end. One end has two-phase conditions, the other single
phase. The script sets up the model, runs it, and produces a plot of
gas saturation along the model at the final time."""

from t2data import *
from t2listing import *
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import time
import scipy

# --- set up the model ---------------------------------
# load experiment data
ct_data = np.loadtxt('ct_imbibe_data.csv', delimiter=',')
# load streamtube permeability data
stream_km2 = np.loadtxt('streamtube_perm.csv')

# core cross-sectional area (m^2)
core_area = 3.1415*(2.54/100)**2/1.39 # the 1.39 accounts for the area of the core that gets cropped off during image processing
length = 0.101
height = 5/100
nblks = 1200
dx = [length / nblks] * nblks
# equivalent area
dy = [core_area/height]
# zblks = round(height/dx[0])
zblks = 20
dz = [height / zblks] * zblks

# set inial condiitions
p0 = 103.0e3 # gas pressure
Sg0 = 10.999 # initial gas saturation Sg+10
T0 = 21 # temperature

# relative permeability exponent
krw_exp = [7.8]

# Directory name
directory_name = 'krw_7_8_1200cells'  
# Let's add that to the path of the current directory
workdir = os.path.join('.', directory_name)

geo = mulgrid().rectangular(dx, dy, dz)
geo.write('geom.dat')

# Create TOUGH2 input data file:
dat = t2data()
dat.title = 'Horizontal 2D cap het'
dat.grid = t2grid().fromgeo(geo)

#TIMES Block
out_times = list(ct_data[1:-4,0]*60) # seconds
# describe exactly when to print output, note that outtimes must be a list!
dat.output_times = {'num_times_specified':len(out_times), 'time': out_times} 

dat.parameter.update(
	{'max_timesteps': 9999,
	  'tstop': 264*60,
	  'const_timestep': 120.,
	  'print_interval': 9999,
	  'print_level':2,
	  'gravity': 9.81,
	  'relative_error':10e-5,
	  'default_incons': [p0, Sg0, T0]})

dat.start = True

# Set MOPs (amount and type of printout):
# MOP(1)  = 1 *** ALLOWS TO GENERATE A SHORT PRINTOUT FOR EACH NEWTON-RAPHSON ITERATION
dat.parameter['option'][1] = 1
# MOP(16) = 5 *** PERMITS TO CHOOSE TIME STEP SELECTION OPTION
dat.parameter['option'][16] = 5
# Set relative permeability and capillarity functions based on core measurements: 
dat.relative_permeability = {'type': 2, 'parameters': [krw_exp, 0.0, 0., 0., 0.]}
dat.capillarity = {'type': 7, 'parameters': [0.345, 0.03, 3.45e-5, 5.0e6, 0.76]}

# Add a second rocktype ('dfalt' is the first rocktype)/
# To query rocktypes use 'dat.grid.rocktype'
dum1 = rocktype('ends', nad=0, conductivity = 0.0, porosity = 0.99, permeability = [1.0e-17]*3)
dat.grid.add_rocktype(dum1)
# Dummy slices with no capillary pressure, all phases fully mobile
dum2 = rocktype('ends2', nad=2, conductivity = 0.0, porosity = 0.99, permeability = [1.0e-17]*3, \
				specific_heat=0.0, relative_permeability = {'type': 5, 'parameters': [0.]}, \
					capillarity = {'type': 0, 'parameters': [0.]})
dat.grid.add_rocktype(dum2)

d={}
core_average_perm = np.mean(stream_km2)
pc_param = dat.capillarity['parameters']
rp_params = dat.relative_permeability['parameters']
core_average_entry_pressure =  1/pc_param[2]
for r in range(0, len(stream_km2)):
	# r1 = rocktype(('rock'+str(r)), nad=0, conductivity = 0.0, porosity = 0.2, permeability = [stream_km2[r-1]]*3, specific_heat=0.0)
	# J-Leverett scaling
	p_entry = core_average_entry_pressure*(core_average_perm/stream_km2[r])**(1/2)
	d["string{0}".format(r)] = rocktype(('roc'+str(r)), nad=2, conductivity = 0.0, \
		porosity = 0.2, permeability = [stream_km2[r]]*3, specific_heat=0.0, \
		relative_permeability = {'type': 2, 'parameters': [rp_params[0], 0.0, 0., 0., 0.]}, \
		capillarity = {'type': 7, 'parameters': [pc_param[0], pc_param[1], 1/p_entry, pc_param[3], pc_param[4]]})
		# original PC    
		# capillarity = {'type': 7, 'parameters': [0.345, 0.03, 1/p_entry, 5.0e6, 0.76]})
		
	dat.grid.add_rocktype(d["string{0}".format(r)])
	
# Loop to assign permeability modifiers this should enable automatic J/Leverett scaling
m = np.random.normal(1, 0.3, (zblks, nblks))
m = m.flatten()
# remove negative values
m[m<0.01]=1
m[m>2.0]=2.0
i = 0
# define grid edges
z_array = -np.cumsum(dz)

for blk in dat.grid.blocklist[0:]:
	# set permeability modifier for each cell
	blk.pmx = m[i]
	i += 1
	# set rocktype, find index that corresponds to streamtube
	rock_index =  next(x[0] for x in enumerate(z_array) if x[1] < blk.centre[2])
	blk.rocktype = d["string{0}".format(rock_index)]
	
# add extra large cells to inlet and outlet 
bvol = 500.0 # should be at least 10x the core volume

inc = dat.grid.incons()
# inc.variable =  initial_cond
for blk in dat.grid.blocklist[:]:
	# explicitly set initial conditions (not necessary if defined in dat.parameter)
	inc[blk.name].variable = [p0, Sg0, T0]

	# Set inlet conditions
	if blk.centre[0] < dx[0]:
		# print(blk.centre)
		inc[blk.name].variable = [p0, 10.01, T0]
		# inc[blk.name].variable = [p0, Sg0, T0]
		# set rocktype to 'ends'
		blk.rocktype = dum2
		# set volume large
		blk.volume = bvol
	elif (blk.centre[0] > dx[0]) and (blk.centre[0] < dx[0]*2):
		inc[blk.name].variable = [p0, 10.35, T0]
	# set outlet conditions
	elif blk.centre[0] > length - dx[0]:
		# print(blk.centre)
		# set rocktype to 'ends'
		blk.rocktype = dum1
		# set volume large
		blk.volume = bvol

# write the incon file
if os.path.isdir(workdir) is False:
	os.mkdir(workdir)
incon_filename = (workdir+'\\INCON')
inc.write(incon_filename)
infilename = (workdir+'\\horiz2D_cap_het.dat')
dat.write(infilename)
output_filename = (infilename[:-3]+'listing')

# --- run the model ------------------------------------
# EOS 1: water (all phases), water with tracer
# Desktop path
simulator="D:\\Dropbox\\Research\\Simulation\\TOUGH2\\TOUGH2_Source\\T2_EOS3_Source_v2.1\\cr2_eos3.exe"

startTime = time.time()
dat.run(save_filename = infilename, incon_filename = incon_filename, simulator=simulator, silent=False, output_filename=output_filename)
executionTime = (time.time() - startTime)
print('Model run time in seconds: ' + str(executionTime))

# --- post-process the output ---------------------------
# Plot saturation profile
lst = t2listing(output_filename)
print('model print times [min] = ' + str(lst.times/60))