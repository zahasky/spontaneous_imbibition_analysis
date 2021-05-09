# -*- coding: utf-8 -*-
"""
horizontal_1D_2phase.py
Revised from https://github.com/acroucher/PyTOUGH/wiki/Horizontal-1D-model
Date: 4/13/2021
@author: Christopher Zahasky 

Horizontal 1D model with temperatures and pressures specified at
each end. One end has two-phase conditions, the other single
phase. The script sets up the model, runs it, and produces a plot of
gas saturation along the model at the final time."""

from t2data import *
from t2listing import *
import os
import matplotlib.pyplot as plt
import time
import numpy as np

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Arial']})
plt.rcParams['font.size'] = 14
# hfont = {'fontname':'Arial'}
# # fontsize
# fs = 16

# --- set up the model ---------------------------------
# set inial condiitions
p0 = 103.0e3 # gas pressure
Sg0 = 10.999999 # initial gas saturation Sg+10
T0 = 21 # temperature
    
    
NB = np.arange(50, 1610, 50)
# NB = np.array((20, 1000))
front_location = np.zeros(NB.size)

for i in range(0, NB.size):
    # model dimensions set nearly identical to model
    core_area = 3.1415*(2.54/100)**2
    # length = 0.10+(0.1/125)
    nblks = NB[i]
    length = 0.10 + (0.2/nblks)
    dx = [length / nblks] * nblks
    # equaivalent 1d area as core
    dy = dz = [core_area**(1/2)]
    # setup grid
    geo = mulgrid().rectangular(dx, dy, dz)
    geo.write('geom.dat')
    
    # Create TOUGH2 input data file:
    dat = t2data()
    dat.title = 'Imbibition fit'
    dat.grid = t2grid().fromgeo(geo)
    # Set EOS conditions (MULTI block) 
    # Any configuration of this seems to produce the following error:
    # " TEMPERATURE = 0.000000E+00  OUT OF RANGE IN SAT " during model initialization
    # Possibly input writing formating error?
    # dat.multi = {'eos':'eos3', 'num_components':2, 'num_equations':3, 'num_inc':3, 'num_phases':2, 'num_secondary_parameters':6}
    
    # PARAM block
    dat.parameter.update(
        {'max_timesteps': 8000,
          'tstop': 120*60,
          'const_timestep': 120.,
          'print_interval': 100,
          'gravity': 9.81,
          'default_incons': [p0, Sg0, T0]})
    
    dat.start = True
    
    # Set MOPs (amount and type of printout):
    # MOP(1)  = 1 *** ALLOWS TO GENERATE A SHORT PRINTOUT FOR EACH NEWTON-RAPHSON ITERATION
    dat.parameter['option'][1] = 1
    # MOP(16) = 5 *** PERMITS TO CHOOSE TIME STEP SELECTION OPTION
    dat.parameter['option'][16] = 5
    ### NOTE ###
    # when using different EOS options you may need to explicitly set MOP 19 and 20
    # dat.parameter['option'][19] = 0 # default
    
    # Set relative permeability and capillarity functions based on core measurements:
    
    dat.relative_permeability = {'type': 2, 'parameters': [7.8, 0.0, 0., 0., 0.]}
    dat.capillarity = {'type': 7, 'parameters': [0.345, 0.03, 3.45e-5, 5.0e6, 0.76]}
    # dat.capillarity = {'type': 7, 'parameters': [0.290, 0.01, 3.79e-5, 5.0e6, 0.65]}
    # pretty good fit
    # dat.relative_permeability = {'type': 2, 'parameters': [7.4, 0.0, 0., 0., 0.]}
    # dat.capillarity = {'type': 7, 'parameters': [0.345, 0.03, 3.45e-5, 1.0e7, 0.76]}
    # dat.capillarity = {'type': 7, 'parameters': [0.47, 0.1, 2.6e-5, 5.0e6, 0.64]}
    
    # Set relative permeability and capillarity functions based on imbibition calculations measurements:
    # dat.relative_permeability = {'type': 2, 'parameters': [8.7, 0.0, 0., 0., 0.]}
    # dat.capillarity = {'type': 7, 'parameters': [0.35, 0.15, 6.6e-5, 5.0e6, 1]}
    
    # dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.0, 0., 0., 0.]}
    
    # original/working
    # dat.relative_permeability = {'type': 3, 'parameters': [0.3, 0.1, 0., 0., 0.]}
    # dat.capillarity = {'type': 1, 'parameters': [0., 0., 1., 0., 0.]}
    
    # Add a second rocktype ('dfalt' is the first rocktype)/
    # To query rocktypes use 'dat.grid.rocktype'
    r1 = rocktype('rock1', nad=0, conductivity = 0.0, porosity = 0.2, permeability = [2.3e-14]*3, specific_heat=0.0)
    dat.grid.add_rocktype(r1)
    r2 = rocktype('ends', nad=0, conductivity = 0.0, porosity = 0.99, permeability = [1.0e-17]*3, specific_heat=0.0)
    dat.grid.add_rocktype(r2)
    
    # Loop to assign permeability modifiers this should enable automatic J/Leverett scaling
    for blk in dat.grid.blocklist[0:]:
        # set permeability modifier for each cell
        blk.pmx = 1.0
        # set rocktype
        blk.rocktype = r1
        
    # add extra large cells to inlet and outlet 
    bvol = 500.0
    # inlet
    blk_in = dat.grid.blocklist[0]
    blk_in.rocktype = r2
    blk_in.volume = bvol
    # outlet
    blk_out = dat.grid.blocklist[-1]
    blk_out.rocktype = r2
    blk_out.volume = bvol

    inc = dat.grid.incons()
    # Loop through all cells and set initial conditions
    for blk in dat.grid.blocklist[1:]:
        inc[blk.name].variable = [p0, Sg0, T0]
    
    # inlet
    # inc['bdy01'].variable = [p0, 10.01, T0]
    inc[blk_in.name].variable = [p0, 10.01, T0]
    
    # outlet
    # inc['bdy02'].variable = [p0, Sg0, T0]
    
    inc.write('INCON')
    dat.write('imb_fit_eos3.dat')
    
    # --- run the model ------------------------------------
    # EOS 1: water (all phases), water with tracer
    # simulator="C:\\Users\\zahas\\Dropbox\\Research\\Simulation\\TOUGH2\\TOUGH2_Source\\T2_EOS1_Source_v2.1\\Executable_EOS1\\xt2_eos1"
    # EOS 3: water, gas
    # Laptop path
    # simulator="C:\\Users\\zahas\\Dropbox\\Research\\Simulation\\TOUGH2\\TOUGH2_Source\\T2_EOS3_Source_v2.1\\Executable_EOS3\\xt2_eos3"
    # Desktop path
    simulator="D:\\Dropbox\\Research\\Simulation\\TOUGH2\\TOUGH2_Source\\T2_EOS3_Source_v2.1\\Executable_EOS3\\xt2_eos3"
    
    startTime = time.time()
    dat.run(simulator=simulator, silent=False)
    executionTime = (time.time() - startTime)
    print('Model run time in seconds: ' + str(executionTime))

    # Plot saturation profile
    lst = t2listing('imb_fit_eos3.listing')
    
    lst.last()
    # omit boundary blocks from the plot results:
    x = [blk.centre[0] for blk in dat.grid.blocklist[:nblks]]
    # to see all output type "lst.element"
    # to loop through the lst data it may be useful to use lst.index=i. You can 
    # then call the time via 'lst.time' and the following can be used to extract 
    # the corresponding data at those times
    # gas saturation
    sg = lst.element['SG']
    # water saturation
    sw = lst.element['SL']
    # capillary pressure
    pcap = lst.element['PCAP']
    # pressure gas
    pressure = lst.element['P']
    
    # find imbibe front
    a = sw>0.1
    b = [i for i, a in enumerate(a) if a]
    # print('imbibition front = ' + str(x[b[-1]]))
    print('Model: ' + str(NB[i]))
    print('last time: ' + str(lst.time/60))
    front_location[i] = x[b[-1]]


fig0, ax00 =  plt.subplots(1, 1, figsize=(6, 4), dpi=200)
ax00.plot(NB, front_location*100, 'ko-')
ax00.set(xlabel='Number of grid cells', ylabel='Imbibition front location [cm]')
fig0.savefig("grid_ref_study.pdf", bbox_inches = 'tight')  

data = np.array([front_location*100, NB])
np.savetxt('grid_ref_study_1600cell.csv', 
            np.transpose(data) , delimiter=',', fmt='%.3e')

# fig0, ax00 =  plt.subplots(1, 1, figsize=(6, 4), dpi=200)
# ax00.plot(NB, (front_location*100)/front_location, 'ko-')
# ax00.set(xlabel='Number of grid cells', ylabel='Imbibition front location [cm]')

# --- post-process the output ---------------------------

# Plot characteristic curves
# Capillary pressure (assuming type 7, van Genucten)
pc_param = dat.capillarity['parameters']
Sw_pc = np.linspace((pc_param[1]+0.05), pc_param[4], num=100)
Se = (Sw_pc - pc_param[1])/(pc_param[4] - pc_param[1])
Pc_vg = (-1/pc_param[2])*(Se**(-1/pc_param[0])-1)**(1-pc_param[0])

# Plot relative permeability
if dat.relative_permeability['type'] == 3:
    print('Relative permeability are Corey Curves')
    rp_param = dat.relative_permeability['parameters']
    Sw_rp = np.linspace((rp_param[0]+0.05), (1- rp_param[0] - rp_param[1]), num=100)
    Se = (Sw_rp - rp_param[0])/(1 - rp_param[0] - rp_param[1])
    krl = Se**4
    krg = (1-Se)**2*(1-Se**2)
elif dat.relative_permeability['type'] == 2:
    print('Water relative permeability is a modified Corey Curve')
    rp_param = dat.relative_permeability['parameters']
    Sw_rp = np.linspace((pc_param[1]+0.05), pc_param[4], num=100)
    krl = Sw_rp**rp_param[0]
    krg = np.ones(Sw_rp.size)
    

# Print all plot times
print('model print times [min] = ' + str(lst.times/60))

# Plot characterisitic curves
fig0, (ax01, ax02) =  plt.subplots(1, 2, figsize=(10, 4), dpi=200)
# rel perm
ax01.plot(Sw_rp, krl, 'k-')
ax01.plot(Sw_rp, krg, 'r-')
ax01.plot(Sw_rp, Sw_rp**8.7, 'g-')
ax01.set(xlabel='Water Saturation [-]', ylabel='Relative Permeability [-]')
# capillary pressure
ax02.plot(Sw_pc, -Pc_vg/1000, 'k-')
ax02.set(xlabel='Water Saturation [-]', ylabel='Capillary Pressure [kPa]')


fig, (ax11, ax12) = plt.subplots(2,  sharex=True, figsize=(8, 6), dpi=200)
ax11.set(ylabel='Water saturation [-]', title = 'Grid blocks = ' + str(nblks), xlim = (0, 0.001))
ax12.set(xlabel='Distance from inlet [cm]', ylabel='Pressure [kPa]', xlim = (0, 0.001))

# Import and plot experimental data
exp_data = np.loadtxt('scaled_jan_si_ct_profile_data.csv', delimiter=',')
ax11.plot(exp_data[:,0], exp_data[:,1], 'k.')

# set colormap
colors = plt.cm.viridis(np.linspace(0, 1, len(lst.times)))

for i in range(0, len(lst.times)):
    lst.index=i
    # gas saturation
    sg = lst.element['SG']
    # water saturation
    sw = lst.element['SL']
    # capillary pressure
    pcap = lst.element['PCAP']
    # pressure gas
    pressure = lst.element['P']
    # Plot points on pc figure
    ax02.plot(sw, -pcap/1000, 'o', color=colors[i])
    
    # Plot water saturation and pressure on second figure
    if lst.time/60 > 50:
        # ax1.plot(x, 1- sg[:nblks], 'ob-')
        ax11.plot(np.array(x)/((lst.time)**(1/2)), sw[:nblks], '-', color=colors[i], label= '%.1f min' %(lst.time/60))
        # plot gas phase pressure
        ax12.plot(np.array(x)/((lst.time)**(1/2)), (pcap-pressure)/1000, 'o-', color=colors[i])
        # ax2.plot(x, pcap/1000, 'or-')
    
ax11.legend()
