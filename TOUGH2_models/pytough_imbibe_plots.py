"""
pytough_imbibe_plots.py
Revised from https://github.com/acroucher/PyTOUGH/wiki/Horizontal-1D-model
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
import matplotlib as mpl
from matplotlib import cm
import numpy as np
# import time
import scipy
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Arial']})
plt.rcParams['font.size'] = 16


# --- set up the model ---------------------------------
# load experiment data
ct_data = np.loadtxt('twoD_capillary_hetero\\ct_imbibe_data.csv', delimiter=',')
# load streamtube permeability data
stream_km2 = np.loadtxt('twoD_capillary_hetero\\streamtube_perm.csv')

 # Let's add that to the path of the current directory
output_filename = os.path.join('.', 'twoD_capillary_hetero\\krw8_2_1200cells\\horiz2D_cap_het.listing')

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

# Define grid for plotting purposes    
geo = mulgrid().rectangular(dx, dy, dz)
# Create TOUGH2 input data file:
dat = t2data()
dat.grid = t2grid().fromgeo(geo)

## EXTRACT GRID INFORMATION
# omit boundary blocks from the plot results:
x = [blk.centre[0] for blk in dat.grid.blocklist[:nblks]]
x.append(x[-1]+dx[0])
y = [blk.centre[1] for blk in dat.grid.blocklist[:nblks*zblks:nblks]]
y.append(y[-1]+dy[0])
z = [blk.centre[2] for blk in dat.grid.blocklist[:nblks*zblks:nblks]]
z.append(z[-1]-dz[0])
# convert to cm
xcm = np.array(x)*100
xcm[0] = 0
zcm = np.array(z)*100
# define grid for 2D maps
X, Y = np.meshgrid(xcm, zcm-zcm[0])
# Quiver grid
Xq, Yq = np.meshgrid(xcm[:-1], zcm[:-1])
# water density (from output)
rho_w = 998.1

## Permeability with flow overlayn
perm_array = np.transpose(np.tile(stream_km2, (nblks, 1)))
array_data = perm_array/9.869233E-16
S = np.array(array_data)
perm_map = np.reshape(S, (len(zcm)-1, len(xcm)-1))

## LOAD DATA
# Plot saturation profile
lst = t2listing(output_filename)
print('Last model print time [min] = ' + str(lst.times/60))

# CT time index
# plot_indices = [12, 20,34, 52]
plot_indices = [8, 18,32, 52]
PV = [0.16, 0.25, 0.40, 0.57]

fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(15, 12), dpi=500)
fs = plt.rcParams['font.size']
axs[0, 0].set_title("Water Saturation")
axs[0, 1].set_title("Flow rate [ml/min]")
axs[3, 1].set_xlabel("Distance from inlet [cm]", fontsize=fs+3)
axs[3, 0].set_xlabel("Distance from inlet [cm]", fontsize=fs+3)
# fig.subplots_adjust(right=1.1)
fig.subplots_adjust(wspace=0.3)

# counter
n=0

# plot_indices = [12]

for i in plot_indices:
    lst.index = i
    
    # water saturation
    axs[n, 0].set_ylabel('PV:' + str(PV[n]), fontsize=fs+3)
    axs[n, 0].set_yticklabels([])
    axs[n, 1].set_yticklabels([])
    
    # extract water saturation
    sw = lst.element['SL']
    S = np.array(sw)
    map_data = np.reshape(S, (len(zcm)-1, len(xcm)-1))
    # Use 'pcolor' function to plot 2d map of water saturation
    im1 = axs[n, 0].pcolor(X, Y, map_data, vmin = 0, vmax = 0.75, cmap='Blues', snap = True)
    # axs[n, 0].clim(0, 0.75) 
    axs[n, 0].set_aspect('equal')  
    
    ##### Show fluxes ####
    # Note that KDATA has to be set to 2. In pyTouhg this is done with: dat.parameter.update({'print_level':2, ...)
    # To printout this full table use the command:
    # lst.connection
    
    # Liquid velocity
    # v = lst.connection['VEL(LIQ.)']*(100*60) # convert to cm/min
    # Liquid flow rate
    v = lst.connection['FLO(LIQ.)']/rho_w*(100**3)*60
    
    # flux_matrix(geo) this is a sparse matrix used to convert the connection flow values to block average fluxes
    a = dat.grid.flux_matrix(geo)
    # plt.spy(a)
    grid = t2grid().fromgeo(geo)
    grid.calculate_connection_centres(geo)
    blkflow = a * v
    blkflow3 = blkflow.reshape((geo.num_underground_blocks, 3))
    # flow along core
    grid_area_x = dz[0]*dy[0]
    # blkflowx = blkflow3[:,0].reshape(zblks, nblks) # in cm/min
    blkflowx = blkflow3[:,0].reshape(zblks, nblks)*grid_area_x # in mL/min
    # plot_2d(blkflowx,  xcm, zcm, 'Block average horizontal flow [mL/min]', 'Greens')
    # flow side-to-side
    # blkflowy = blkflow3[:,1].reshape(zblks, nblks)
    # plot_2d(blkflowy,  xcm, zcm, 'Flow', 'PiYG')
    # flow up down (positive flow is down)
    grid_area_z = dx[0]*dy[0]
    blkflowz = blkflow3[:,2].reshape(zblks, nblks)*grid_area_z # in mL/min
    # blkflowz = blkflow3[:,2].reshape(zblks, nblks) # in cm/min
    # plot_2d(blkflowz,  xcm, zcm, 'Block average vertical flow [mL/min]', 'PRGn')
    # clim = np.max(np.abs(blkflowz))*0.4
    # plt.clim(-clim, clim) 
    
    krw_map = perm_map*(map_data**7.8)
    # im2 = axs[n, 1].pcolor(X, Y, krw_map, vmin = 0, vmax = 3, cmap='gray_r', alpha=1, snap = True)
    im2 = axs[n, 1].pcolor(X, Y, perm_map, vmin = 10, vmax = 32, cmap='gray_r', alpha=1, snap = True)
    # axs[n, 1].clim(10, 32) 
    axs[n, 1].set_aspect('equal')  
    
    # Use 'pcolor' function to plot 2d map of concentration
    # Note that we are flipping map_data and the yaxis to so that y increases downward
    gaps = 30
    stop = -30
    clim = 4E-6*30
    im3 = axs[n, 1].quiver(Xq[:,:stop:gaps], Yq[:,:stop:gaps], blkflowx[:,:stop:gaps], 
               blkflowz[:,:stop:gaps]*30*10, blkflowz[:,:stop:gaps]*30, 
               clim=[-clim, clim], cmap='RdYlGn_r', width=0.005)
    # update counter
    n +=1

# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

# divider = make_axes_locatable(axs[0,0])
# cax = divider.append_axes('right', size='5%', pad=0.05)
fraction = 0.05
aspect = 40
fig.colorbar(im1, ax=axs[0:4, 0], fraction=fraction, pad=0.05, aspect=aspect)

fig.colorbar(im2, ax=axs[0:4, 1], location='left', fraction=fraction, pad=0.001, aspect=aspect, label='milliDarcy')
fig.colorbar(im3, ax=axs[0:4, 1], fraction=fraction, pad=0.05, aspect=aspect, format='%.0e')

fig.savefig("model_maps.pdf")   

##### Flow rate analysis by layer
# lst.index = 34
# v = lst.connection['FLO(LIQ.)']/rho_w*(100**3)*60

# # flux_matrix(geo) this is a sparse matrix used to convert the connection flow values to block average fluxes
# a = dat.grid.flux_matrix(geo)
# # plt.spy(a)
# grid = t2grid().fromgeo(geo)
# grid.calculate_connection_centres(geo)
# blkflow = a * v
# blkflow3 = blkflow.reshape((geo.num_underground_blocks, 3))
# # flow along core
# # blkflowx = blkflow3[:,0].reshape(zblks, nblks) # in cm/min
# blkflowx = blkflow3[:,0].reshape(zblks, nblks)*grid_area_x # in mL/min
# blkflowz = blkflow3[:,2].reshape(zblks, nblks)*grid_area_z # in mL/min

# # color=cm.gray_r(np.linspace(0,1,20))
# norm = mpl.colors.Normalize(vmin=10, vmax=32)
# cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.gray_r)
# cmap.set_array([])

# fig = plt.figure(figsize=(10, 4), dpi=200)
# for i in range(0,20):
#     plt.plot(xcm[:-1], blkflowx[i,:], c=cmap.to_rgba(perm_map[i,1]))

# fig.colorbar(cmap, label='milliDarcy')
# plt.ylabel('Axial-parallel flow rate [ml/min]')
# plt.xlabel('Distance from inlet [cm]'


##### Extract imbibition rate data and save it
# imbibe_rate = np.zeros((len(lst.times),1))
# Q = np.zeros((zblks, 1))

# for n in range(0, len(lst.times)):
#     lst.index = n
#     q = lst.connection['FLO(LIQ.)']
#     for i in range(0, zblks):
#         Q[i] = q[i*nblks*2-i]
#     # imbibition rate ml/min
#     imbibe_rate[n] = -Q.sum()/rho_w*(100**3)*60

    
# time_min = lst.times/60
# plt.figure(figsize=(5, 4), dpi=200)
# plt.plot(time_min, imbibe_rate, '-ok')
# plt.plot(ct_data[:,0], ct_data[:,1], '.r')
# plt.xscale('log')
# plt.yscale('log')
# plt.ylabel('Imibe rate')
# plt.xlabel('Time [min]')

# data = np.array([time_min, imbibe_rate.flatten()])

# np.savetxt('total_imbibe_rate_data_relperm_7_8_1200cells.csv', 
#            np.transpose(data) , delimiter=',', fmt='%.3e')


### PLOT nondimensioal comparison
# fig, (ax11) = plt.subplots(1,  sharex=True, figsize=(8, 4), dpi=200)
# ax11.set(ylabel='Water saturation [-]', xlabel='\lambda', 
#          title = 'Grid blocks: ' + str(nblks), xlim = (0, 0.001))

# # Import and plot experimental data
# exp_data = np.loadtxt('scaled_jan_si_ct_profile_data.csv', delimiter=',')
# ax11.plot(exp_data[:,0], exp_data[:,1], 'k.')

# # set colormap
# colors = plt.cm.viridis(np.linspace(0, 1, len(lst.times)))
# x = [blk.centre[0] for blk in dat.grid.blocklist[:nblks]]
# for i in range(0, len(lst.times)):
#     lst.index=i
#     # water saturation
#     sw = lst.element['SL']
#     S = np.array(sw)
#     SW = np.reshape(S, (len(zcm)-1, len(xcm)-1))

#     # Plot water saturation and pressure on second figure
#     if lst.time/60 > 50:
#         # ax1.plot(x, 1- sg[:nblks], 'ob-')
#         ax11.plot(np.array(x)/((lst.time)**(1/2)), SW[10,:], '-', color=colors[i], label= '%.1f min' %(lst.time/60))


# plt.savefig("overlay.pdf")   
# ax11.legend()
# Import core shape mask
# core_mask = np.loadtxt('core_mask.csv', delimiter=',')

dat.relative_permeability = {'type': 2, 'parameters': [7.8, 0.0, 0., 0., 0.]}
# New fitting based on 2D model results
dat.capillarity = {'type': 7, 'parameters': [0.345, 0.03, 3.45e-5, 5.0e6, 0.76]}
                
# Plot characteristic curves
# Capillary pressure (assuming type 7, van Genucten)
pc_param = dat.capillarity['parameters']
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
    # Sw_rp = np.linspace((pc_param[1]+0.05), pc_param[4], num=100)
    Sw_rp = np.linspace((pc_param[1]+0.05), 1, num=100)
    krl = Sw_rp**rp_param[0]
    krg = np.ones(Sw_rp.size)
    
# Capillary pressure (assuming type 7, van Genucten)
Sw_pc = np.linspace((pc_param[1]+0.05), pc_param[4], num=100)
Se = (Sw_pc - pc_param[1])/(pc_param[4] - pc_param[1])
Pc_vg = (-1/pc_param[2])*(Se**(-1/pc_param[0])-1)**(1-pc_param[0])
    
# Plot characterisitic curves
plt.rcParams['font.size'] = 15
fig0, (ax01, ax02) =  plt.subplots(1, 2, figsize=(10, 4), dpi=300)
fig0.subplots_adjust(wspace=0.35)
# fig0.subplots_adjust(top=0.9)
# rel perm
ax01.plot(Sw_rp, krl, 'k-', label='model water')
# ax01.plot(Sw_rp, krg, 'k--', label='gas')

swexp = [1.000, 0.984, 0.970, 0.956, 0.951, 0.936, 0.924]
krwexp = [1.00, 0.846, 0.750, 0.655, 0.600, 0.528, 0.477]
ax01.plot(swexp, krwexp, 'ko', alpha = 0.3, label='experimental water')

# ax01.plot(Sw_rp, krl, 'k-')
ax01.set(xlabel='Water Saturation [-]', ylabel='Relative Permeability [-]', xlim =(0,1))
ax01.legend()

cp_param = dat.capillarity['parameters']
core_average_perm = np.mean(stream_km2)
core_average_entry_pressure =  1/cp_param[2]

for r in range(0, len(stream_km2)):
    # r1 = rocktype(('rock'+str(r)), nad=0, conductivity = 0.0, porosity = 0.2, permeability = [stream_km2[r-1]]*3, specific_heat=0.0)
    # J-Leverett scaling
    p_entry = core_average_entry_pressure*(core_average_perm/stream_km2[r])**(1/2)
    Pc_vg_layer = (-p_entry)*(Se**(-1/pc_param[0])-1)**(1-pc_param[0])
    ax02.plot(Sw_pc, -Pc_vg_layer/1000, 'k-', alpha=0.1)
    
# capillary pressure
# ax02.plot(Sw_pc, -Pc_vg/1000, 'k-')
ax02.plot(Sw_pc, -Pc_vg/1000, 'r-')

ax02.set(xlabel='Water Saturation [-]', ylabel='Capillary Pressure [kPa]', ylim = (5,5000), xlim =(0,1))
ax02.set_yscale('log')

fig0.savefig("char_curves.pdf", bbox_inches = 'tight')  