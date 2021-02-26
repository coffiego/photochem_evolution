#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 14:15:57 2020

@author: shungokoyama
"""



import h5py
import matplotlib.pyplot as plt
import numpy as np

species_color ={
   'HOCO': "#dead91",
   'H2O2' : "#bebddb",
   'HO2': "#9e9ac8",
   'O3': "#83d0c1",
   'OH' :"#b593e2",
   'H2O' :"#1f78b4",
   'O1D': "#269e56",
   'CO2pl': "#eed2d4",
   'H' :"#ef3b2c",
   'H2' : "#fc9aa1",
   'O2': "#41ae76",
   'O': "#00702d",
   'CO': "#fd8d3c",
   'Ar': "#808080",
   'N2': "#cccccc",
   'CO2' : "#fdd0a2"
   }

def flux_pressure_comp(file_for_Hesc,file_for_species,Oescflux_after=2.63e8):
    #file = file_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    fluxvals_mat = []
    h5file = h5py.File(file_for_Hesc,"r")
    time_array = h5file["fluxes/times"].value #get time array here
    #time_array = time_array/3.15e7
    
    
    ##ax0 H escape flux
    fluxvals_mat = h5file["fluxes/fluxvals"].value
    H2O2dep=h5file["depositions/H2O2"].value
    O3dep=h5file["depositions/O3"].value
    #HO2dep=h5file["depositions/HO2"].value
    #Ototalflux_redox=H2O2dep + O3dep*3 + HO2dep*1.5 + Oescflux_after
    #Ototalflux = h5file["outfluxes/O"].value
    CO2esc = h5file["outfluxes/CO2esc"].value
    CO2outgassing = h5file["outfluxes/CO2out"].value
    Hesc = h5file["outfluxes/Hesc"].value
    #H2outgassing = h5file["outfluxes/H2out"].value
    #Oreturn = h5file["outfluxes/Oreturn"].value
    Othermal = h5file["outfluxes/Othermal"].value

    Oesc = np.zeros(len(time_array)+1) #change later
    Oesc[:]=Oescflux_after
    
    
    array_species_overtime_CO2 = species_pressure_overtime(file_for_species,"CO2",len(time_array-1))
    array_species_overtime_O2 = species_pressure_overtime(file_for_species,"O2",len(time_array-1))
    array_species_overtime_CO = species_pressure_overtime(file_for_species,"CO",len(time_array-1))
    array_species_overtime_H2 = species_pressure_overtime(file_for_species,"H2",len(time_array-1))
    output = [time_array, fluxvals_mat[0:-1], O3dep[0:-1],H2O2dep[0:-1], CO2esc[0:-1],CO2outgassing[0:-1],Hesc[0:-1],Othermal[0:-1],
            array_species_overtime_CO2,array_species_overtime_O2,array_species_overtime_CO,array_species_overtime_H2]

    return output
"""
lowest_13euv = flux_pressure_comp("Tinf480_lowestCO2_euv13_1382K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_lowestCO2_euv13_1382K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_135euv = flux_pressure_comp("Tinf480_lowestCO2_euv13.5_1526K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_lowestCO2_euv13.5_1526K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_14euv = flux_pressure_comp("Tinf480_lowestCO2_euv14_1670K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_lowestCO2_euv14_1670K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_145euv = flux_pressure_comp("Tinf480_lowestCO2_euv14.5_1826K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_lowestCO2_euv14.5_1826K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_15euv = flux_pressure_comp("Tinf480_lowestCO2_euv15_2026K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_lowestCO2_euv15_2026K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_155euv = flux_pressure_comp("Tinf480_lowestCO2_euv15.5_2252K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_lowestCO2_euv15.5_2252K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
"""

lowest_10euv = flux_pressure_comp("Tinf480_lowestCO2_euv10/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_lowestCO2_euv10/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_11euv = flux_pressure_comp("Tinf480_lowestCO2_euv11/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_lowestCO2_euv11/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_12euv = flux_pressure_comp("Tinf480_lowestCO2_euv12/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_lowestCO2_euv12/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_125euv = flux_pressure_comp("Tinf480_lowestCO2_euv12.5_1254K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_lowestCO2_euv12.5_1254K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_13euv = flux_pressure_comp("Tinf480_lowestCO2_euv13_1382K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_lowestCO2_euv13_1382K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_135euv = flux_pressure_comp("Tinf480_lowestCO2_euv13.5_1526K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_lowestCO2_euv13.5_1526K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_14euv = flux_pressure_comp("Tinf480_lowestCO2_euv14_1670K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_lowestCO2_euv14_1670K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_145euv = flux_pressure_comp("Tinf480_lowestCO2_euv14.5_1826K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_lowestCO2_euv14.5_1826K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_15euv = flux_pressure_comp("Tinf480_lowestCO2_euv15_2026K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_lowestCO2_euv15_2026K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
lowest_155euv = flux_pressure_comp("Tinf480_lowestCO2_euv15.5_2252K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_lowestCO2_euv15.5_2252K/CO2_1mbar_Tsurf_210_Tian_BD.h5")

lowest = [lowest_155euv,lowest_15euv,lowest_145euv,lowest_14euv,lowest_135euv,lowest_13euv]
#lowest = [lowest_155euv,lowest_15euv,lowest_145euv,lowest_14euv,lowest_135euv,lowest_13euv,lowest_125euv,lowest_12euv,lowest_11euv,lowest_10euv]


"""
middle_13euv = flux_pressure_comp("Tinf480_middleCO2_euv13_1382K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_middleCO2_euv13_1382K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_135euv = flux_pressure_comp("Tinf480_middleCO2_euv13.5_1526K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_middleCO2_euv13.5_1526K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_14euv = flux_pressure_comp("Tinf480_middleCO2_euv14_1670K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_middleCO2_euv14_1670K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_145euv = flux_pressure_comp("Tinf480_middleCO2_euv14.5_1826K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_middleCO2_euv14.5_1826K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_15euv = flux_pressure_comp("Tinf480_middleCO2_euv15_2026K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_middleCO2_euv15_2026K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_155euv = flux_pressure_comp("Tinf480_middleCO2_euv15.5_2252K_flexo_1mbar/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_middleCO2_euv15.5_2252K_flexo_1mbar/CO2_1mbar_Tsurf_210_Tian_BD.h5")
"""
middle_10euv = flux_pressure_comp("Tinf480_middleCO2_euv10/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_middleCO2_euv10/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_11euv = flux_pressure_comp("Tinf480_middleCO2_euv11/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_middleCO2_euv11/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_12euv = flux_pressure_comp("Tinf480_middleCO2_euv12/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_middleCO2_euv12/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_125euv = flux_pressure_comp("Tinf480_middleCO2_euv12.5_1254K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_middleCO2_euv12.5_1254K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_13euv = flux_pressure_comp("Tinf480_middleCO2_euv13_1382K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                   "Tinf480_middleCO2_euv13_1382K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_135euv = flux_pressure_comp("Tinf480_middleCO2_euv13.5_1526K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_middleCO2_euv13.5_1526K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_14euv = flux_pressure_comp("Tinf480_middleCO2_euv14_1670K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_middleCO2_euv14_1670K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_145euv = flux_pressure_comp("Tinf480_middleCO2_euv14.5_1826K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_middleCO2_euv14.5_1826K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_15euv = flux_pressure_comp("Tinf480_middleCO2_euv15_2026K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                  "Tinf480_middleCO2_euv15_2026K/CO2_1mbar_Tsurf_210_Tian_BD.h5")
middle_155euv = flux_pressure_comp("Tinf480_middleCO2_euv15.5_2252K/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5",
                                    "Tinf480_middleCO2_euv15.5_2252K/CO2_1mbar_Tsurf_210_Tian_BD.h5")

middle = [middle_155euv,middle_15euv,middle_145euv,middle_14euv,middle_135euv,middle_13euv]
#middle = [middle_155euv,middle_15euv,middle_145euv,middle_14euv,middle_135euv,middle_13euv,middle_125euv,middle_12euv,middle_11euv,middle_10euv]

#non-thermal escape
#euvlist = np.array([15.5,15.5, 15, 14.5, 14,13.5,13,12.5,12, 11,10])
euvlist = np.array([15.5,15.5, 15, 14.5, 14,13.5,13])
O_nonthermal = (1+0.16*euvlist)*1.e8 #3<EUV<20
C_nonthermal = (-2.7+0.43*euvlist)*1.e8 #10<EUV<20

#### Lowest ###################################################################
ts=[]
flux_val=[]
O3dep=[]
H2O2dep=[]
Cesc=[]
CO2out = []
Hesc =[]
Oesc=[]
CO2_num =[]
O2_num =[]
CO_num =[]
H2_num = []

all_array = [ts,flux_val,O3dep,H2O2dep,Cesc,CO2out,Hesc,Oesc,CO2_num,O2_num,CO_num,H2_num]

#input all arrays except time array
for euv in lowest:
    for i,ar in enumerate(all_array):
        ar.extend(euv[i])

#time array is needed to be rearanged in continueous time series.
ts=[0.]
for euv in lowest:
    ts.extend(euv[0]+[ts[-1]]*len(euv[0]))
ts.pop(0)


#get every 3994 grids
from operator import itemgetter
scatter_index=[]
for i in range(len(lowest)+1):
    scatter_index.append(3994*i-1)
scatter_index[0]=2000


scatter_array = [ts,flux_val,O3dep,H2O2dep,Cesc,CO2out,Hesc,Oesc,CO2_num,O2_num,CO_num,H2_num]
for i,ar in enumerate(scatter_array):
    scatter_array[i] = itemgetter(*scatter_index)(ar)

"""plot all points

ts = np.array(ts)
ts = ts/3.15e13
fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
fig.subplots_adjust(hspace=0)

axs[0].plot(ts, Hesc, color='dodgerblue',label='H escape')
axs[0].plot(ts,Oesc, color='g',label='O thermal escape flux')
axs[0].plot(ts,Cesc,label='C escape flux',color='orange')
axs[0].plot(ts,CO2out,label='CO2 outgassing flux',color='orange',linestyle='-.')
axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
axs[0].set_ylim((5e5,2e12))
axs[0].set_yscale("log")
axs[0].legend(prop={'size':8.0}, loc= "upper right")

axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
axs[1].set_xlabel('Time [Myr]')
axs[1].set_yscale("log")
#axs[1].set_xscale("log")
axs[1].set_ylim((1e-5,5e3))
axs[1].plot(ts,CO2_num,label="CO2",color=species_color["CO2"])
axs[1].plot(ts,O2_num,label="O2",color=species_color["O2"])
axs[1].plot(ts,CO_num,label="CO",color=species_color["CO"])
axs[1].plot(ts,H2_num,label="H2",color=species_color["H2"])
axs[1].legend(prop={'size':8.0},loc='upper left')
"""

#"""scatter plot
[ts,flux_val,O3dep,H2O2dep,Cesc,CO2out,Hesc,Oesc,CO2_num,O2_num,CO_num,H2_num] = scatter_array
ts = np.array(ts)
ts = ts/3.15e13

"""
#### scatter ##################################################################
fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
fig.subplots_adjust(hspace=0)
axs[0].scatter(ts, Hesc, color='dodgerblue',label='H escape')
axs[0].scatter(ts,Oesc, color='g',label='O thermal escape flux')
axs[0].scatter(ts,Cesc,label='C escape flux',color='orange')
axs[0].scatter(ts,CO2out,label='CO2 outgassing flux',color='orange',marker='x')
axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
axs[0].set_ylim((5e5,2e12))
axs[0].set_yscale("log")
axs[0].legend(prop={'size':8.0}, loc= "upper right")

axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
axs[1].set_xlabel('Time [Myr]')
axs[1].set_yscale("log")
#axs[1].set_xscale("log")
axs[1].set_ylim((1e-6,5e3))
axs[1].scatter(ts,CO2_num,label="CO2",color=species_color["CO2"])
axs[1].scatter(ts,O2_num,label="O2",color=species_color["O2"])
axs[1].scatter(ts,CO_num,label="CO",color=species_color["CO"])
axs[1].scatter(ts,H2_num,label="H2",color=species_color["H2"])
axs[1].legend(prop={'size':8.0},loc='upper left')
###############################################################################
"""

#O escape adjustment
#Oesc include only thermal escape, so non-thermal escape is added
Oesc = Oesc + O_nonthermal

#### marker style line ########################################################
fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
fig.subplots_adjust(hspace=0)
axs[0].plot(ts, Hesc, color='dodgerblue',label='H escape',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,Oesc, color='g',label='O escape',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,O3dep, color='#83d0c1',label='O3 deposition',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,H2O2dep, color='#bebddb',label='H2O2 deposition',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,Cesc,label='C escape',color='orange',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,CO2out,label='CO2 outgassing',color=species_color["CO2"],marker='x',fillstyle="none",linestyle="--")
axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
axs[0].set_ylim((5e5,2e12))
axs[0].set_yscale("log")
axs[0].legend(prop={'size':6.0}, loc= "upper right",ncol=2)

axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
axs[1].set_xlabel('Time [Myr]')
axs[1].set_yscale("log")
#axs[1].set_xscale("log")
axs[1].set_ylim((1e-7,5e3))
axs[1].plot(ts,CO2_num,label="CO2",color=species_color["CO2"],marker="o",fillstyle="none",linestyle="--")
axs[1].plot(ts,O2_num,label="O2",color='g',marker="o",fillstyle="none",linestyle="--")
axs[1].plot(ts,CO_num,label="CO",color=species_color["CO"],marker="o",fillstyle="none",linestyle="--")
#axs[1].plot(ts,H2_num,label="H2",color=species_color["H2"],marker="o",fillstyle="none",linestyle="--")
axs[1].legend(prop={'size':6.0})
###############################################################################
#"""

### redox balance #############################################################
C_O = []
H_O = []
C_H_O = []
for i in range(len(Hesc)):
    Closs = 4*Cesc[i]
    Hloss = Hesc[i]
    Oloss = 2*(Oesc[i]+3*O3dep[i]+H2O2dep[i])
    C_O.append(Closs/Oloss)
    H_O.append(Hloss/Oloss)
    C_H_O.append((Hloss+Closs)/Oloss)
    print(C_H_O[-1])

plt.figure(dpi=150)
plt.plot(ts, H_O, color='dodgerblue',label=r"$\Phi_H / 2\Phi_O$",marker="x",fillstyle="none",linestyle="--",markersize=10)
plt.plot(ts, C_O, color='orange',label=r"$2\Phi_C / \Phi_O$",marker="^",fillstyle="none",linestyle="--",markersize=10)
plt.plot(ts, C_H_O, color='k',label=r"$(\Phi_H+4\Phi_C) / 2\Phi_O$",marker="o",fillstyle="none",linestyle="--")
plt.xlabel('Time [Myr]')
plt.ylabel('Redox Balance')
plt.ylim((-0.1,1.5))
plt.legend(prop={'size':7.0})


#### middle ###################################################################

ts=[]
flux_val=[]
O3dep=[]
H2O2dep=[]
Cesc=[]
CO2out = []
Hesc =[]
Oesc=[]
CO2_num =[]
O2_num =[]
CO_num =[]
H2_num = []

all_array = [ts,flux_val,O3dep,H2O2dep,Cesc,CO2out,Hesc,Oesc,CO2_num,O2_num,CO_num,H2_num]

#all arrays except time array
for euv in middle:
    for i,ar in enumerate(all_array):
        ar.extend(euv[i])
#time array
ts=[0.]
for euv in middle:
    ts.extend(euv[0]+[ts[-1]]*len(euv[0]))
ts.pop(0)

#get every 3994 grids
from operator import itemgetter
scatter_index=[]
for i in range(len(middle)+1):
    scatter_index.append(3994*i-1)
scatter_index[0]=2000


scatter_array = [ts,flux_val,O3dep,H2O2dep,Cesc,CO2out,Hesc,Oesc,CO2_num,O2_num,CO_num,H2_num]
for i,ar in enumerate(scatter_array):
    scatter_array[i] = itemgetter(*scatter_index)(ar)

""" plot all points
ts = np.array(ts)
ts = ts/3.15e13
fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
fig.subplots_adjust(hspace=0)

axs[0].plot(ts, Hesc, color='dodgerblue',label='H escape')
axs[0].plot(ts,Oesc, color='g',label='O thermal escape flux')
axs[0].plot(ts,Cesc,label='C escape flux',color='orange')
axs[0].plot(ts,CO2out,label='CO2 outgassing flux',color='orange',linestyle='-.')
axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
axs[0].set_ylim((5e5,2e12))
axs[0].set_yscale("log")
axs[0].legend(prop={'size':8.0}, loc= "upper right")

axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
axs[1].set_xlabel('Time [Myr]')
axs[1].set_yscale("log")
#axs[1].set_xscale("log")
axs[1].set_ylim((1e-5,5e3))
axs[1].plot(ts,CO2_num,label="CO2",color=species_color["CO2"])
axs[1].plot(ts,O2_num,label="O2",color=species_color["O2"])
axs[1].plot(ts,CO_num,label="CO",color=species_color["CO"])
axs[1].plot(ts,H2_num,label="H2",color=species_color["H2"])
axs[1].legend(prop={'size':8.0},loc='upper left')
"""

#"""scatter plot
[ts,flux_val,O3dep,H2O2dep,Cesc,CO2out,Hesc,Oesc,CO2_num,O2_num,CO_num,H2_num] = scatter_array
ts = np.array(ts)
ts = ts/3.15e13

"""
#### scatter ##################################################################
fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
fig.subplots_adjust(hspace=0)
axs[0].scatter(ts, Hesc, color='dodgerblue',label='H escape')
axs[0].scatter(ts,Oesc, color='g',label='O thermal escape')
axs[0].scatter(ts,Cesc,label='C escape',color='orange')
axs[0].scatter(ts,CO2out,label='CO2 outgassing',color='orange',marker='x')
axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
axs[0].set_ylim((5e5,2e12))
axs[0].set_yscale("log")
axs[0].legend(prop={'size':8.0}, loc= "upper right")

axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
axs[1].set_xlabel('Time [Myr]')
axs[1].set_yscale("log")
#axs[1].set_xscale("log")
axs[1].set_ylim((1e-6,5e3))
axs[1].scatter(ts,CO2_num,label="CO2",color=species_color["CO2"])
axs[1].scatter(ts,O2_num,label="O2",color=species_color["O2"])
axs[1].scatter(ts,CO_num,label="CO",color=species_color["CO"])
axs[1].scatter(ts,H2_num,label="H2",color=species_color["H2"])
axs[1].legend(prop={'size':8.0},loc='upper left')
###############################################################################
"""

#O escape adjustment
#Oesc include only thermal escape, so non-thermal escape is added
Oesc = Oesc + O_nonthermal

#### marker style line ########################################################
fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
fig.subplots_adjust(hspace=0)
axs[0].plot(ts, Hesc, color='dodgerblue',label='H escape',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,Oesc, color='g',label='O escape',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,O3dep, color='#83d0c1',label='O3 deposition',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,H2O2dep, color='#bebddb',label='H2O2 deposition',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,Cesc,label='C escape',color='orange',marker="o",fillstyle="none",linestyle="--")
axs[0].plot(ts,CO2out,label='CO2 outgassing',color=species_color["CO2"],marker='x',fillstyle="none",linestyle="--")
axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
axs[0].set_ylim((5e5,2e12))
axs[0].set_yscale("log")
axs[0].legend(prop={'size':6.0}, loc= "upper left",ncol=2)

axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
axs[1].set_xlabel('Time [Myr]')
axs[1].set_yscale("log")
#axs[1].set_xscale("log")
axs[1].set_ylim((1e-7,5e3))
axs[1].plot(ts,CO2_num,label="CO2",color=species_color["CO2"],marker="o",fillstyle="none",linestyle="--")
axs[1].plot(ts,O2_num,label="O2",color='g',marker="o",fillstyle="none",linestyle="--")
axs[1].plot(ts,CO_num,label="CO",color=species_color["CO"],marker="o",fillstyle="none",linestyle="--")
#axs[1].plot(ts,H2_num,label="H2",color=species_color["H2"],marker="o",fillstyle="none",linestyle="--")
axs[1].legend(prop={'size':6.0},loc='upper left')
###############################################################################
#"""

### redox balance
C_O = []
H_O = []
C_H_O = []
for i in range(len(Hesc)):
    Closs = 4*Cesc[i]
    Hloss = Hesc[i]
    Oloss = 2*(Oesc[i]+3*O3dep[i]+H2O2dep[i])
    C_O.append(Closs/Oloss)
    H_O.append(Hloss/Oloss)
    C_H_O.append((Hloss+Closs)/Oloss)
    print(C_H_O[-1])

plt.figure(dpi=150)
plt.plot(ts, H_O, color='dodgerblue',label=r"$\Phi_H / 2\Phi_O$",marker="x",fillstyle="none",linestyle="--",markersize=10)
plt.plot(ts, C_O, color='orange',label=r"$2\Phi_C / \Phi_O$",marker="^",fillstyle="none",linestyle="--",markersize=10)
plt.plot(ts, C_H_O, color='k',label=r"$(\Phi_H+4\Phi_C) / 2\Phi_O$",marker="o",fillstyle="none",linestyle="--")
plt.xlabel('Time [Myr]')
plt.ylabel('Redox Balance')
plt.ylim((-0.1,1.5))
plt.legend(prop={'size':7.0})