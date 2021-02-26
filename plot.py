#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 15:08:57 2018

@author: shungokoyama
"""

#this file is made to plot the results of chaffin code
#ex)
#1 compare_converged('converged_Tsurf_209_CO2_6.3mbar_Oesc_1e7_depv_0.02.h5',all_species,same_figure=False,file_1="converged_standardwater_chaffin.h5")
#  compare_converged_mixing('converged_standardwater_CO2_1bar.h5',all_species, same_figure=False,file_1="converged_standardwater_chaffin.h5")
#2 plot_H_escape_history('H_esc_flux_history.h5')
#3 plot_n_current('O_escape_factor_2_to_normal.h5', 1299, ['CO2','H2O','H','H2','O2','OH','HO2','O3'],same_figure=True)
#4 plot_mixing('ppm_80_alt_60.h5', 999, ['CO2','H2O'])
#  plot_mixing_compare('CO2_1bar_Oesc_normal_to_2.h5',1,99, all_species, same_figure=True)
#5 calc_CO2_pressure_file('ppm_80_alt_80.h5', 999)
#6 calc_CO2_pressure_array([1,23,3,,,,])
#  reaction_rates_compare('CO2_1bar_Oesc_normal_to_2.h5',1,1000,[])
#  flux_history_compare('CO2_1bar_Oesc_normal_to_2.h5',1,900,['O','H','H2','CO'])
#...

#################################################################################
import h5py
import matplotlib.pyplot as plt
import numpy as np
import copy as cp
from matplotlib import animation as animation
import time

#initial_condition.h5のファイル形式に合わせたdic
chaffin_dic = {'H2O':4, 'N2':5, 'H2O2':6, 'OH':7, 'O3':8, 'O1D':9, 'CO2':11,
               'O2':15,'CO':19, 'HO2':24, 'H':25, 'O':26, 'CO2pl':27, 'HOCO':29, 'Ar':30,'H2':33}
#converged_standardwater.h5の出力形式に合わせたdic
my_dic_con = {'H2O':20, 'N2':30, 'H2O2':7, 'OH':21, 'O3':33, 'O1D':18, 'CO2':22,
               'O2':1,'CO':28, 'HO2':5, 'H':26, 'O':0, 'CO2pl':25, 'HOCO':15, 'Ar':11,'H2':19}
#one_year_response_to_80ppm_at_60km.h5に合わせたdic
my_dic_one = {'H2O':20, 'N2':30, 'H2O2':7, 'OH':21, 'O3':33, 'O1D':18, 'CO2':22,
               'O2':1,'CO':28, 'HO2':5, 'H':26, 'O':32, 'CO2pl':25, 'HOCO':15, 'Ar':11,'H2':19}
#ppm_80_alt_60.h5に合わせたdic
my_dic_ppm = {'H2O':20, 'N2':30, 'H2O2':7, 'OH':21, 'O3':0, 'O1D':18, 'CO2':22,
               'O2':1,'CO':28, 'HO2':5, 'H':26, 'O':32, 'CO2pl':25, 'HOCO':15, 'Ar':11,'H2':19}
all_species = ['H2O','N2','H2O2','OH','O3','O1D','CO2','O2','CO','HO2','H','O','CO2pl','HOCO','Ar','H2']
all_species_except_ArN2=['H2O','H2O2','OH','O3','O1D','CO2','O2','CO','HO2','H','O','CO2pl','HOCO','H2']
all_species_except_ArN2HOCO=['H2O','H2O2','OH','O3','O1D','CO2','O2','CO','HO2','H','O','CO2pl','H2']
all_species_num = [0,1,5,7,11,15,18,19,20,21,22,25,26,28,30,33]
all_species_num_ppm = [0,1,5,7,11,15,18,19,20,21,22,25,26,28,30,32,33]
#all_species_num
all_species_num_chaffin = [4,5,6,7,8,9,11,15,19,24,25,26,27,29,30,33]
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
reactions_color=["#dead91","#bebddb", "#9e9ac8","#83d0c1","#b593e2","#1f78b4", "#269e56", "#eed2d4","#ef3b2c",
"#fc9aa1","#41ae76","#00702d","#fd8d3c","#808080","#cccccc","#fdd0a2",]
color_list=['dodgerblue','yellowgreen','green','orange','saddlebrown','tomato','grey','black','darkgoldenrod']

###############################################################################
#this function creates 2 figures of species density profile with respect ot altitude
#   input: file name of your original converged case
#   input: array of species you want to display ['H2O', 'CO2',,,,]
#   input: whether plot in same figure or not, default is false, so plot in a different figure
#compare_converged("./converged_Tsurf_270_CO2_1000mbar_Oesc_1.2e8_depv_0.02_nofixedCO2_noHOCO_Tmod_1billion.h5",all_species_except_ArN2HOCO,same_figure=False,file_1="converged_standardwater_chaffin.h5")
def compare_converged(file_2,species_displayed,same_figure=False,file_1="converged_standardwater_chaffin.h5"):
    
    file = file_1
    file2 = file_2
    
    #one_year fileとppm fileとconverged fileによって対応するspeciesリストが変わるためその調整をする
    if 'one' in file2:
        my_dic = my_dic_one
    elif 'ppm' in file2:
        my_dic = my_dic_ppm
    else:
        my_dic = my_dic_con
    
    if 'chaffin' in file:
        my_dic2 = chaffin_dic
    else:
        my_dic2 = my_dic_con
    
    

    h5file = h5py.File(file,"r") #read a h5 file
    h5file2 = h5py.File(file2,"r") #read a h5 file

    #1
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    initial_mat = h5file["n_current/n_current_mat"].value
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list

    #2
    alt_list2 = h5file2["n_current/alt"].value #altitude list
    alt_list2 = alt_list2/1e5 #cm to km
    initial_mat2 = h5file2["n_current/n_current_mat"].value
    species2 = h5file2["n_current/species"].value #species type: object containing string element
    species_list2 = []
    i = 0
    for i in range(len(species2)):
        species_list2.append(species2[i]) #add string element to the list

    #chaffin converged
    plt.figure(1,dpi=150)
    for j in species_displayed:
        #plt.plot(initial_mat[chaffin_dic[j],:],alt_list[1:-1],label = species_list[chaffin_dic[j]],color=species_color[j],linestyle=':')
        plt.plot(initial_mat[my_dic2[j],:],alt_list[1:-1],color=species_color[j],linestyle=':')
    plt.xscale("log")
    #plt.legend()
    plt.title('Chaffin converged Standard water case')
    plt.xlabel(r"Number Density$[cm^{-3}]$")
    plt.ylabel("Altitude[km]")
    plt.tick_params(labelsize=13)


    #自分で作った平衡状態
    if same_figure==False:
        plt.figure(2,dpi=150)
    for j in species_displayed:
        plt.plot(initial_mat2[my_dic[j],:],alt_list2[1:-1],label = species_list2[my_dic[j]],color=species_color[j])
    plt.xscale("log")
    plt.legend(loc="upper right", fontsize=9.5)
    #plt.title('your converged')
    plt.xlabel(r"Number Density$[cm^{-3}]$",size=15)
    plt.ylabel("Altitude[km]",size=15)
    plt.tick_params(labelsize=13)
    plt.xlim((1e-2,1e20))
    
    print("CO: "+str(sum(initial_mat2[28,:])))
    print("CO mixing ratio is: " + str(sum(initial_mat2[28,:])/(sum(initial_mat2[22,:])+sum(initial_mat2[11,:])+sum(initial_mat2[30,:])+sum(initial_mat2[28,:])))) #CO, CO2
    print("CO2 pressure:" + str(calc_CO2_pressure_array(initial_mat2[22,:])))
    print("HOCO dep flux:" + str((initial_mat2[15,0]*0.02)))
    print("H2O col: "+str(sum(initial_mat2[my_dic['H2O'],:])*2e5))
    #print("H2O2 deposition flux: " + str(initial_mat2[7,1]*0.02))
    #print("O3 depostion flux: " + str(initial_mat2[33,1]*0.02))
    #print("H escape flux: " + str(966.45*initial_mat2[26,98]+2*3.1163*initial_mat2[19,98]))
    #return sum(initial_mat2[28,:])/sum(initial_mat2[22,:]) #CO, CO2
    #return initial_mat2[my_dic['O'],:], initial_mat2[my_dic['O2'],:],initial_mat2[my_dic['O3'],:]

###############################################################################
#this function creates 2 figures of species density profile with respect to altitude
#   input: file name of your original converged case
#   input: array of species you want to display ['H2O', 'CO2',,,,]
#   input: whether plot in same figure or not, default is false, so plot in a different figure
def compare_converged_mixing(file_2,species_displayed, same_figure=False,file_1="converged_standardwater_chaffin.h5"):

    file = file_1
    file2 = file_2

    #one_year fileとppm fileとconverged fileによって対応するspeciesリストが変わるためその調整をする
    if 'one' in file2:
        my_dic = my_dic_one
    elif 'ppm' in file2:
        my_dic = my_dic_ppm
    else:
        my_dic = my_dic_con
    
    
    #this controls which species dictionary should be used for file_1, default = chaffin converged one
    if 'chaffin' in file:
        my_dic2 = chaffin_dic
        all_species_num_1 = all_species_num_chaffin#this is needed to calc mixing ratio
    else:
        my_dic2 = my_dic_con
        all_species_num_1 = all_species_num
    

    h5file = h5py.File(file,"r") #read a h5 file
    h5file2 = h5py.File(file2,"r") #read a h5 file

    #1
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    initial_mat = h5file["n_current/n_current_mat"].value
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list

    #2
    alt_list2 = h5file2["n_current/alt"].value #altitude list
    alt_list2 = alt_list2/1e5 #cm to km
    initial_mat2 = h5file2["n_current/n_current_mat"].value
    species2 = h5file2["n_current/species"].value #species type: object containing string element
    species_list2 = []
    i = 0
    for i in range(len(species2)):
        species_list2.append(species2[i]) #add string element to the list

    #number density to mixing ratio
    initial_temp = cp.deepcopy(initial_mat)
    for i in range(len(initial_mat[0,:])):
        for j in all_species_num_1:
            initial_mat[j,i] = initial_temp[j,i]/sum(initial_temp[all_species_num_1,i])
    #number density to mixing ratio
    initial_temp2 = cp.deepcopy(initial_mat2)
    for i in range(len(initial_mat2[0,:])):
        for j in all_species_num:
            initial_mat2[j,i] = initial_temp2[j,i]/sum(initial_temp2[all_species_num,i])


    #chaffin のconverged_standardwater.h5形式
    plt.figure(1,dpi=150)
    for j in species_displayed:
        #plt.plot(initial_mat[chaffin_dic[j],:],alt_list[1:-1],label = species_list[chaffin_dic[j]],color=species_color[j],linestyle=':')
        plt.plot(initial_mat[my_dic2[j],:],alt_list[1:-1],color=species_color[j],linestyle=':')
    plt.xscale("log")
    #plt.legend()
    plt.title('Chaffin converged Standard water case')
    plt.xlabel("mixing ratio")
    plt.ylabel("Altitude[km]")


    #自分で作った平衡状態
    if same_figure==False:
        plt.figure(2,dpi=150)
    for j in species_displayed:
        plt.plot(initial_mat2[my_dic[j],:],alt_list2[1:-1],label = species_list2[my_dic[j]],color=species_color[j])
    plt.xscale("log")
    plt.legend()
    #plt.title('CO atmo/ 100mbar,Oesc1.2e8,Tsurf209K,nofixedCO2')
    plt.xlabel("mixing ratio")
    plt.ylabel("Altitude[km]")

    #h2O=initial_mat[20,:]
    return initial_mat2[28,:], initial_mat2[22,:] #CO, CO2

##############################################################################
#this function plots histrory of all outfluxes
#plot_outflux("Hesc_12e8_Tmod/Hesc_depfluxes_CO2_6.3mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5")
#HO 
def plot_outflux(readfile,title="C escape from 4.1Ga"):
    file = readfile
    h5file = h5py.File(file,"r")
    time = h5file["fluxes/times"].value
    Hescflux = h5file["fluxes/fluxvals"].value
    Htotalflux=h5file["outfluxes/H"].value
    Ototalflux=h5file["outfluxes/O"].value
    H2O2dep=h5file["depositions/H2O2"].value
    O3dep=h5file["depositions/O3"].value
    HO2dep=h5file["depositions/HO2"].value
    CO2esc = h5file["outfluxes/CO2esc"].value
    CO2outgassing = h5file["outfluxes/CO2out"].value
    Oreturn = h5file["outfluxes/Oreturn"].value
    
    H_O_ratio = np.array(Htotalflux)/np.array(Ototalflux)
    
    time=time/3.15e7
    plt.figure(1,dpi=150)
    plt.plot(time, Hescflux[0:-1],'lightblue',label="H escape flux")
    #plt.plot(time, Htotalflux[0:-1], label="total H outflux")
    #plt.plot(time, Ototalflux[0:-1], label="total O outflux")
    #plt.plot(time, H2O2dep[0:-1],label="H2O2 deposition flux")
    plt.plot(time, O3dep[0:-1],label="O3 deposition flux")
    plt.plot(time, CO2esc[0:-1],label="C escape flux")
    plt.plot(time, CO2outgassing[0:-1],label="CO2 outgassing")
    plt.plot(time, Oreturn[0:-1],label="O return flux")
    plt.ylim((1e10,3e11))
    plt.xlim((1e6,1e8))
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Time[year]",fontsize=15)
    plt.ylabel(r"Flux$[cm^{-2}s^{-1}]$",fontsize=15)
    plt.legend(loc="upper left", fontsize=10)
    plt.tick_params(labelsize=13)
    plt.title(title)
    
    return H_O_ratio

###############################################################################
#plot_outflux_redox("Hesc_12e8_Tmod/Hesc_depfluxes_CO2_6.3mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5")
# This function calcs outfluxes as redox parameters not by total fluxes.
# H2O=0, H2O2= one oxygen, O3 = 3 oxygens, etc.
def plot_outflux_redox(readfile,title="time-response of outflux"):
    file = readfile
    h5file = h5py.File(file,"r")
    time = h5file["fluxes/times"].value
    Hescflux = h5file["fluxes/fluxvals"].value
    #Htotalflux=h5file["outfluxes/H"].value
    #Ototalflux=h5file["outfluxes/O"].value
    H2O2dep=h5file["depositions/H2O2"].value
    O3dep=h5file["depositions/O3"].value
    HO2dep=h5file["depositions/HO2"].value
    
    Ototalflux_redox=H2O2dep + O3dep*3 + HO2dep*1.5 + 2.4e8
    
    time=time/3.15e7
    plt.figure(1,dpi=150)
    plt.plot(time, Hescflux[0:-1],'lightblue',label="H escape flux")
    #plt.plot(time, Htotalflux[0:-1], label="total H outflux")
    plt.plot(time, Ototalflux_redox[0:-1], label="total O outflux")
    plt.plot(time, H2O2dep[0:-1],label="H2O2 deposition flux")
    plt.plot(time, O3dep[0:-1],label="O3 deposition flux")
    plt.ylim((3e7,5e11))
    plt.xlim((1e-1,1e9))
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Time[year]",fontsize=15)
    plt.ylabel(r"Flux$[cm^{-2}s^{-1}]$",fontsize=15)
    plt.legend(loc="upper left", fontsize=10)
    plt.tick_params(labelsize=13)
    #plt.title( title+"Oesc1.2e8->2.4e8")
    
    #return H_O_ratio



###############################################################################
#This function creates the figure of H escape flux history over time
#   input: file name of H escape flux history
#   output: H escape flux over time
def plot_H_escape_history(file_2):
    file = file_2
    h5file = h5py.File(file,"r")
    time = h5file["fluxes/times"].value
    fluxvals_mat = h5file["fluxes/fluxvals"].value
    #fluxvals_mat: (1302,6,4) = (flux, Altitude, water)
    # x elements corrspond to H escape flux over time
    # y elements corresponds to Altitude [20,40,60,80,100,120km]
    # z elemets correspond to water parcel added [20,40,60,80ppm]
    #ex) [:,2,2] -> [flux vector at 60km with 60ppm water added]


    #make a plot of H escape flux history
    plt.figure(1,dpi=100)
    '''
    for i in range(0,4):
        plt.plot(time, fluxvals_mat[2:-1, 0, i],'lightblue')
        plt.plot(time, fluxvals_mat[2:-1, 1, i], 'steelblue')
        plt.plot(time, fluxvals_mat[2:-1, 2, i], 'silver')
        plt.plot(time, fluxvals_mat[2:-1, 3, i], 'gray')
        plt.plot(time, fluxvals_mat[2:-1, 4, i], 'darksalmon')
        plt.plot(time, fluxvals_mat[2:-1, 5, i], 'red')
    '''
    plt.plot(time, fluxvals_mat[2:-1, 0, 0],'lightblue')
    #plt.plot(time, fluxvals_mat[2:-1, 1, 0],'steelblue')
    plt.ylim((9e7,1e9))
    plt.xlim((1,1e19))
    plt.yscale("log")
    plt.xscale("log")
    plt.title('CO2:1bar,   O escape normal -> 2')
    plt.xlabel("time(s)")
    plt.ylabel(r"H Escape Flux$[cm^{-2}s^{-1}]$")

    plt.show()
    #return fluxvals_mat

###############################################################################
#This function creates the figure of H escape flux history over time
#   input: an array of file names of H escape flux history
#        [r'$0.5\Phi_O$',r'$\Phi_O$',,]    
#   output: Each H escape flux over time
#ex)plot_H_escape_histories(['H_esc_flux_history_O_escape_factor_normal_to_0.5.h5','H_esc_flux_history_nochange.h5','H_esc_flux_history_O_escape_factor_normal_to_1.5.h5','H_esc_flux_history_O_escape_factor_normal_to_2.h5','H_esc_flux_history_O_escape_factor_normal_to_3.h5','H_esc_flux_history_O_escape_factor_normal_to_5.h5','H_esc_flux_history_O_escape_factor_normal_to_10.h5'],[r'$0.5\Phi_O$',r'$\Phi_O$',r'$1.5\Phi_O$',r'$2\Phi_O$',r'$3\Phi_O$',r'$5\Phi_O$',r'$10\Phi_O$'])
def plot_H_escape_histories(file_array,label_array):
    files = file_array
    #h5files = []
    fluxvals_mat = []
    h5file = h5py.File(files[0],"r")
    time = h5file["fluxes/times"].value #get time array here
    plt.figure(1,dpi=150)
    i=0
    for file in files:
        h5file = h5py.File(file,"r")
        fluxvals_mat = h5file["fluxes/fluxvals"].value
        plt.plot(time, fluxvals_mat[2:-1, 0, 0],label=label_array[i])
        i+=1

    plt.legend(prop={'size':15})
    plt.ylim((1e8,1e10))
    plt.xlim((1,1e17))
    plt.yscale("log")
    plt.xscale("log")
    plt.title('H escape flux change')
    plt.xlabel("time(s)")
    plt.ylabel(r"H Escape Flux$[cm^-2s^-1]$")

    plt.show()

###############################################################################
#This function creates the figure of H escape flux histories and some specific species mixing ratio over time
#vertically adjusted each other
#plot_Hesc_species("Hesc_12e8_Tmod/Hesc_depfluxes_CO2_500mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5","CO2_12e8_Tmod/CO2_500mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5",['CO','O2'],2.4e8)
def plot_Hesc_species(file_for_Hesc,file_for_species,speciesname_array,Oescflux_after=4.2e8):
    #file = file_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    fluxvals_mat = []
    h5file = h5py.File(file_for_Hesc,"r")
    time_array = h5file["fluxes/times"].value #get time array here
    time_array = time_array/3.15e7
    
    #create 2 subplots vertically adjusted
    fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
    fig.subplots_adjust(hspace=0)
    
    ##ax0 H escape flux
    fluxvals_mat = h5file["fluxes/fluxvals"].value
    H2O2dep=h5file["depositions/H2O2"].value
    O3dep=h5file["depositions/O3"].value
    HO2dep=h5file["depositions/HO2"].value
    #Ototalflux_redox=H2O2dep + O3dep*3 + HO2dep*1.5 + Oescflux_after
    Ototalflux = h5file["outfluxes/O"].value
    CO2esc = h5file["outfluxes/CO2esc"].value
    CO2outgassing = h5file["outfluxes/CO2out"].value
    Oreturn = h5file["outfluxes/Oreturn"].value
    H2outgassing = h5file["outfluxes/H2out"].value
    Oesc = np.zeros(len(time_array)+1)
    Oesc[:]=Oescflux_after
    
    axs[0].plot(time_array, fluxvals_mat[0:-1], color='dodgerblue',label='H escape')
    #axs[0].plot(time_array,Ototalflux[0:-1], color='g',label='total O outflux')
    axs[0].plot(time_array,Oesc[0:-1], color='g',label='O escape')
    #axs[0].plot(time_array,H2O2dep[0:-1], color='mediumpurple',label='H2O2 deposition')

    #axs[0].plot(time_array,CO2esc[0:-1],label='C escape flux',color='orange')
    axs[0].plot(time_array,CO2outgassing[0:-1],label='CO outgassing',color='orange',linestyle='-.')
    axs[0].plot(time_array,H2outgassing[0:-1],label='H2 outgassing',color='dodgerblue',linestyle='-.')

    #axs[0].plot(time_array,Oreturn[0:-1],label='O return flux', color='g',linestyle='-.')
    #axs[0].plot(time_array,O3dep[0:-1],label='O3 deposition flux',color='limegreen')

    #axs[0].set_xlabel('time(s)')
    axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
    axs[0].set_ylim((5e6,1e14))
    #axs[0].set_ylim((1e8,5e12))
    axs[0].set_yscale("log")
    #axs[0].tick_params(axis='y')
    axs[0].legend(prop={'size':8.0},loc='upper right')
    
    #ax2 species mixing ratio
    axs[1].set_ylabel('Mixing ratio')  # we already handled the x-label with ax1
    axs[1].set_xlabel('Time [year]')
    axs[1].set_yscale("log")
    axs[1].set_xscale("log")
    axs[1].set_ylim((1e-5,1.2))
    #axs[1].set_ylim((2e-4,1.2e-3))
    for speciesname in speciesname_array:
        array_species_overtime = species_mixing_overtime(file_for_species,speciesname,len(time_array-1))
        axs[1].plot(time_array,array_species_overtime,label=speciesname,color=species_color[speciesname])
        #axs[1].tick_params(axis='y')
        print(array_species_overtime[0]-array_species_overtime[len(time_array)-1])
    
    ##setting for both ax1 and ax2
    axs[1].legend(prop={'size':9.5})
    #axs[1].legend(prop={'size':9.5},loc='upper right')
    axs[1].set_xlim((1e-1,1e7))
    #plt.xscale("log")
    #plt.title('H escape flux [H20:9.5pr]')
    plt.show()


###############################################################################
#This function creates the figure species flux histories and some specific species partial pressure over time
#vertically adjusted each other
#flux_pressure("output/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5","output/CO2_1mbar_Tsurf_210_Tian_BD.h5",["O2","CO2"])
#flux_pressure("after3.95Ga_lowCO2/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5","after3.95Ga_lowCO2/CO2_1mbar_Tsurf_210_Tian_BD.h5",["O2","CO2"])
def flux_pressure(file_for_Hesc,file_for_species,speciesname_array,Oescflux_after=2.63e8):
    #file = file_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    fluxvals_mat = []
    h5file = h5py.File(file_for_Hesc,"r")
    time_array = h5file["fluxes/times"].value #get time array here
    time_array = time_array/3.15e7
    
    #create 2 subplots vertically adjusted
    fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
    fig.subplots_adjust(hspace=0)
    
    ##ax0 H escape flux
    fluxvals_mat = h5file["fluxes/fluxvals"].value
    #H2O2dep=h5file["depositions/H2O2"].value
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
    
    #axs[0].plot(time_array, fluxvals_mat[0:-1], color='dodgerblue',label='H escape')
    axs[0].plot(time_array, Hesc[0:-1], color='dodgerblue',label='H escape')

    #axs[0].plot(time_array,Ototalflux[0:-1], color='g',label='total O outflux')
    axs[0].plot(time_array,Othermal[0:-1], color='g',label='O thermal escape flux')
    axs[0].plot(time_array,CO2esc[0:-1],label='C escape flux',color='orange')
    axs[0].plot(time_array,CO2outgassing[0:-1],label='CO2 outgassing flux',color='orange',linestyle='-.')
    #axs[0].plot(time_array,H2outgassing[0:-1],label='H2 outgassing flux',color='dodgerblue',linestyle='-.')

    #axs[0].plot(time_array,Oreturn[0:-1],label='O return flux', color='g',linestyle='-.')
    #axs[0].plot(time_array,O3dep[0:-1],label='O3 deposition flux',color='limegreen')

    #axs[0].set_xlabel('time(s)')
    axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
    #axs[0].set_ylim((5e6,1e11))
    axs[0].set_ylim((5e5,2e12))
    axs[0].set_yscale("log")
    #axs[0].tick_params(axis='y')
    axs[0].legend(prop={'size':8.0}, loc= "upper left")
    
    #ax2 species mixing ratio
    axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
    axs[1].set_xlabel('Time [year]')
    axs[1].set_yscale("log")
    axs[1].set_xscale("log")
    axs[1].set_ylim((1e-5,5e3))
    #axs[1].set_ylim((2e-4,1.2e-3))
    for speciesname in speciesname_array:
        array_species_overtime = species_pressure_overtime(file_for_species,speciesname,len(time_array-1))
        axs[1].plot(time_array,array_species_overtime,label=speciesname,color=species_color[speciesname])
        #axs[1].tick_params(axis='y')
    
    
    ##setting for both ax1 and ax2
    #axs[1].legend(prop={'size':9.5})
    axs[1].legend(prop={'size':9.5},loc='upper left')
    #axs[1].set_xlim((1e-1,6e8))
    axs[1].set_xlim((1e2,1e8))

    #axs[1].set_xlim((1e7,1.1E7))
    #plt.xscale("log")
    #plt.title('H escape flux [H20:9.5pr]')
    plt.show()



def flux_pressure_temp(file_for_Hesc,file_for_species,speciesname_array,Oescflux_after=2.63e8):
    #file = file_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    fluxvals_mat = []
    h5file = h5py.File(file_for_Hesc,"r")
    time_array = h5file["fluxes/times"].value #get time array here
    time_array = time_array/3.15e7/1e6
    
    #create 2 subplots vertically adjusted
    fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
    fig.subplots_adjust(hspace=0)
    
    ##ax0 H escape flux
    fluxvals_mat = h5file["fluxes/fluxvals"].value
    #H2O2dep=h5file["depositions/H2O2"].value
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
    
    #axs[0].plot(time_array, fluxvals_mat[0:-1], color='dodgerblue',label='H escape')
    #axs[0].plot(time_array, Hesc[0:-1], color='dodgerblue',label='H escape')

    #axs[0].plot(time_array,Ototalflux[0:-1], color='g',label='total O outflux')
    #axs[0].plot(time_array,Othermal[0:-1], color='g',label='O thermal escape flux')
    axs[0].plot(time_array,CO2esc[0:-1],label='C escape flux',color='orange')
    axs[0].plot(time_array,CO2outgassing[0:-1],label='CO2 outgassing flux',color='orange',linestyle='-.')
    #axs[0].plot(time_array,H2outgassing[0:-1],label='H2 outgassing flux',color='dodgerblue',linestyle='-.')

    #axs[0].plot(time_array,Oreturn[0:-1],label='O return flux', color='g',linestyle='-.')
    #axs[0].plot(time_array,O3dep[0:-1],label='O3 deposition flux',color='limegreen')

    #axs[0].set_xlabel('time(s)')
    axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
    #axs[0].set_ylim((5e6,1e11))
    axs[0].set_ylim((5e8,2e12))
    axs[0].set_yscale("log")
    #axs[0].tick_params(axis='y')
    axs[0].legend(prop={'size':8.0}, loc= "upper left")
    
    #ax2 species mixing ratio
    axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
    axs[1].set_xlabel('Time [Myr]')
    axs[1].set_yscale("log")
    #axs[1].set_xscale("log")
    axs[1].set_ylim((1e-1,5e3))
    #axs[1].set_ylim((2e-4,1.2e-3))
    for speciesname in speciesname_array:
        array_species_overtime = species_pressure_overtime(file_for_species,speciesname,len(time_array-1))
        axs[1].plot(time_array,array_species_overtime,label=speciesname,color=species_color[speciesname])
        #axs[1].tick_params(axis='y')
    
    
    ##setting for both ax1 and ax2
    #axs[1].legend(prop={'size':9.5})
    axs[1].legend(prop={'size':9.5},loc='upper left')
    #axs[1].set_xlim((1e-1,6e8))
    axs[1].set_xlim((2e1,1.5e2))

    #axs[1].set_xlim((1e7,1.1E7))
    #plt.xscale("log")
    #plt.title('H escape flux [H20:9.5pr]')
    plt.show()


###############################################################################
#This function creates the figure species flux histories and some specific species partial pressure over time
#vertically adjusted each other
#plot flux during 4.1~3.5Ga continuously
#flux_pressure("after3.95Ga_lowCO2/Hesc_depfluxes_CO2_1mbar_Tsurf_210_Tian_BD.h5","after3.95Ga_lowCO2/CO2_1mbar_Tsurf_210_Tian_BD.h5",["O2","CO2"])
def evol_flux_pressure(file_for_Hesc, file_for_species, file_for_Hesc2, file_for_species2, speciesname_array,Oescflux_after=2.63e8):
    #file = file_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    fluxvals_mat = []
    h5file = h5py.File(file_for_Hesc,"r")
    time_array1 = h5file["fluxes/times"].value #get time array here
    time_array1 = 4.1 - np.array(time_array1)/3.15e16
    
    h5file2 = h5py.File(file_for_Hesc2,"r")
    time_array2 = h5file2["fluxes/times"].value #get time array here
    time_array2 = time_array1[-1] - np.array(time_array2)/3.15e16
    
    time_array = np.append(time_array1,time_array2)
    #create 2 subplots vertically adjusted
    fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
    fig.subplots_adjust(hspace=0)
    
    ##first file
    fluxvals_mat = np.append(h5file["fluxes/fluxvals"].value, h5file2["fluxes/fluxvals"].value)
    H2O2dep = np.append(h5file["depositions/H2O2"].value, h5file2["depositions/H2O2"].value)
    O3dep = np.append(h5file["depositions/O3"].value, h5file2["depositions/O3"].value)
    HO2dep = np.append(h5file["depositions/HO2"].value, h5file2["depositions/HO2"].value)
    #Ototalflux_redox=H2O2dep + O3dep*3 + HO2dep*1.5 + Oescflux_after
    Ototalflux = np.append(h5file["outfluxes/O"].value, h5file2["outfluxes/O"].value)
    CO2esc = np.append(h5file["outfluxes/CO2esc"].value, h5file2["outfluxes/CO2esc"].value)
    CO2outgassing = np.append(h5file["outfluxes/CO2out"].value, h5file2["outfluxes/CO2out"].value)
    Oreturn = np.append(h5file["outfluxes/Oreturn"].value, h5file2["outfluxes/Oreturn"].value)
    Oesc = np.zeros(len(time_array)+2)#change later
    Oesc[:]=Oescflux_after
    
    
    axs[0].plot(time_array,CO2outgassing[1:-1],label='CO2 outgassing flux',color='k',linestyle='--')
    #axs[0].plot(time_array, fluxvals_mat[1:-1], color='dodgerblue',label='H escape')
    #axs[0].plot(time_array,Ototalflux[0:-1], color='g',label='total O outflux')
    axs[0].plot(time_array,CO2esc[1:-1],label='C escape flux',color='orange')
    #axs[0].plot(time_array,Oreturn[0:-1],label='O return flux', color='g',linestyle='-.')
    axs[0].plot(time_array,O3dep[1:-1],label='O3 deposition flux',color='g')
    #axs[0].plot(time_array,Oesc[1:-1], color='g',label='O escape flux')

    #axs[0].set_xlabel('time(s)')
    axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
    axs[0].set_ylim((1e9,1e11))
    #axs[0].set_ylim((1e8,5e12))
    axs[0].set_yscale("log")
    #axs[0].tick_params(axis='y')
    axs[0].legend(prop={'size':8.0})
    
    #ax2 species mixing ratio
    axs[1].set_ylabel('Partial Pressure [mbar]')  # we already handled the x-label with ax1
    axs[1].set_xlabel('Time [Ga]')
    axs[1].set_yscale("log")
    #axs[1].set_xscale("log")
    axs[1].set_ylim((1e-2,5e3))
    #axs[1].set_ylim((2e-4,1.2e-3))
    for speciesname in speciesname_array:
        array_species_overtime1 = species_pressure_overtime(file_for_species,speciesname,len(time_array1-1))
        array_species_overtime2 = species_pressure_overtime(file_for_species2,speciesname,len(time_array2-1))
        array_species_overtime = np.append(array_species_overtime1,array_species_overtime2)
        axs[1].plot(time_array,array_species_overtime,label=speciesname,color=species_color[speciesname])
        #axs[1].tick_params(axis='y')
    
    
    ##setting for both ax1 and ax2
    axs[1].legend(prop={'size':8.0})
    #axs[1].legend(prop={'size':9.5},loc='upper left')
    axs[1].set_xlim((4.1,3.5))
    #plt.xscale("log")
    #plt.title('H escape flux [H20:9.5pr]')
    plt.show()
    
    return array_species_overtime, time_array


def CO2outgassing(t,start,initial_flux):
    t1_Gyr=(t/3.15e16) + start #Gyr
    f_outgas= initial_flux*np.exp(-2.624*t1_Gyr)
    return f_outgas

def CO2outgassing_with_escape(file_for_Hesc):
    #file = file_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    start_time = 0.479677377
    
    h5file = h5py.File(file_for_Hesc,"r")
    time_array = h5file["fluxes/times"].value #get time array here
    time_array = time_array
    
    lowestCO2 = [CO2outgassing(t,start_time,2.8e10) for t in time_array]
    lowerCO2 = [CO2outgassing(t,start_time,7.0e10) for t in time_array]
    higherCO2 = [CO2outgassing(t,start_time,1.4e11) for t in time_array]
    
    
    CO2esc = h5file["outfluxes/CO2esc"].value
    #CO2outgassing = h5file["outfluxes/CO2out"].value
    
    plt.figure(1,dpi=150)
    plt.plot(time_array,CO2esc[0:-1],label='C escape flux',color='orange')
    plt.plot(time_array,higherCO2,label='higher CO2 outgassing',color='k')
    plt.plot(time_array,lowerCO2,label='lower CO2 outgassing',color='k',linestyle='--')
    plt.plot(time_array,lowestCO2,label='lowest CO2 outgassing',color='k',linestyle='-.')

    #axs[0].set_xlabel('time(s)')
    plt.ylabel(r"Flux$[cm^{-2}s^{-1}]$")
    plt.xlabel("Time [years]")
    plt.ylim((1e9,1e11))
    plt.yscale("log")
    plt.legend(prop={'size':10.0})
    
    #axs[1].legend(prop={'size':9.5},loc='upper left')
    plt.xlim((1e-1,6e8))
    plt.xscale("log")
    #plt.title('H escape flux [H20:9.5pr]')
    plt.show()

###############################################################################
#This function creates the figure of H escape flux histories and some specific species mixing ratio over time
#   input: an array of file names of H escape flux history,label_array,speciesname
#        [r'$0.5\Phi_O$',r'$\Phi_O$',,]    
#   output: Each H escape flux over time
#plot_H_escape_histories_with_species(["H_esc_flux_history_CO2_6.3mbar_H2O_9.5pr_Oesc_normal_to_2.h5","H_esc_flux_history_CO2_20mbar_H2O_9.5pr_Oesc_normal_to_2.h5","H_esc_flux_history_CO2_50mbar_H2O_9.5pr_Oesc_normal_to_2.h5"],["CO2_6.3mbar_H2O_9.5pr_Oesc_normal_to_2.h5","CO2_20mbar_H2O_9.5pr_Oesc_normal_to_2.h5","CO2_50mbar_H2O_9.5pr_Oesc_normal_to_2.h5"],'CO',["CO2:6.3mbar","CO2:20mbar","CO2:50mbar"])
def plot_H_escape_histories_with_species(file_array_for_Hesc,file_array_for_species,speciesname,label_array):
    files = file_array_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    fluxvals_mat = []
    h5file = h5py.File(files[0],"r")
    time_array = h5file["fluxes/times"].value #get time array here
    
    fig, ax1 = plt.subplots(dpi=100)
    ##ax1 H escape flux
    ax1.set_xlabel('time(s)')
    ax1.set_ylabel(r"H Escape Flux$[cm^{-2}s^{-1}]$")
    ax1.set_ylim((1e8,2e9))
    ax1.set_yscale("log")
    ax1.tick_params(axis='y')
    i=0
    for file in files:
        h5file = h5py.File(file,"r")
        fluxvals_mat = h5file["fluxes/fluxvals"].value
        ax1.plot(time_array, fluxvals_mat[0:-1],label=label_array[i],color=color_list[i])
        ax1.tick_params(axis='y')
        i+=1
    
    #ax2 CO mixing ratio
    ax2 = ax1.twinx()
    ax2.set_ylabel('CO mixing ratio')  # we already handled the x-label with ax1
    ax2.set_yscale("log")
    ax2.set_ylim((1e-5,1))
    i=0
    for file in file_array_for_species:
        array_species_overtime = species_mixing_overtime(file,speciesname,len(time_array-1))
        ax2.plot(time_array,array_species_overtime,label="CO",color=color_list[i],linestyle="--")
        ax2.tick_params(axis='y')
        i+=1
        print(array_species_overtime[0]-array_species_overtime[len(time_array)-1])
    
    ##setting for both ax1 and ax2
    ax1.legend(prop={'size':8},loc="upper left")
    plt.xlim((1e8,1e17))
    plt.xscale("log")
    plt.title('H escape flux [H20:9.5pr]')
    plt.show()



#############################################################################
#This function produces an array of specified species mixing ratio over time
#input: filename
#output: an array of specified species mixing ratio
def species_mixing_overtime(filename,speciesname, max_iter_num):
    my_dic = my_dic_ppm
    
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    #time_list = h5file["n_current/timelist"].value #time list
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list
    
    
    array_overtime=[]
    for iter_num in range(1,max_iter_num+1):
        iteration = h5file["n_current/iter_"+str(iter_num)].value #iteration n_current
        sum_all_species=0
        iteration_temp = cp.deepcopy(iteration)
        sum_specified_species = sum(iteration_temp[my_dic[speciesname],:])
        for j in all_species_num:
            sum_all_species += sum(iteration_temp[j,:])
        array_overtime.append(sum_specified_species/sum_all_species)
    
    return array_overtime

#############################################################################
#This function produces an array of specified species partial pressure over time
#input: filename
#output: an array of specified species mixing ratio
spmole = {'H2O':18, 'N2':28, 'H2O2':34, 'OH':17, 'O3':48, 'O1D':16, 'CO2':44,
               'O2':32,'CO':28, 'HO2':33, 'H':1, 'O':16, 'CO2pl':44, 'HOCO':45, 'Ar':40,'H2':2}
def species_pressure_overtime(filename,speciesname, max_iter_num):
    my_dic = my_dic_ppm
    
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    #time_list = h5file["n_current/timelist"].value #time list
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list
    
    array_overtime=[]
    for iter_num in range(1,max_iter_num+1):
        iteration = h5file["n_current/iter_"+str(iter_num)].value #iteration n_current
        iteration_temp = cp.deepcopy(iteration)
        p_specified_species = calc_pressure_array(iteration_temp[my_dic[speciesname],:],spmole[speciesname])
        #for j in all_species_num:
        #    sum_all_species += sum(iteration_temp[j,:])
        array_overtime.append(p_specified_species)
    
    return array_overtime


#############################################################################
#calc pOx change
#pOx after - pOx before
#calc_pox_timechange("Hesc_12e8_Tmod/Hesc_depfluxes_CO2_500mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5","CO2_12e8_Tmod/CO2_500mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5",2.4e8)
def calc_pox_timechange(file_for_Hesc,file_for_species,Oescflux_after):
    fluxvals_mat = []
    h5file = h5py.File(file_for_Hesc,"r")
    time_array = h5file["fluxes/times"].value #get time array here
    time_array = time_array/3.15e7
    
    #create 2 subplots vertically adjusted
    fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
    fig.subplots_adjust(hspace=0)
    
    ##ax0 H escape flux
    fluxvals_mat = h5file["fluxes/fluxvals"].value
    H2O2dep=h5file["depositions/H2O2"].value
    O3dep=h5file["depositions/O3"].value
    HO2dep=h5file["depositions/HO2"].value
    Ototalflux_redox=H2O2dep + O3dep*3 + HO2dep*1.5 + Oescflux_after
    
    axs[0].plot(time_array, fluxvals_mat[0:-1], color='k',label='H escape')
    axs[0].plot(time_array,Ototalflux_redox[0:-1], color='k',linestyle=':',label='net O outflux')
    #axs[0].set_xlabel('time(s)')
    axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
    axs[0].set_ylim((1e8,2e9))
    #axs[0].set_ylim((1e8,5e12))
    axs[0].set_yscale("log")
    #axs[0].tick_params(axis='y')
    axs[0].legend(prop={'size':9.5})
    
    #ax2 species mixing ratio
    axs[1].set_ylabel('Mixing ratio')  # we already handled the x-label with ax1
    axs[1].set_xlabel('Time [year]')
    axs[1].set_yscale("log")
    axs[1].set_xscale("log")
    #axs[1].set_ylim((1e-4,3.5e-4))
    #axs[1].set_ylim((2e-4,1.2e-3))
    
    #calc pox timechange
    species_timechange = {}
    species_timechange['O2'] = species_mixing_overtime(file_for_species,'O2',len(time_array-1))
    species_timechange['CO'] = species_mixing_overtime(file_for_species,'CO',len(time_array-1))
    species_timechange['H2'] = species_mixing_overtime(file_for_species,'H2',len(time_array-1))
    pox_timechange = [0]*len(species_timechange['O2'])
    
    for i in range(len(species_timechange['O2'])):
        pox_timechange[i] = 2*species_timechange['O2'][i] - species_timechange['CO'][i] - species_timechange['H2'][i]
    #|pOx|変換
    abspox_timechange = [abs(i) for i in pox_timechange]
    print("pOx change: ", pox_timechange[-1]-pox_timechange[0]) #最後-最初
    print("pOx(final) / pOx(initial): ", pox_timechange[-1]/pox_timechange[0])
    print(pox_timechange[0],pox_timechange[-1])
    
    #axs[1].plot(time_array,pox_timechange,label='fOx')
    axs[1].plot(time_array,abspox_timechange,label='|pOx|')
    #axs[1].tick_params(axis='y')
    
    
    ##setting for both ax1 and ax2
    axs[1].legend(prop={'size':9.5})
    #axs[1].legend(prop={'size':9.5},loc='upper left')
    axs[1].set_xlim((1e2,1e9))
    #plt.xscale("log")
    #plt.title('H escape flux [H20:9.5pr]')
    plt.show()

    
    
    


###############################################################################
#This function creates initial and ?iteration number density profile with respect to altitude
    #input: file name
    #input: iteration number you wanna see
    #input: species list array you want to display ['H2O', 'CO',,,]
    #input: whether plot in same figure or not, default is false, so plot in a different figure
#plot_n_current("Tinf480_middleCO2_euv15.5_2252K/CO2_1mbar_Tsurf_210_Tian_BD.h5", 3999,all_species_except_ArN2HOCO, same_figure=False)
#plot_n_current("CO2_12e8_Tmod/CO2_500mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5", 2000,all_species_except_ArN2HOCO, same_figure=True)
def plot_n_current(file_ppm_alt,iteration_number, species_displayed,same_figure=False):

    file = file_ppm_alt
    h5file = h5py.File(file,"r") #read a h5 file
    time_list = h5file["n_current/timelist"].value #time list
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    species_list = []
    for i in range(len(species)):
        print(species[i])
        species_list.append(species[i]) #add string element to the list
    initial_mat = h5file["n_current/init"].value #initial n_current
    iteration = h5file["n_current/iter_"+str(iteration_number)].value #iteration n_current


    #iteration number to time unit
    time_second = time_list[iteration_number-1]
    proper_time = time_transfer(time_second)

    #one_year fileとppm fileとconverged fileによって対応するspeciesリストが変わるためその調整をする
    my_dic = my_dic_ppm
    print(calc_pressure_array(iteration[22,:],44))
    plt.figure(1,dpi=150)
    for j in species_displayed:
        plt.plot(iteration[my_dic[j],:],alt_list[1:-1],label = species_list[my_dic[j]],color=species_color[j])
    plt.xscale("log")
    plt.legend()
    #plt.title(proper_time+' species profile')
    plt.xlabel(r"Number Density$[cm^{-3}]$",size=15)
    plt.ylabel("Altitude[km]",size=15)
    plt.xlim(1e-2,1e20)

    if same_figure == False:
        plt.figure(2,dpi=150)
    for j in species_displayed:
        plt.plot(initial_mat[my_dic[j],:],alt_list[1:-1],label = species_list[my_dic[j]],color=species_color[j],linestyle=':')
    plt.xscale("log")
    #plt.legend()
    #plt.title('initial species profile')
    plt.xlabel(r"Number Density$[cm^{-3}]$",size=15)
    plt.ylabel("Altitude[km]",size=15)
    plt.xlim(1e-2,1e20)
    plt.tick_params(labelsize=13)
    #return species_list
    return iteration
    

################################################################################
#plot_mixing("CO2_12e8_Tmod/CO2_1000mbar_Tsurf_270_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5", 1, ['H2O'], same_figure=False)
#this function plot mixing ratio profile for each species
    #input: file name
    #input: iteration number you wanna see
    #input: species list array you want to display ['H2O', 'CO',,,]
    #input: whether plot in same figure or not, default is false, so plot in a different figure
def plot_mixing(file_ppm_alt, iteration_number, species_displayed, same_figure=True):
    file = file_ppm_alt
    h5file = h5py.File(file,"r") #read a h5 file
    time_list = h5file["n_current/timelist"].value #time list
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list
    initial_mat = h5file["n_current/init"].value #initial n_current
    iteration = h5file["n_current/iter_"+str(iteration_number)].value #iteration n_current

    all_species_num = all_species_num_ppm 
    #number density to mixing ratio
    iteration_temp = cp.deepcopy(iteration)
    for i in range(len(iteration[0,:])):
        for j in all_species_num:
            iteration[j,i] = iteration_temp[j,i]/sum(iteration_temp[all_species_num,i])

    #number density to mixing ratio
    initial_temp = cp.deepcopy(initial_mat)
    for i in range(len(initial_mat[0,:])):
        for j in all_species_num:
            initial_mat[j,i] = initial_temp[j,i]/sum(initial_temp[all_species_num,i])
    #print(iteration[32,:])
    #iteration number to time unit
    time_second = time_list[iteration_number-1]
    proper_time = time_transfer(time_second)


    #one_year fileとppm fileとconverged fileによって対応するspeciesリストが変わるためその調整をする
    my_dic = my_dic_ppm

    plt.figure(1,dpi=150)
    for j in species_displayed:
        plt.plot(iteration[my_dic[j],:],alt_list[1:-1],label = j,color=species_color[j])
    plt.xscale("log")
    plt.legend()
    plt.title('initial ->' + proper_time+' species profile')
    plt.xlabel("volume mixing ratio")
    plt.ylabel("Altitude[km]")

    if same_figure == False:
        plt.figure(2,dpi=150)
        plt.title('initial species profile')
    for j in species_displayed:
        plt.plot(initial_mat[my_dic[j],:],alt_list[1:-1],label = j,color=species_color[j],linestyle=':')
    plt.xscale("log")
    #plt.legend()
    plt.xlabel("volume mixing ratio")
    plt.ylabel("Altitude[km]")
    plt.xlim((1e-10,1))

    #return time_list


###############################################################################
#this function plot mixing ratio profile for each species
    #input: file name
    #input: iteration number you wanna see
    #input: species list array you want to display ['H2O', 'CO',,,]
    #input: whether plot in same figure or not, default is false, so plot in a different figure
def plot_mixing_compare(file_ppm_alt, iteration_number1,iteration_number2, species_displayed, same_figure=True):
    file = file_ppm_alt
    h5file = h5py.File(file,"r") #read a h5 file
    time_list = h5file["n_current/timelist"].value #time list
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list

    iteration1 = h5file["n_current/iter_"+str(iteration_number1)].value #iteration n_current
    iteration2 = h5file["n_current/iter_"+str(iteration_number2)].value #iteration n_current

    all_species_num = all_species_num_ppm
    #number density to mixing ratio
    iteration_temp = cp.deepcopy(iteration1)
    for i in range(len(iteration1[0,:])):
        for j in all_species_num:
            iteration1[j,i] = iteration_temp[j,i]/sum(iteration_temp[all_species_num,i])
    #number density to mixing ratio
    iteration_temp2 = cp.deepcopy(iteration2)
    for i in range(len(iteration2[0,:])):
        for j in all_species_num:
            iteration2[j,i] = iteration_temp2[j,i]/sum(iteration_temp2[all_species_num,i])

    #iteration number to time unit
    time_second1 = time_list[iteration_number1-1]
    proper_time1 = time_transfer(time_second1)

    time_second2 = time_list[iteration_number2-1]
    proper_time2 = time_transfer(time_second2)


    #one_year fileとppm fileとconverged fileによって対応するspeciesリストが変わるためその調整をする
    my_dic = my_dic_ppm

    plt.figure(1,dpi=150)
    for j in species_displayed:
        plt.plot(iteration2[my_dic[j],:],alt_list[1:-1],label = species_list[my_dic[j]],color=species_color[j])

    plt.legend()
    for j in species_displayed:
        plt.plot(iteration1[my_dic[j],:],alt_list[1:-1],label = species_list[my_dic[j]],color=species_color[j],linestyle=':')
    plt.xscale("log")
    plt.xscale("log")
    #plt.xlim((1e-9,1e-1))
    plt.title(proper_time1+'->'+proper_time2)
    plt.xlabel("volume mixing ratio[]")
    plt.ylabel("Altitude[km]")

    #return species_list
################################### time change ##########################################################################
#This function creates animation of time shift of abundances of all species
#input: filename, species_displayed
#output: animation of plot
#import matplotlib.animation as animation
def timechange_mix(filename, species_displayed,savefilename):
    my_dic = my_dic_ppm #one_year fileとppm fileとconverged fileによって対応するspeciesリストが変わるためその調整をする
    
    #read data from HDF file
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    time_list = h5file["n_current/timelist"].value #time list
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list
    time_array = range(1,int(len(time_list)))
    
    fig = plt.figure(dpi=100)
    #images=[]
    #plt.legend()
    #Here it iterates over all time
    ani = animation.FuncAnimation(fig, update_mix, fargs=(file,time_array,time_list,alt_list,species_displayed,species_list,my_dic),interval = 100, frames = 50)
    #fig.show()
    ani.save(savefilename, writer='imagemagick', fps=4)
    
"""
    for t in range(1,int(len(time_list)),50):
        print(t)
        iteration = h5file["n_current/iter_"+str(t)].value #iteration n_current
        
        #def temporary array to calc mixing ratio
        iteration_temp = cp.deepcopy(iteration)#def temporary array to calc mixing ratio
        for i in range(len(iteration[0,:])):
            for j in all_species_num:
                iteration[j,i] = iteration_temp[j,i]/sum(iteration_temp[all_species_num,i])
        
        #transfer second to proper time unit
        time_second = time_list[t-1]
        proper_time = time_transfer(time_second)
        
        #plot here
        for j in species_displayed:
            image = plt.plot(iteration[my_dic[j],:],alt_list[1:-1],label = species_list[my_dic[j]],color=species_color[j])
        #plt.legend()
        plt.xscale("log")
        plt.xscale("log")
        plt.title(proper_time)
        plt.xlabel("volume mixing ratio")
        plt.ylabel("Altitude[km]")
        images.append(image)
        #time.sleep(5)
        #plt.cla()
    ani = animation.ArtistAnimation(fig, images, interval=1000)
    plt.show()
    """
#This is for number density
#timechange_num("CO2_12e8_Tmod/CO2_1000mbar_Tsurf_270_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5",all_species,"CO2_1000mbar_Tsurf_270_Tmod.gif")
def timechange_num(filename, species_displayed,savefilename):
    my_dic = my_dic_ppm #one_year fileとppm fileとconverged fileによって対応するspeciesリストが変わるためその調整をする
    
    #read data from HDF file
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    time_list = h5file["n_current/timelist"].value #time list
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    species = h5file["n_current/species"].value #species type: object containing string element
    species_list = []
    for i in range(len(species)):
        species_list.append(species[i]) #add string element to the list
    time_array = range(1,int(len(time_list)))
    
    fig = plt.figure(dpi=100)
    #images=[]
    #plt.legend()
    #Here it iterates over all time
    ani = animation.FuncAnimation(fig, update_num, fargs=(file,time_array,time_list,alt_list,species_displayed,species_list,my_dic),interval = 100, frames = 50)
    #fig.show()
    ani.save(savefilename, writer='imagemagick', fps=4)
    
#This function helpes timechange of mixing ratio function above
def update_mix(i,filename, time_array,time_list,alt_list,species_displayed,species_list,my_dic):
    plt.clf()
    t = time_array[i*25]
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    iteration = h5file["n_current/iter_"+str(t)].value #iteration n_current
    
    #def temporary array to calc mixing ratio
    iteration_temp = cp.deepcopy(iteration)#def temporary array to calc mixing ratio
    for i in range(len(iteration[0,:])):
        for j in all_species_num:
            iteration[j,i] = iteration_temp[j,i]/sum(iteration_temp[all_species_num,i])
    
    #transfer second to proper time unit
    time_second = time_list[t-1]
    proper_time = time_transfer(time_second)
    
    #plot here
    for j in species_displayed:
        plt.plot(iteration[my_dic[j],:],alt_list[1:-1],label = species_list[my_dic[j]],color=species_color[j])
    plt.xscale("log")
    plt.xscale("log")
    plt.xlim((1e-18,1))
    plt.title(proper_time)
    plt.xlabel("volume mixing ratio")
    plt.ylabel("Altitude[km]")
    plt.legend()
    plt.draw()
#This function helpes timechange of number density function above
def update_num(i,filename, time_array,time_list,alt_list,species_displayed,species_list,my_dic):
    plt.clf()
    t = time_array[i*50]
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    iteration = h5file["n_current/iter_"+str(t)].value #iteration n_current
    
    
    #transfer second to proper time unit
    time_second = time_list[t-1]
    proper_time = time_transfer(time_second)
    
    #plot here
    for j in species_displayed:
        plt.plot(iteration[my_dic[j],:],alt_list[1:-1],label = species_list[my_dic[j]],color=species_color[j])
    plt.xscale("log")
    plt.xscale("log")
    plt.xlim((1,1e20))
    plt.title(proper_time)
    plt.xlabel(r"number density$[cm^{-3}]$")
    plt.ylabel("Altitude[km]")
    plt.legend()
    plt.draw()
################################################### timechange #####################################################



#this function calc CO2 pressure from a file like'ppm_alt'
    #input string file name
    #input iteration number you want to display
    #output: initial condition CO2 pressure(mbar), iteration?'s CO2 pressure(mbar)
def calc_CO2_pressure_file(file_any,iteration_number):
    file = file_any
    h5file = h5py.File(file,"r") #read a h5 file
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    initial_mat = h5file["n_current/init"].value #initial n_current
    iteration = h5file["n_current/iter_"+str(iteration_number)].value #iteration n_current

    G = 6.674e-11
    M = 0.64171e24
    p_initial = 0
    p_iter = 0
    alt = [(j+1)*2*1000 for j in range(99)]
    for i in range(99):
        p_initial += initial_mat[22,i]*2e5*44/(6e23)*1e-3*G*M/((3389500+alt[i])**2) #N/cm^2

    for i in range(99):
        p_iter += iteration[22,i]*2e5*44/(6e23)*1e-3*G*M/((3389500+alt[i])**2) #N/cm^2

    return p_initial*100, p_iter*100 #mbar

#this function calc CO2 pressure from CO2 number density array
    #input array of CO2 num density with respect to altitude(200km)
    #output:  CO2 pressure(mbar)
def calc_CO2_pressure_array(density_array):
    G = 6.674e-11
    M = 0.64171e24
    p = 0
    alt = [(j+1)*2*1000 for j in range(len(density_array))]
    for i,val in enumerate(density_array):
        p += val*2e5*44/(6e23)*1e-3*G*M/((3389500+alt[i])**2) #N/cm^2
    return p*100 #mbar

def calc_O2_pressure_array(density_array):
    G = 6.674e-11
    M = 0.64171e24
    p = 0
    alt = [(j+1)*2*1000 for j in range(len(density_array))]
    for i,val in enumerate(density_array):
        p += val*2e5*32/(6e23)*1e-3*G*M/((3389500+alt[i])**2) #N/cm^2
    return p*100 #mbar

def calc_CO_pressure_array(density_array):
    G = 6.674e-11
    M = 0.64171e24
    p = 0
    alt = [(j+1)*2*1000 for j in range(len(density_array))]
    for i,val in enumerate(density_array):
        p += val*2e5*28/(6e23)*1e-3*G*M/((3389500+alt[i])**2) #N/cm^2
    return p*100 #mbar

def calc_pressure_array(density_array, molarmass):
    G = 6.674e-11
    M = 0.64171e24
    p = 0
    alt = [(j+1)*2*1000 for j in range(len(density_array))]
    for i,val in enumerate(density_array):
        p += val*2e5*molarmass/(6e23)*1e-3*G*M/((3389500+alt[i])**2) #N/cm^2
    return p*100 #mbar

#function second to proper time unit(second, hour,year,million year, Gyr)
def time_transfer(second):
    output = second
    if second<1:
        output = round(output,3)
        return str(output)+'seconds'
    elif second <3600:#1 hour = 3600s
        output = round(output)
        return str(output)+'seconds'
    elif second<86400:# =1 day
        output = int(round(second/3600))
        return str(output)+'hour'
    elif second<31536000:#=1 year
        output = int(round(second/86400))
        return str(output)+'days'
    elif second<31536000e2:#=1 million year
        output = int(round(second/31536000))
        return str(output)+'years'
    elif second<31536000e4:#=1 million year
        output = int(round(second/31536000,-2))
        return str(output)+'years'
    elif second<31536000e5:#=1 million year
        output = int(round(second/31536000,-3))
        return str(output)+'years'
    elif second<31536000e6:#=1 million year
        output = int(round(second/31536000,-4))
        return str(output)+'years'
    elif second<31536000e9:#=1 billion year
        output = int(round(second/(31536000e6)))
        return str(output)+'million years'
    elif second<31536000e12:
        output = int(round(second/(31536000e9),1))
        return str(output)+'Gyrs'
    else:
        print('too large or minus')


########################################################################


flux_dic = {'CO2':0, 'O2':1, 'O3':2, 'H2':3, 'OH':4, 'HO2':5, 'H2O':6, 'H2O2':7, 'O':8, 'CO':9,
                         'O1D':10, 'H':11, 'N2':12, 'Ar':13, 'CO2pl':14, 'HOCO':15}


def flux_history(filename,iteration_num, species_displayed):
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    time_list = h5file["n_current/timelist"].value #time list
    fluxes = h5file["fluxes/flux_history"].value #16*99*1300
    plt.figure(dpi=150)
    for sp in species_displayed:
        flux_list = np.array(fluxes[flux_dic[sp],:,iteration_num])
        flux_list_positive = cp.deepcopy(flux_list)
        flux_list_negative = cp.deepcopy(flux_list)
        flux_list_positive[flux_list_positive<0] = 0
        flux_list_negative[flux_list_negative>=0] = 0
        flux_list_negative = -flux_list_negative
        plt.plot(flux_list_positive,alt_list[1:-1],label=sp, color=species_color[sp])
        plt.plot(flux_list_negative,alt_list[1:-1],label=sp, color=species_color[sp],linestyle='-.')
        #plt.plot(fluxes[flux_dic[sp],:,iteration_num],alt_list[1:-1],label=sp, color=species_color[sp])

    #iteration number to time unit
    time_second = time_list[iteration_num-1]
    proper_time = time_transfer(time_second)


    plt.xscale("log")
    plt.legend()
    plt.title('Flux'+proper_time)
    plt.xlim(1e5,1e12)
    plt.xlabel("Flux[$cm^{-2}s^{-1}$]")
    plt.ylabel("Altitude[km]")

#flux_history_compare("CO2_12e8_Tmod/CO2_300mbar_Tsurf_250_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5",1,1000,['H2'])
def flux_history_compare(filename,iteration_num1,iteration_num2,species_displayed):
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    time_list = h5file["n_current/timelist"].value #time list
    fluxes = h5file["fluxes/flux_history"].value #16*99*1300

    plt.figure(dpi=150)
    for sp in species_displayed:
        flux_list = np.array(fluxes[flux_dic[sp],:,iteration_num2])
        flux_list_positive = cp.deepcopy(flux_list)
        flux_list_negative = cp.deepcopy(flux_list)
        flux_list_positive[flux_list_positive<0] = 0
        flux_list_negative[flux_list_negative>=0] = 0
        flux_list_negative = -flux_list_negative
        plt.plot(flux_list_positive,alt_list[1:-1],label=sp+' up', color=species_color[sp])
        plt.plot(flux_list_negative,alt_list[1:-1],label=sp+' down', color=species_color[sp],linestyle='-.')
    plt.legend()
    for sp in species_displayed:
        flux_list = np.array(fluxes[flux_dic[sp],:,iteration_num1])
        flux_list_positive = cp.deepcopy(flux_list)
        flux_list_negative = cp.deepcopy(flux_list)
        flux_list_positive[flux_list_positive<0] = 0
        flux_list_negative[flux_list_negative>=0] = 0
        flux_list_negative = -flux_list_negative
        plt.plot(flux_list_positive,alt_list[1:-1],label=sp, color=species_color[sp],alpha=0.3)
        plt.plot(flux_list_negative,alt_list[1:-1],label=sp, color=species_color[sp],linestyle='-.',alpha = 0.3)

    #iteration number to time unit
    time_second1 = time_list[iteration_num1-1]
    proper_time1 = time_transfer(time_second1)
    time_second2 = time_list[iteration_num2-1]
    proper_time2 = time_transfer(time_second2)

    plt.xscale("log")
    #plt.legend()
    plt.title('Flux'+proper_time1+' to '+proper_time2)
    plt.xlim(1e5,1e12)
    plt.xlabel("Flux[$cm^{-2}s^{-1}$]")
    plt.ylabel("Altitude[km]")



#################################################################################
#this function compare flux history between fluxes of each file

def flux_history_compare_files(filename1,filename2,iteration_num1,iteration_num2,species_displayed):
    file1 = filename1
    h5file1 = h5py.File(file1,"r") #read a h5 file
    alt_list = h5file1["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    time_list = h5file1["n_current/timelist"].value #time list
    fluxes1 = h5file1["fluxes/flux_history"].value #16*99*1300

    file2 = filename2
    h5file2 = h5py.File(file2,"r") #read a h5 file
    fluxes2 = h5file2["fluxes/flux_history"].value #16*99*1300


    plt.figure(dpi=150)
    for sp in species_displayed:
        flux_list = np.array(fluxes2[flux_dic[sp],:,iteration_num2])
        flux_list_positive = cp.deepcopy(flux_list)
        flux_list_negative = cp.deepcopy(flux_list)
        flux_list_positive[flux_list_positive<0] = 0
        flux_list_negative[flux_list_negative>=0] = 0
        flux_list_negative = -flux_list_negative
        plt.plot(flux_list_positive,alt_list[1:-1],label=sp+' up', color=species_color[sp],linewidth=1.5)
        plt.plot(flux_list_negative,alt_list[1:-1],label=sp+' down', color=species_color[sp],linestyle='-.',linewidth=1.5)
    plt.legend()
    for sp in species_displayed:
        flux_list = np.array(fluxes1[flux_dic[sp],:,iteration_num1])
        flux_list_positive = cp.deepcopy(flux_list)
        flux_list_negative = cp.deepcopy(flux_list)
        flux_list_positive[flux_list_positive<0] = 0
        flux_list_negative[flux_list_negative>=0] = 0
        flux_list_negative = -flux_list_negative
        plt.plot(flux_list_positive,alt_list[1:-1],label=sp, color=species_color[sp],linewidth=0.8)
        plt.plot(flux_list_negative,alt_list[1:-1],label=sp, color=species_color[sp],linestyle='-.',linewidth=0.8)

    #iteration number to time unit
    time_second1 = time_list[iteration_num1-1]
    proper_time1 = time_transfer(time_second1)
    time_second2 = time_list[iteration_num2-1]
    proper_time2 = time_transfer(time_second2)

    plt.xscale("log")
    #plt.legend()
    plt.title('Flux'+' initial'+' to '+proper_time2)
    plt.xlim(1e5,1e12)
    plt.xlabel("Flux[$cm^{-2}s^{-1}$]")
    plt.ylabel("Altitude[km]")


reactions = [
#photodissociation
r'$CO2 + h\nu -> CO + O$',
r'$CO2 + h\nu -> CO + O1D$',
r'$O2 + h\nu -> O + O$',
r'$O2 + h\nu -> O + O1D$',
r'$O3 + h\nu -> O2 + O$',
r'$O3 + h\nu -> O2 + O1D$',#6
r'$O3 + h\nu -> O + O + O$',
r'$H2 + h\nu -> H + H$',
r'$OH + h\nu -> O + H$',
r'$OH + h\nu -> O1D + H$',#10
r'$HO2 + h\nu -> OH + O$',
r'$H2O + h\nu -> H + OH$',#12
r'$H2O + h\nu -> H2 + O1D$',
r'$H2O + h\nu -> H + H + O$',
r'$H2O2 + h\nu -> OH + OH$',#15
r'$H2O2 + h\nu -> HO2 + H$',#16
r'$H2O2 + h\nu -> H2O + O1D$',#17

#recombination of O
r'$O + O + M -> O2 + M$',
r'$O + O2 + N2 -> O3 + N2$',#19
r'$O + O2 + CO2 -> O3 + CO2$',#20
r'$O + O3 + M -> O2 + O2$',#21
r'$O + CO + M -> CO2 + M$',#22

#O1D attack
r'$O1D + O2 -> O + O2$',
r'$O1D + O3 -> O2 + O2$',#24
r'$O1D + O3 -> O + O + O2$',
r'$O1D + H2 -> H + OH$',
r'$O1D + CO2 -> O + CO2$',
r'$O1D + H2O -> OH + OH$',#28

#loss of H2
r'$H2 + O -> OH + H$',#29
r'$OH + H2 -> H2O + H$',#30

#recombination of H
r'H + H + CO$_{2}$ $\rightarrow$ H$_{2}$ + CO$_{2}$',#31
r'$H + OH + CO2 -> H2O + CO2$',
r'$H + HO2 -> OH + OH$',
r'$H + HO2 -> H2O + O1D$',
r"H + HO$_{2}$ $\rightarrow$ H$_{2}$ + O$_{2}$",#35
r'$H + H2O2 -> HO2 + H2$',
r'$H + H2O2 -> H2O + OH$',

#Interconversion of odd H
r'H + O$_{2}$ $\rightarrow$ HO$_{2}$',#38
r'$H + O3 -> OH + O2$',#39
r'$O + OH -> O2 + H$',
r'$O + HO2 -> OH + O2$',#41
r'$O + H2O2 -> OH + HO2$',
r'$OH + OH -> H2O + O$',
r'$OH + OH -> H2O2$',
r'$OH + O3 -> HO2 + O2$',#45
r'$OH + HO2 -> H2O + O2$',
r'$OH + H2O2 -> H2O + HO2$',
r'$HO2 + O3 -> OH + O2 + O2$',
r'$HO2 + HO2 -> H2O2 + O2$',
r'$HO2 + HO2 + M -> H2O2 + O2 + M$',#50

#CO2 recombination due to odd H(with HOCO intermidiate)
r'$CO + OH -> CO2 + H$',#51
r'$OH + CO -> HOCO$',
r'$HOCO + O2 -> HO2 + CO2$',

#CO2+ attack on moelcular hydrogen
r'$CO2^+ + H2 -> CO2 + H + H$'
]

O_array_right=[23,25,27,43]
O_array=[18,19,20,21,22,29,40,41,42]
O3_del_array=[5,6,7,21,24,25,39,45,48]
O3_pro_array=[20]
H_array=[8,9,10,12,14,16,26,29,30,31,32,33,34,35,36,37,38,39,40,51,54]
H_array_delete=[31,32,33,34,35,36,37,38,39]
H_array_produce=[8,9,10,12,14,16,26,29,30,40,51,54]
H_pro_dominant=[40,51,54]
H2_array_delete=[8,26,29,30,54]
H2_array_produce=[13,31,35,36]
dominant=[19,20,21,40,39,35,54]
CO2_left=[1,2,27]
CO_array=[1,2,22,51,52]
O2_pro = [5,6,18,21,23,24,25,35,39,40,41,45,46,48,49,50]
O2_pro_dominant=[5,6,18,35,39,40,41,46,49]
O2_del=[3,4,19,20,23,38]
##############################################################################
    
#reaction_rates_compare("CO2_12e8_Tmod/CO2_300mbar_Tsurf_250_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5",1,2000,H2_array_produce)
#this function plots the reaction rate with respect to altitude at specific time(iteration_num)
#input: filename; string
#       iteration_num;int
#       reaction_displayed; array filled with int number[1~54]
def reaction_rates_compare(filename,iteration_num1,iteration_num2, reaction_displayed):
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    time_list = h5file["n_current/timelist"].value #time list
    reaction_rates = h5file["rates/reaction_rates_history"].value #54*99*1300
    plt.figure(dpi=150)
    count=0
    for i in reaction_displayed:
        plt.plot(reaction_rates[i-1,1:len(alt_list[2:-1])+1,iteration_num2],alt_list[2:-1],label=reactions[i-1],color=color_list[count%8])
        count += 1
    plt.legend(prop={'size':10})
    count=0
    for i in reaction_displayed:
        plt.plot(reaction_rates[i-1,1:len(alt_list[2:-1])+1,iteration_num1],alt_list[2:-1],label=reactions[i-1],linestyle=':' ,color=color_list[count%8])
        count+=1
    #iteration number to time unit
    time_second1 = time_list[iteration_num1-1]
    proper_time1 = time_transfer(time_second1)
    time_second2 = time_list[iteration_num2-1]
    proper_time2 = time_transfer(time_second2)

    plt.xscale("log")
    #plt.legend()
    plt.title('reaction rates '+proper_time1+' -> '+proper_time2)
    plt.xlabel(r"Reaction rates$[cm^{-3}s^{-1}]$")
    plt.ylabel("Altitude[km]")
    plt.xlim((1e-3,1e7))


#this function plots the reaction rate of 2 files with respect to altitude at specific time(iteration_num)
#input: filename1; string: initial reaction rates file(nochange)
#       filename2; string: next reaction rates file(after changed)
#       iteration_num1,2: iteration number you want to compare in each file
#       reaction_displayed; array filled with int number[1~54]
def reaction_rates_compare_files(filename1,filename2,iteration_num1, iteration_num2,reaction_displayed,title_memo=""):
    #file1
    file1 = filename1
    h5file1 = h5py.File(file1,"r") #read a h5 file
    alt_list = h5file1["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    time_list = h5file1["n_current/timelist"].value #time list
    reaction_rates1 = h5file1["rates/reaction_rates_history"].value #54*99*1300
    #file2
    file2 = filename2
    h5file2 = h5py.File(file2,"r") #read a h5 file
    reaction_rates2 = h5file2["rates/reaction_rates_history"].value #54*99*1300

    plt.figure(dpi=150)
    count=0
    for i in reaction_displayed:
        plt.plot(reaction_rates2[i-1,:,iteration_num2],alt_list[1:-1],label=reactions[i-1],color=reactions_color[count%16])
        count += 1
    plt.legend(prop={'size':8})
    count=0
    for i in reaction_displayed:
        plt.plot(reaction_rates1[i-1,:,iteration_num1],alt_list[1:-1],label=reactions[i-1],linestyle=':' ,color=reactions_color[count%16])
        count+=1
    #iteration number to time unit
    time_second1 = time_list[iteration_num1-1]
    proper_time1 = time_transfer(time_second1)
    time_second2 = time_list[iteration_num2-1]
    proper_time2 = time_transfer(time_second2)

    plt.xscale("log")
    #plt.legend()
    plt.title(title_memo)
    plt.xlabel(r"reaction rates$[cm^{-3}s^{-1}]$")
    plt.ylabel("Altitude[km]")
    plt.xlim((1e-16,1e7))
    
def reaction_rates_columm_history(filename, reaction_displayed):
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    time_list = h5file["n_current/timelist"].value/3.15e7 #time list
    reaction_rates = h5file["rates/reaction_rates_history"].value #54*99*1300
    plt.figure(dpi=200)
    col_reaction_hist=[]
    for i in reaction_displayed:
        for j in range(len(time_list)):
            col_reaction_hist.append(sum(reaction_rates[i-1,:,j]*2e5))
        plt.plot(time_list,col_reaction_hist,label=reactions[i-1])
        col_reaction_hist=[]
    plt.legend(prop={'size':4})
    #iteration number to time unit
    plt.xscale("log")
    plt.yscale("log")
    #plt.legend()
    #plt.title('reaction rates '+proper_time1+' -> '+proper_time2)
    plt.ylabel(r"columm reaction rates$[cm^{-2}s^{-1}]$")
    plt.xlabel("time[years]")
    #plt.ylabel("Altitude[km]")
    #plt.xlim((1e-5,1e6))


reactions_H2CO = [
#photodissociation
r'$CO2 + h\nu -> CO + O$',
r'$CO2 + h\nu -> CO + O1D$',
r'$O2 + h\nu -> O + O$',
r'$O2 + h\nu -> O + O1D$',
r'$O3 + h\nu -> O2 + O$',
r'$O3 + h\nu -> O2 + O1D$',#6
r'$O3 + h\nu -> O + O + O$',
r'$H2 + h\nu -> H + H$',
r'$OH + h\nu -> O + H$',
r'$OH + h\nu -> O1D + H$',#10
r'$HO2 + h\nu -> OH + O$',
r'$H2O + h\nu -> H + OH$',#12
r'$H2O + h\nu -> H2 + O1D$',
r'$H2O + h\nu -> H + H + O$',
r'$H2O2 + h\nu -> OH + OH$',#15
r'$H2O2 + h\nu -> HO2 + H$',#16
r'$H2O2 + h\nu -> H2O + O1D$',#17

# H2CO related reaction!!
r'$HCO + h\nu -> H+CO$',
r'$H2CO + h\nu -> H+HCO$',
r'$H2CO + h\nu-> H2+CO$',

r'$H+CO+M -> HCO+M$',
r'$H+HCO -> H2+CO$',
r'$HCO+HCO -> H2CO+CO$',
r'$OH+HCO -> H2O+CO$',
r'$O+HCO -> H+CO2$',
r'$O+HCO -> OH+CO$',
r'$O2+HCO -> HO2+CO$',
r'$HO2+HCO -> H2O2+CO$',
r'$H+H2CO -> H2+HCO$',
r'$OH+H2CO -> H2O+HCO$',

#recombination of O
r'$O + O + M -> O2 + M$',
r'$O + O2 + N2 -> O3 + N2$',#19
r'$O + O2 + CO2 -> O3 + CO2$',#20
r'$O + O3 + M -> O2 + O2$',#21
r'$O + CO + M -> CO2 + M$',#22

#O1D attack
r'$O1D + O2 -> O + O2$',
r'$O1D + O3 -> O2 + O2$',#24
r'$O1D + O3 -> O + O + O2$',
r'$O1D + H2 -> H + OH$',
r'$O1D + CO2 -> O + CO2$',
r'$O1D + H2O -> OH + OH$',#28

#loss of H2
r'$H2 + O -> OH + H$',#29
r'$OH + H2 -> H2O + H$',#30

#recombination of H
r'H + H + CO$_{2}$ $\rightarrow$ H$_{2}$ + CO$_{2}$',#31
r'$H + OH + CO2 -> H2O + CO2$',
r'$H + HO2 -> OH + OH$',
r'$H + HO2 -> H2O + O1D$',
r"H + HO$_{2}$ $\rightarrow$ H$_{2}$ + O$_{2}$",#35
r'$H + H2O2 -> HO2 + H2$',
r'$H + H2O2 -> H2O + OH$',

#Interconversion of odd H
r'H + O$_{2}$ $\rightarrow$ HO$_{2}$',#38
r'$H + O3 -> OH + O2$',#39
r'$O + OH -> O2 + H$',
r'$O + HO2 -> OH + O2$',#41
r'$O + H2O2 -> OH + HO2$',
r'$OH + OH -> H2O + O$',
r'$OH + OH -> H2O2$',
r'$OH + O3 -> HO2 + O2$',#45
r'$OH + HO2 -> H2O + O2$',
r'$OH + H2O2 -> H2O + HO2$',
r'$HO2 + O3 -> OH + O2 + O2$',
r'$HO2 + HO2 -> H2O2 + O2$',
r'$HO2 + HO2 + M -> H2O2 + O2 + M$',#50

#CO2 recombination due to odd H(with HOCO intermidiate)
r'$CO + OH -> CO2 + H$',#51
r'$OH + CO -> HOCO$',
r'$HOCO + O2 -> HO2 + CO2$',

#CO2+ attack on moelcular hydrogen
r'$CO2^+ + H2 -> CO2 + H + H$'
]


#this function plots the reaction rate with respect to altitude
#input: filename; string
#       reaction_displayed; array filled with int number[1~54]
def reaction_rates_conv(filename,reaction_displayed):
    file = filename
    h5file = h5py.File(file,"r") #read a h5 file
    alt_list = h5file["n_current/alt"].value #altitude list
    alt_list = alt_list/1e5 #cm to km
    reaction_rates = h5file["rates/reaction_rates"].value #54*99*1300

    plt.figure(dpi=150)
    for i in reaction_displayed:
        plt.plot(reaction_rates[i-1,1:len(alt_list[2:-1])+1],alt_list[2:-1],label=reactions_H2CO[i-1])
        print(reactions[i-1] + ":" + str(sum(reaction_rates[i-1,:])*2e5))
    
    plt.xscale("log")
    plt.legend()
    plt.xlabel(r"reaction rates$[cm^{-3}s^{-1}]$")
    plt.ylabel("Altitude[km]")
    plt.xlim((1e-8,1e7))
    plt.ylim((0,60))

