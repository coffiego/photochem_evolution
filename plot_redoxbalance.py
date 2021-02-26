#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:52:06 2021

@author: shungokoyama
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df_redox = pd.read_csv("/Users/shungokoyama/Programming/result/photochemistry_c/redox_CHO.csv",delimiter=',')
CO = df_redox['CO']
HO = df_redox['HO']
HCO = df_redox['HCO']
feuv = df_redox['feuv']
CO = [co for euv, co in sorted(zip(feuv, CO))]
HO = [ho for euv, ho in sorted(zip(feuv, HO))]
HCO = [hco for euv, hco in sorted(zip(feuv, HCO))]
feuv = sorted(feuv)


plt.figure(dpi=150)
CO = np.array(CO,dtype=np.float)
HO = np.array(HO,dtype=np.float)
HCO = np.array(HCO,dtype=np.float)
feuv = np.array(feuv,dtype=np.float)

plt.plot(feuv, HO, color='dodgerblue',label=r"$\Phi_H / 2\Phi_O$",marker="x",fillstyle="none",linestyle="--",markersize=10)
plt.plot(feuv, CO, color='orange',label=r"$2\Phi_C / \Phi_O$",marker="^",fillstyle="none",linestyle="--",markersize=10)
plt.plot(feuv, HCO, color='k',label=r"$(\Phi_H+4\Phi_C) / 2\Phi_O$",marker="o",fillstyle="none",linestyle="--")
plt.xlabel('EUV',size=15)
plt.ylabel('Redox Balance',size=15)
plt.ylim((-0.1,1.1))
plt.xlim((16,2))
plt.tick_params(labelsize=13)
plt.legend(prop={'size':10.0})

