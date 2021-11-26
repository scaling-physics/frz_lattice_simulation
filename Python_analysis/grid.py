#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:33:48 2021

@author: hannohennighausen
"""

import hexalattice.hexalattice as hex
import numpy as np
import matplotlib.pyplot as plt
from pathlib import  Path

p = Path('/home/hannohennighausen/Documents/frz_lattice_model')

#%%

grid = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_J_1_alpha_0.txt', skiprows=99899)
#%%

def color(a,b,c):
    a1=a/255
    b1=b/255
    c1=c/255
    d=np.array((a1,b1,c1))
    return d



colors=np.zeros([50*50,3])
for row in range(100):
    for x in range(2500):
        if grid[row,x]==0:
            colors[x]=color(0,0,0)
        elif grid[row,x]==1:
            colors[x]=[195/255,192/255,192/255]
        elif grid[row,x]==2:
             colors[x]=[0, 1, 0]
        elif grid[row,x]==3:
             colors[x]=[1,0,0]
        elif grid[row,x]==4:
         colors[x]=[0,0,1]

    
    plt.figure(figsize=(20,15))
    hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors,do_plot=True)
    centers_x = hex_centers[:, 0]
    centers_x = hex_centers[:, 1]
    plt.savefig(f'1_grid_19899+{row}.svg')
    plt.show()