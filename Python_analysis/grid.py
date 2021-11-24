#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:33:48 2021

@author: hannohennighausen
"""

import hexalattice.hexalattice as hex
import numpy as np
import matplotlib.pyplot as plt

grid = np.loadtxt('grid.txt', skiprows=1)
#%%





colors=np.zeros([50*50,3])
for row in range(100):
    for x in range(2500):
        if grid[19899+row,x]==0:
            colors[x]=[0.,0.,0.]
        elif grid[19899+row,x]==1:
            colors[x]=[0.3, 0.3, 0.3]
        else:
            colors[x]=[0.8,0.8,0.8]
    
    plt.figure(figsize=(20,15))
    hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors,do_plot=True)
    centers_x = hex_centers[:, 0]
    centers_x = hex_centers[:, 1]
    plt.savefig(f'grid_19899+{row}.svg')