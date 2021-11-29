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


#grid = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_J_-1000_alpha_0.txt', skiprows=99899)
#grid = np.loadtxt('C:/Users/hanno/OneDrive/Documents/HannoTablet/Physics/Marburg/Murray_Praktikum/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=99899)

grid = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=1, max_rows=10)
grid1 = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=1000, max_rows=10)
grid2 = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=100000, max_rows=10)

#%%
grid_full = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=1, max_rows=950000)
#%%

def color(a,b,c):
    a1=a/255
    b1=b/255
    c1=c/255
    d=np.array((a1,b1,c1))
    return d



colors=np.zeros([50*50,3])
for row in range(10):
    for x in range(2500):
        if grid2[row,x]==0:
            colors[x]=color(0,0,0)
        elif grid2[row,x]==1:
            colors[x]=[195/255,192/255,192/255]
        elif grid2[row,x]==2:
             colors[x]=[0, 1, 0]
        elif grid2[row,x]==3:
             colors[x]=[1,0,0]
        elif grid2[row,x]==4:
         colors[x]=[0,0,1]

    
    plt.figure(figsize=(20,15))
    hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors,do_plot=True)
    centers_x = hex_centers[:, 0]
    centers_x = hex_centers[:, 1]
    #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
    
    plt.show()
#%%
grid_full = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=1, max_rows=950000)
def num_particles(grid,i):
    num =0
    j=0
    while j<2500:
        if grid[i,j]>0:
            num+=1
        j+=1
    return num

X=[]
    
    
for i in range(5*10**5):
    if i%1000==0:
        X=np.append(X,num_particles(grid_full, i))
        print(num_particles(grid_full, i))

MC = np.arange(95)

plt.figure()
plt.plot(MC,X)
plt.show()