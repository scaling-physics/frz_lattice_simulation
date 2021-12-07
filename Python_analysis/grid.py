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

# grid = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=1, max_rows=10)
# grid1 = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=1000, max_rows=10)
# grid2 = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=9500000, max_rows=10)

# grid = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid.txt', skiprows=1, max_rows=10)
# grid1 = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid.txt', skiprows=1000, max_rows=10)
# grid2 = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid.txt', skiprows=9500000, max_rows=10)
# #%%
# grid_full = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_rand_J_1_alpha_0.txt', skiprows=1, max_rows=950000)
# #%%

# def color(a,b,c):
#     a1=a/255
#     b1=b/255
#     c1=c/255
#     d=np.array((a1,b1,c1))
#     return d



# colors=np.zeros([50*50,3])
# for row in range(10):
#     for x in range(2500):
#         if grid2[row,x]==0:
#             colors[x]=color(0,0,0)
#         elif grid2[row,x]==1:
#             colors[x]=[195/255,192/255,192/255]
#         elif grid2[row,x]==2:
#              colors[x]=[0, 1, 0]
#         elif grid2[row,x]==3:
#              colors[x]=[1,0,0]
#         elif grid2[row,x]==4:
#          colors[x]=[0,0,1]

    
#     plt.figure(figsize=(20,15))
#     hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors,do_plot=True)
#     centers_x = hex_centers[:, 0]
#     centers_x = hex_centers[:, 1]
#     #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
    
#     plt.show()

# #%%
# def num_particles(grid,i):
#     num =0
#     j=0
#     while j<2500:
#         if grid[i,j]>0:
#             num+=1
#         j+=1
#     return num

# num_part=[]
# def diffuse_part(grid,i):
#     num =0
#     j=0
#     while j<2500:
#         if grid[i,j]==1:
#             num+=1
#         j+=1
#     return num
# dif_part=[]
# particles_in_clusters=[]
    
# for i in range(950000):
#     q=0
#     if i%10==0:
#         num_part=np.append(num_part,num_particles(grid_full, i))
#         dif_part=np.append(dif_part,diffuse_part(grid_full, i))
#         particles_in_clusters=np.append(particles_in_clusters,(num_part[q]-dif_part[q]))
#         print(num_particles(grid_full, i))
#         q+=1

# MC = np.arange(95000)
# #%%
# plt.figure()
# plt.plot(MC,num_part)
# plt.plot(MC,dif_part)
# plt.plot(MC,particles_in_clusters)
# plt.show()
#%%
# grid_full = np.loadtxt('/home/hannohennighausen/Documents/frz_lattice_model/grid_J_large_half_0.txt', skiprows=1)
# for x in range(2500):
#     if grid_full[x]==0:
#         colors[x]=color(0,0,0)
#     elif grid_full[x]==1:
#         colors[x]=[195/255,192/255,192/255]
#     elif grid_full[x]==2:
#          colors[x]=[0, 1, 0]
#     elif grid_full[x]==3:
#          colors[x]=[1,0,0]
#     elif grid_full[x]==4:
#      colors[x]=[0,0,1]

    
# plt.figure(figsize=(20,15))
# hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors,do_plot=True)
# centers_x = hex_centers[:, 0]
# centers_x = hex_centers[:, 1]
# #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')

# plt.show()
#%%
J1=[1,2]
alpha1=[0,0.2,0.5]
density=0.2
#%%
counter=0
for J in J1:
    for alpha in alpha1:
        labels = np.loadtxt(f'/home/hannohennighausen/Documents/frz_lattice_model/outputlabels_{J}_{alpha}_{density}_0.txt',dtype=int,skiprows=1)
        
        labels_name = np.zeros(50*50)
        
        num_of_clusters=[]
        part_in_clusters=[]
        part_dif=[]
        avg_cluster_size=[]
        for i in range(1000):
            
            num_of_clusters=np.append(num_of_clusters,np.max(labels[i*3]))
            q=0
            p=0
            for j in range(500):
                if labels[i*3+2,j]>1:
                    
                    q+=1
                if labels[i*3+2,j]==1:
                    p+=1
            part_in_clusters=np.append(part_in_clusters,q)
            part_dif=np.append(part_dif,p)
        
        
        plt.figure()
        plt.plot(np.arange(1000),num_of_clusters,label = 'num of clusters')
        plt.plot(np.arange(1000),part_in_clusters,label = 'particles in clusters')
        plt.plot(np.arange(1000),part_dif,label = 'diffuse particles')
        plt.xlim((0,1000))
        plt.vlines([100,500,950], 0, 400, colors='r',linestyles='dashed')
        plt.xlabel('MC_steps *10')
        plt.title(f'J={J}, alpha={alpha},density={density}')
        plt.legend()
        plt.savefig(f'{counter}_{J}_{alpha}_{density}.svg')
        plt.show()
        
        
        steps=[100,500,950]
        for step in steps:
            cluster_size = np.bincount(labels[step*3])
            cluster_size_distribution = np.bincount(cluster_size[1:])    
            plt.figure()
            plt.bar(np.arange(cluster_size_distribution.size), cluster_size_distribution)
            plt.suptitle('cluster size distribution')
            plt.title(f'J={J}, alpha={alpha},density={density}')
            plt.savefig(f'{counter}_size_dist_at_{step}_{J}_{alpha}_{density}.svg')
            plt.show()
            
            colors_labels = np.zeros([50*50,3])
            coloring_labels = np.zeros([50*50,3])
            i=0
            j=0
            Cluster_num = np.max(labels[step*3])+20
            for x in labels[step*3+1]:
                if labels[step*3,j]==0:
                    coloring_labels[x]=[195/255,192/255,192/255] 
                else:
                    coloring_labels[x]=[1,0.5+labels[0,j]/Cluster_num,0]
                j+=1


            for x in labels[step*3+1]:

                # else:
                #     coloring_labels[x]=[195/255,192/255,192/255]
                    # labels_name[x]=labels[0,i]
                if labels[step*3+2,i]==1:
                    colors_labels[x]=[195/255,192/255,192/255]  
                elif labels[step*3+2,i]==2:
                    colors_labels[x]=[0, 1, 0]
                elif labels[step*3+2,i]==3:
                    colors_labels[x]=[1,0,0]
                elif labels[step*3+2,i]==4:
                    colors_labels[x]=[0,0,1]
                i+=1
                

            plt.figure(figsize=(20,15))
            hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors_labels,do_plot=True)
            centers_x = hex_centers[:, 0]
            centers_x = hex_centers[:, 1]
            #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
            plt.savefig(f'{counter}_labels_{step}_{J}_{alpha}_{density}.svg')

            plt.show()

            plt.figure(figsize=(20,15))
            hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=coloring_labels,do_plot=True)
            centers_x = hex_centers[:, 0]
            centers_x = hex_centers[:, 1]
            #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
            plt.savefig(f'{counter}_ori_{step}_{J}_{alpha}_{density}.svg')

            plt.show()
        counter+=1
#%%

i=0
j=0
Cluster_num = np.max(labels[999*3])+20
for x in labels[999*3+1]:
    if labels[999*3,j]==0:
        coloring_labels[x]=[195/255,192/255,192/255] 
    else:
        coloring_labels[x]=[1,0.5+labels[0,j]/Cluster_num,0]
    j+=1


for x in labels[999*3+1]:

    # else:
    #     coloring_labels[x]=[195/255,192/255,192/255]
        # labels_name[x]=labels[0,i]
    if labels[999*3+2,i]==1:
        colors_labels[x]=[195/255,192/255,192/255]  
    elif labels[999*3+2,i]==2:
        colors_labels[x]=[0, 1, 0]
    elif labels[999*3+2,i]==3:
        colors_labels[x]=[1,0,0]
    elif labels[999*3+2,i]==4:
        colors_labels[x]=[0,0,1]
    i+=1
    

plt.figure(figsize=(20,15))
hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors_labels,do_plot=True)
centers_x = hex_centers[:, 0]
centers_x = hex_centers[:, 1]
#plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')

plt.show()

plt.figure(figsize=(20,15))
hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=coloring_labels,do_plot=True)
centers_x = hex_centers[:, 0]
centers_x = hex_centers[:, 1]
#plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')

plt.show()