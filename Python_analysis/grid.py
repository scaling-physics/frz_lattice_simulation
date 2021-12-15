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
import imageio as iio

#%%

color_scale=[[0,0,1],[0,1,0],[1,0,0],[1,0,0.6],[1,0.6,0],[0.6,1,0],[0,1,0.6],[0,0.6,1],[0.6,0,1],[1,1,0.6],[1,0.6,1],[0.6,1,1],[0,0,0.6],[0,0.6,0],[0.6,0,0],
             [0,0,0.2],[0,0.2,0],[0.2,0,0],[1,0,0.2],[1,0.2,0],[0.2,1,0],[0,1,0.2],[0,0.2,1],[0.2,0,1],[1,1,0.2],[1,0.2,1],[0.2,1,1],
             [0,0.4,0.2],[0,0.2,0.4],[0.2,0.4,0],[0.4,0,0.2],[0.4,0.2,0],[0.2,0.4,0],[0.4,0.4,0.2],[0.4,0.2,0.4],[0.2,0.4,0.4],[0.4,1,0.2],[0.4,0.2,1],[0.2,1,0.4]]
alpha=0.5
density=0.2
for J in [0.5,1,2,2.5,3,3.5,4]:
    for alpha in [0,0.2,0.5]:
        labels = np.loadtxt(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/batch/outputlabels_{J}_{alpha}_{density}_1.txt',dtype=int,skiprows=1)
        Path(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}').mkdir( exist_ok=True)
        saved_steps = int(labels.shape[0]/3)
        
        num_of_clusters=[]
        part_in_clusters=[]
        part_dif=[]
        avg_cluster_size=[]
        for i in range(saved_steps):
            
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
        plt.plot(np.arange(saved_steps),num_of_clusters,label = 'num of clusters')
        plt.plot(np.arange(saved_steps),part_in_clusters,label = 'particles in clusters')
        plt.plot(np.arange(saved_steps),part_dif,label = 'diffuse particles')
        plt.xlim((0,saved_steps))
        plt.vlines([100,500,950], 0, 400, colors='r',linestyles='dashed')
        plt.xlabel('MC_steps *10')
        plt.title(f'J={J}, alpha={alpha},density={density}')
        plt.legend()
        plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image/clusternum_{J}_{alpha}_{density}.svg')
        plt.show()
        for step in range(1,50):
            
            if step in [1,10,49]:
                cluster_size = np.bincount(labels[step*3])
                cluster_size_distribution = np.bincount(cluster_size[1:])
                plt.figure()
                plt.bar(np.arange(cluster_size_distribution.size), cluster_size_distribution)
                plt.ylim((0,55))
                plt.suptitle('cluster size distribution')
                plt.title(f'J={J}, alpha={alpha},density={density}')
                plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/hist_{step}_{J}_{alpha}_{density}.svg')
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
                    coloring_labels[x]=color_scale[labels[step*3,j]%39]
                j+=1
        
            plt.figure(figsize=(20,15))
            hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=coloring_labels,do_plot=True)
            centers_x = hex_centers[:, 0]
            centers_x = hex_centers[:, 1]
            plt.suptitle(f'{step}')
            plt.title(f'J={J} alpha={alpha}')
            #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
            plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/label_{step}.jpg')
            plt.show()
            
            # for x in labels[step*3+1]:
        
            #     # else:
            #     #     coloring_labels[x]=[195/255,192/255,192/255]
            #         # labels_name[x]=labels[0,i]
            #     if labels[step*3+2,i]==1:
            #         colors_labels[x]=[195/255,192/255,192/255]  
            #     elif labels[step*3+2,i]==2:
            #         colors_labels[x]=[0, 1, 0]
            #     elif labels[step*3+2,i]==3:
            #         colors_labels[x]=[1,0,0]
            #     elif labels[step*3+2,i]==4:
            #         colors_labels[x]=[0,0,1]
            #     i+=1
                
        
            # plt.figure(figsize=(20,15))
            # hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors_labels,do_plot=True)
            # centers_x = hex_centers[:, 0]
            # centers_x = hex_centers[:, 1]
            # #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
            # # plt.savefig(f'{counter}_labels_{step}_{J}_{alpha}_{density}.svg')
        
            # plt.show()
        
        images=[]
        for file in range(1,49):
            im = iio.imread(f"/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/label_{file}.jpg")
            images.append(im)
        iio.mimsave(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/movies/movie_{J}_{alpha}.gif', images, duration=0.8)
        
    
    #%%
J=4
alpha=0.5
density=0.2

labels = np.loadtxt(f'/home/hannohennighausen/Documents/frz_lattice_model/testlabels_{J}_{alpha}_{density}_11.txt',dtype=int,skiprows=1)
color_scale=[[0,0,1],[0,1,0],[1,0,0],[1,0,0.6],[1,0.6,0],[0.6,1,0],[0,1,0.6],[0,0.6,1],[0.6,0,1],[1,1,0.6],[1,0.6,1],[0.6,1,1],[0,0,0.6],[0,0.6,0],[0.6,0,0],
             [0,0,0.2],[0,0.2,0],[0.2,0,0],[1,0,0.2],[1,0.2,0],[0.2,1,0],[0,1,0.2],[0,0.2,1],[0.2,0,1],[1,1,0.2],[1,0.2,1],[0.2,1,1],
             [0,0.4,0.2],[0,0.2,0.4],[0.2,0.4,0],[0.4,0,0.2],[0.4,0.2,0],[0.2,0.4,0],[0.4,0.4,0.2],[0.4,0.2,0.4],[0.2,0.4,0.4],[0.4,1,0.2],[0.4,0.2,1],[0.2,1,0.4]]
for step in range(1,50):
    cluster_size = np.bincount(labels[step*3])
    cluster_size_distribution = np.bincount(cluster_size[1:])
    # plt.figure()
    # plt.bar(np.arange(cluster_size_distribution.size), cluster_size_distribution)
    # plt.ylim((0,55))
    # plt.suptitle('cluster size distribution')
    # plt.title(f'J={J}, alpha={alpha},density={density}')
    # plt.savefig(f'{counter}_size_dist_at_{step}_{J}_{alpha}_{density}.svg')
    # plt.show()
    
    colors_labels = np.zeros([50*50,3])
    coloring_labels = np.zeros([50*50,3])
    i=0
    j=0
    Cluster_num = np.max(labels[step*3])+20
    for x in labels[step*3+1]:
        if labels[step*3,j]==0:
            coloring_labels[x]=[195/255,192/255,192/255] 
        else:
            coloring_labels[x]=color_scale[labels[step*3,j]%39]
        j+=1


    # for x in labels[step*3+1]:

    #     # else:
    #     #     coloring_labels[x]=[195/255,192/255,192/255]
    #         # labels_name[x]=labels[0,i]
    #     if labels[step*3+2,i]==1:
    #         colors_labels[x]=[195/255,192/255,192/255]  
    #     elif labels[step*3+2,i]==2:
    #         colors_labels[x]=[0, 1, 0]
    #     elif labels[step*3+2,i]==3:
    #         colors_labels[x]=[1,0,0]
    #     elif labels[step*3+2,i]==4:
    #         colors_labels[x]=[0,0,1]
    #     i+=1
        

    # plt.figure(figsize=(20,15))
    # hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors_labels,do_plot=True)
    # centers_x = hex_centers[:, 0]
    # centers_x = hex_centers[:, 1]
    # #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
    # # plt.savefig(f'{counter}_labels_{step}_{J}_{alpha}_{density}.svg')

    # plt.show()

    plt.figure(figsize=(20,15))
    hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=coloring_labels,do_plot=True)
    centers_x = hex_centers[:, 0]
    centers_x = hex_centers[:, 1]
    plt.suptitle(f'{step}')
    plt.title(f'J={J} alpha={alpha}')
    #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
    # plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/label_{step}.jpg')
    plt.show()
    #%%
images=[]
for file in range(1,49):
    im = iio.imread(f"/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/label_{file}.jpg")
    images.append(im)
iio.mimsave(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/movie_{J}_{alpha}.gif', images, duration=0.5)