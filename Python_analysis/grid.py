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
import pandas as pd
import seaborn as sns

ALPHA=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]
J_r=[0.5,1,1.25,1.5,1.75,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.25,3.5,3.75,4]
#%%
dnum_of_clusters=[]
dpart_in_clusters=[]
dpart_dif=[]
X=[]
for run in range(1,20):
            labels = np.loadtxt(f'/scratch2/hannohennighausen/batch_4_0.5/outputlabels_2_0.5_0.2_{run}.txt',dtype=int,skiprows=1)
            # print(J,alpha,run)
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
            
            dnum_of_clusters=np.append(dnum_of_clusters,num_of_clusters)
            dpart_in_clusters=np.append(dpart_in_clusters ,part_in_clusters)
            dpart_dif=np.append(dpart_dif ,part_dif)
            X=np.append(X,np.arange(0,50,1))
dtraj_long = {"X":X,"part_in_cluster":dpart_in_clusters}
traj_long = pd.DataFrame(data=dtraj_long)

plt.figure()
sns.lineplot(data=traj_long,x='X',y='part_in_cluster',color='g', ci='sd')
plt.legend()
plt.show()
#%%
#%%
density=0.2
for J in J_r:
    for alpha in ALPHA:
        dnum_of_clusters=[]
        dpart_in_clusters=[]
        dpart_dif=[]
        X=[]
        for run in range(1,20):
            labels = np.loadtxt(f'/scratch2/hannohennighausen/Parameter_sweep/outputlabels_{J}_{alpha}_0.2_{run}.txt',dtype=int,skiprows=1)
            print(J,alpha,run)
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
            
            dnum_of_clusters=np.append(dnum_of_clusters,num_of_clusters)
            dpart_in_clusters=np.append(dpart_in_clusters ,part_in_clusters)
            dpart_dif=np.append(dpart_dif ,part_dif)
            X=np.append(X,np.arange(0,50,1))
        dtraj_long = {"X":X,"part_in_cluster":dpart_in_clusters,"part_dif":dpart_dif,"num_of_clusters":dnum_of_clusters}
        traj_long = pd.DataFrame(data=dtraj_long)
        
        
        fig, axs = plt.subplots(3, sharex=True, figsize=(10,8))
        fig.suptitle(f'AVG 20 J={J}, alpha={alpha},density={density}')
        sns.lineplot(data=traj_long,x='X',y='part_in_cluster',color='orange', ci='sd', ax=axs[0])
        sns.lineplot(data=traj_long,x='X',y='part_dif',color='g', ci='sd',ax=axs[1])     
        sns.lineplot(data=traj_long,x='X',y='num_of_clusters',color='b', ci='sd',ax=axs[2])
        plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_par_sweep/avg_{J}_{alpha}_{density}.svg')
        plt.show()
#%%

color_scale=[[0,0,1],[0,1,0],[1,0,0],[1,0,0.6],[1,0.6,0],[0.6,1,0],[0,1,0.6],[0,0.6,1],[0.6,0,1],[1,1,0.6],[1,0.6,1],[0.6,1,1],[0,0,0.6],[0,0.6,0],[0.6,0,0],
             [0,0,0.2],[0,0.2,0],[0.2,0,0],[1,0,0.2],[1,0.2,0],[0.2,1,0],[0,1,0.2],[0,0.2,1],[0.2,0,1],[1,1,0.2],[1,0.2,1],[0.2,1,1],
             [0,0.4,0.2],[0,0.2,0.4],[0.2,0.4,0],[0.4,0,0.2],[0.4,0.2,0],[0.2,0.4,0],[0.4,0.4,0.2],[0.4,0.2,0.4],[0.2,0.4,0.4],[0.4,1,0.2],[0.4,0.2,1],[0.2,1,0.4]]
alpha=0.5
density=0.2
run=8
index = 12
off = 0.03
FrzB = 210
for run in [20]:
    for J in [4]:
        for titration in [1]:
            # labels = np.loadtxt(f'/scratch2/hannohennighausen/Parameter_sweep/rectangular_outputlabels_{J}_{alpha}_0.2_11.txt',dtype=int,skiprows=1)
            labels = np.loadtxt(f'/home/subraman/Documents/Github_projects/frz_lattice_simulation/{index}_dFrzB_labels_J_{J}_alpha_{alpha}_FrzB_{FrzB}_off_{off}.txt',dtype=int,skiprows=1)
            Frz = np.loadtxt(f'/home/subraman/Documents/Github_projects/frz_lattice_simulation/{index}_dFrzB_labels_J_{J}_alpha_{alpha}_FrzB_{FrzB}_off_{off}.txt',dtype=int,skiprows=1)
            # f"/home/hannohennighausen/Documents/frz_lattice_model/
            Path(f'/home/subraman/Documents/Github_projects/frz_lattice_simulation/{index}_dFrzB_labels_J_{J}_alpha_{alpha}_FrzB_{FrzB}_off_{off}').mkdir( exist_ok=True)
            Nx,Ny=np.loadtxt(f'/home/subraman/Documents/Github_projects/frz_lattice_simulation/{index}_dFrzB_labels_J_{J}_alpha_{alpha}_FrzB_{FrzB}_off_{off}.txt',dtype=int,max_rows=1)
            saved_steps = int(labels.shape[0]/3)
  
            for step in range(0,saved_steps):
                
                # if step in [1,10,49]:
                #     cluster_size = np.bincount(labels[step*3])
                #     cluster_size_distribution = np.bincount(cluster_size[1:])
                #     plt.figure()
                #     plt.bar(np.arange(cluster_size_distribution.size), cluster_size_distribution)
                #     plt.ylim((0,55))
                #     plt.suptitle('cluster size distribution')
                #     plt.title(f'J={J}, alpha={alpha},density={density}')
                #     plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/hist_{step}_{J}_{alpha}_{density}.svg')
                #     plt.show()
                
                colors_labels = np.zeros([int(Nx)*int(Ny),3])
                coloring_labels = np.zeros([int(Nx)*int(Ny),3])
                colors_Frz=np.zeros([int(Nx)*int(Ny),3])
                i=0
                j=0
                ii=0
                
                for x in labels[step*3+1]:
                    if labels[step*3,j]==0:
                        coloring_labels[x]=[195/255,192/255,192/255] 
                    else:
                        coloring_labels[x]=color_scale[labels[step*3,j]%39]
                    j+=1
                fig, axs = plt.subplots(2,sharex=(True))
                # plt.figure(figsize=(20,15))
                hex_centers, _ = hex.create_hex_grid(nx= int(Nx),ny=int(Ny), face_color=coloring_labels,do_plot=True,h_ax=axs[0])
                centers_x = hex_centers[:, 0]
                centers_x = hex_centers[:, 1]
                # plt.suptitle(f'{step}')
                fig.suptitle(f'{(step+1)} $\cdot 10^4$ Monte Carlo Steps')
                #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
                # plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/FrzB_test_{run}_image_{J}_{alpha}_{density}/label_{step}.jpg')
                # plt.show()
                
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
                    
            
                # # plt.figure(figsize=(20,15))
                # hex_centers, _ = hex.create_hex_grid(nx= int(Nx),ny=int(Ny), face_color=colors_labels,do_plot=True,h_ax=axs[2])
                # centers_x = hex_centers[:, 0]
                # centers_x = hex_centers[:, 1]
                # #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
                # # plt.savefig(f'{counter}_labels_{step}_{J}_{alpha}_{density}.svg')
            
                # # plt.show()
                
                for x in labels[step*3+1]:
            
                    # else:
                    #     coloring_labels[x]=[195/255,192/255,192/255]
                        # labels_name[x]=labels[0,i]
                    if Frz[step,ii]==0:
                        colors_Frz[x]=[195/255,192/255,192/255]  
                    elif Frz[step,ii]==1:
                        colors_Frz[x]=[0, 1, 0]
                    elif Frz[step,ii]==2:
                        colors_Frz[x]=[0, 0, 1]
                    elif Frz[step,ii]==3:
                        colors_Frz[x]=[1, 0, 0]
                    ii+=1
                    
            
                # plt.figure(figsize=(20,15))
                
                hex_centers, axs[1] = hex.create_hex_grid(nx= int(Nx),ny=int(Ny), face_color=colors_Frz,do_plot=True,h_ax=axs[1])
                centers_x = hex_centers[:, 0]
                centers_x = hex_centers[:, 1]
                
                # plt.title(f'{step}')
                plt.savefig(f'/home/subraman/Documents/Github_projects/frz_lattice_simulation/{index}_dFrzB_labels_J_{J}_alpha_{alpha}_FrzB_{FrzB}_off_{off}/Frz_B_{step}.jpg')
                #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
                # plt.savefig(f'{counter}_labels_{step}_{J}_{alpha}_{density}.svg')
            
                plt.show()
            
            images=[]
            for file in range(1,49):
                im = iio.imread(f"/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/weird_FrzB_{run}_image_{J}_{alpha}_{titration}/Frz_B_{file}.jpg")
                images.append(im)
            iio.mimsave(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/weird_FrzB_{run}_image_{J}_{alpha}_{titration}/movie_{J}_{alpha}_{run}.gif', images, duration=0.8)
      #%%
FrzB_tot = []
for i in range(50):
    unique,counts=np.unique(Frz[i],return_counts=True)
    
    S = 0
    if len(counts)==4:
        S= counts[1]+counts[2]+2*counts[3]
    elif len(counts)==3:
        S=counts[1]+counts[2]    
    FrzB_tot = np.append(FrzB_tot,S)
    print(dict(zip(unique, counts)), S)

    
pT = 420
xT = 610
k_off = 0.3
Nxy = Nx*Ny

p = 0.5*(xT+2*pT-k_off)
q= 2*pT*xT
pq = p+np.sqrt(p**2-q)
print(pq,p-np.sqrt(p**2-q))    
#%%
# J=4
# alpha=0.5
# density=0.2

# labels = np.loadtxt(f'/home/hannohennighausen/Documents/frz_lattice_model/testlabels_{J}_{alpha}_{density}_11.txt',dtype=int,skiprows=1)
# color_scale=[[0,0,1],[0,1,0],[1,0,0],[1,0,0.6],[1,0.6,0],[0.6,1,0],[0,1,0.6],[0,0.6,1],[0.6,0,1],[1,1,0.6],[1,0.6,1],[0.6,1,1],[0,0,0.6],[0,0.6,0],[0.6,0,0],
#              [0,0,0.2],[0,0.2,0],[0.2,0,0],[1,0,0.2],[1,0.2,0],[0.2,1,0],[0,1,0.2],[0,0.2,1],[0.2,0,1],[1,1,0.2],[1,0.2,1],[0.2,1,1],
#              [0,0.4,0.2],[0,0.2,0.4],[0.2,0.4,0],[0.4,0,0.2],[0.4,0.2,0],[0.2,0.4,0],[0.4,0.4,0.2],[0.4,0.2,0.4],[0.2,0.4,0.4],[0.4,1,0.2],[0.4,0.2,1],[0.2,1,0.4]]
# for step in range(1,50):
#     cluster_size = np.bincount(labels[step*3])
#     cluster_size_distribution = np.bincount(cluster_size[1:])
#     # plt.figure()
#     # plt.bar(np.arange(cluster_size_distribution.size), cluster_size_distribution)
#     # plt.ylim((0,55))
#     # plt.suptitle('cluster size distribution')
#     # plt.title(f'J={J}, alpha={alpha},density={density}')
#     # plt.savefig(f'{counter}_size_dist_at_{step}_{J}_{alpha}_{density}.svg')
#     # plt.show()
    
#     colors_labels = np.zeros([50*50,3])
#     coloring_labels = np.zeros([50*50,3])
#     i=0
#     j=0
#     Cluster_num = np.max(labels[step*3])+20
#     for x in labels[step*3+1]:
#         if labels[step*3,j]==0:
#             coloring_labels[x]=[195/255,192/255,192/255] 
#         else:
#             coloring_labels[x]=color_scale[labels[step*3,j]%39]
#         j+=1


#     # for x in labels[step*3+1]:

#     #     # else:
#   #     #     coloring_labels[x]=[195/255,192/255,192/255]
#     #         # labels_name[x]=labels[0,i]
#     #     if labels[step*3+2,i]==1:
#     #         colors_labels[x]=[195/255,192/255,192/255]  
#     #     elif labels[step*3+2,i]==2:
#     #         colors_labels[x]=[0, 1, 0]
#     #     elif labels[step*3+2,i]==3:
#     #         colors_labels[x]=[1,0,0]
#     #     elif labels[step*3+2,i]==4:
#     #         colors_labels[x]=[0,0,1]
#     #     i+=1
        

#     # plt.figure(figsize=(20,15))
#     # hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=colors_labels,do_plot=True)
#     # centers_x = hex_centers[:, 0]
#     # centers_x = hex_centers[:, 1]
#     # #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
#     # # plt.savefig(f'{counter}_labels_{step}_{J}_{alpha}_{density}.svg')

#     # plt.show()

#     plt.figure(figsize=(20,15))
#     hex_centers, _ = hex.create_hex_grid(nx= 50,ny=50, face_color=coloring_labels,do_plot=True)
#     centers_x = hex_centers[:, 0]
#     centers_x = hex_centers[:, 1]
#     plt.suptitle(f'{step}')
#     plt.title(f'J={J} alpha={alpha}')
#     #plt.savefig(f'D:/Hanno/Physics/Marburg/Murray/frz_lattice_model/rand_grid_19899+{row}.png')
#     # plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/label_{step}.jpg')
#     plt.show()
    #%%
# images=[]
# for file in range(1,49):
#     im = iio.imread(f"/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/label_{file}.jpg")
#     images.append(im)
# iio.mimsave(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_{J}_{alpha}/movie_{J}_{alpha}.gif', images, duration=0.5)

#%%
ALPHA1=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
mesh=np.meshgrid(ALPHA,J_r)
coarsening=np.zeros((len(J_r),len(ALPHA1)))
num_dif=np.zeros((len(J_r),len(ALPHA1)))
n=0

for J in J_r:
    m=0
    for alpha in ALPHA1:
        dnum_of_clusters=[]
        dnum_dif=[]
        for run in range(1,21):
            # print(J,alpha,run)
            labels = np.loadtxt(f'/scratch2/hannohennighausen/Parameter_sweep/outputlabels_{J}_{alpha}_0.2_{run}.txt',dtype=int,skiprows=1)
            u=np.bincount(labels[49*3])
            # bonds=np.loadtxt(f'/scratch2/hannohennighausen/Parameter_sweep/output_bonds{J}_{alpha}_0.2_{run}.txt',skiprows=49,max_rows=1)
            
            ordering=u[1:].max()
            dnum_of_clusters=np.append(dnum_of_clusters,ordering)
            dnum_dif=np.append(dnum_dif,u[0])
        coarsening[n,m]=np.mean(dnum_of_clusters)/500
        num_dif[n,m]=np.mean(dnum_dif)
        m+=1
    n+=1
plt.figure()
ax = sns.heatmap(coarsening,xticklabels=ALPHA1,yticklabels=J_r)
ax.invert_yaxis()
ax.set_xlabel("alpha")
ax.set_ylabel("J")
ax.set_title('(#part in largest cluster)/(#tot part)')
plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_par_sweep/heatmap_numpartincluster_higheralpha.jpg')
plt.show()    

plt.figure()
ax = sns.heatmap(num_dif,xticklabels=ALPHA1,yticklabels=J_r)
ax.invert_yaxis()
ax.set_xlabel("alpha")
ax.set_ylabel("J")
ax.set_title('dif part at the end')
plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_par_sweep/heatmap_numdif.jpg')
plt.show()    

#%%
ALPHA1=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
def hamiltonian(J,alpha,num_bonds,num_part):
    return -J*num_bonds+alpha*J*num_part
Energy=np.zeros((len(J_r),len(ALPHA1)))
n=0
for J in J_r:
    m=0
    for alpha in ALPHA1:
        dEnergy=[]
        for run in range(1,21):
            print(J,alpha,run)
            bonds=pd.read_csv(f'/scratch2/hannohennighausen/Parameter_sweep/output_bonds{J}_{alpha}_0.2_{run}.txt',skiprows=50,sep='\t',header=None,keep_default_na=False)
            labels = np.loadtxt(f'/scratch2/hannohennighausen/Parameter_sweep/outputlabels_{J}_{alpha}_0.2_{run}.txt',dtype=int,skiprows=1)
            u=np.bincount(labels[49*3])
            # bonds=np.loadtxt(f'/scratch2/hannohennighausen/Parameter_sweep/output_bonds{J}_{alpha}_0.2_{run}.txt',skiprows=49,max_rows=1)
            
            ordering=hamiltonian(J,alpha,bonds.sum(axis=1),u[1:].sum())
            dEnergy=np.append(dEnergy,ordering)
        Energy[n,m]=np.mean(dEnergy)
        m+=1
    n+=1

plt.figure()
ax = sns.heatmap(Energy,xticklabels=ALPHA1,yticklabels=J_r)
ax.invert_yaxis()
ax.set_xlabel("alpha")
ax.set_ylabel("J")
ax.set_title('Energy of system')
plt.savefig(f'/home/hannohennighausen/Documents/frz_lattice_model/Python_analysis/image_par_sweep/heatmap_Energy_system.jpg')
plt.show()    

#%%

print(labels[0])