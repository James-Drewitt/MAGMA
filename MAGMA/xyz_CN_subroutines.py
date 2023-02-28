#!/usr/bin/env python3
#
## AUTHOR: DR JAMES DREWITT, 01/11/2020
##
## james.drewitt@bristol.ac.uk
##
## Last update: 28/02/2023

import numpy as np
import os
from pathlib import Path

def d_pbc(coord1, coord2, L):

    '''
    Calculates the distance between two coordinates within periodic boundary conditions
    coord1 -> coordinate 1
    coord2 -> coordinate 2
    L-> configuration box length
    '''

    if coord1 < 0:          # First check for negative coordinates
        coord1 = L + coord1
    if coord2 < 0:
        coord2 = L + coord2 
    d1 = abs(coord2 - coord1) # find minimum distance within periodic boundary conditions
    d2 = abs(coord2 - coord1 + L)
    d3 = abs(coord2 - coord1 - L)

    d = min(abs(d1), abs(d2), abs(d3))

    if abs(d1) == d:
        d = d1
    elif abs(d2) == d:
        d = d2
    elif abs(d3) == d:
        d = d3
        
    return d


def alpha_beta(L, rcut, coord_a, b, coord_b, data, n_beta_tot):

    n_beta = [] # initialise coordination number list
    
    for i in range(len(coord_a)):
            pair=[[0,0,0,0,0],[0,0,0,0,0]] # initialise np array for current alpha-beta pair
            n=0 # initialise counter
            D=0 # initialise distance
            a_atom=coord_a[i][0]

            for j in range(len(coord_b)):
                x1, y1, z1, x2, y2, z2 = coord_a[i][1], coord_a[i][2], coord_a[i][3], coord_b[j][1], coord_b[j][2], coord_b[j][3] # extract coordinates of aplha-beta atoms

                dx = abs(d_pbc(x1, x2, L)) # calculate distance from central atom in each dimension
                dy = abs(d_pbc(y1, y2, L))
                dz = abs(d_pbc(z1, z2, L))
                d=np.sqrt( np.square(dx) + np.square(dy) + np.square(dz) )#compute distance

                if d <= rcut:
                    n+=1 # count beta atoms within distance rcut from alpha
                    D+=d

                    b_atom = coord_b[j][0]
                    pair.append( [b_atom, x2, y2, z2, d] ) # populate array with alpha-beta pair and distance

            if n == 0:
                D = "--"
            else:
                D=D/n # average alpha-beta distance
            n_beta.append(n)
            n_beta_tot.append(n)
            string1="n("+a_atom+"-"+b+") ="
            string2="r("+a_atom+"-"+b+") ="
            p_coord=[string1, int(n), " ", string2, D] # array containing partial coordination and distance
            alpha_atom=[a_atom, x1, y1, z1, " "] # array containing atom alpha configuration
            pair[0]=p_coord
            pair[1]=alpha_atom
            for m in range(len(pair)):
                data.append(pair[m])
            data.append( [" ", " ", " ", " ", " "] )

    return data, n_beta, n_beta_tot


def calc_n_data(n_beta, traj, n_data):

    if len(n_beta) >= 1:
        
        # Add the minimum and maximum values in "n_beta" and the number of elements between these
        # values to the "n_data" list. 
        
        n_min = min(n_beta)
        n_max = max(n_beta)
        num_n = (n_max - n_min) + 1
        n_data += [[traj+1, " " , num_n]]
        
        # List comprehension to add a list containing the partial coordination value 'k' 
        # the count of 'k' in 'n_beta', and its percent fraction
        
        n_data += (
            [[k, n_beta.count(k), (n_beta.count(k)/len(n_beta))*100] 
             for k in range(n_min, n_max+1)]
            )
        # Add empty strings to n_data
        
        n_data += [[" ", " ", " "]]
    
    else:
        
        # Executes if n_beta is empty
        
        n_data += [traj+1, " ", 0]
        n_data += [" ", " ", " "]

    return n_data


def av_cn(alpha, beta, n_beta_tot, n_traj_T, n_data):

    # calculate average coordination across all trajectories
    
    n_min_tot = min(n_beta_tot) # find minimum partial coordination
    n_max_tot = max(n_beta_tot) # find maxium parial coordination
    n_tot=n_max_tot-n_min_tot+1 # number of partial coordinations

    ini = [0]*n_tot #initialise array for partial coordinations
    ini2 = ini[:]
    ini3 = ini[:]
    cn_tot = [ini,ini2,ini3]

    for a in range(n_tot): # populate first column of "cn_tot" with each partial coordination
        cn_tot[0][a]=a+n_min_tot
        
    n=0 # initialise iterator
    N=0 # initialise total average coordination
    for b in range(n_traj_T):
         n += 2
         cn_num = n_data[n][2] # number of partial coordinations in current trajectory
         for p in range(cn_num):
             n += 1
             for q in range(n_tot):
                 if n_data[n][0] == cn_tot[0][q]:
                     cn_tot[1][q] += n_data[n][1] # populate each row with number of units of specified coordination
    
    total=0
    for val in cn_tot[1]:
        total += val
    
    for r in range(n_tot):
        cn_tot[2][r] = ( cn_tot[1][r] / total ) * 100 # fraction (per cent) of each partial coordination
        cn_tot[2][r] = round(cn_tot[2][r],3)
        N += cn_tot[0][r] * (cn_tot[2][r] / 100) # calculate average coordination
    
    print(f" average {alpha}-{beta} coordination number = {N}")

    return cn_tot, N 

def bond_lifetime(data):

    n_traj = data[0][0]
    
    bond_list_list = []
    bond_list_tot = []

    n = 0
    
    if n_traj >=2:
        

        for i in range(n_traj):
    
            bond_list = []
    
            n += 1
    
            n_a = data[n][2]
    
            num_b = 0
    
            for j in range(n_a):
                n += 1
                cn_num = data[n][1] # number of partial coordinations in trajectory
                a_atom = data[n+1][0]
                
                if cn_num > 0:
                    num_b += cn_num
                    for k in range(2,cn_num+2):
                        bond_name = a_atom+"-"+data[n+k][0]
                        bond_list.append(bond_name)
                        bond_list_list.append(bond_name)
    
                n += cn_num + 2
    
            bond_list_tot.append(num_b)
            bond_list_tot.append(bond_list)
            
        bond_set_list = list(set(bond_list_list)) # generate list of unique bonds only
    
        tot_bonds = len(bond_set_list)
    
        life = np.zeros((tot_bonds, n_traj+1), dtype = object)
        life[:,0] = bond_set_list[:]
    
        n = 0 # initialise iterator
    
        for i in range(n_traj): # loop for number of trajectories 
    
            n_bonds = bond_list_tot[n] # get number of bonds in current trajecty
    
            b_list = []
            n += 1
    
            for j in range(n_bonds):
                bond = bond_list_tot[n][j]
                b_list.append(bond)
            for k in range(tot_bonds):
                if life[k,0] in b_list:
                    life[k,i+1] = 1
    
            n += 1
    
        life_time = []
    
        for i in range(tot_bonds):
            L = 0
            for j in range(n_traj):
                if life[i, j+1] == 1:
                    L += 1
                    if j == n_traj:
                        life_time.append(L)
                else:
                    if L != 0:
                        life_time.append(L)
                    L = 0 # reset L count
                    
    else:
        life_time = 1,1

    mean_lifetime = np.mean(life_time)
    median_lifetime = np.median(life_time)
    max_lifetime = max(life_time)
    min_lifetime = min(life_time)

    return life_time, mean_lifetime, median_lifetime, max_lifetime, min_lifetime
        

def save_files(alpha, beta, data, n_data, cn_tot, N, save_config, working_dir):

    CWD=os.getcwd()
    datadir=Path(CWD+"/"+working_dir)
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    str_n="n("+alpha+"-"+beta+")"

    np_cn_tot = np.array(cn_tot)
    np_cn_tot = np.transpose(np_cn_tot)
    mt_row=np.empty([1, 3], dtype=object)
    for i in range(3):
        mt_row[0,i]=" "

    mt_row2=np.copy(mt_row)
    mt_row2[0,0]="average CN:"
    mt_row2[0,1]=N

    head=np.copy(mt_row)
    head[0,0]=str_n
    head[0,1]="number"
    head[0,2]="fraction (%)"

    head = np.append(head, mt_row, axis=0)
    np_cn_tot = np.append(head, np_cn_tot, axis=0)
    np_cn_tot = np.append(np_cn_tot, mt_row, axis=0)
    np_cn_tot = np.append(np_cn_tot, mt_row2, axis=0)

    av_CN = Path(CWD+"/"+working_dir+"/"+alpha+"-"+beta+"-CN-av.dat")
    print(f" ... saving {av_CN} ...")
    np.savetxt(av_CN , np_cn_tot , delimiter=" " , fmt="%s" )
    
    if save_config ==1:
        np_data = np.array(data)
        CN_config = Path(CWD+"/"+working_dir+"/"+alpha+"-"+beta+"-CN-config.dat")
        print(f"\n ... saving {CN_config} ...")
        np.savetxt(CN_config , np_data , delimiter=" " , fmt="%s")

        np_n_data = np.array(n_data)
        partial_n = Path(CWD+"/"+working_dir+"/"+alpha+"-"+beta+"-CN-traj.dat")
        print(f" ... saving {partial_n} ...")
        np.savetxt(partial_n , np_n_data , delimiter=" " , fmt="%s" )
