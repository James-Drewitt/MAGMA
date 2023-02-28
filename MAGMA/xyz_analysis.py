#!/usr/bin/env python3
#
# Dr James Drewitt, 23/06/2020 # last update: 28/02/2023
#
import numpy as np
import time
from xyz_CN_subroutines import alpha_beta, calc_n_data, av_cn, bond_lifetime, save_files
import os
from pathlib import Path

def CN(filename,T_step,alpha,beta,rcut,L, save_config, xyz_numT, working_dir):

    start = time.time() # initiate runtime timer
    
    # Read xyz file
    print(f"\n reading file {filename} ....")
    #
    with open(filename, 'r') as f:
        data_list = f.readlines()
    
    # determine the total number of atoms in the system:
    n_atoms = int(data_list[0])

    # determine the number of trajectories:
    n_traj = int(len(data_list) / (n_atoms+2))
    #
    print(f"\n xyz trajectory file contains {n_atoms} atoms and {n_traj} trajectories")
    
    # truncate xyz trajectory if required
    
    if n_traj > xyz_numT:
        trunc_traj = int(n_traj - xyz_numT)
        trunc_lines = int(trunc_traj *(n_atoms+2))
        print(f"\n truncating the first {trunc_traj} trajectories ({trunc_lines} lines)")
        del data_list[0:trunc_lines]

    n_traj = int(len(data_list) / (n_atoms+2))
    n_traj_T = int(n_traj/T_step)
    print(f"\n *** Calculating {alpha}-{beta} coordination ***")
    print(f" sampling every {T_step} trajectories, total trajectories in analysis = {n_traj_T}")
    
    # define strings for output
    
    str_n="n("+alpha+"-"+beta+")"
    str_r="r("+alpha+"-"+beta+") [rcut= "+str(rcut)+"])"

    # initialise lists for alpha-beta calculation
    #
    # data contains the coordinates of alpha-beta atoms for specified coordinations
    # data2 is the same as data2 but for the beta-alpha pairs
    # n_data provides the fractional coordinations for each trajectory
    
    data = [[n_traj_T, "  Partial coordination", str_n, "and mean distance", str_r]] # number of trajectories sampled header
    n_data = [[n_traj_T, " ", " "]]
    n_data.append([alpha, beta, "fractional coordination"])
    n_beta_tot = []

    data2 = data.copy() # copy lists for beta-alpha calculation
    n_data2 = n_data.copy() 
    n_beta_tot2 = []
    
    for traj in range(0, n_traj, T_step): # iterate over all trajectories with step T_step
           
        count_a = 0 # (re)initilise counter
        count_b = 0
        
        # each trajectory in an xyz file contains two header lines followed by {natoms} lines.
        # iterate over atoms in current trajectory "traj" using a list comprehension
        # if atom is "alpha" then append a tuple containing the atom symbol and its cartesian coordinates
        # to the list "alpha_atoms_and_coords"
        
        alpha_atoms_and_coords = (
                
            [(a, float(x), float(y), float(z)) 
             for line in data_list[2+traj*(n_atoms+2):n_atoms+2+traj*(n_atoms+2)]
             for a, x, y, z in [line.split()] if a == alpha]
            )
        
        # add labels to the atomic coordinates in "alpha_atoms_and_coords" and count the number of alpha atoms
        
        coord_a = [(f"{a}{i+1}", x, y, z) for i, (a, x, y, z) in enumerate(alpha_atoms_and_coords)]
        count_a = len(coord_a)
        
        # if atom is "beta" then append a tuple containing the atom symbol and its cartesian coordinates
        # to the list "beta_atoms_and_coords"
        
        beta_atoms_and_coords = (
                
            [(b, float(x), float(y), float(z)) 
             for line in data_list[2+traj*(n_atoms+2):n_atoms+2+traj*(n_atoms+2)]
             for b, x, y, z in [line.split()] if b == beta]
            )
        
        # add labels to the atomic coordinates in "alpha_atoms_and_coords" and count the number of alpha atoms
        
        coord_b = [(f"{b}{i+1}", x, y, z) for i, (b, x, y, z) in enumerate(beta_atoms_and_coords)]
        count_b = len(coord_b)
        
        data += [[traj+1, " ", count_a, " ", " "]] # provide current trajectory number for output
        data2 += [[traj+1, " ", count_b, " ", " "]] # provide current trajectory number for output
        
        data, n_beta, n_beta_tot = alpha_beta(L, rcut, coord_a, beta, coord_b, data, n_beta_tot)
        n_data = calc_n_data(n_beta, traj, n_data) # alpha-beta

        data2, n_beta2, n_beta_tot2 = alpha_beta(L, rcut, coord_b, alpha, coord_a, data2, n_beta_tot2)
        n_data2 = calc_n_data(n_beta2, traj, n_data2) #beta-alpha
 
    cn_tot, N = av_cn(alpha, beta, n_beta_tot, n_traj_T, n_data)
    cn_tot2 , N2 = av_cn(beta, alpha, n_beta_tot2, n_traj_T, n_data2)

    end = time.time() # end runtime timer
    elapsed = round(end - start , 4)
    print(f"\n runtime for coordination number calculations = {elapsed} s")

    save_files(alpha, beta, data, n_data, cn_tot, N, save_config, working_dir) # save alpha-beta
    save_files(beta, alpha, data2, n_data2, cn_tot2, N2, save_config, working_dir) # save beta-alpha

    print(f"\n *** Compute {alpha}-{beta} bond lifetimes ***")

    start = time.time() # initiate runtime timer
    
    b_lifetime, mean_b_lifetime, median_b_lifetime, max_b_lifetime, min_b_lifetime = bond_lifetime(data)

    end = time.time() # end runtime timer
    elapsed = round(end - start , 4)

    b_lifetime = b_lifetime*T_step
    mean_b_lifetime = mean_b_lifetime*T_step
    median_b_lifetime = median_b_lifetime*T_step
    max_b_lifetime = max_b_lifetime*T_step
    min_b_lifetime = min_b_lifetime*T_step

    str1 = " mean lifetime of " + alpha + "-" + beta + " bonds = " + str(mean_b_lifetime) + " timesteps"
    str2 = " median lifetime = " + str(median_b_lifetime) + " timesteps"
    str3 = " min lifetime = " + str(min_b_lifetime) + " timesteps"
    str4 = " max lifetime = " + str(max_b_lifetime) + " timesteps"
    print(str1)
    print(str2)
    print(str3)
    print(str4)
    print(f"\n runtime for lifetime calculations = {elapsed} s")

    life_array = np.zeros(4, dtype=object)

    life_array[0] = str1
    life_array[1] = str2
    life_array[2] = str3
    life_array[3] = str4

    CWD=os.getcwd()
    fname = Path(CWD+"/"+working_dir+"/"+alpha + "-" + beta + "_bond_lifetime.dat")
        
    np.savetxt(fname, life_array, delimiter =" ", fmt ="%s")

    return data, data2
    
    
