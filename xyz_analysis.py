# James Drewitt, 23/06/2020 # last update: 26/08/2020
#
import numpy as np
import time
from xyz_CN_subroutines import *


def CN(file,T_step,alpha,beta,rcut,L, save_config):

    start = time.time() # initiate runtime timer
    
    filename=file
    #Read xyz file
    print(f"\n reading file {filename} ....")
    #
    with open(filename, 'r') as f:
        data_list = f.readlines()
    
    #determine number of atoms:
    n_atoms = int(data_list[0])

    #determine number of trajectories:
    n_traj = int(len(data_list) / (n_atoms+2))
    n_traj_T = int(n_traj/T_step)
    #
    print(f"\n xyz trajectory file contains {n_atoms} atoms and {n_traj} trajectories")
    print(f"\n *** Calculating {alpha}-{beta} coordination ***")
    print(f" sampling every {T_step} trajectories, total trajectories in analysis = {n_traj_T}")
    #calculate CN(alpha-beta)
    str_n="n("+alpha+"-"+beta+")"
    str_r="r("+alpha+"-"+beta+") [rcut= "+str(rcut)+"])"
    #str_rcut="[rcut= "+str(rcut)+"])"
    str_xa="x_"+alpha;str_ya="y_"+alpha;str_za="z_"+alpha
    str_xb="x_"+beta;str_yb="y_"+beta;str_zb="z_"+beta

    # initialise lists for alpha-beta calculation
    data = [[n_traj_T, "  Partial coordination", str_n, "and mean distance", str_r]] # number of trajectories sampled header
    n_data = [[n_traj_T, " ", " "]]
    n_data.append([alpha, beta, "fractional coordination"])
    n_beta_tot = []

    data2 = data.copy() # copy lists for beta-alpha calculation
    n_data2 = n_data.copy() 
    n_beta_tot2 = []

    count_a = 0 # initilise counter
    count_b = 0
    
    for traj in range(0, n_traj, T_step): # iterate over all trajectories with step T_step
           
        count_a = 0 # reinitilise counter
        count_b = 0
        
        for i in range(2+traj*(n_atoms+2), n_atoms+2+traj*(n_atoms+2)): # iterate over each atom in current trajectory

            coord=data_list[i].split() # get atomic coordinates of current atom
            
            if coord[0]==alpha: # execute if current atom is desired alpha

                x = float(coord[1])
                y = float(coord[2])
                z = float(coord[3])
                a = coord[0]
                a_label = count_a + 1
                a_atom = a + str(a_label) # label current alpha atom
                
                if count_a == 0:
                    coord_a = [[a_atom, x, y, z]]
                else:
                    coord_a.append( [a_atom, x, y, z] )

                count_a += 1 # iterate counter
                
            elif coord[0]==beta: # execute if current atom is desired alpha

                x = float(coord[1])
                y = float(coord[2])
                z = float(coord[3])
                b = coord[0]
                b_label = count_b + 1
                b_atom = b + str(b_label) # label current beta atom
                
                if count_b == 0:
                    coord_b = [[b_atom, x, y, z]]
                else:
                    coord_b.append( [b_atom, x, y, z] )

                count_b += 1 # iterate counter

        data.append([traj+1, " ", count_a, " ", " "])# provide current trajectory number for output
        data2.append([traj+1, " ", count_b, " ", " "])# provide current trajectory number for output
        
        data, n_data, n_beta, n_beta_tot = alpha_beta(L, rcut, coord_a, b, coord_b, data, n_data, n_beta_tot)
        n_data = calc_n_data(n_beta, traj, n_data) # alpha-beta

        data2, n_data2, n_beta2, n_beta_tot2 = alpha_beta(L, rcut, coord_b, a, coord_a, data2, n_data2, n_beta_tot2)
        n_data2 = calc_n_data(n_beta2, traj, n_data2) #beta-alpha

    print(f"\n There are {count_a} {alpha} atoms, {count_b} {beta} atoms") 
    
    cn_tot, N = av_cn(alpha, beta, n_beta, n_traj_T, n_data)
    cn_tot2 , N2 = av_cn(beta, alpha, n_beta2, n_traj_T, n_data2)

    end = time.time() # end runtime timer
    elapsed = round(end - start , 4)
    print(f"\n runtime for coordination number calculations = {elapsed} s")

    save_files(alpha, beta, data, n_data, cn_tot, N, save_config) # save alpha-beta
    save_files(beta, alpha, data2, n_data2, cn_tot2, N2, save_config) # save beta-alpha

    if T_step ==1:

        print(f"\n *** Compute {alpha}-{beta} bond lifetimes ***")

        start = time.time() # initiate runtime timer
    
        b_lifetime, mean_b_lifetime, median_b_lifetime, max_b_lifetime, min_b_lifetime = bond_lifetime(data)

        end = time.time() # end runtime timer
        elapsed = round(end - start , 4)

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

        fname = alpha + "-" + beta + "_bond_lifetime.dat"
        
        np.savetxt(fname, life_array, delimiter =" ", fmt ="%s")

    return data, data2
    
    
