#!/usr/bin/env python3
#
# Dr James Drewitt, 27/08/2020. Last update: 02/11/2022
#
import time
import numpy as np
from xyz_CN_subroutines import calc_n_data, av_cn
import os
from pathlib import Path

def get_n_list(beta, alpha, data, data2, p_data, p_data2, p_CN, n_traj, rcut):
    
    n = 0    
    x = 0 # initialise iterators
    y = 0

    n_data2 = [[n_traj, " ", " "]]
    n_data2.append([beta, alpha, "fractional coordination"])

    n_beta_tot = [] #initialise list containing number of beta atoms in specified coordination shell for all trajectories

    a_list_list = []
    a_list_tot = [n_traj] # create array to list alpha atoms in all trajectories for lifetime calculation

    bond_distance = []
    
    for i in range(n_traj):

        n_beta = []
        a_list = []
        n += 1
        n_a = data[n][2]

        num_a = 0
        p_data.append([i+1, " ", " ", " ", " "])
        index_a = len(p_data)-1
        
        for j in range(n_a):
            n += 1
            cn_num = data[n][1] # number of partial coordinations in current trajectory
            if cn_num == p_CN:
                num_a += 1
                a_list.append( data[n+1][0] ) # populate list with all alpha atoms in trajectory with desired coordination
                a_list_list.append( data[n+1][0] )
                for k in range(cn_num+2):
                    p_data.append([ data[n+k][0], data[n+k][1], data[n+k][2], data[n+k][3], data[n+k][4] ])
                for k in range(cn_num):
                    bond_distance.append(data[n+k+2][4])
                p_data.append([" "] * 5 )
            n += cn_num + 2

        a_list_tot.append(num_a) # save number of alpha atoms in each trajectory with desired coordination
        a_list_tot.append(a_list) # append the alpha atom labels to this list
    
        p_data[index_a][2]=num_a

        num_b = 0
        x += 1
        n_b = data2[x][2]
        p_data2.append([i+1, " ", " ", " ", " "])
        index_b = len(p_data2)-1
        for j in range(n_b):
            x += 1
            b_cn_num = data2[x][1] # number of partial coordinations in current trajectory
            at = data2[x+2][0]
            if b_cn_num > 0:
                y+=1
                test = 0
                for z in range(1,b_cn_num+2):
                    at = data2[x+z][0]
                    if at in a_list:
                        test += 1
                if test > 0:
                    for z in range(b_cn_num+2):
                        p_data2.append([ data2[x+z][0], data2[x+z][1], data2[x+z][2], data2[x+z][3], data2[x+z][4] ])
                    n_beta.append(b_cn_num)
                    n_beta_tot.append(b_cn_num)
                    num_b += 1
                        
            x += b_cn_num + 2

            if p_data2[len(p_data2)-1] != [" "] * 5:
                if num_a == 0:
                   p_data2[index_b][2]=0
                elif p_data2[len(p_data2)-1] == [i+1, " ", " ", " ", " "]:
                    continue
                else:
                    p_data2.append([" "] * 5 )
                    p_data2[index_b][2]=num_b
        
        n_data2 = calc_n_data(n_beta, i, n_data2) # partial coordinations beta-alpha

    a_set_list = list(set(a_list_list)) # extract unique entries from a_list_list

    max_bond = []

    for i in range(len(bond_distance)):
        if (i+1)%p_CN==0:
            bonds=bond_distance[i+1-p_CN:i]
            max_bond.append(max(bonds))

    #print(f"max bond = {max_bond}\n")
    
    bins = np.linspace(0, rcut, int(rcut/0.01 +1))
    np_b_dist_hist, np_hist1 = np.histogram(bond_distance,bins) # generate histogram of bond angles

    np.insert(np_b_dist_hist, 0, 0.01)
    return p_data, p_data2, n_data2, n_beta_tot, a_set_list, a_list_tot, bond_distance, np_b_dist_hist

def lifetime(a_set_list, a_list_tot, T_step):
    
    n_traj = a_list_tot[0]
    tot_a = len(a_set_list)

    life = np.zeros((tot_a,n_traj+1), dtype=object)
    life[:,0] = a_set_list[:]

    n = 0 # initialise iterator
    
    if n_traj >= 2:

        for i in range(n_traj): # loop for number of trajectories
            n += 1 
            num_a = a_list_tot[n] # get number of central (alpha) atoms in current trajectory
            a_list = []
            n += 1
            for j in range(num_a):
                unit = a_list_tot[n][j]
                a_list.append(unit)
            for k in range(tot_a):
                if life[k,0] in a_list:
                    life[k,i+1] = 1
    
        life_time = []
    
        for i in range(tot_a):
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

    mean_lifetime = np.mean(life_time)*T_step
    median_lifetime = np.median(life_time)*T_step
    max_lifetime = max(life_time)*T_step
    min_lifetime = min(life_time)*T_step

    return life_time, mean_lifetime, median_lifetime, max_lifetime, min_lifetime

def get_Qn(beta, alpha, data, data2, p_CN, n_traj, rcut):
    
    print("\n *** Calculate Si[4] Q-speciation ***")
    
    n, x, y = 0, 0, 0 # initialise iterators
    
    O_index = []
    O_CN = []
    
    cBO = 0 #initialise array for number of bridging oxygen

    n_Q0, n_Q1, n_Q2, n_Q3, n_Q4 = 0, 0, 0, 0, 0
    
    tot_num_a = 0
    
    for i in range(n_traj):
        
        cBO_total = []

        n += 1
        n_a = data[n][2] # number of Si atoms

        num_a = 0
        
        for j in range(n_a):
            n += 1
            cn_num = data[n][1] # number of partial coordinations in current trajectory
            if cn_num == p_CN:
                num_a += 1
                tot_num_a += 1
                for k in range(2,cn_num+2):
                    data[n+k][0]=data[n+k][0].replace('O','') # strip out index for all oxygen bonded to each Si
                    O_index.append(data[n+k][0]) # array contains oxygen indices for each Si tetrahedron        
            n += cn_num + 2
    
        num_b = 0
        x += 1
        n_b = data2[x][2]


        for j in range(n_b):
            x += 1
            b_cn_num = data2[x][1] # number of O-Si partial coordinations in current trajectory
            y+=1
            O_CN.append(b_cn_num) #extract partial coordination of each oxygen in current trajectory
            num_b += 1
            x += b_cn_num + 2
    
        for i in range(0,len(O_index)-3,4): #add up number of bridging oxygen's per silicon
                  cBO += O_CN[int(O_index[i])-1]-1
                  cBO += O_CN[int(O_index[i+1])-1]-1
                  cBO += O_CN[int(O_index[i+2])-1]-1
                  cBO += O_CN[int(O_index[i+3])-1]-1
                  cBO_total.append(cBO)
                  cBO = 0
                 
    n_Q0 += cBO_total.count(0) / (tot_num_a)
    n_Q1 += cBO_total.count(1) / (tot_num_a)
    n_Q2 += cBO_total.count(2) / (tot_num_a)
    n_Q3 += cBO_total.count(3) / (tot_num_a)
    n_Q4 += cBO_total.count(4) / (tot_num_a)
             
    
    Q_list = ([ "Q0 "+str(round(n_Q0*100,2))+"%", "Q1 "+str(round(n_Q1*100,2))+"%", \
             "Q2 "+str(round(n_Q2*100,2))+"%", "Q3 "+str(round(n_Q3*100,2))+"%", \
             "Q4 "+str(round(n_Q4*100,2))+"%"])              
              
    print(Q_list)
        
    return Q_list


def nCN(p_CN, data, alpha, data2, beta, T_step, save_config, working_dir, rcut, Qn):

    start = time.time() # initiate runtime timer
    print(f"\n *** Consider partial coordination {alpha}-{beta} = {p_CN} ***")

    n_traj = data[0][0]

    str_n="n("+alpha+"-"+beta+") = "+str(p_CN)
    p_data = [[ n_traj , " Configuration for ", str_n, " ", " " ]]
    p_data2 = [[ n_traj , " Configuration for ", str_n, " ", " " ]]
    
    p_data, p_data2, n_data2, n_beta_tot, a_set_list, a_list_tot, bond_distance, np_b_dist_hist = get_n_list(beta, alpha, data, data2, p_data, p_data2, p_CN, n_traj, rcut)

    cn_tot2 , N2 = av_cn(beta, alpha, n_beta_tot, n_traj, n_data2) # average beta-alpha CN

    if Qn == 1:
        Q_list = get_Qn(beta, alpha, data, data2, p_CN, n_traj, rcut)

    end = time.time() # end runtime timer
    elapsed = round(end - start , 4)
    print(f"\n runtime for partial coordination calculations = {elapsed} s")

    start = time.time() # initiate runtime timer
    print(f"\n *** Compute {alpha}{beta}{p_CN} lifetime ***")  

    life_time, mean_lifetime, median_lifetime, max_lifetime, min_lifetime = lifetime(a_set_list, a_list_tot, T_step)

    end = time.time() # end runtime timer
    elapsed = round(end - start , 4)
    print(f" mean lifetime of {alpha}{beta}{p_CN} units = {mean_lifetime} timesteps")
    print(f" median lifetime = {median_lifetime} timesteps")
    print(f" min lifetime = {min_lifetime} timesteps")
    print(f" max lifetime = {max_lifetime} timesteps")
    print(f"\n runtime for lifetime calculations = {elapsed} s")
    
######## DATA OUTPUT ####################

    CWD=os.getcwd()
  
    np_cn_tot2 = np.array(cn_tot2)
    np_cn_tot2 = np.transpose(np_cn_tot2)
    mt_row=np.empty([1, 3], dtype=object)
    for i in range(3):
        mt_row[0,i]=" "

    mt_row2=np.copy(mt_row)
    mt_row2[0,0]="average CN:"
    mt_row2[0,1]=N2

    head=np.copy(mt_row)
    head[0,0]=str_n
    head[0,1]="number"
    head[0,2]="fraction (%)"

    head = np.append(head, mt_row, axis=0)
    np_cn_tot2 = np.append(head, np_cn_tot2, axis=0)
    np_cn_tot2 = np.append(np_cn_tot2, mt_row, axis=0)
    np_cn_tot2 = np.append(np_cn_tot2, mt_row2, axis=0)
    
    if Qn==1:
        np_Q_list=np.array(Q_list)
        Qn_path = Path(CWD+"/"+working_dir+"/Qn.dat")
        np.savetxt(Qn_path, np_Q_list, delimiter=" " , fmt="%s" )

    if T_step == 1:

        foot1=np.copy(mt_row)
        foot1[0,0] = "mean lifetime:"
        foot1[0,1] = mean_lifetime
        foot1[0,2] = "timesteps"

        foot2=np.copy(mt_row)
        foot2[0,0] = "median lifetime:"
        foot2[0,1] = median_lifetime
        foot2[0,2] = "timesteps"

        foot3=np.copy(mt_row)
        foot3[0,0] = "min lifetime:"
        foot3[0,1] = min_lifetime
        foot3[0,2] = "timesteps"

        foot4=np.copy(mt_row)
        foot4[0,0] = "max lifetime:"
        foot4[0,1] = max_lifetime
        foot4[0,2] = "timesteps"
    
        np_cn_tot2 = np.append(np_cn_tot2, mt_row, axis=0)
        np_cn_tot2 = np.append(np_cn_tot2, foot1, axis=0)
        np_cn_tot2 = np.append(np_cn_tot2, foot2, axis=0)
        np_cn_tot2 = np.append(np_cn_tot2, foot3, axis=0)
        np_cn_tot2 = np.append(np_cn_tot2, foot4, axis=0)

    av_CN = Path(CWD+"/"+working_dir+"/"+beta+"-"+alpha+str(p_CN)+"-av.dat")

    print(f" ... saving {av_CN} ...")
    np.savetxt(av_CN , np_cn_tot2 , delimiter=" " , fmt="%s" )

    b_dist_path = Path(CWD+"/"+working_dir+"/"+beta+"-"+alpha+str(p_CN)+"-bond_dist.dat")
    np.savetxt(b_dist_path, np_b_dist_hist, delimiter=" " , fmt="%s" )

    if save_config == 1:

        np_bond_distance = np.array(bond_distance)
        bond_distance_path = Path(CWD+"/"+working_dir+"/"+beta+"-"+alpha+str(p_CN)+"-bond_distance.dat")
        np.savetxt(bond_distance_path , np_bond_distance , delimiter=" " , fmt="%s" )
        
        filename = Path(CWD+"/"+working_dir+"/"+alpha+str(p_CN)+"-"+beta+"-config.dat")
        np_p_data = np.array(p_data)
        print(f"\n ... saving {filename} ...")
        np.savetxt(filename, np_p_data, delimiter=" ", fmt ="%s")

        filename = Path(CWD+"/"+working_dir+"/"+beta+"-"+alpha+str(p_CN)+"-config.dat")
        np_p_data2 = np.array(p_data2)
        print(f" ... saving {filename} ...")
        np.savetxt(filename, np_p_data2, delimiter=" ", fmt ="%s")

        filename = Path(CWD+"/"+working_dir+"/"+beta+"-"+alpha+str(p_CN)+"-traj.dat")
        np_n_data2 = np.array(n_data2)
        print(f" ... saving {filename} ...")
        np.savetxt(filename, np_n_data2, delimiter=" ", fmt ="%s")

    return p_data, p_data2
