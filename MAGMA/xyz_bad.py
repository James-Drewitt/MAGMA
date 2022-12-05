#!/usr/bin/env python3
#
# Dr James Drewitt, 17/08/2020. Last update: "02/11/2022"
#
import numpy as np
import time
from math import acos, degrees
import os
from pathlib import Path

def pbc(coord1, coord2, L):

    '''
    Calculates the distance between two coordinates within periodic boundary conditions
    coord1 -> coordinate 1
    coord2 -> coordinate 2
    L-> configuration box length
    '''

    if coord1 < 0:          #First check for negative coordinates
        coord1 = L + coord1
    if coord2 < 0:
        coord2 = L + coord2 
    d1 = coord2 - coord1 #find minimum distance within periodic boundary conditions
    d2 = coord2 - coord1 + L
    d3 = coord2 - coord1 - L
    d = min(abs(d1), abs(d2), abs(d3))

    if abs(d1) == d:
        d = d1
    elif abs(d2) == d:
        d = d2
    elif abs(d3) == d:
        d = d3
        
    return d

def bad_calc(data, alpha, beta, L):

    n_traj = data[0][0]
    n=0 # initialise iterator
    BAD_list = [[0]*4]
    BADs = [0]
    for i in range(n_traj):
        n += 1
        num_a = data[n][2]
        for j in range(num_a):
            n += 1
            cn_num = data[n][1] # number of partial coordinations in current trajectory
            if cn_num == 0:
                n += 2
            elif cn_num == 1:
                n += 3
            elif cn_num >= 2:
                n += 1
                coord_list = [[0] * 4 for i in range(cn_num+1)] # initialise list of beta coordinates
                for k in range(cn_num+1):
                    coord_list[k][0] = data[n][0] # current atom        
                    coord_list[k][1] = data[n][1] # x
                    coord_list[k][2] = data[n][2] # y
                    coord_list[k][3] = data[n][3] # z
                    n += 1

                for a in range(1, cn_num):
                    for b in range(a+1, cn_num+1):
                        x0 = coord_list[0][1]
                        x1 = coord_list[a][1]
                        x2 = coord_list[b][1]
                        x1 = pbc(x0, x1, L)
                        x2 = pbc(x2, x0, L)
                        rx_dp = x1*x2
                        

                        y0 = coord_list[0][2]
                        y1 = coord_list[a][2]
                        y2 = coord_list[b][2]
                        y1 = pbc(y0, y1, L)
                        y2 = pbc(y2, y0, L)
                        ry_dp = y1*y2

                        z0 = coord_list[0][3]
                        z1 = coord_list[a][3]
                        z2 = coord_list[b][3]
                        z1 = pbc(z0, z1, L)
                        z2 = pbc(z2, z0, L)
                        rz_dp = z1*z2

                        r_dp = rx_dp + ry_dp + rz_dp # dot product ra.rb

                        ra_mag = np.sqrt( np.square(x1) + np.square(y1) + np.square(z1) ) # magnitude |ra|
                        rb_mag = np.sqrt( np.square(x2) + np.square(y2) + np.square(z2) ) # magnitude |rb|

                        cos_theta = -r_dp / (ra_mag * rb_mag) # ra.rb / |ra||rb|

                        if (cos_theta > 1.0): # handle rounding errors
                            cos_theta = 1.0
                        elif (cos_theta < -1.0):
                            cos_theta = -1.0

                        theta = degrees(acos(cos_theta)) #Bond angle in degrees

                        BAD_list.append([ coord_list[a][0], coord_list[0][0], coord_list[b][0], theta])
                        BADs.append(theta)
    BAD_list=BAD_list[1:][:] # remove empty first line
    BADs=BADs[1:][:] # remove empty first line
    bins = np.linspace(0, 180, 181)
    np_hist, np_hist1 = np.histogram(BADs,bins) # generate histogram of bond angles

    return(BAD_list, np_hist)

def bad(data, alpha, data2, beta, n_CN, p_CN, L, save_config, working_dir):

    start = time.time() # initiate runtime timer
    print("\n *** Calculating bond angle distributions ***")

    BAD_list, hist = bad_calc(data, alpha, beta, L) # calc beta-alpha-beta angles
    BAD_list2, hist2 = bad_calc(data2, beta, alpha, L) # calc alpha-beta-alpha angles

    end = time.time() # end runtime timer
    elapsed = round(end - start , 4)
    print(f"\n runtime for bond angle calculations = {elapsed} s")

    save_files(1, alpha, beta, BAD_list, hist, n_CN, p_CN, save_config, working_dir)
    save_files(2, beta, alpha, BAD_list2, hist2, n_CN, p_CN, save_config, working_dir)

def save_files(meth, alpha, beta, BAD_list, np_hist, n_CN, p_CN, save_config, working_dir):

    CWD=os.getcwd()
    
    np_BAD_list = np.array(BAD_list)

    if n_CN == 1:
        alpha2 = alpha
        beta2 = beta
        if meth == 1:
            alpha2 = alpha+str(p_CN)
        if meth == 2:
            beta2 = beta+str(p_CN)
        BAD_file = Path(CWD+"/"+working_dir+"/"+beta2+"-"+alpha2+"-"+beta2+"_BAD.dat")
        BAD_hist = Path(CWD+"/"+working_dir+"/"+beta2+"-"+alpha2+"-"+beta2+"_BAD_hist.dat")
    else:
        BAD_file = Path(CWD+"/"+working_dir+"/"+beta+"-"+alpha+"-"+beta+"_BAD.dat")
        BAD_hist = Path(CWD+"/"+working_dir+"/"+beta+"-"+alpha+"-"+beta+"_BAD_hist.dat")

    print(f" ... saving {BAD_hist}...")
    np.savetxt(BAD_hist , np_hist , delimiter=" " , fmt="%s")

    if save_config == 1:
        print(f"\n ... saving {BAD_file} ...")
        np.savetxt(BAD_file , np_BAD_list , delimiter=" " , fmt="%s")
    
