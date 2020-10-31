
'''
Input analysis parameters here
'''

xyz = "movie.xyz" # xyz file containing MD trajectories
T_step = 10 # number of sampling steps between each trajectory
alpha = "Mg" # alpha
beta = "O" # beta
r_cut = 2.93 # cut-off distance (MgO 2.93)
L = 12.733 # box length

BAD = 1 # compute BAD (1 = True; 0 = False)

n_CN = 1 # consider environment n-centred units only (1 = True; 0 = False)
p_CN = 4 # specify partial coordination to consider

save_config = 0 # flag to save configuration files (1 = True; 0 = False)

