
'''
Input analysis parameters here
'''

CONT_path = "../CONTCAR"
xyz = "../rings/data/out.xyz" # xyz file containing MD trajectories
xyz_numT = 50 # total number of trajectories to consider
T_step = 1 # number of sampling steps between each trajectory
alpha = "C" # alpha
beta = "O" # beta
r_cut_offs = [2.5, 2.6, 2.7, 2.8] # cut-off distance

BAD = 1 # compute BAD (1 = True; 0 = False)

n_CN = 1 # consider specific n-centred units (1 = True; 0 = False)
p_CN = [3,4] # specify partial coordination to consider

save_config = 0 # flag to save configuration files (1 = True; 0 = False)

