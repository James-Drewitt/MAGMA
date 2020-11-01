
'''
Input analysis parameters here
'''

xyz = "movie.xyz" # xyz file containing MD trajectories
T_step = 1 # number of sampling steps between each trajectory
xyz_numT = 1000 # total number of trajectories to consider
alpha = "Mg" # alpha
beta = "O" # beta
r_cut = 2.3 # cut-off distance
L = 14.0572118759 # box length

BAD = 1 # compute BAD (1 = True; 0 = False)

n_CN = 1 # consider environment n-centred units only (1 = True; 0 = False)
p_CN = 4 # specify partial coordination to consider

save_config = 0 # flag to save configuration files (1 = True; 0 = False)

