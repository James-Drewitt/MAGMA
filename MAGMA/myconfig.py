CONT_path = "../CONTCAR"
xyz = "../out.xyz" # xyz file containing MD trajectories
xyz_numT = 100 # total number of trajectories to consider
T_step = 1 # number of sampling steps between each trajectory
alpha = "Si" # alpha
beta = "O" # beta
r_cut_offs = [2.3] # cut-off distance
BAD = 1 # compute BAD (1 = True; 0 = False)
n_CN = 1 # consider specific n-centred units (1 = True; 0 = False)
p_CN = [4] # specify partial coordination to consider
save_config = 0 # flag to save configuration files (1 = True; 0 = False)
Qn = 1 # Q speciation analysis (1 = True; 0 = False)
carbonate = 0 #