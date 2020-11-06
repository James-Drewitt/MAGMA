import numpy as np
import os
import sys
from pathlib import Path

print(' Go Magma');sys.stdout.flush()
data=np.loadtxt("input_magma.txt", dtype=np.str)
for i in range(data.shape[0]):
        task=data[i]
        print(task)
        CWD=os.getcwd()
        taskdir=Path(CWD+"/"+task+"/MAGMA")
        
        # Create target directory if not exist
        if not os.path.exists(taskdir):
            os.makedirs(taskdir)
            print("Creating Directory " , taskdir);sys.stdout.flush()
        else:    
            print("Creating Directory " , taskdir ,  " ... already exists");sys.stdout.flush()

        os.chdir(taskdir)#move into directory to run magma
        print(' Generating config file...');sys.stdout.flush()
        
        inp=r'CONT_path = "../CONTCAR"'
        inp=inp+'\n'+r'xyz = "../rings/data/out.xyz" # xyz file containing MD trajectories'
        inp=inp+'\n'+r'xyz_numT = 25000 # total number of trajectories to consider'
        inp=inp+'\n'+r'T_step = 1 # number of sampling steps between each trajectory'
        inp=inp+'\n'+r'alpha = "C" # alpha'
        inp=inp+'\n'+r'beta = "O" # beta'
        inp=inp+'\n'+r'r_cut_offs = [2.3] # cut-off distance'
        inp=inp+'\n'+r'BAD = 1 # compute BAD (1 = True; 0 = False)'
        inp=inp+'\n'+r'n_CN = 1 # consider specific n-centred units (1 = True; 0 = False)'
        inp=inp+'\n'+r'p_CN = [3,4] # specify partial coordination to consider'
        inp=inp+'\n'+r'save_config = 0 # flag to save configuration files (1 = True; 0 = False)'

        config_out = Path('C:/Users/jd14363/OneDrive - University of Bristol/0_Research/00_BLUECRYSTAL_4/MAGMA_CODE/myconfig.py')
        with open(config_out, 'w') as f:
                f.write(inp)

        cmd = 'run_magma.py'
        os.system(cmd)
        
        os.chdir(CWD)#move back to original directory        

        
