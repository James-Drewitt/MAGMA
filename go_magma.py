import numpy as np
import os
import sys
from pathlib import Path

print('\n Go Magma');sys.stdout.flush()

data=np.loadtxt("input_magma.txt", dtype=str)

try:
    shape = data.shape[0]
except IndexError:
    shape = 1

for i in range(shape):
        if shape == 1:
            task = str(data)
        else:
            task = data[i]

        CWD=os.getcwd()
        taskdir=Path(CWD+"/"+task+"/MAGMA")
        
        # Create target directory if not exist
        if not os.path.exists(taskdir):
            os.makedirs(taskdir)
            print("\n Creating Directory " , taskdir, "\n");sys.stdout.flush()
        else:    
            print("\n Creating Directory " , taskdir ,  " ... already exists\n");sys.stdout.flush()

        os.chdir(taskdir)#move into directory to run magma
        print(' Generating config file...\n' );sys.stdout.flush()
        
        inp=r'CONT_path = "../CONTCAR"'
        inp=inp+'\n'+r'xyz = "../out.xyz" # xyz file containing MD trajectories'
        inp=inp+'\n'+r'xyz_numT = 1000 # total number of trajectories to consider'
        inp=inp+'\n'+r'T_step = 1 # number of sampling steps between each trajectory'
        inp=inp+'\n'+r'alpha = "Si" # alpha'
        inp=inp+'\n'+r'beta = "O" # beta'
        inp=inp+'\n'+r'r_cut_offs = [2.3] # cut-off distance'
        inp=inp+'\n'+r'BAD = 1 # compute BAD (1 = True; 0 = False)'
        inp=inp+'\n'+r'n_CN = 1 # consider specific n-centred units (1 = True; 0 = False)'
        inp=inp+'\n'+r'p_CN = [4] # specify partial coordination to consider'
        inp=inp+'\n'+r'save_config = 0 # flag to save configuration files (1 = True; 0 = False)'
        inp=inp+'\n'+r'Qn = 1 # Q speciation analysis (1 = True; 0 = False)'

        config_out = Path('C:/Users/jd14363/OneDrive - University of Bristol/0_Research/00_BLUECRYSTAL_4/MAGMA_CODE/myconfig.py')
        with open(config_out, 'w') as f:
                f.write(inp)

        cmd = 'run_magma.py'
        os.system(cmd)
        
        os.chdir(CWD)#move back to original directory        
