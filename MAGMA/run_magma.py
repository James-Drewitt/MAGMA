##
## AUTHOR: DR JAMES DREWITT
##
## james.drewitt@bristol.ac.uk

Mod_date = "01/11/2020"

from wel import welcome
from xyz_analysis import CN
from xyz_bad import bad
from xyz_ncn import nCN
import os
from pathlib import Path
from myconfig import * # import analysis parameters

welcome(Mod_date)

'''
Calculate coordination numbers

outputs:

data   -> list of alpha atoms with all beta atom coordinates within coordination shell
n_data -> list of alpha-centred partial coordinations for each trajectory
cn_tot -> average alpha-centred partial coordinations
n_a  -> number of alpha atoms

data2   -> list of beta atoms with all alpha atom coordinates within coordination shell
n_data2 -> list of beta-centred partial coordinations for each trajectory
cn_tot2 -> average beta-centred partial coordinations
n_b  -> number of beta atoms
'''

f = open(CONT_path)
lines=f.readlines()
f.close()
boxlen, *rest = lines[2].rstrip().split('0.0')
L=round(float(boxlen),9)
V=round(L**3,4)

print(f"\n box length = {L} Ang; volume = {V} Ang**3")

for i in range(len(r_cut_offs)):
    r_cut = r_cut_offs[i]
    working_dir = alpha+"-"+beta+"_"+str(r_cut)
    data, data2, = CN(xyz, T_step, alpha, beta, r_cut, L, save_config, xyz_numT, working_dir)

#Calculate bond angle distributions

    if BAD == 1:
        bad(data, alpha, data2, beta, 0, p_CN, L, save_config, working_dir)

#Calculate bond lengths, beta-alpha CN, and BAD for specific partial coordinations

    if n_CN == 1:
        for j in range(len(p_CN)):
            p_data, p_data2 = nCN(p_CN[j], data, alpha, data2, beta, T_step, save_config, working_dir, r_cut)
            if BAD == 1:
                bad(p_data, alpha, p_data2, beta, 1, p_CN[j], L, save_config, working_dir)
