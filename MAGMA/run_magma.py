#!/usr/bin/env python3
#
##
## AUTHOR: DR JAMES DREWITT, 06/11/2020
##
## james.drewitt@bristol.ac.uk

Mod_date = "19/07/2026"

from wel import welcome
from xyz_analysis import CN
from xyz_bad import bad
from xyz_ncn import nCN
from xyz_carbonate import carbonate_units
import os
from types import ModuleType

config_text = os.environ.get("MAGMA_CONFIG_TEXT")
if config_text is None:
    raise RuntimeError("Run go_magma.py to start MAGMA.")
config = ModuleType("magma_runtime_config")
exec(compile(config_text, "<MAGMA runtime configuration>", "exec"), config.__dict__)

CONT_path = config.CONT_path
xyz = config.xyz
xyz_numT = config.xyz_numT
T_step = config.T_step
alpha = config.alpha
beta = config.beta
r_cut_offs = config.r_cut_offs
save_detailed_analysis_data = config.save_detailed_analysis_data

analysis = config.analysis
BAD = bool(analysis.get("bond_angle_distribution", False))
p_CN = list(analysis.get("partial_coordination", []))
Qn = bool(analysis.get("q_speciation", False))
carbonate_options = analysis.get("carbonate", {})

carbonate = bool(carbonate_options.get("enabled", False))
carbonate_short_cutoff = carbonate_options.get("short_cutoff", 2.0)
carbonate_CC_cutoff = carbonate_options.get("CC_cutoff", 2.0)
carbonate_long_min = carbonate_options.get("long_min", 2.4)
carbonate_long_max = carbonate_options.get("long_max", 2.9)
carbonate_angle_tolerance = carbonate_options.get("angle_tolerance", 20.0)
carbonate_axial_tolerance = carbonate_options.get("axial_tolerance", 20.0)

if Qn and beta != "O":
    raise ValueError("Q-speciation requires beta='O'.")
if Qn and 4 not in p_CN:
    raise ValueError("Q-speciation requires 4 in analysis['partial_coordination'].")

welcome(Mod_date)

'''
Calculate coordination numbers

outputs:

data    -> list of alpha atoms with all beta atom coordinates within coordination shell
n_data  -> list of alpha-centred partial coordinations for each trajectory
cn_tot  -> average alpha-centred partial coordinations
n_a     -> number of alpha atoms

data2   -> list of beta atoms with all alpha atom coordinates within coordination shell
n_data2 -> list of beta-centred partial coordinations for each trajectory
cn_tot2 -> average beta-centred partial coordinations
n_b     -> number of beta atoms
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
    data, data2, = CN(xyz, T_step, alpha, beta, r_cut, L, save_detailed_analysis_data, xyz_numT, working_dir)

#Calculate bond angle distributions

    if BAD == 1:
        bad(data, alpha, data2, beta, 0, p_CN, L, save_detailed_analysis_data, working_dir)

#Calculate bond lengths, beta-alpha CN, and BAD for specific partial coordinations

    if p_CN:
        for partial_cn in p_CN:
            p_data, p_data2 = nCN(
                partial_cn, data, alpha, data2, beta, T_step, save_detailed_analysis_data,
                working_dir, r_cut, Qn and partial_cn == 4,
            )
            if BAD == 1:
                bad(p_data, alpha, p_data2, beta, 1, partial_cn, L, save_detailed_analysis_data, working_dir)

    if carbonate == 1:
        carbonate_units(
            xyz, T_step, alpha, beta, L, xyz_numT, working_dir,
            carbonate_short_cutoff, carbonate_CC_cutoff, carbonate_long_min, carbonate_long_max,
            carbonate_angle_tolerance, carbonate_axial_tolerance,
        )
