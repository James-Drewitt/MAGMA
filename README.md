#### Melt And Glass Molecular Analysis
     __  __   _   ___ __  __   _ 
    |  \/  | /_\ / __|  \/  | /_\             
    | |\/| |/ _ \ (_ | |\/| |/ _ \           
    |_|  |_/_/ \_\___|_|  |_/_/ \_\         
           
    Dr James W E Drewitt, June 2020      
      james.drewitt@bristol.ac.uk

    Last modified: 19/07/2026
       
Post processing script for partial coordination analysis of melt and glass 
molecular dynamics trajectories.

## Serial configuration

Run MAGMA from the project root. To choose a trajectory folder interactively 
in serial mode, use:

> python go_magma.py --select

Add '--confirm' to review and edit the configuration before the run:

> python go_magma.py --select --confirm

Choose the folder containing `CONTCAR` and the XYZ trajectory.

## Batch configuration

To run in batch mode, put `input_magma.txt` in the folder from which you run 
the command, list task folders (one per line), then run 

> python go_magma.py

All relative task and configuration paths are resolved from the current 
launch directory.

By default, a batch run uses a shared 'magma_config.txt' in the location 
specified in `input_magma.txt`.

```
config_folder = .
data/001
data/002
data/003
```

The configuration used during the batch run will, therefore, be 
`config_folder/magma_config.txt`. This is the simplest setup when all datasets 
have the same atom types, trajectory filename, and analysis settings.

For mixed datasets, or when every dataset requires its own settings, place an 
independent `magma_config.txt` file in the data directory and overide the 
shared config using the command:

> python go_magma.py --config-override

The `--config-override` (alias: `-o`) will use the each data folder's 
`magma_config.txt` file. If this is absent, MAGMA falls back to the shared 
configuration within the project root.

MAGMA will generate a default configuration when no configuration file exists. 
To create a shared configuration file for a batch run, copy 
`magma_config_template.txt` to the project root, for a shared batch
configuration, or into a data folder for a dataset-specific configuration and
rename it to `magma_config.txt`.

### Analysis Options

MAGMA will always calculate the alpha--beta coordination number and bond
lifetimes for every value in `r_cut_offs`. Optional analysis settings include:

`bond_angle_distribution [default = true]` calculates the beta--alpha--beta and 
alpha--beta--alpha bond-angle distributions for the selected alpha--beta pair. 

`partial_coordination [default = 4]` Values may be comma-sperated, e.g., 3, 4.
This identifies units with the specified number of beta atoms around each alpha
atom, and write their coordination, bond-distance, lifetime, and (if enabled) 
bond angles. For example, use  `4` for SiO4 tetrahedra or `3` for planar CO3 
groups. Use an empty value (`partial_coordination =`) to skip this analysis.

`q_speciation [default = false]`. Set to true to calculate Q0--Q4 populations 
for four-fold alphaO4 units. It may be used with any `alpha` atom type, but 
requires `beta = O` and that `partial_coordination` contains `4`.

`carbonate_enabled [default = false]` enables the carbonate-specific CO3+1 and
C2O5 analyses. For this analysis, use `alpha = C` and `beta = O`. The C--O 
cutoff, homopolar C--C cutoff, and angular tolerances are configurable.

CO3+1 is a trigonal-planar carbonate with three short C--O bonds and one long,
approximately axial C--O bond.

C2O5 is defined as either a pair of planar CO3 groups sharing one short C--O--C 
bridge, or a short homopolar C--C pair with three and two distinct short C--O 
neighbours (O3--C=C--O2).

