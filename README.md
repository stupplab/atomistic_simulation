This repository contains all you need top run the atomistic simulations of PAs.

`script.py` is the main script that is run first. Look inside to see how it works and how the parameters are defined.
It will create a folder simulation for the simulation instance of the defined parameters and automatically interacts with the command line input using subprocess Popen. 

Essential folder/files (basically means that you don't need to change them):
- charmm36-jul2020.ff
- quaternion.py
- utils.py
- top_all36_prot_forPAs.rtf
- PA_generic.pdb (for C16 PAs and peptides)

Custom files (You probably need to tune them based on your use case)
- script.py
- *.mdp
- jobscript.sh
- run.py / rerun.py

# Note

Remember to add the custom made PA residues `C16` `C12` `12C` `SPI` to `residuetypes.dat`. On macos, this should be located as 
`/usr/local/gromacs/share/gromacs/top/residuetypes.dat`\
In this file add lines
```
C16 Protein
C12 Protein
12C Protein
SPI Protein
```


For a new molecule (`drug.mol2`) generate the force field using the following site \
https://cgenff.umaryland.edu/initguess/
It will generate `drug.str` file containing the force field parameters.
To generate the GROMACS comapatible parameters files, download `cgenff_charmm2gmx.py` from https://cgenff.umaryland.edu/commonFiles/utility.php \
Used as `./cgenff_charmm2gmx.py <RESNAME> drug.mol2 drug.str charmm36.ff` 

