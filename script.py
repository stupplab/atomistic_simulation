"""Basic Gromacs script for atomistic simulation

All in all you need the following files to run the conplete simulation:
    - itp and pdb files for the molecule, 
    - topol.top file that included the itp files of forcefield, our molecule, water, ions
    - mdp files - ions,em,nvt,npt,md, and 
    - forcefield - ex. charmm36-jul2020.ff
    - Custom python scripts: script.py, run.py, rerun.py


For a new molecule, use SwissParam.ch, it takes mol2 input and 
generate CHARMM - GROMACS compatible files
mol2 files are usually available online or can be created using software Avogadro

This script is tested using gromacs 2020.2
"""



import numpy as np
import os, subprocess

import utils

############################# Simulation parameters #############################

molname = 'PA'
PA_seq             = 'C16VVAAEE' #+ 'Z' # X:C12, Z:12C
num_PA             = 22*9
# residuecharge_PA1  = [['E',0], ['E',-1]]
Lx                 = 12
Ly                 = 12
Lz                 = 12
topfile            = 'topol.top'
PA_generic='PA_generic.pdb'

################################################################################



#################################### Setup ####################################

# Create a directory path for this example and go inside it
path = os.getcwd()+'/simulation'
if not os.path.exists(path):
    os.mkdir(path)
os.chdir(path)

utils.generic_to_specific_PA(PA_seq, PA_generic, molname)

utils.make_aa_pdb(molname)


#-------------------------------------------------------------------------------
# Copy relevant files
os.system(f'cp -r ../charmm36-jul2020.ff ./')
os.system(f'cp -r ../*.mdp ./')
os.system(f'cp -r ../*run.py ./')
#-------------------------------------------------------------------------------


#------------------------------------------------------------------------------- 

# Generate gro, itp and top files using the pdb file. 
# Remember this assumes that the force field identifies all the residues in the pdb
cmd = f'gmx pdb2gmx -f {molname}_aa.pdb -o {molname}.gro -ter -p {topfile} -i {molname}.itp -water spce -ff charmm36-jul2020'
subprocess.run(cmd, shell=True).check_returncode()


# Convert the pdb to gro
# cmd = f'gmx editconf -f {molname1}.pdb -o {molname1}.gro'
# subprocess.run(cmd, shell=True).check_returncode()


# Randomly insert molecules in a box, generates <molname>_box.gro file
cmd = f'gmx insert-molecules -box {Lx} {Ly} {Lz} -nmol {num_PA} -ci {molname}.gro -radius 0.2 -o {molname}_box.gro'
subprocess.run(cmd, shell=True).check_returncode()
input(f'>>> Manually correct the number of {molname} molecules in the top file to how many are added (See above). \n \
        Then press ENTER to continue...')



# Change the random positions of the PA to an initial configuration
num_atoms = utils.get_num_atoms_fromGRO('PA.gro')
num_molecules = utils.get_num_molecules_fromtop('topol.top', molname='Protein')
utils.init_fiber_config(
    'PA_box.gro', 
    num_atoms, num_molecules,
    Lx,Ly,Lz,
    start_from_nth_atom=0,
    invert=False,
    C_indices=[0,2])



# Solvate the box, generates <molname>_water.gro
cmd = f'gmx solvate -cp {molname}_box.gro -cs spc216.gro -p topol.top -o {molname}_water.gro'
subprocess.run(cmd, shell=True).check_returncode()
input(f'>>> Check the number of solvent molecules added in the top file is correct (See above). \n \
        Then press ENTER to continue...')

# Add ions
cmd = f'gmx grompp -f ions.mdp -c {molname}_water.gro -p topol.top -o ions.tpr'
subprocess.run(cmd, shell=True).check_returncode()
# Make changes in the below cmd to add extra ions apart from what needed for neutralization)
cmd = f'gmx genion -s ions.tpr -o {molname}_water.gro -p topol.top -pname NA -nname CL -neutral'
subprocess.run(cmd, shell=True).check_returncode()
input(f'>>> Check if the number of ions added in the top file is correct (See above). \n \
        Then press ENTER to continue...')



# Indexing to group multiple molecules. Add -n index.ndx to all the run commands
# USually Sol and Ions will be grouped together. 
# Also check tc-grps in the mdp files nvt,npt,md for correctness of residue groups included
# cmd = f'gmx make_ndx -f em.gro -o index.ndx'
# subprocess.run(cmd, shell=True).check_returncode()
# input(f'>>> In <nvt/npt/md>.mdp verify that tc-grps has the correct coupling groups as formed in index.ndx (See above). \n \
#         Then press ENTER to continue...')

# Remove backup files
os.system('rm \#*') 

#-------------------------------------------------------------------------------


