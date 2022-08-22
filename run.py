# run.py

import os, subprocess


gmx = 'gmx'
index = '' #'-n index.ndx'

# Energy minimization run
cmd = f'gmx grompp -f em.mdp -c *_water.gro -p topol.top -o em.tpr'
subprocess.run(cmd, shell=True).check_returncode()
cmd = f'gmx mdrun -v -deffnm em'
subprocess.run(cmd, shell=True).check_returncode()

# NVT run
cmd = f'{gmx} grompp -f nvt.mdp -c em.gro -r em.gro {index} -p topol.top -o nvt.tpr'
subprocess.run(cmd, shell=True).check_returncode()
cmd = f'{gmx} mdrun -deffnm nvt -v'
subprocess.run(cmd, shell=True).check_returncode()

# NPT run
cmd = f'{gmx} grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro {index} -p topol.top -o npt.tpr'
subprocess.run(cmd, shell=True).check_returncode()
cmd = f'{gmx} mdrun -deffnm npt -v'
subprocess.run(cmd, shell=True).check_returncode()

# Remove backup files
os.system('rm \#*') 

# MD run
cmd = f'gmx grompp -f md.mdp -c em.gro -p topol.top -o md.tpr'
subprocess.run(cmd, shell=True).check_returncode()
cmd = f'gmx mdrun -deffnm md -v -cpt 5'
#cmd = 'mpirun -np 4 gmx_mpi mdrun -deffnm PA_water_eq -v -cpt 5 -rdd 2.0 -nb gpu -pme gpu -npme 1 -ntomp 4 &> out.log'
subprocess.run(cmd, shell=True).check_returncode()


