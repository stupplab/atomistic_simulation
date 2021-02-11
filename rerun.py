# rerun.py

import os, subprocess

# Remove backup files
os.system('rm \#*') 

gmx = 'gmx'

# MD re-run
cmd = f'{gmx} grompp -f md.mdp -c npt.gro -t md.cpt -p topol.top -o md.tpr'
subprocess.run(cmd, shell=True).check_returncode()
cmd = f'{gmx} mdrun -deffnm md -v -cpt 5'
#cmd='mpirun -np 2 gmx_mpi mdrun -deffnm PA_water_eq -v -cpt 5 -cpi PA_water_eq.cpt -append -nb gpu -pme gpu -npme 1 -ntomp 4 &> out.log'
subprocess.run(cmd, shell=True).check_returncode()