#SBATCH -A p30009               # Allocation
#SBATCH -p gengpu               # Queue
#SBATCH --gres gpu:2
#SBATCH -t 10:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 2                    # Number of Cores (Processors)
#SBATCH -J "C16V2A2E2-atomistic"            # Name of job

# unload any modules that carried over from your command line session
#module purge


# load modules you need to use

np=2

# run
gmx grompp -f em.mdp -c *_water.gro -p topol.top -o em.tpr
mpirun -np $np gmx_mpi mdrun -v -deffnm em -ntomp 4
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
mpirun -np $np gmx_mpi mdrun -v -deffnm nvt -ntomp 4
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr
mpirun -np $np gmx_mpi mdrun -v -deffnm npt -ntomp 4
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
mpirun -np $np gmx_mpi mdrun -deffnm md -v -cpt 5 -rdd 2.0 -nb gpu -pme gpu -npme 1 -ntomp 4 &> out.log

# re-run
#mpirun -np $np gmx_mpi mdrun -deffnm md -v -cpt 5 -cpi md.cpt -append -nb gpu -pme gpu -npme 1 -ntomp 4 &> out.log
