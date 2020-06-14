#PBS -N 2D
#PBS -l nodes=node25:ppn=28
#PBS -q default
#PBS -e err

cd $PBS_O_WORKDIR

# For intel mpi
source /public/apps/intel/bin/compilervars.sh intel64
source /public/apps/intel/composer_xe_2015.0.090/mkl/bin/mklvars.sh intel64
source /public/apps/MPI/impi/5.0.1.035/bin64/mpivars.sh
source ~/miniconda3/bin/activate

elastool >log
