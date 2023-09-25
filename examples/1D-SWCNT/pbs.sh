#PBS -N elastool
#PBS -l nodes=1:ppn=28
#PBS -q default
#PBS -e err

cd $PBS_O_WORKDIR
source ~/miniconda3/bin/activate
source /public/apps/intel/bin/compilervars.sh intel64
source /public/apps/intel/composer_xe_2015.0.090/mkl/bin/mklvars.sh intel64
source /public/apps/MPI/impi/5.0.1.035/bin64/mpivars.sh

elastool >log
