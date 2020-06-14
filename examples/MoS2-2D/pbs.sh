#PBS -N elastool
#PBS -l nodes=1:ppn=28
#PBS -q default
#PBS -e err

cd $PBS_O_WORKDIR

elastool >log
