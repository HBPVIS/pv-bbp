#!/bin/bash

# This function writes a slurm script. 
# We can call it with different parameter 
# settings to create different experiments

function write_script
{
TASKS=$[$NPERNODE * $NODES]
JOB_NAME=$(printf 'cct-%s-%04d' ${CIRCUIT} ${TASKS})
DIR_NAME=$(printf '%s-N%03d-n%02d' ${JOB_NAME} ${NODES} ${NPERNODE})
TASKS_PER_GPU=$[$TASKS / 2]

if [ -f $DIR_NAME/timing-full.txt ] ; then
	echo "$DIR_NAME/timing-full.txt already exists, skipping..."
	return 0
else
	echo "Creating job $DIR_NAME"
fi

mkdir -p $DIR_NAME

cat << _EOF_ > ${DIR_NAME}/slurm-exp.bash
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=${QUEUE}
#SBATCH --nodes=${NODES}
#SBATCH --ntasks-per-node=${NPERNODE}
#SBATCH --distribution=cyclic
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --gres=gpu:2

module load turbovnc
module load virtualgl
module load cmake
module load gcc/4.8.1
module load boost/1.53.0
module load gcc/4.8.1

export LD_LIBRARY_PATH=/scratch/castor/biddisco/build/pv4/lib:/apps/castor/mvapich2/1.9-gcc-4.8.0/mvapich2-gnu/lib:/scratch/castor/biddisco/apps/hdf5_1_8_cmake/lib/:/scratch/castor/biddisco/apps/bbp/lib:$LD_LIBRARY_PATH

export PV_PLUGIN_PATH=/scratch/castor/biddisco/build/plugins/bin
export OMP_NUM_THREADS=1
export HYDRA_LAUNCHER_EXTRA_ARGS=--gres=gpu:2 --startx
#export CUDA_PCI_DEVICES="2:0:0,3:0:0,132:0:0"
export CUDA_VISIBLE_DEVICES="0,1"

/apps/castor/mvapich2/1.9-gcc-4.8.0/mvapich2-gnu/bin/mpiexec -binding rr -ppn 1 -n ${TASKS_PER_GPU} -env DISPLAY :0.0 ${EXECUTABLE} ${PYSCRIPT} -t ${CIRCUIT} -n ${NEURONS} -p ${BASEDIR}/${DIR_NAME}/ : -n ${TASKS_PER_GPU} -env DISPLAY :0.1 ${EXECUTABLE} ${PYSCRIPT} -t ${CIRCUIT} -n ${NEURONS} -p ${BASEDIR}/${DIR_NAME}/ 

# --param1 ${I} --param2 ${J} --param3 ${K} 

_EOF_

chmod 775 ${DIR_NAME}/slurm-exp.bash

echo "cd ${DIR_NAME}; sbatch slurm-exp.bash; cd \$BASEDIR" >> run_jobs.bash

}

# get the path to this generate script, works for most cases
pushd `dirname $0` > /dev/null
BASEDIR=`pwd`
popd > /dev/null
echo "Generating jobs using base directory $BASEDIR"

# Create another script to submit all generated jobs to the scheduler
echo "#!/bin/bash" > run_jobs.bash
echo "BASEDIR=$BASEDIR" >> run_jobs.bash
chmod 775 run_jobs.bash

# set some vars which are fixed in this test
QUEUE=production
EXECUTABLE=/scratch/castor/biddisco/build/pv4/bin/pvbatch
PYSCRIPT=/project/csvis/biddisco/src/plugins/pv-bbp/scripts/5k-castor.py
CIRCUIT=

# Loop through all the parameter combinations generating jobs for each
for NEURONS in 5000
do
	if [ $NEURONS -eq 1000 ]; then
	  CIRCUIT=1K
	elif [ $NEURONS -eq 2000 ]; then
	  CIRCUIT=2K
	elif [ $NEURONS -eq 5000 ]; then
	  CIRCUIT=5K
	fi	
	
	NPERNODE=6
	for NODES in 17
	do
		write_script
	done

done
