#!/bin/bash

# This function writes a slurm script. 
# We can call it with different parameter 
# settings to create different experiments

function write_script
{
TASKS=$[$NPERNODE * $NODES]
JOB_NAME=$(printf 'cct-%s-%04d' ${CIRCUIT} ${TASKS})
DIR_NAME=$(printf '%s-N%03d-n%02d' ${JOB_NAME} ${NODES} ${NPERNODE})

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
####SBATCH --partition=${QUEUE}
#SBATCH --nodes=${NODES}
#SBATCH --ntasks-per-node=${NPERNODE}
#SBATCH --distribution=cyclic
#SBATCH --time=48:00:00

#SBATCH --partition=long 
#SBATCH --res=apache



module switch gcc/4.6.2 gcc/4.7.3

export LD_LIBRARY_PATH=${LIB_PATH}
export PV_PLUGIN_PATH=${PLUGIN_PATH}

export OMP_NUM_THREADS=1

aprun -n ${TASKS} -N ${NPERNODE} ${EXECUTABLE} ${PYSCRIPT} -c ${CONFIG} -t ${CIRCUIT} -n ${NEURONS} -d ${BASEDIR}/${DIR_NAME}/ -p ${PLUGIN_PATH}

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

#
# settings for paraview on todi
#
QUEUE=night
EXECUTABLE=/project/csvis/biddisco/todi/build/pv4/bin/pvbatch
PYSCRIPT=/project/csvis/biddisco/src/plugins/pv-bbp/scripts/5K-voltage-differentials.py
CONFIG="/project/csvis/biddisco/bbpdata/egpgv/centralV.cfg"
CIRCUIT=""
LIB_PATH=/apps/todi/hdf5/1.8.12-snap3/gnu_47/lib/:/project/csvis/biddisco/todi/build/pv4/lib:/opt/cray/mpt/5.6.4/gni/mpich2-gnu/47/lib:/opt/gcc/mpc/0.8.1/lib:/opt/gcc/mpfr/2.4.2/lib:/opt/gcc/gmp/4.3.2/lib:/opt/gcc/4.7.3/snos/lib64:/apps/todi/boost/1.51.0/gnu_471/lib
PLUGIN_PATH=/project/csvis/biddisco/todi/build/plugins/bin/

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
	for NODES in 128
	do
		write_script
	done

done
