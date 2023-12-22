#!/bin/sh
#An example for serial job.
#SBATCH -J Cresil_EJM
#SBATCH -o /home/lifesci/luosongwen/eccDNA/stimulate_data/test_data/logs/Cresil_EJM_test_%j.log
#SBATCH -e /home/lifesci/luosongwen/eccDNA/stimulate_data/test_data/logs/Cresil_EJM_test_%j.err
#SBATCH --qos=aepmemqos
#SBATCH --exclusive
#SBATCH --time=240:00:00
#SBATCH -N 1 -n 40 -p 2TB-AEP-Mem
echo Running on $SBATCH_PARTITION paratation
echo Start time is `date`
source /home/lifesci/luosongwen/miniconda3/etc/profile.d/conda.sh
conda activate ecc_pipe_old
#echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu core


cd /gpfs/home/lifesci/luosongwen/eccDNA/ecc_pipe
python3 ecc_pipe_master.py --Detect --tool cresil -n 40 --config /home/lifesci/luosongwen/eccDNA/ecc_pipe/config/cresil_config_EJM_test.yaml

echo End time is `date`