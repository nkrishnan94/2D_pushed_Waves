import numpy as np
import glob
import os 
from datetime import datetime

Bcoops = np.array([1,2,2,2.5,3,3.5,4])
MutFs = np.array([10,25,50,75,100,125,150,175])
samps = 10
files= []
time_str=datetime.now()
flogstr ="submit_log_1"+".log"
flog = open(flogstr,"a")
flog.write('date/time:\n')
flog.write(str(time_str)+'\n')
flog.write('B Values:\n')
flog.write(str(Bcoops)+'\n')
flog.write('Initial Mut. :\n')
flog.write(str(MutFs)+'\n')
ID =0

for B,Bcoop in enumerate(Bcoops) :
	for M,MutF in enumerate(MutFs):
		param_str="B"+str(Bcoop)+"_I"+str(MutF)
		f = open("run_het_sims_sde_"+param_str+".slurm","a")
		files.append("run_het_sims_sde_"+param_str+".slurm")
		f.write('#!/bin/bash\n')
		f.write('#SBATCH -J run_sims'+param_str+'\n')
		f.write('#SBATCH -A FUSCO-SL3-CPU\n')
		f.write('#SBATCH --nodes=1\n')
		f.write('#SBATCH --ntasks=1\n')
		f.write('#SBATCH --cpus-per-task=1\n')


		f.write('#SBATCH --time=00:10:00\n')


		f.write('#SBATCH --mem=5980mb\n')
		f.write('#SBATCH --array=1-'+str(samps)+'\n')
		f.write('#SBATCH -p skylake\n')
		f.write('./etc/profile.d/modules.sh\n')
		f.write('module purge \n')
		f.write('module load rhel7/default-peta4\n')
		f.write('echo "This is job" $SLURM_ARRAY_TASK_ID\n')
		#f.write('g++ -o outanc eff_pop_2d_radial_neutral.cpp -std=c++11\n')
		f.write('./out -T 650 -B '+str(Bcoop)+' -I'  +  str(MutF) +  ' -G $SLURM_ARRAY_TASK_ID  \n')
		f.close()

for f in files:
	B = f.split('B')[1].split('_')[0]
	m= f.split('_I')[1].split('_')[0]

	flog.write("B: "+str(B) +" Init. Mut. : " +str(m)+ '\n')
	os.system('sbatch '+str(f)+ ' >> '+flogstr)
os.system('mkdir slurm_scripts\n')
os.system('mv *.slurm slurm_scripts/\n')

		


