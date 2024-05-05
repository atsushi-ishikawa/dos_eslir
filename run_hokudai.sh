#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=n22240a
#PJM -L node=1
#PJM --mpi proc=40
#PJM -L elapse=20:00:00 
#PJM -g n22240
#PJM -j
#------- Program execution -------#
module load intel

PRG=/lustre0/home/n22240/vasp/vasp.5.4.4/bin/vasp_std
OUT=stdout_$$

VASP_VDW=/lustre0/home/n22240/vasp/vasp.5.4.4/vdw_kernel.bindat

cp $VASP_VDW ./

# remove old results
rm stdout* stderr*

# INP variables are given by submit_ads.py: INP1 = elem1, INP2 = ads

if [ "${INP3}" != "" ]; then
	# alloy
	python adsorption.py --elem1=${INP1} --elem2=${INP2} --comp1=${INP3} --ads=${INP4} --calculator="vasp" 1> std_${INP1}_${INP2}_${INP3}_${INP4}.out 2> err_${INP1}_${INP2}_${INP3}_${INP4}.out
else
	# pure metals
	python adsorption.py --elem1=${INP1} --ads=${INP2} --calculator="vasp" 1> std_${INP1}_${INP2}.out 2> err_${INP1}_${INP2}.out
fi

#python adsorption.py ${INP1} ${INP2} ${INP3} 1> std_${INP1}_${INP2}_${INP3}_ads.out 2> err_${INP1}_${INP2}_${INP3}_ads.out # alloy
#python adsorption.py ${INP1} ${INP2} ${INP3} 1> std_${INP1}_${INP2}_${INP3}_ads.out 2> err_${INP1}_${INP2}_${INP3}_ads.out # alloy
#python adsorption.py --elem1=${INP1} --ads=${INP2} --calculator="vasp" 1> std_${INP1}_${INP2}.out 2> err_${INP1}_${INP2}.out

