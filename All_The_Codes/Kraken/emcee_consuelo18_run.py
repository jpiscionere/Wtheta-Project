from __future__ import print_function
import sys
import numpy as np
import emcee
import _chi2_fof
import socket
from datetime import datetime
from subprocess import call
from datetime import datetime
from mpi4py import MPI
from emcee.utils import MPIPool
import pickle
import os



comm=MPI.COMM_WORLD

def lnprob(p):
	lum_sample=18
	box, siglogM,logM0, logM1, alpha, gamma, fgal = p
	if not 0<= box < 10:
		return -np.inf
	if not 0 < alpha < 2:
		return -np.inf
	if not 10<=logM1 < 16:
		return -np.inf
      	if not 0 < gamma <= 3:
		return -np.inf
	if not 0 < fgal < 2:
		return -np.inf
	if not 5 < logM0 < 16:
		return -np.inf
	if not 0 < siglogM < 1:
		return -np.inf 
	value =_chi2_fof.chi2_fof(lum_sample,*p)
	return -0.5*value


resume=1

status_file="status_file_consuelo18.out"
pickle_file="status_file_consuelo18.pkl"
ndim=7
nwalkers=500
nburn_in=5
box=np.array([np.random.randint(0,9) for col in range(nwalkers)])
box=box.reshape(nwalkers,1)
start_positions=np.array([0.19,9.81, 12.6050679001 ,  1.15246263627 ,  1.62330725248 ,  0.220777562472])
pos=[ start_positions * (1 + 0.01 * np.random.randn(len(start_positions)))  for col in range(nwalkers)]
p0=np.append(box,pos,axis=1)


files_copied=0
pool=MPIPool(comm=comm,debug=False)
if not pool.is_master():
	pool.wait()
	sys.exit(0)

sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,pool=pool)

if resume == 0:
        f=open(status_file,"w")
        f.write("## Beginning the calculations....{}\n".format(datetime.time(datetime.now())))
        f.close()

        ### run the burn-in sequence
        pos,prob,rstate=sampler.run_mcmc(p0,nburn_in)

        f=open(status_file,"a")
        f.write("## Done with the burn-in sequence.... {}\n".format(datetime.time(datetime.now())))
        f.write("## Mean acceptance fraction: {0:.3f}\n".format(np.mean(sampler.acceptance_fraction)))


else:
        
	try:
		pkl_file = open(pickle_file, 'rb')
        	pos = pickle.load(pkl_file)
        	prob = pickle.load(pkl_file)
        	rstate = pickle.load(pkl_file)
        	pkl_file.close()

	except:
		pool.close()
		exit
	if len(prob) != nwalkers:
                print("ERROR: Nwalkers from pickle file does not agree with Nwalkers in script")
               	sys.exit(0) 

        f=open(status_file,"a")
        f.write("## Resuming: the calculation....{}\n".format(datetime.time(datetime.now())))
	f.close()

sampler.reset()
iter=0
for pos, prob, rstate in sampler.sample(pos,prob,rstate,iterations=10000,storechain=False):


	pos_matrix=pos.reshape(nwalkers,ndim)
    	prob_array=prob.reshape(nwalkers,1)
    	array=np.append(pos_matrix,prob_array,axis=1)

	execstring="cp -p {0:s} {0:s}.bak".format(pickle_file)
	ret=call(execstring,shell=True)
	pkl_file = open(pickle_file, 'wb')
        pickle.dump(pos, pkl_file,-1)
        pickle.dump(prob, pkl_file,-1)
        pickle.dump(rstate, pkl_file,-1)
        pkl_file.close()



	f = open(status_file,"a")
	f.write("Mean acceptance fraction: {0:.3f}\n".format(np.mean(sampler.acceptance_fraction)))


# Write the current position to a file, one line per walker
	f.write("\n".join(["\t".join([str(q) for q in p]) for p in array]))
	f.write("\n")
	f.write("Done with this iteration.\n")
    	f.close()
	iter = iter + 1


pool.close()

f = open(status_file,"a")
f.write("Final Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
f.write("You've got to see the show, it's a dynamo. You've got to see the show, it's rock and roll {}\n".format(datetime.time(datetime.now())))
f.close()
