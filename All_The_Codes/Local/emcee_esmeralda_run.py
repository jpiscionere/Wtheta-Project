import numpy as np
import emcee
import _chi2_fof
import matplotlib as pl
import psutil
from datetime import datetime

def lnprob(p):
	lum_sample=20
        machine=socket.gethostname()
        if(machine=="bender"):
                bender=1
        else:   
                bender=0
        box, siglogM, logM0, logM1, alpha, gamma, fgal = p
        if not 0<= box < 9:
                return -np.inf
        if not 0<= siglogM <= 2:
                return -np.inf
        if not 0 < alpha < 2:
                return -np.inf
        if not 10<=logM1 < 16:
                return -np.inf
        if not 8 <= logM0 <= logM1:
                return -np.inf


        value =_chi2_fof.chi2_fof(lum_sample,bender,*p)
        return -0.5*value

ndim=7

nwalkers=500

box=np.array([np.random.randint(0,8) for col in range(nwalkers)])

box=box.reshape(nwalkers,1)
        



#start_positions=np.array([0.0972008997931,12.7430056861,13.0177787909,1.04953050812,1.7682103135,0.204769578295])
start_positions=np.array([0.0976826863347, 12.0372625837 ,  13.2711595479 ,  1.17625300352 ,  1.82346906053 ,  0.210540118623])

pos=[ start_positions * (1 + 0.01 * np.random.randn(len(start_positions)))  for col in range(nwalkers)]
#p0= start_positions[None,:] * (1 + np.random.randn(len(start_positions) * nwalkers).reshape((nwalkers,len(start_positions))))

p0=np.append(box,pos,axis=1)

f=open("status_file.out","a")
f.write("Welcome back my friend to the show that never ends....{}\n".format(datetime.time(datetime.now())))
f.close()


sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=32)
pos,prob,state=sampler.run_mcmc(p0,2)
#	f=open("status_file.out","a")
#	f.write("Iteration....{}\n".format(datetime.time(datetime.now())))
#	f.close()
#	print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
#        f = open("status_file.out","a")
#        f.write("sampler.acceptance_fraction {}\n".format(np.mean(sampler.acceptance_fraction)))
#	f.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
#        for i in pos.tolist():
#                f.write("{}\n".format(i))
#        f.write("\n")
#	f.close()



print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))





f=open("status_file.out","a")
f.write("We're fof glad you could attend, come inside come inside {}\n".format(datetime.time(datetime.now())))
f.write("Mean acceptance fraction: {0:.3f}\n".format(np.mean(sampler.acceptance_fraction)))
f.close()

sampler.reset()

for pos, prob, rstate in sampler.sample(pos,prob,state,iterations=1000,storechain=False):
#for pos, prob, rstate in sampler.sample(p0,iterations=1000):
#	print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))	
	psutil.virtual_memory()
	pos_matrix=pos.reshape(nwalkers,ndim)
    	prob_array=prob.reshape(nwalkers,1)
    	array=np.append(pos_matrix,prob_array,axis=1)

	f = open("status_file.out","a")
	f.write("Mean acceptance fraction: {0:.3f}\n".format(np.mean(sampler.acceptance_fraction)))
#	f.write("Autocorrelation time: {}\n".format(sampler.acor))
# Write the current position to a file, one line per walker
	f.write("\n".join(["\t".join([str(q) for q in p]) for p in array]))
	f.write("\n")
	f.write("Done with this iteration.\n")
    	f.close()




f = open("chain_fof.out","a")
f.write("Final Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
f.write("You've got to see the show, it's a dynamo. You've got to see the show, it's rock and roll {}\n".format(datetime.time(datetime.now())))
f.close()


import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
    
pl.figure()
for k in range(nwalkers):
        pl.plot(sampler.chain[k, :, 0])
pl.xlabel("time")
pl.savefig("chain_fof_time.png")

pl.figure(figsize=(8,8))
x, y = sampler.flatchain[:,4], sampler.flatchain[:,5]
pl.plot(x, y, "ok", ms=1, alpha=0.1)
pl.savefig("chain_2d.png")

fig = pl.figure()
ax = fig.add_subplot(111, projection="3d")

for k in range(nwalkers):
        x, y = sampler.chain[k,:,4], sampler.chain[k,:,5]
        z = sampler.lnprobability[k,:]
ax.scatter(x, y, z, marker="o", c="k", alpha=0.5, s=10)
pl.savefig("chain_3d.png")
