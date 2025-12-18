import Solverv2 as Solver
import LoadSave_Libv2 as lslib
import time
import os
import sys
#---------------------
from numpy import *

#======================
N=int(2**512) # rapidity sampling
Lambda=160	# rapidity cutoff
#======================
L=4 #worldsheet lenght


#==================
#== SAVING DIRECTORY
sim_tag=f'L{L}'
path=os.getcwd()
os.chdir('{}/data/'.format(path))
os.makedirs('{}'.format(sim_tag),exist_ok=True)
os.chdir('{}'.format(sim_tag))
#==================
sys.stdout = lslib.Logger('output.txt'.format(sim_tag))

print('\nString length L=',L)
print('Lambda: ', Lambda)
print('Discretization N=', N)


#==============================
kset=arange(-Lambda,Lambda,(2*Lambda)/N) #rapidity set from -Lambda to Lambda

#================
#== SOLVING
#================
t0=time.time()

# Solution of the model
f0,f,gamma1,gamma2,en,list_energy=Solver.full_solution(L,Lambda,N)
print('Elapsed time:',time.time()-t0)


#================
#== SAVING
#================
lslib.master_save(f0,f,gamma1,gamma2,en)
lslib.single_save(list_energy,tag='list_energy_plot')
lslib.single_save(kset,tag='kset')
# ---------------------------

print('--All data are saved--')

# --------------------------
print('--end--')
