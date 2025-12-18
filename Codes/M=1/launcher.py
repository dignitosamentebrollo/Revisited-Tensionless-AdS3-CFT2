import os
import Solver as Solver
import LoadSave_Lib as lslib
import time
import sys
#---------------------
from numpy import *
#======================
#======================
N=int(2**11) # rapidity sampling
Lambda=160	# rapidity cutoff
#======================
Lset=array([16]) #values of string length to compute


for L in Lset:
    #==================
    #== SAVING DIRECTORY
    sim_tag=f'L{L}'
    path=os.getcwd()
    os.chdir('{}/data/'.format(path))
    os.makedirs('{}'.format(sim_tag),exist_ok=True)
    os.chdir('{}'.format(sim_tag))
    # os.chdir('../')
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
    f0,f,gamma,en,list_gamma,list_energy,list_c=Solver.full_solution(L,Lambda,N)
    print('Elapsed time:',time.time()-t0)


    #================
    #== SAVING
    #================
    lslib.master_save(f0,f,gamma,en)
    lslib.single_save(list_gamma,tag='list_gamma_plot')
    lslib.single_save(list_energy,tag='list_energy_plot')
    lslib.single_save(list_c,tag='list_c_plot')
    lslib.single_save(kset,tag='kset')
    # lslib.single_filesave(en,tag='en',folder_name='{}/{}'.format(foldername,sim_tag))
    # lslib.single_filesave(gamma,tag='gamma',folder_name='{}/{}'.format(foldername,sim_tag))

    # ---------------------------
    print('--All data are saved--')

    #================
    #== PLOTTING
    #================
    nset=arange(1,L//2+1,1)
    import matplotlib.pyplot as plt
    plt.axes()
    plt.plot(linspace(0,L//2,2**10),4*sin((pi*linspace(0,L//2,2**10))/L),'c',linestyle='dashed',label=r'Free')
    plt.plot(nset,en[1::].real,'k.',label=r'Exact')

    # x axis
    plt.xlabel('Mode',fontsize=12)
    if L<100:
        a=2
    elif L>=500:
        a=10
    else: a=5
    plt.xticks(range(0,L//2+1,a))

    #y axis
    plt.ylabel('Energy',fontsize=11)
    plt.title(f'L={L}')

    plt.legend()
    plt.tight_layout()
    # os.chdir('{}/data/{}'.format(path,sim_tag))
    plt.savefig(f'EnergyL{L}.pdf')
    plt.close()
    os.chdir('../')
    os.chdir('../')

# --------------------------
print('--end--')
