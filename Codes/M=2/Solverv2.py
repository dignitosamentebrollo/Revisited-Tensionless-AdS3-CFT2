import LoadSave_Libv2 as lslib
import TBAv2 as tba
#---------------------
from numpy import *

    
#--------------------------------------------------------

def full_solution(L,Lambda,N):
    f0,f,gamma1,gamma2,en=lslib.Bethe_Set(N,L) #creates the arrays
    list_energy_plot = []

    #------Solution of the tba--------
    for n1 in range(1,L//2+1): #iteration over the modes of the excitations
        for n2 in [2**i for i in range(int(log2(L)))]:
            f0[n1,n2],f[n1,n2],gamma1[n1,n2],gamma2[n1,n2],en[n1,n2],energy_plot=tba.tba_solver(L,n1,n2,Lambda,N)
            list_energy_plot.append(energy_plot)
    print('--Done--')
    return f0,f,gamma1,gamma2,en,vstack(list_energy_plot)
