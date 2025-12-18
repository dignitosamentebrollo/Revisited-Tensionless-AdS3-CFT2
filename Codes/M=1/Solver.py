import LoadSave_Lib as lslib
import TBA as tba
#---------------------
from numpy import *

def full_solution(L,Lambda,N):
	y0,y,gamma,en=lslib.Bethe_Set(N,L) #creates the arrays
	#----------------------------------------
	list_gamma_plot = []
	list_energy_plot = []
	list_c_plot = []
	#---------------------------------

	#------Solution of the tba--------
	for n in range(1,L//2+1): #iteration over the modes of the excitations
		y0[n],y[n],gamma[n],en[n],gamma_plot,energy_plot,c_plot=tba.tba_solver(L,n,Lambda,N)
		list_gamma_plot.append(gamma_plot)
		list_energy_plot.append(energy_plot)
		list_c_plot.append(c_plot)
	print('--Done--')
	return y0,y,gamma,en,vstack(list_gamma_plot),vstack(list_energy_plot),vstack(list_c_plot)

