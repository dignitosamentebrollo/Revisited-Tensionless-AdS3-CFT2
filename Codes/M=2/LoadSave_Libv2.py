from numpy import *
import sys

class Logger:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()
    
    def close(self):
        self.log.close()
        
#============================
def Bethe_Set(N,L):
	f0=array([array([zeros(N)]*(L//2+1))]*(L//2+1),dtype=complex)
	f=array([array([zeros(N)]*(L//2+1))]*(L//2+1),dtype=complex)
	gamma1=array([zeros(L//2+1)]*(L//2+1),dtype=complex)
	gamma2=array([zeros(L//2+1)]*(L//2+1),dtype=complex)
	en=array([zeros(L//2+1)]*(L//2+1),dtype=complex)

	return  f0,f,gamma1,gamma2,en
#==============================================================
# Note: NEED TO CREATE BEFOREHAND A FOLDER CALLED 'data' in the current Directory
#==============================================================
def single_save(array, tag='test'):
	with open('{}.npy'.format(tag), 'wb') as f:
		save(f,array)

def single_filesave(array, tag='test'):
	with open('{}.dat'.format(tag), 'wb') as f:
		savetxt(f,array)

def single_load(tag='test'):
	with open('{}.npy'.format(tag), 'rb') as f:
		array=load(f)
	return array
#===============================================================
def master_save(f0,f,gamma1,gamma2,en):
	single_save(f0,tag='f0')
	single_save(f,tag='f')
	single_save(gamma1,tag='gamma1')
	single_save(gamma2,tag='gamma2')
	single_save(en,tag='energy')

#-------------------------------

