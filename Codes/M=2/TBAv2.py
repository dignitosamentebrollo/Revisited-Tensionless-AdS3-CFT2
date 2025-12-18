# Solver for TBA Eqs of tensionless AdS3
#################
import cmath as c
from scipy import integrate
from scipy import fft
from scipy import optimize
from numpy import *
import matplotlib.pyplot as plt
from functools import partial
# -----------------

#-----------------
def exp2(x,A=500,vec=True): # exp function that avoids overflows
    if vec==False: #the exp function is from numpy, so I have to pass vectors
        x=[x] 
    N=len(x)
    res=zeros(N,dtype=complex)
    for j in range(N):
        if abs(x[j]) >= A: # Note: exp(-500)=7.12e-218~0; exp(500)=1.40e+217~ +infty << our approximation
            res[j]= exp(sign(x[j])*A)
        else:
            res[j]= exp(x[j])
    if vec==False:
        return res[0].real
    return res.real

def log2(x,A=500,vec=True): # log function that avoids overflows
    if vec==False:	#the log function is from numpy, so I have to pass vectors
        x=[x] 
    N=len(x)
    res=zeros(N,dtype=complex)
    for j in range(N):
        if abs(x[j]) >= exp(-A): # Note: exp(-500)=7.12e-218~0; exp(500)=1.40e+217~ +infty << our approximation
            res[j]= log(complex(x[j]))
        elif x[j]==0:
            res[j]= -A
        else:
            res[j]= log(complex(sign(x[j])*exp(-A)))
    if vec==False:
        return res[0]
    return res

def sinh2(x,A=400,B=exp(-250),vec=True): # sinh function that avoids overflows
    if vec==False: #the sinh function is from numpy, so I have to pass vectors
        x=[x] 
    N=len(x)
    res=zeros(N,dtype=complex)
    for j in range(N):
        if abs(x[j]) >= A: # Note: sinh(400)=2.61e+163~ +infty << our approximation
            res[j]= sinh(sign(x[j])*A)
        elif x[j]==0:
            res[j]= B
        elif abs(x[j]) <= B: # Note: sinh(B)~B
            res[j]= sign(x[j])*B
        else:
            res[j]= sinh(x[j])
    if vec==False:
        return res[0]
    return res


#----------------	
def shift(vec): # it assumes N even # Note: twice a shift == original order
    N=len(vec)
    n=int(N//2)
    new_vec=zeros(N,dtype=complex)
    for j in range(n):
        new_vec[j]=vec[j+n]
        new_vec[j+n]=vec[j]
    return new_vec
#-------------------
# The following function compute the conv. defined as \int dk' f(k-k')g(k')
# In the paper the convs. as defined as \int dk' f(k'-k)g(k')
# So in the next part of the code the kernel K(k) are defined as K(-k)
from scipy.fft import fft,ifft
def fast_conv(fK,g,Lambda,N):
    fg=fft(g)
    conv=ifft(fK*fg)
    return shift(real(conv))*(2*Lambda/(N)) 



# The following function compute the conv. for the exact Bethe eqns
# We have to cure the case 1/sinh(0) which is bad defined
# but can be extend continuosly the integrand as explained in the paper
def int_conv(x0,gamma2,g,conv1,kset,lim=2j/pi):
    def der(gamma1,gamma2,ind):
        return exp(conv1[ind])*(0.5*tanh(-gamma1))*(tanh((gamma2-gamma1)/2)*tanh((-gamma2-gamma1)/2)) #derivate of Y in the pole
    div=False
    for i,x in enumerate(kset):
        if abs(x0-x)<=exp(-250):
            div=True
            ind=i
            #print('Wow')
    integrand=g*sstar(x0-kset)
    if div==True:
        integrand[ind]=der(x0,gamma2,ind)*lim
    conv = integrate.trapz(integrand,kset)
    return conv

# ---------------------
###############################
# Kernels and S-matrix
def Sstar(x):
    return -1j*tanh(x/2)

def s(x):
    return 1/(2*pi*cosh(x))

def sstar(x): #sign inverted as commented before
    return -1/(2j*pi*sinh2(x))
#-----------------------------
###############################
# Energies and momenta
def pmirr(x): #divided by h
    return -2/(sinh2(x))

def E(x): #divided by h
    return 2/(cosh(x))

def Emirr(x):
    return -log2(((1-exp2(x))/(1+exp2(x)))**2)

def p(x):
    return -1j*log2(((exp2(x,vec=False)-1j)/(exp2(x,vec=False)+1j))**2,vec=False)

###############################
def error(a,b):
    N=len(a)
    #---------------
    # Abs error 
    abs_err=max([abs(a[j]-b[j]) for j in range(N)])
    #---------------
    # Rel error
    # rel_err=max([abs((a[j]-b[j])/b[j]) for j in range(N)])

    return abs_err
#-------------------------------



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Thermodynamic-Bethe-Ansatz Equations
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def exact_bethe(L,f,n,gamma1,gamma2,conv1,kset):
    def f_pol(x):
        w = -1j*L*p(x) - 4*int_conv(gamma1,gamma2,log2(1-f),conv1,kset) - 1j*pi*2*(n)
        return w.imag #the equation is purely immaginary
    pol_eqn=lambda x: f_pol(x)
    sol=optimize.root_scalar(pol_eqn,bracket=[0,7.5],xtol=1e-14)
    return sol.root


def asymptotic_bethe(L,n):
    def f_pol(x):
        w=-1j*L*p(x) - 1j*pi*2*(n)
        return w.imag #the equation is purely immaginary
    pol_eqn=lambda x: f_pol(x)
    sol=optimize.root_scalar(pol_eqn,bracket=[0,7.5],xtol=1e-14)
    return sol.root


def tba_solver(L, n1, n2, Lambda, N, prec=1e-14, min_prec=1e-12, iterations=int(1e03), damping_param=0.6):
    
    #-- Rapidity space
    kset=arange(-Lambda,Lambda,(2*Lambda)/N)
    N=len(kset)
    
    #-- Integral Kernels dicretized
    s_set=s(kset)
    fs = fft(s_set)
    
    #-------------------
    # 	INITIALIZATION 
    #-------------------
    a_1 = -L*Emirr(kset)
    gamma1 = asymptotic_bethe(L, n1)  # initialize =>> asymptotic bethe roots
    gamma2 = asymptotic_bethe(L, n2)
    en=4/cosh(gamma1)+4/cosh(gamma2) #asymptotic energy

    f0 = zeros(N)  # initialize =>> Y0 asymptotics
    #-------------------


    # -------------------
    # -- Starting output
    print('Mode \u03BD1,\u03BD3=', n1,n2)
    print('Asymptotic Bethe Roots:')
    print(gamma1, gamma2)
    print('Asymptotic Energy:',4/cosh(gamma1)+4/cosh(gamma2))
    
    asy_gamma1=gamma1
    asy_gamma2=gamma2
    asy_en=4/cosh(gamma1)+4/cosh(gamma2)


    # -------------------
    # 		RUN
    # -------------------
    err = 1
    cont = 0
    
    # Collecting energy at any iteration
    energy_plot = zeros(iterations)
    energy_plot[0] = en
    # -------------------------------

    while err > prec and cont < iterations:
        
        # 		Evo 
        #-------------------
        conv1 = fast_conv(fs,log2(1+f0), Lambda, N)
        f = exp2(conv1)*Sstar(-gamma1-kset)*Sstar(gamma1-kset)*Sstar(-gamma2-kset)*Sstar(gamma2-kset)
        
        conv2 = 4*fast_conv(fs,log2(1 - f), Lambda, N)
        corr = a_1+conv2
        f0_tmp = exp(corr)
        
        gamma1_tmp = exact_bethe(L,f,n1,gamma1,gamma2,conv1,kset)
        gamma2_tmp = exact_bethe(L,f,n2,gamma2,gamma1,conv1,kset)
        entmp = energy(f0, gamma1_tmp, gamma2_tmp, Lambda, N)
        energy_plot[cont]=entmp #saving
        
        
        # 		Error 
        #-------------------
        err = max((abs(en-entmp),abs(gamma1-gamma1_tmp),abs(gamma2-gamma2_tmp)))
        
        if cont%100==-1: #Disabled
            print('Iteration:',cont)
            print('Error',err)
        
        # 	    Update 
        #-------------------
        en = entmp
        f0 = (1-damping_param)*f0_tmp+damping_param*f0
        gamma1 = gamma1_tmp
        gamma2 = gamma2_tmp
        cont += 1
    # -----------------------------
    
    
    
    # 	 Final outputs 
    #-------------------
    
    if cont == iterations:
        print('Max iterations reached, error: ', err)
        if err > min_prec :
            print('Failed convergence, error:', err)

    print('The final gamma are:')
    print(gamma1,gamma2)
    print('The energy in this state is', energy(f0, gamma1, gamma2, Lambda, N))

    print('Iterations:',cont)
    print('Error: ', err)

    return f0,f,gamma1,gamma2,en,array(energy_plot)


def energy(f0, gamma1, gamma2, Lambda, N):
    kset = arange(-Lambda, Lambda, (2*Lambda)/N)
    dp = [2*cosh(j)/(sinh2(j, vec=False)**2) for j in kset]
    dp[N//2] = 0  # in the pole the integrand is zero
    # no asymptotic problems due to exp suppression
    integ = -(1/(2*pi))*integrate.simps(dp*(log2(1+f0)), kset)
    return (integ + 4/(cosh(gamma1))+4/(cosh(gamma2))).real
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
