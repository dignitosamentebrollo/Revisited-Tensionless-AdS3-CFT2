#Solver for TBA and Bethe roots Eqs of tensionless AdS3
#################
from numpy import *
import LoadSave_Lib as lslib
import matplotlib.pyplot as plt
from functools import partial
import datetime
from scipy.optimize import curve_fit
from scipy import integrate
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

def sinh2(x,A=500,B=exp(-250),vec=True): # sinh function that avoids overflows
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

def cosh2(x,A=500,vec=True): # cosh function that avoids overflows
    if vec==False: #the cosh function is from numpy, so I have to pass vectors
        x=[x] 
    N=len(x)
    res=zeros(N,dtype=complex)
    for j in range(N):
        if abs(x[j]) >= A: # Note: cosh(400)=2.61e+173~ +infty << our approximation
            res[j]= cosh(A)
        else:
            res[j]= cosh(x[j])
    if vec==False:
        return res[0]
    return res

def restrict_interval(kset, a, b):
    # Select values in [a, b]
    mask = (kset >= a) & (kset <= b)
    return mask

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
def int_conv(x0,g,conv1,kset,lim=2j/pi):
    def der(gamma,ind):
        return exp(conv1[ind])*(0.5*tanh(-gamma)) #derivate of Y in the pole
    div=False
    for i,x in enumerate(kset):
        if abs(x0-x)<=exp(-250):
            div=True
            ind=i
            #print('Wow')
    integrand=g*sstar(x0-kset)
    if div==True:
        integrand[ind]=der(x0,ind)*lim
    conv = integrate.trapz(integrand,kset)
    return conv

#---------------------
###############################	
# Kernels and S-matrix
def Sstar(x):
    return -1j*tanh(x/2)

def s(x):
    return 1/(2*pi*cosh2(x))

def sstar(x): #sign inverted as commented before
    return -1/(2j*pi*sinh2(x))
#-----------------------------
###############################
# Energies and momenta
def pmirr(x): #divided by h
    return -2/(sinh2(x))

def E(x): #divided by h
    return 2/(cosh2(x))

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



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Thermodynamic-Bethe-Ansatz Equations
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# The program is more stable if the integral convolutions
# of the exact Bethe are computed in the previous gamma
# of the iteration process (as commented in the paper)
from scipy import optimize
def exact_bethe(L,f,n,gamma,conv1,kset):
    def f_pol(x):
        w=-1j*L*p(x) - 4*int_conv(gamma,log2(1-f),conv1,kset) - 1j*pi*(2*n)
        return w.imag #the equation is purely immaginary
    pol_eqn=lambda x: f_pol(x)
    sol=optimize.root_scalar(pol_eqn,bracket=[0,7.5],xtol=1e-14)
    return sol.root

def asymptotic_bethe(L,n):
    def f_pol(x):
        w=-1j*L*p(x) - 1j*pi*(2*n)
        return w.imag #the equation is purely immaginary
    pol_eqn=lambda x: f_pol(x)
    sol=optimize.root_scalar(pol_eqn,bracket=[0,7.5],xtol=1e-14)
    return sol.root

def A0(x,c): #Asymptotic Y0 function
    return (4/(pi**4))*((abs(x)-c)**4)

def A(x,c): #Asymptotic Y function
    return 0.5 - (2/(pi**2))*((abs(x)-c)**2)


def asyconv1(c,Lambda,N):
    # Extended rapidity space for the convolution
    ksetBIG = arange(-4*Lambda,4*Lambda,(2*Lambda)/N)
    s_vals = s(ksetBIG)
    fs_vals = fft(s_vals)
    log_vals = log(1 + A0(ksetBIG,c))
    conv_result_full = fast_conv(fs_vals,log_vals,Lambda,N)
    
    # Extract the physical region [-Lambda, Lambda]
    mask = (ksetBIG >= -Lambda) & (ksetBIG < Lambda)
    conv_result = conv_result_full[mask]
    return conv_result.real

def asyconv2(c,Lambda,N):
    # Extended rapidity space for the convolution
    ksetBIG = arange(-4*Lambda,4*Lambda,(2*Lambda)/N)
    s_vals = s(ksetBIG)
    fs_vals = fft(s_vals)
    log_vals = log(1 - A(ksetBIG,c))
    conv_result_full = fast_conv(fs_vals,log_vals,Lambda,N)
    
    # Extract the physical region [-Lambda, Lambda]
    mask = (ksetBIG >= -Lambda) & (ksetBIG < Lambda)
    conv_result = conv_result_full[mask]
    return 4*conv_result.real


def InitializeAsymptotics(c,kset,Lambda,N):
    a0 = A0(kset,c)
    a = A(kset,c)
    loga0 = log2(a0)
    loga = log2(a)
    convA1 = asyconv1(c,Lambda,N)
    convA2 = asyconv2(c,Lambda,N)
    return a0,a,exp2(convA1-loga),exp2(convA2-loga0)
    # return a0,a,loga0,loga,convA1,convA2

def Updatec(c,f0,kset):
    red = restrict_interval(kset,40,60) #from \gamma=40 the driving term are already (almost) zero
    ksetFIT = kset[red]
    parameters, covariance = curve_fit(A0, ksetFIT, (f0)[red],p0=c)
    c_tmp = parameters[0]
    return c_tmp

    
def tba_solver(L,n,Lambda,N,prec=1e-13,min_prec=1e-12,iterations=int(10e04),damping_param=0.1):
    #------------------------

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
    gamma=asymptotic_bethe(L,n)    #initialize =>> asymptotic bethe roots
    en=4/(cosh2(gamma,vec=False)).real #asymptotic energy

    y0 = exp2(a_1) #starting values of y0
    y = -1*Sstar(-gamma-kset)*Sstar(gamma-kset) #starting values of y
    
    c = 5 #starting value of the unkonwn parameter for the asymptotic functions
    a0, a, Acorr1, Acorr2 = InitializeAsymptotics(c,kset,Lambda,N)
    f0 = y0*a0 #Y0
    f = y*a #Y
    #-------------------
    
    #-- Starting output
    print('Mode \u03BD=', n)
    print('Asymptotic gamma:',gamma)


    #-------------------
    # 		RUN 
    #-------------------
    err=1
    cont=0

    # Collecting data at any iteration
    gamma_plot = zeros(iterations)
    energy_plot = zeros(iterations)
    c_plot = zeros(iterations)
    gamma_plot[0] = gamma 
    energy_plot[0] = en
    c_plot[0] = c
    # -------------------------------
    
    
    #    Iterations 
    #-------------------
    while err > prec and cont< iterations :
        
        # 		Evo 
        #-------------------
        conv1 = fast_conv(fs,log2( (1+f0)/(1+a0) ),Lambda,N) 
        y = Acorr1*exp2(conv1)*Sstar(-gamma-kset)*Sstar(gamma-kset)
        f = y*a
        
        conv2 = 4 * fast_conv(fs,log2( (1-f)/(1-a) ),Lambda,N)
        corr = a_1 + conv2
        y0 = Acorr2*exp2(corr)
        y0[N//2]=0
        f0_tmp = y0*a0
        
        gamma_tmp = exact_bethe(L,f,n,gamma,conv1,kset)
        gamma_plot[cont]=gamma_tmp #saving
        entmp = energy(f0,gamma_tmp,Lambda,N)
        energy_plot[cont]=entmp #saving
        
        # 		Error 
        #-------------------
        err = max((abs(en-entmp),abs(gamma-gamma_tmp)))
        
        if cont%10000==-1: #Disabled
            print('Iteration:',cont)
            print('Error',err)
        
        # 	    Update 
        #-------------------
        c_tmp = Updatec(c,f0_tmp,kset)
        c_plot[cont] = c_tmp
        a0, a, Acorr1, Acorr2 = InitializeAsymptotics(c_tmp,kset,Lambda,N)
        
        y0 = f0_tmp/a0
        y = f/a 
        c = c_tmp
        gamma = gamma_tmp
        en = entmp  
        f0 = (1-damping_param)*f0_tmp + damping_param*f0
        cont += 1 
    
    
    # 	 Final outputs 
    #-------------------
    
    if cont == iterations:
        print('Max iterations reached, error: ', err)
        if err > min_prec :
            print('Failed convergence, error:', err)
            
    
    print('The constant for the asymptotic is:', c)
    print('The final gamma is:', gamma)
    print('The energy in this state is', energy(a0*y0,gamma,Lambda,N))

    print('Iterations:',cont)
    print('Error: ', err)

    return f0,f,gamma,en,array(gamma_plot),array(energy_plot),array(c_plot)


def asy_tba_solver(L,n,Lambda=50,N=int(2**10)):
    gamma = asymptotic_bethe(L,n)    #initialize =>> asymptotic bethe roots
    en = 4/(cosh2(gamma,vec=False)).real

    return gamma, en
    
    
def energy(f0,gamma,Lambda,N):
    kset=arange(-Lambda,Lambda,(2*Lambda)/N)
    dp=[(2*cosh2(j,vec=False)/sinh2(j,vec=False))/ sinh2(j,vec=False) for j in kset] 
    dp[N//2]=0 #in the pole the integrand is zero
    integ=-(1/(2*pi))*integrate.simps(dp*(log2(1+f0)),kset)	#no asymptotic problems due to exp suppression
    return (integ + 4/(cosh2(gamma,vec=False))).real
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
