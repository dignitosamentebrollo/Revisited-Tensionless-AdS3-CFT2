from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset



N=int(2**11) # rapidity sampling
Lambda = 160 # rapidity cutoff

##================
#== ENERGY
##================
# L=4 #string length
# nset4=arange(0,L//2+1,1)
# sim_tag=f'L{L}'
# with open('data/{}/energy.npy'.format(sim_tag), 'rb') as w:
# 	en4=load(w)
# with open('data/{}/asy_en.npy'.format(sim_tag), 'rb') as w:
# 	asy_en4=load(w)

L=8 #string length
nset8=arange(0,L//2+1,1)
sim_tag=f'L{L}'
with open('M=1/data/{}/energy.npy'.format(sim_tag), 'rb') as w:
	en8=load(w)

L=16 #string length
nset16=arange(0,L//2+1,1)
sim_tag=f'L{L}'
with open('M=1/data/{}/energy.npy'.format(sim_tag), 'rb') as w:
	en16=load(w)

L=32 #string length
nset32=arange(0,L//2+1,1)
sim_tag=f'L{L}'
with open('M=1/data/{}/energy.npy'.format(sim_tag), 'rb') as w:
	en32=load(w)

L=64 #string length
nset64=arange(0,L//2+1,1)
sim_tag=f'L{L}'
with open('M=1/data/{}/energy.npy'.format(sim_tag), 'rb') as w:
	en64=load(w)

L=128 #string length
nset128=arange(0,L//2+1,1)
sim_tag=f'L{L}'
with open('M=1/data/{}/energy.npy'.format(sim_tag), 'rb') as w:
	en128=load(w)



# setting plot widht and height
f=plt.figure()
f.set_figwidth(17.8)
f.set_figheight(4)


L=8
plt.subplot(1,4,1)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.87)
plt.plot(linspace(0,L//2,2**8),4*sin((pi*linspace(0,L//2,2**8))/L),'k',linestyle='dashed',linewidth=2,label=r'Free')
plt.plot(nset8,en8,'r.',label=r'Exact',markersize=14)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,1))
plt.ylabel(r'$E_{(1)}$',fontsize=24,labelpad=-50)
plt.yticks(range(0,5))
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=20)
plt.tight_layout()

L=16
plt.subplot(1,4,2)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.87)
plt.plot(linspace(0,L//2,2**8),4*sin((pi*(linspace(0,L//2,2**8)))/L),'k',linestyle='dashed',linewidth=2,label=r'Free')
plt.plot(nset16,en16,'r.',label=r'Exact',markersize=14)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,2))
ax = plt.gca()
ax.get_yaxis().set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.tight_layout()

L=32
plt.subplot(1,4,3)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.87)
plt.plot(linspace(0,L//2,2**8),4*sin((pi*(linspace(0,L//2,2**8)))/L),'k',linestyle='dashed',linewidth=2,label=r'Free')
plt.plot(nset32,en32,'r.',label=r'Exact',markersize=14)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,L//8))
ax = plt.gca()
ax.get_yaxis().set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.tight_layout()

L=128
plt.subplot(1,4,4)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.87)
plt.plot(linspace(0,L//2,2**8),4*sin((pi*linspace(0,L//2,2**8))/L),'k',linestyle='dashed',label=r'Free')
plt.plot(concatenate((ones(1),nset128[::4])),concatenate((en128[1]*ones(1),en128[::4])),'r.',label=r'Exact',markersize=14)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,L//8))
ax = plt.gca()
ax.get_yaxis().set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.legend(prop={'size': 22},loc=(0.5,0.12))
plt.tight_layout()


plt.savefig('plt1.pdf')
plt.show()

##================
#== DEV
##================
f=plt.figure()
f.set_figwidth(17.8/2)
f.set_figheight(4)

# plt.title(r'Finite-size correction',fontsize=22,x=0.5,y=0.99)


plt.plot((arange(0,32//2+1,1)/32)[1::],(en32-4*sin(pi*(arange(0,32//2+1,1)/32)))[1::],'r.',label=r'L=32',markersize=16)
plt.plot((arange(0,64//2+1,1)/64)[2::2],(en64-4*sin(pi*(arange(0,64//2+1,1)/64)))[2::2],'r+',label=r'L=64',markersize=16)
plt.plot((arange(0,128//2+1,1)/128)[::4][1::],(en128-4*sin(pi*(arange(0,128//2+1,1)/128)))[2::4],'b.',label=r'L=128',markersize=16)
plt.plot((arange(0,256//2+1,1)/256)[1::],zeros(256//2),'k',linestyle='dashed',lw=1)

ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel(r'$\nu_1/L$',fontsize=24,labelpad=-50)
plt.xticks(arange(0,0.6,0.1))
plt.ylabel(r'$E_{(1)}-E_{(1)}^\text{free}$',fontsize=24,labelpad=-105)


plt.legend(prop={'size': 20})
plt.tight_layout()
plt.savefig('plt2.pdf')
plt.show()


##================
#== M=2
##================
L=16
f=plt.figure()
f.set_figwidth(17.8/2)
f.set_figheight(4)

#plt.subplot(2,1,2)
plt.title(r'$L=16$',fontsize=24,x=0.5,y=0.87)

n1=1
en=[1.624344205, 2.391045315, 3.095815577, 3.712694022, 4.218544252, 4.594243942, 4.825535141, 4.903622828]
# asy_en=[0.784137123, 1.553207270, 2.277655509, 2.929641698, 3.484110375, 3.919753619, 4.219829904, 4.372807468]
plt.plot(linspace(0,L//2,2**8),4*sin(pi*n1/L)+4*sin((pi*linspace(0,L//2,2**8))/L),color='tab:gray',linestyle='dashed',linewidth=2)
# plt.plot(nset16[1::],asy_en,'kx',markersize=12)
plt.plot(nset16[1::],en,'k.',label=r'$\nu_1=1$',markersize=16)

'''n1=2
en=[2.3712135226950624, 3.1068257934781567, 3.7861102552033246, \
4.38642718090383, 4.887883790865961, 5.275011580638568, \
5.537485375349091, 5.669769716797533]
asy_en=[2.2941048379003437, 3.0258098880427124, 3.7038162386502496, \
4.3041875835119825, 4.806238992523134, 5.194048969723625, \
5.457053965955099, 5.589617450258511]
plt.plot(linspace(0,L//2,2**8),4*sin(pi*n1/L)+4*sin((pi*linspace(0,L//2,2**8))/L),color='tab:gray',linestyle='dashed',linewidth=2)
plt.plot(nset16[1::],asy_en,'bx',markersize=12)
plt.plot(nset16[1::],en,'b.',label=r'$\nu_1=2$',markersize=14)'''

n1=4
en=[3.712694022, 4.446274543, 5.123922147, 5.718456076, 6.206607959, 6.569428836, 6.792882659, 6.868338195]
# asy_en=[2.929641698, 3.698711846, 4.423160084, 5.075146273, 5.629614950, 6.065258194, 6.365334480, 6.518312043]
plt.plot(linspace(0,L//2,2**8),4*sin(pi*n1/L)+4*sin((pi*linspace(0,L//2,2**8))/L),color='tab:gray',linestyle='dashed',linewidth=2)
# plt.plot(nset16[1::],asy_en,'gx',markersize=12)
plt.plot(nset16[1::],en,'g.',label=r'$\nu_1=4$',markersize=16)

n1=8
en=[4.903622828, 5.620003409, 6.284331818, 6.868338195, 7.348389597, 7.705428556, 7.925404976, 8.000000000]
# asy_en=[4.372807468, 5.141877616, 5.866325854, 6.518312043, 7.072780720, 7.508423964, 7.808500250, 7.961477813]
plt.plot(linspace(0,L//2,2**8),4*sin(pi*n1/L)+4*sin((pi*linspace(0,L//2,2**8))/L),color='tab:gray',linestyle='dashed',linewidth=2)
# plt.plot(nset16[1::],asy_en,'rx',markersize=12)
plt.plot(nset16[1::],en,'r.',label=r'$\nu_1=8$',markersize=16)

plt.xlabel(r'$\nu_3$',fontsize=24,labelpad=-50)
plt.xticks(range(0,16//2+1))
plt.ylabel(r'$E_{(1)}$',fontsize=24,labelpad=-48)
plt.yticks(range(0,10,2))
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=20)


plt.legend(prop={'size': 20})
plt.tight_layout()

plt.savefig('plt3.pdf')
plt.show()


##================
#== OLD vs NEW
##================

f=plt.figure()
f.set_figwidth(17.8)
f.set_figheight(4)

olden16 = [0,0.465420935, 1.260651724, 1.994543892, 2.642471977, 3.183015730, 3.599708149, 3.881788814, 4.023771543]
olden8 = [0,0.963587386, 2.436226614, 3.480392136, 4.008807950]

L=8
plt.subplot(1,2,1)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.875)
plt.plot(linspace(0,L//2,2**8),4*sin((pi*(linspace(0,L//2,2**8)))/L),'k',linestyle='dashed',linewidth=2,label=r'$E^\text{free}$')
plt.plot(nset8,en8,'r.',label=r'$E$',markersize=14)
# plt.plot(nset16[1::],asy_en16[1::],'bx',label=r'Bethe-Yang',markersize=12)
plt.plot(nset8,olden8,'bx',label=r'$\mathcal{E}$',markersize=16,markeredgewidth=1.6)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,2))
plt.ylabel(r'$E_{(1)}$',fontsize=24,labelpad=-50)
plt.yticks(range(0,5))
ax = plt.gca()
# ax.get_yaxis().set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.legend(prop={'size': 20},loc=(0.775,0.05))

L=16
plt.subplot(1,2,2)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.875)
plt.plot(linspace(0,L//2,2**8),4*sin((pi*(linspace(0,L//2,2**8)))/L),'k',linestyle='dashed',linewidth=2,label=r'$E^\text{free}$')
plt.plot(nset16,en16,'r.',label=r'$E$',markersize=14)
# plt.plot(nset16[1::],asy_en16[1::],'bx',label=r'Bethe-Yang',markersize=12)
plt.plot(nset16,olden16,'bx',label=r'$\mathcal{E}$',markersize=16,markeredgewidth=1.6)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,2))
ax = plt.gca()
ax.get_yaxis().set_visible(False) 
ax.tick_params(axis='both', which='major', labelsize=20)



plt.tight_layout()
plt.savefig('plt4.pdf')
plt.show()


##================
#== OLD vs NEW v2
##================

f=plt.figure()
f.set_figwidth(17.8)
f.set_figheight(4)

olden16 = [0,0.465420935, 1.260651724, 1.994543892, 2.642471977, 3.183015730, 3.599708149, 3.881788814, 4.023771543]
olden8 = [0,0.963587386, 2.436226614, 3.480392136, 4.008807950]

L=8
plt.subplot(1,2,1)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.875)
# plt.plot(linspace(0,L//2,2**8),4*sin((pi*(linspace(0,L//2,2**8)))/L),'k',linestyle='dashed',linewidth=2,label=r'Free')
plt.plot(nset8[1::],en8[1::] - 4*sin((pi*nset8[1::])/L),'r.',label=r'$E - E^\text{free}$',markersize=16)
# plt.plot(nset16[1::],asy_en16[1::],'bx',label=r'Bethe-Yang',markersize=12)
plt.plot(nset8[1::],olden8[1::]- 4*sin((pi*nset8[1::])/L),'bx',label=r'$\mathcal{E} - E^\text{free}$',markersize=16,markeredgewidth=1.6)
plt.plot((arange(-1,L//2+1,1))[1::],zeros(L//2+1),'k',linestyle='dashed',lw=0.8)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,1))
plt.ylabel(r'$E_{(1)} - E^\text{free}$',fontsize=24,labelpad=-90)
plt.ylim(-0.65,0.1)
# plt.yticks(range(0,5))
ax = plt.gca()
# ax.get_yaxis().set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.legend(prop={'size': 20},loc=(0.7,0.05))

L=16
plt.subplot(1,2,2)
plt.title(fr'$L=${L}',fontsize=24,x=0.5,y=0.875)
# plt.plot(linspace(0,L//2,2**8),4*sin((pi*(linspace(0,L//2,2**8)))/L),'k',linestyle='dashed',linewidth=2,label=r'Free')
plt.plot(nset16[1::],en16[1::]- 4*sin((pi*nset16[1::])/L),'r.',label=r'$E- E^\text{free}$',markersize=16)
# plt.plot(nset16[1::],asy_en16[1::],'bx',label=r'Bethe-Yang',markersize=12)
plt.plot(nset16[1::],olden16[1::]- 4*sin((pi*nset16[1::])/L),'bx',label=r'$\mathcal{E}- E^\text{free}$',markersize=16,markeredgewidth=1.6)
plt.plot((arange(-1,L//2+1,1))[1::],zeros(L//2+1),'k',linestyle='dashed',lw=0.8)
plt.xlabel(r'$\nu_1$',fontsize=24,labelpad=-50)
plt.xticks(range(0,L//2+1,2))
plt.ylim(-0.35,0.05)
ax = plt.gca()
# ax.get_yaxis().set_visible(False) 
ax.tick_params(axis='both', which='major', labelsize=20)



plt.tight_layout()
plt.savefig('plt5.pdf')
plt.show()

##================
#== Diverging Y functions
##================

L=16 #string length
nset=arange(0,L//2+1,1)
sim_tag=f'L{L}'
with open('M=1/data/{}/f0.npy'.format(sim_tag), 'rb') as w:
	f0=load(w)
with open('M=1/data/{}/f.npy'.format(sim_tag), 'rb') as w:
	f=load(w)
with open('M=1/data/{}/kset.npy'.format(sim_tag), 'rb') as w:
	kset=load(w)
with open('M=1/data/{}/gamma.npy'.format(sim_tag), 'rb') as w:
	gamma=load(w)


n=1
s=0.45
bd=int(N/2)-int(N/2*(10/40))

print(kset[bd:N-bd+1])

# setting plot widht and height
fi=plt.figure()
fi.set_figwidth(17.8) 
fi.set_figheight(6) 

plt.subplot(1, 2, 1)
#plt.title(fr'$Y_0$\quad$L=${L}\quad$n=${n}',fontsize=50,x=0.5,y=0.92)
plt.plot(kset[bd:N-bd+1],f0[n].real[bd:N-bd+1],'k-',linewidth=5*s)
plt.plot(gamma[n]*ones(1),zeros(1),'rx',markersize=30*s,markeredgewidth=2.5)
plt.plot(-gamma[n]*ones(1),zeros(1),'rx',markersize=30*s,markeredgewidth=2.5)
plt.plot(kset[bd:N-bd+1],zeros(N-2*bd+1),'k-',linestyle='dashed',linewidth=4*s,alpha=0.5)
# plt.xlabel(r'$\gamma$',fontsize=54*s)
#plt.xticks(range(0,L//2+1,1))
plt.ylabel(r'$Y_0(\gamma)$',fontsize=64*s,labelpad=-55)
#plt.yticks(range(0,5))
ax = plt.gca()
ax.tick_params(axis='both', which='major', length=12*s, labelsize=50*s)
ax.ticklabel_format(axis='y',style='sci',scilimits=(0, 0),useMathText=True)
offset = ax.yaxis.get_offset_text()
offset.set_size(40*s)
offset.set_y(1.02)
#plt.legend(prop={'size': 46},loc='lower right')

axins1 = inset_axes(ax,width="50%",height="45%",loc='upper center',borderpad=2)
axins1.plot(kset[bd:N-bd+1],f0[n].real[bd:N-bd+1],'k-',linewidth=4*s)
axins1.plot(gamma[n]*ones(1),zeros(1),'rx',markersize=20*s,markeredgewidth=1.5)
axins1.plot(-gamma[n]*ones(1),zeros(1),'rx',markersize=20*s,markeredgewidth=1.5)

axins1.set_xlim(-4.2, 4.2)
axins1.set_ylim(-0.2, 1.5)

axins1.axhline(0, color='k', linestyle='dashed', linewidth=3*s, alpha=0.5)
axins1.tick_params(axis='both',which='major',labelsize=40*s,length=6*s)


plt.subplot(1, 2, 2)
#plt.title(fr'$Y_0$\quad$L=${L}\quad$n=${n}',fontsize=50,x=0.5,y=0.92)
plt.plot(kset[bd:N-bd+1],f[n].real[bd:N-bd+1],'k-',linewidth=5*s)
plt.plot(gamma[n]*ones(1),zeros(1),'rx',markersize=30*s,markeredgewidth=2.5)
plt.plot(-gamma[n]*ones(1),zeros(1),'rx',markersize=30*s,markeredgewidth=2.5)
plt.plot(kset[bd:N-bd+1],zeros(N-2*bd+1),'k-',linestyle='dashed',linewidth=4*s,alpha=0.5)
# plt.xlabel(r'$\gamma$',fontsize=54*s)
#plt.xticks(range(0,L//2+1,1))
plt.ylabel(r'$Y(\gamma)$',fontsize=64*s,labelpad=-105)
plt.yticks([-250,-150,-50,0])
ax = plt.gca()
ax.tick_params(axis='both', which='major', length=12*s, labelsize=50*s)
#plt.legend(prop={'size': 46},loc='lower right')

axins2 = inset_axes(ax,width="50%",height="45%",loc='lower center',borderpad=2)
axins2.plot(kset[bd:N-bd+1],f[n].real[bd:N-bd+1],'k-',linewidth=4*s)
axins2.plot(gamma[n]*ones(1),zeros(1),'rx',markersize=20*s,markeredgewidth=1.5)
axins2.plot(-gamma[n]*ones(1),zeros(1),'rx',markersize=20*s,markeredgewidth=1.5)

axins2.set_xlim(-4.2, 4.2)
axins2.set_ylim(-1.2, 0.8)

axins2.axhline(0, color='k', linestyle='dashed', linewidth=3*s, alpha=0.5)
axins2.tick_params(axis='both',which='major',labelsize=40*s,length=6*s)

plt.tight_layout()
plt.savefig(f'YfuncL{L}.pdf')

plt.show()
