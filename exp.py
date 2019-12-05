import numpy as np
import matplotlib.pyplot as plt



filename= "235U.exp"
filename2 = "239Pu.new.xsec_3d_no_j_var" 
filename3 = "232Th.new.xsec_3d_no_j_var"

x,dx,y,dy = np.loadtxt(filename,usecols=(0,1,2,3),unpack=True)
ener, xsec = np.loadtxt(filename2,usecols=(0,4),unpack=True)
a, b = np.loadtxt(filename3,usecols=(0,4),unpack=True)
for i in range(len(x)):
    y[i]=y[i]*1000 # a mb
    dy[i]=dy[i]*1000
fig, ax = plt.subplots()

ratio=[]

for i in range(len(ener)):
    ratio.append(2*(xsec[i]-b[i])/(xsec[i]+b[i]))
#ax.plot(ener,xsec,'--b',zorder=4)
#ax.plot(a,b,"-r",zorder=3)
#ax.errorbar(x,y,xerr=dx,yerr=dy,fmt='.g',zorder=1)
ax.plot(ener,ratio,'--r')
#----------------------------------
ax.set_xlim([3*10**-2,25])
ax.set_ylim([-0.07,0.07])
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_title(r"R($^{239}$Pu, $^{232}$Th)")
ax.set_xlabel("Energy (MeV)")
ax.set_ylabel(r"R($^{239}$Pu, $^{232}$Th)")
ax.legend(["FRESCO V3.3"])
ax.grid()
#-------------------------------
#fig.savefig("R-239Pu-232Th.pdf",format="pdf")
plt.show()