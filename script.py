import numpy as np
import matplotlib.pyplot as plt



def graph(name, x, y):
    title = name.replace('.xsec','')
    fig,ax = plt.subplots()
    ax.plot(x, y, '--r')
    ax.set_title('n + ' + title + '\n Total cross section')
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Cross section (mb)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(["FRESCO 3.3"])#,loc='center right')
    ax.grid()
    #fig.savefig(title + ".pdf", format='pdf')
    #ax.set(xlabel='Time (s)', ylabel='Voltage (mV)',
    #    title='Plot example')
    plt.show()
def graph2(name1, x, y, xx, yy):
    title = name1.replace('.xsec','')
    fig,ax = plt.subplots()
    ax.plot(x, y, '--r',xx, yy, '-.b')
    ax.set_title('n + ' + title + '\n Total cross section')
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Cross section (mb)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.legend(["FRESCO 3.3", "FRESCO frx6j"])#,loc='center right')
    ax.grid()
    #ax.set(xlabel='Time (s)', ylabel='Voltage (mV)',
    #    title='Plot example')
    plt.show()
    #fig.savefig(title + ".pdf", format='pdf')
def graph3(name1, x, y, xx, yy, xxx, yyy):
    title = name1.replace('.xsec','')
    fig,ax = plt.subplots()
    ax.plot(x, y, '--r',xx, yy, '-.b', xxx, yyy, '+k')
    ax.set_title('n + ' + title + '\n Total cross section')
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Cross section (mb)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(["FRESCO 3.3 abs=0.001 j=jmax", "FRESCO frxy6j abs=0.001", "FRESCO 3.3 abs=-1.000 j=j_cal"])#,loc='center right')
    ax.grid()
    #ax.set(xlabel='Time (s)', ylabel='Voltage (mV)',
    #    title='Plot example')
    plt.show()
    #fig.savefig(title + ".pdf", format='pdf')
def graph4(name1, x, y, xx, yy, xxx, yyy):
    title = name1.replace('.new.xsec','')
    for i in range(len(xxx)):
        xxx[i]=xxx[i]*1000 # a mb
        yyy[i]=yyy[i]*1000
    
    fig,ax = plt.subplots()
    ax.plot(x, y, '--k',zorder=2)
    ax.errorbar(xx,xxx,xerr=yy,yerr=yyy,fmt='.c',zorder=1)
    ax.set_title('n + ' + title + '\n Total cross section')
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Cross section (mb)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([np.amin(xx),30])
    ax.legend(["FRESCO V3.3","Experimental data"])#,loc='center right')
    ax.grid()
    #ax.set(xlabel='Time (s)', ylabel='Voltage (mV)',
    #    title='Plot example')
    plt.show()
    #fig.savefig(title + ".pdf", format='pdf')
#fig, ax = plt.subplots()
filename = "prueba.txt" #reading input's names

with open(filename) as file:
     inputs = [l.strip() for l in file]
     #outfile = inputs.replace('.inp','.xsec')
#print(inputs)

for input in inputs:
    outfile = input.replace('.inp','.new.xsec_3d_no_j_var')
    outfile2 = input.replace('.inp','.exp')
    #outfile2 = input.replace('.inp','.xsec4')
    #outfile3 = input.replace('.inp','.xsec3')
    energies, total_xsec = np.loadtxt(outfile,usecols=(0,4),unpack=True)
    x,dx,y,dy = np.loadtxt(outfile2,usecols=(0,1,2,3),unpack=True)
    #energies2, total_xsec2 = np.loadtxt(outfile2,usecols=(0,4),unpack=True)
    #energies3, total_xsec3 = np.loadtxt(outfile3,usecols=(0,4),unpack=True)
    #graph(outfile, energies, total_xsec)
    #graph2(outfile, energies, total_xsec, energies2, total_xsec2)
    #graph3(outfile, energies, total_xsec, energies2, total_xsec2, energies3, total_xsec3)
    graph4(outfile,energies,total_xsec,x,dx,y,dy)
    #print(energies, total_xsec)
    #print(outfile)
#np.rea(filename,)
