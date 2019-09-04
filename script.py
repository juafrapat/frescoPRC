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
    fig.savefig(title + ".pdf", format='pdf')
    #ax.set(xlabel='Time (s)', ylabel='Voltage (mV)',
    #    title='Plot example')
    #plt.show()
def graph2(name1, x, y, xx, yy):
    title = name1.replace('.xsec','')
    fig,ax = plt.subplots()
    ax.plot(x, y, '--r',xx, yy, '-.b')
    ax.set_title('n + ' + title + '\n Total cross section')
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Cross section (mb)")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(["FRESCO frxy6j", "FRESCO 3.3"])#,loc='center right')
    ax.grid()
    #ax.set(xlabel='Time (s)', ylabel='Voltage (mV)',
    #    title='Plot example')
    #plt.show()
    fig.savefig(title + ".pdf", format='pdf')
#fig, ax = plt.subplots()
filename = "prueba.txt" #reading input's names

with open(filename) as file:
     inputs = [l.strip() for l in file]
     #outfile = inputs.replace('.inp','.xsec')
#print(inputs)

for input in inputs:
    outfile = input.replace('.inp','.xsec')
    outfile2 = input.replace('.inp','.xsec2')
    energies, total_xsec = np.loadtxt(outfile,usecols=(0,4),unpack=True)
    energies2, total_xsec2 = np.loadtxt(outfile2,usecols=(0,4),unpack=True)
    #graph(outfile, energies, total_xsec)
    graph2(outfile, energies, total_xsec, energies2, total_xsec2)


    #print(energies, total_xsec)
    #print(outfile)
#np.rea(filename,)
