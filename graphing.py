'''
This program loops over all the .txt files that have * "__" * in the filename, turns the first column into a loggamma array, second into ne_ss array, then plots them all onto the same graph.
'''


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import glob
import itertools
import pylab


fig = plt.figure()
ax = plt.gca()
#marker = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))




#for secondary code
#muons

'''
for file in glob.glob('*qe_nu_mu*.dat'):
    data = np.loadtxt(file)
    data = np.transpose(data)
    eps = data[0]
    logeps = np.log10(eps)
    muons = data[1]
    logmu = np.log10(muons)

    
    color = next(ax._get_lines.color_cycle) 
    plt.plot(logeps, logmu, 'o',  color=color)
    print file
    
fig.suptitle('Plot of Muons. \n Varying Pi from 1.9 to 2.05, and 2.55 to 2.625, stepsize=0.025, 11 files, | B |=4.e-4, \n Ecut=1e+15.')
plt.xlabel('log10(eps)   [ev]')
plt.ylabel('log10(qe_nu_mu)   [ergs/cm^2/s]')
plt.show()


fig.savefig('qe_nu_mu.jpg')  #saving
'''






#neutrinos
'''
for file in glob.glob('*qe_nu_e*.dat'):
    data = np.loadtxt(file)
    data = np.transpose(data)
    eps = data[0]
    logeps = np.log10(eps)
    neutrinos = data[1]
    logelec = np.log10(neutrinos)

    
    color = next(ax._get_lines.color_cycle) 
    plt.plot(logeps, logelec, 'o',  color=color)
    print file
    
fig.suptitle('Plot of Neutrinos. \n Varying Pi from 1.9 to 2.05,and 2.55 to 2.625, stepsize=0.025, 11 files, | B |=4.e-4, \n Ecut=1e+15.')
plt.xlabel('log10(eps)   [ev]')
plt.ylabel('log10(qe_nu)   [ergs/cm^2/s]')
plt.show()


fig.savefig('qe_nu.jpg')  #saving
'''








#electrons and positrons
'''
for file in glob.glob('*qe_e_minus*.dat'):
    data = np.loadtxt(file)
    data = np.transpose(data)
    eps = data[0]
    logeps = np.log10(eps)
    electrons = data[1]
    logelec = np.log10(electrons)
    
    plt.plot(logeps, logelec, 'o',  color='b')
    
for file in glob.glob('*qe_e_plus*.dat'):
    data = np.loadtxt(file)
    data = np.transpose(data)
    eps = data[0]
    logeps = np.log10(eps)
    positrons = data[1]
    logelec = np.log10(positrons)

    
    plt.plot(logeps, logelec, 'o',  color='r')
    
fig.suptitle('Plot of Electrons and Positrons. \n Varying Pi from 1.9 to 2.05, and 2.55 to 2.625, stepsize=0.025, 18 files, | B |=4.e-4, \n Ecut=1e+15.')
plt.xlabel('log10(eps)   [ev]')
plt.ylabel('log10(qe_e)   [ergs/cm^2/s]')
plt.show()


fig.savefig('qe_e.jpg')  #saving
'''






#gammas
'''
for file in glob.glob('*pp_gam*.dat'):
    data = np.loadtxt(file)
    data = np.transpose(data)
    eps = data[0]
    logeps = np.log10(eps)
    gammas = data[1]
    loggammas = np.log10(gammas)

    
    color = next(ax._get_lines.color_cycle) 
    plt.plot(logeps, loggammas, 'o',  color=color)
    
fig.suptitle('Plot of Gamma Values. \n Varying Pi from 1.9 to 2.05, and 2.55 to 2.625, stepsize=0.025, 11 files, | B |=4.e-4, \n Ecut=1e+15.')
plt.xlabel('log10(eps)   [ev]')
plt.ylabel('log10(pp_gam)   [ergs/cm^2/s]')
plt.show()


fig.savefig('pp_gam.jpg')  #saving
'''


#for normalized proton code
'''
for file in glob.glob('*normalized_pf*.txt'):

    
    data = np.loadtxt(file)
    data = np.transpose(data)
    loggamma = data[1]  #loggamma - .1
    y = data[2]
    ne_ss = y
    ne_ss = ne_ss - ne_ss[0] #remove the y-intercept by looking at the first ne_ss value   

    color = next(ax._get_lines.color_cycle) 
    label1 = file
 
    plt.plot(loggamma, ne_ss, 'o',  color=color) 


#overlay lines of best fit obtained from pi_and_pf_values.txt

x = np.arange(0,4)
data = np.loadtxt('shorter_normalized_pi_and_pf_values_data.txt')
data = np.transpose(data)
pf = data[1]

for i in pf:
   slope = 1 - i #conversion from a to pf
   line = slope*x
   plt.plot(line)


#plotting

fig.suptitle('Varying Pi from 1.9 to 2.05 and 2.55 to 2.625, stepsize=0.025, 11 files, | B |=4.e-4, \n Ecut=1e+15.')
plt.xlabel('Gamma')
plt.ylabel('Proton Distribution')

#handles, labels = ax.get_legend_handles_labels()  #legend
#plt.legend(handles, labels, loc=3)

plt.show()


#fig.savefig('shorter_normalized_pf_best_fit.jpg')  #saving
'''
