''' This is a function that allows you to input an arbitrary x value, and it will interpolate the corresponding y value based on an interpolation between the model points. It graphs your value alongside the model values, and outputs the chisquare'''


import math
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.stats import chisquare

fig = plt.figure()

chisquares = []

#open model manaully
'''
pi = raw_input('Pi value: ')
if len(pi) == 1:
    zeros = '.000'
if len(pi) == 2:
    zeros = '000'
if len(pi) == 3:
    zeros = '00'
if len(pi) == 4:
    zeros = '0'
if len(pi) == 5:
    zeros = ''

string = 'pp_gam_spec_ss_finite_B_400.0_plindex_'+ pi + zeros + '_Ecut_15.0.dat'
'''
for modelfile in glob.glob('*pp_gam*.dat'):
    openfile = open(modelfile, 'r')
    lines = openfile.readlines()
    

    x_model = []
    y_model = []
    

    for line in lines:              #read in lines
        p = line.split()
        x_model.append(float(p[0]))
        y_model.append(float(p[1]))

    log_x = np.array(np.log10(x_model)) #make lists arrays, also take log of both
    log_y = np.array(np.log10(y_model))

    plt.plot(x_model, y_model, 'bo')  #pl0t model

#import data file

    for datafile in glob.glob('*_1data*.txt'):
        data = np.loadtxt(datafile)
        data = np.transpose(data)
   
        x_values = data[0] * 1000000
        y_values = data[3]*1000000 * 1.602e-12  #convert from Mev to ev, then from ev to ergs
        error = data[4]*1000000 * 1.602e-12
    
        plt.plot (x_values, y_values, 'r8')   #plot data

        xnew = np.linspace(6, 15, num=100, endpoint=True) #num determines how many times the linear line bend
        chi_list = []
  
        number_of_points = len(data[0])
        for i in range(number_of_points):

            f = interp1d(np.log10(x_model), np.array(y_model))              #interpolate spaces between model points
    
            print 'log_x=', np.log10(x_values[i]), 'y=', y_values[i], 'error=', error[i]
    
            if datafile == 'fermi_1data.txt' and i <= 3:
                numerator = y_values[i] - f(np.log10(x_values[i]))
                denominator = error[i]

                math = (numerator/denominator)**2 
                chi_list.append(math)
                
            if datafile == 'hess_1data.txt' and i<=4:
                numerator = y_values[i] - f(np.log10(x_values[i]))
                denominator = error[i]

                math = (numerator/denominator)**2 
                chi_list.append(math)
                
  
          
        chisq = sum(chi_list)
        
        if datafile == 'hess_1data.txt':
            print "chisq: ", chisq
             
            chisquares.append(chisq) 
            
        plt.plot (10**xnew, f(xnew), 'g-') # plot interpolation
        yerr=error    
        plt.errorbar(x_values, y_values, yerr=yerr, fmt = ' ', color="green")
           
#plotting data, model, and interpolation line
  
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Plot of Gammas and Fermi+Hess Data. chisq = ' + str(chisq) + '\n model =' + modelfile) #+ ', | B |=4.e-4, Ecut=1e+15.')
    plt.xlabel('log10(eps)   [ev]')
    plt.ylabel('log10(pp_gam)   [ergs/cm^2/s]')
    plt.show()
  
    

print 'Compare all chisquares: ', chisquares





#manual inputs
'''
for i in range(int(raw_input('How many data points? '))):
    interp = raw_input('interpoloation type? l, c: ')
    if interp == 'l':
        f = interp1d(x, np.array(y_model))              #interpolate spaces between model points
    if interp == 'c':
        f = interp1d(x, np.array(y_model), kind='cubic')
    if interp == 's':
        #tck = interpolate.splrep(x, y, s=0)
        #f = interpolate.splev(xnew, tck, der=0)
    
    x_data = float(raw_input('Enter a value for x_data: '))
    xholder.append(x_data)
    y_data = float(raw_input('Enter a value for y_data: '))
    yholder.append(y_data)

    numerator = y_data - f(x_data)
    denominator = float(raw_input('vertical error bar size: '))

    math = (numerator/denominator)**2
    stuff.append(math)

    
    plt.plot (x, y, 'bo', xnew, np.log10(f(xnew)), 'g-')
    plt.plot (np.array(xholder), np.log10(np.array(yholder)), 'r8')
    
    
    fig.suptitle('Plot of Gammas and Data Points. \n  Pi =' + pi + ', | B |=4.e-4, Ecut=1e+15.')
    plt.xlabel('log10(eps)   [ev]')
    plt.ylabel('log10(pp_gam)   [ergs/cm^2/s]')
    plt.show()

chisq = sum(stuff)
print "chisq: ", chisq
'''