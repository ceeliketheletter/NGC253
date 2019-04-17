'''This function will perform the least-squares function to find 'pf' for each *finding_pf* txt file, then put the filename, pi, and pf values into an output file'''



import numpy as np
import glob


f = open("shorter_normalized_pi_and_pf_values_data_with_names.txt", "a") #file about to be created

#finding pf
for file in glob.glob('*normalized_pf*.txt'):


    data = np.loadtxt(file)
    data = np.transpose(data)
    loggamma = data[0] 
    t = data[1]
    y = data[2]
    b = data[3]
    inputplleindex = data[4]
    pi = inputplleindex[0]


    numerator = sum((b-y)*t)
    denominator = sum (t*t)
    pf = numerator / denominator
    print file, pi, pf

'''
#save file
                                           
    s = str(file)
    s1 = str(pi)
    s2 = str(pf)
    f.write(s)
    f.write("\t")
    f.write(s1)
    f.write("\t")
    f.write(s2)
    f.write("\n")                         
f.close()  
'''










#previous code that used other least squares function that didnt normalize lines so I had to manually normalize

#finding2_a.jpg
'''

for file in glob.glob('*finding2_a*.txt'):


    data = np.loadtxt(file)
    data = np.transpose(data)
    loggamma = data[0] - .1
    ne_ss = data[1]

#this bit adjusts so each line starts in the same place
    if file == 'finding2_a_0327670.txt':
        ne_ss = ne_ss - 3.79768

    if file == 'finding2_a_5327670.txt':
        ne_ss = ne_ss + 2.66001

    if file == 'finding2_a_4327670.txt':
        ne_ss = ne_ss + 1.37465

    if file == 'finding2_a_3327670.txt':
        ne_ss = ne_ss + 0.0866606

    if file == 'finding2_a_2327670.txt':
        ne_ss = ne_ss - 1.20437

    if file == 'finding2_a_1327670.txt':
        ne_ss = ne_ss - 2.49895

#this is the least squares function
    numerator = sum(np.multiply(ne_ss, loggamma))
 
    denominator = sum(np.multiply(loggamma, loggamma))
    
    a = (numerator/denominator)
   
 
    print file, a


    value = file, a
    s = str(value)
    f.write(s)
    f.write("\n")


f.close()

'''


