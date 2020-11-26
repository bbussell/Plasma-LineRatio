# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 20:20:31 2020

@author: beaub
"""

import os
import numpy as np
from spec_integrate_Ne import integrate 

#Program to calculate the electron density from 451/750nm
#ratio
os.chdir(r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\081020')

def e_density(T,file):
    
    nec_750 = 3E12
    nec_451 = 4E11
    nec_419 = 3.2E11 #calculated from reciprocal of transition probability
    nec_367 = 1E11 #estimate based on Zhu 2007
    nec_415 = 5E11
    nec_383 = 1E11 #Zhu 2007
    nec_365 = 1E11
    nec_360 = 1E11
    
    K = 567 - (160*Te)+(16*(Te**2))
    
    #filename = input("What is the name of the file you want to analyse?")
    #file = filename+'_IntegratedIntensity.txt'
    #change directory 
    #os.chdir(r'..\..\OES\Calibrated\081020')
    
    
    lamda, I = np.loadtxt(file+'_IntegratedIntensity.txt', comments='#', delimiter=',', skiprows=2, unpack=True, 
                                            usecols=(0,1))
    
    #os.chdir(r'..\..\..\Plasma Spec Code\Plasma-LineRatio')
    
    I_750 = I[5]
    I_451 = I[9]
    I_419 = I[10]
    I_367 = I[11]
    I_415 = I[12]
    I_383 = I[13]
    I_365 = I[14]
    I_360 = I[15]
    
    # n_e = (1-((I_750/I_365)/(K)))/(((I_750/I_383)/(K*nec_750))-(1/nec_383))
    n_e = (1-((I_750/I_360)/K)) / (((I_750/I_360)/(K*nec_750))-(1/nec_360))
    N_e = abs(n_e)
    print("the electron density for", file, " is: ", "%6.2e" % N_e, "cm-3")
    
    return N_e

#List of datafiles to be used in electron density calculations.
filelist = ['BACK_PP0009',
            'BACK_PP0005',
            'BACK_PP0003']

#integrating over all files in filelist and prints electron density result to file
#labeled with the orignal datafile name.

for i in filelist:
    #Ask user for electon temperature value to be used in electron density calculation
    #calling integrate function to determine intensity of 451 and 750 emission lines
    integrate(i)
    
    print("")
    T_user = input("What is the electron temperature? for " + i)
    print("")

    Te = float(T_user)
    #calling e-density function
    N_e = e_density(Te,i)
    N_e_str = str(N_e)
          
    with open(i+'360_E_density.txt', 'w+') as resultsfile:
        resultsfile.write('Datafile: '+i+'\n')
        resultsfile.write('Electron density (cm^-3) \n')
        resultsfile.write(N_e_str)



