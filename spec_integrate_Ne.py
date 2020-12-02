# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 21:35:25 2020

@author: beaub
"""

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import auc

from scipy.interpolate import InterpolatedUnivariateSpline

#change directory path to one where files are located.

#os.chdir(r'..\..\OES\Calibrated\061020')

#wrapping the integral into a function so that any file
#can be used.

def integrate(fileN):
    #filename = input("What is the name of the file you want to analyse?")
    file = fileN+'.tit'
    
    with open(file, 'r') as fileR:
        print("")
        print("--------------------------------------------------------------")
        print("File: ", fileN)
        print("The parameter used in the experiment was: ")
        print("")
        print(fileR.readline())
        
        int_time = next(fileR)
        it = str(int_time)
        newit = it.replace("Integration time: ",'')
        print(newit)
        IT_stripped = newit.replace("ms", "")
        IT = float(IT_stripped)
        print(IT)
                    
    lamda, Irr = np.loadtxt(file, comments='#', 
                        delimiter=';', skiprows=7, unpack=True, usecols=(0,4))
    
    irr = Irr/IT
        
    f = InterpolatedUnivariateSpline(lamda, irr, k=1)  # k=1 gives linear interpolation

    fig = plt.figure()
    est = fig.add_subplot(1,1,1)
    est.plot(lamda,f(lamda), color='tab:orange',lw=0.5,alpha=0.5)
    
    real=fig.add_subplot(1,1,1)
    
    real.plot(lamda,irr, color='tab:blue',lw=0.5,linestyle='dashed')
    est.set_xlim([690,800])
    real.set_xlim([690,800])
    #est.set_ylim([-1000,6000])

    plt.show()

    PeakA_696 = f.integral(694.78, 698.15)
    PeakA_705 = f.integral(705.44,708.24)
    PeakA_714 = f.integral(713.84,716.64)
    PeakA_727 = f.integral(725.58,728.94)
    PeakA_738 = f.integral(736.19,740.65)
    PeakA_750 = f.integral(747.9,754.02)
    PeakA_763 = f.integral(761.25,766.25)
    PeakA_772 = f.integral(770.14,775.13)
    PeakA_794 = f.integral(792.86,797.29)
    PeakA_451 = f.integral(449.5,451.88)
    PeakA_419 = f.integral(418.3,421.1)
    PeakA_367 = f.integral(366.42,368.75)
    PeakA_415 = f.integral(414.77,417.09)
    PeakA_383 = f.integral(382.17,383.92)
    PeakA_365 = f.integral(364.66,366.42)
    PeakA_360 = f.integral(359.99,361.74)
    PeakA_488 = f.integral(487.1,488.83)
    
    I = [PeakA_696,PeakA_705,PeakA_714, PeakA_727,PeakA_738,PeakA_750,PeakA_763,PeakA_772,PeakA_794,PeakA_451,PeakA_419,PeakA_367,PeakA_415,PeakA_383,PeakA_365,PeakA_360,PeakA_488]
    I_data = {'Wavelength (nm)':[696,706,714,727,738,750,763,772,794,451,419,367,415,383,365,360,488],'Integrated Intensity':I}
    #print("")
    #print("The result of using an integrate function is: ")
    #print("")
    df = pd.DataFrame(I_data,index=None,columns=["Wavelength (nm)","Integrated Intensity"])
    with open(fileN+'_IntegratedIntensity.txt', 'w') as resultsfile:
        resultsfile.write('Datafile: '+fileN+'\n')   
    df.to_csv(fileN+'_IntegratedIntensity.txt', index=False,mode="a")
    return fileN

#integral_result = integrate('LOCAL_RF_IRRA0001')
#print(integral_result)










