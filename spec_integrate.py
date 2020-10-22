# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 21:35:25 2020

@author: beaub
"""

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sympy as sy
from scipy.integrate import quad

from sklearn.metrics import auc

from scipy.interpolate import InterpolatedUnivariateSpline

#change directory path to one where files are located.
os.chdir(r'..\..\OES\Calibrated\071020')

#wrapping the integral into a function so that any file
#can be used.

def integrate():
    filename = input("What is the name of the file you want to analyse?")
    file = filename+'.txt'
    lamda, irr = np.loadtxt(file, comments='#', 
                        delimiter=' ', skiprows=6, unpack=True, usecols=(0,1))
    
    f = InterpolatedUnivariateSpline(lamda, irr, k=1)  # k=1 gives linear interpolation

    fig = plt.figure()
    est = fig.add_subplot(1,1,1)
    est.plot(lamda,f(lamda), color='tab:orange',lw=0.5,alpha=0.5)
    
    real=fig.add_subplot(1,1,1)
    
    real.plot(lamda,irr, color='tab:blue',lw=0.5,linestyle='dashed')
    est.set_xlim([690,800])
    real.set_xlim([690,800])
    est.set_ylim([-1000,6000])

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
    
    I = [PeakA_696,PeakA_705,PeakA_727,PeakA_738,PeakA_750,PeakA_763,PeakA_772,PeakA_794]
    I_data = {'Wavelength (nm)':[696,706,727,738,750,763,772,794],'Integrated Intensity':I}
    print("")
    print("The result of using an integrate function is: ")
    print("")
    df = pd.DataFrame(I_data,index=None,columns=["Wavelength (nm)","Integrated Intensity"])
    with open(filename+'_IntegratedIntensity.txt', 'w') as resultsfile:
        resultsfile.write('Datafile: '+filename+'\n')   
    df.to_csv(filename+'_IntegratedIntensity.txt', index=False,mode="a")
    return df

df = integrate()
print(df)










