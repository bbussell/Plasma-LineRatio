# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 21:35:25 2020

@author: beaub
"""

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon

from sklearn.metrics import auc

from scipy.interpolate import InterpolatedUnivariateSpline as IUS

#change directory path to one where files are located.

#os.chdir(r'..\..\OES\Calibrated\061020')

#wrapping the integral into a function so that any file
#can be used.

def ReadAndNormaliseIntegrationTime(fileN):
    file = fileN+'.tit.'

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
        
    return irr, lamda

def InterpolateNormalisedSpectrum(irr,lamda):
    
    f = IUS(lamda, irr, k=1)  # k=1 gives linear interpolation

    #Plotting the raw and interpolated data on the same figure
    fig, ax = plt.subplots(1)
    fig.suptitle('Raw Spectra and Linear Inerpolation')
    
    ax.plot(lamda,irr,color='blue', lw=1.0, label='Raw irradiance spectra')
    ax.plot(lamda,f(lamda),color='orange',lw=1.0, label='Interpolated spectra',linestyle='dashed')
    
    plt.legend(loc='upper left',frameon=False)
    
    ax.set_xlim(350,450)
    ax.set_ylim(0,0.4)
    
    return f