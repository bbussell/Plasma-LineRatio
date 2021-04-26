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

def ReadAndNormaliseIntegrationTime(fileN):
    FullFilename = fileN+'.tit.'

    with open(FullFilename, 'r') as fileR:
        print("")
        print("--------------------------------------------------------------")
        print("File: ", fileN)
        print("The parameter used in the experiment was: ")
        print("")
        print(fileR.readline())
        
        TopFileLine = next(fileR)
        TopFileLineStr = str(TopFileLine)
        LineMinusText1 = TopFileLineStr.replace("Integration time: ",'')
        print(LineMinusText1)
        LineMinusText2 = LineMinusText1.replace("ms", "")
        IntegrationTime = float(LineMinusText2)
        print(IntegrationTime)
        
        Wavelength, IrradianceWithoutNormalisation = np.loadtxt(FullFilename, comments='#', 
                        delimiter=';', skiprows=7, unpack=True, usecols=(0,4))
    
        NormalisedIrradiance = IrradianceWithoutNormalisation/IntegrationTime
        
    return NormalisedIrradiance, Wavelength

def InterpolateNormalisedSpectrum(NormalisedIrradiance,Wavelength):
    
    InterpolatedSpectrum = IUS(Wavelength, NormalisedIrradiance, k=1)  # k=1 gives linear interpolation

    #Plotting the raw and interpolated data on the same figure
    fig, ax = plt.subplots(1)
    fig.suptitle('Raw Spectra and Linear Inerpolation')
    
    ax.plot(Wavelength,NormalisedIrradiance,color='blue', lw=1.0, label='Raw irradiance spectra')
    ax.plot(Wavelength,InterpolatedSpectrum(Wavelength),color='orange',lw=1.0, label='Interpolated spectra',linestyle='dashed')
    
    plt.legend(loc='upper left',frameon=False)
    
    ax.set_xlim(350,450)
    ax.set_ylim(0,0.4)
    
    return InterpolatedSpectrum

# os.chdir(r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\181220')
# Result = ReadAndNormaliseIntegrationTime('RFPOWER0005')
# NormalisedIrradiance = Result[0]
# Wavelength = Result[1]
# InterpolatedSpectrum = InterpolateNormalisedSpectrum(NormalisedIrradiance,Wavelength)


