# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 19:38:03 2021

@author: beaub
"""

import os 
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from IntTimeAndSpectaFit import ReadAndNormaliseIntegrationTime, InterpolateNormalisedSpectrum

def str_to_class(str):
    return getattr(sys.modules[__name__], str)

def FetchSpectra(filename):
    os.chdir(r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\181220')
    NormalisedSpectraData = ReadAndNormaliseIntegrationTime(filename)
    NormalisedIrradiance = NormalisedSpectraData[0]
    Wavelength = NormalisedSpectraData[1]
#NormalisedSpectraData_df = pd.DataFrame()
#NormalisedSpectraData_df['Wavelength']=pd.Series(Wavelength)
#NormalisedSpectraData_df['Normalised Irradiance']=pd.Series(NormalisedIrradiance)
    InterpolatedSpectrum = InterpolateNormalisedSpectrum(NormalisedIrradiance,Wavelength)
    
    return Wavelength, InterpolatedSpectrum

def Background_Calculation(LowerPeakLimit,UpperPeakLimit,InterpolatedSpectrum):
    RawAreaOfPeak = InterpolatedSpectrum.integral(LowerPeakLimit,UpperPeakLimit)
    Background= 0.5*(InterpolatedSpectrum(LowerPeakLimit)+InterpolatedSpectrum(UpperPeakLimit))*(UpperPeakLimit-LowerPeakLimit)
    print("The background is ", Background)
    print("The raw area is ", RawAreaOfPeak)
    PeakAreaMinusBackground = 0
    # if Background > 0.7*RawAreaOfPeak:
    #     PeakAreaMinusBackground = RawAreaOfPeak
    #     print('The background is too large')
    if Background <= 0 :
        PeakAreaMinusBackground = RawAreaOfPeak 
        print('The background is negative')
    else:
        PeakAreaMinusBackground = RawAreaOfPeak - Background
        print('The background is positive and has been corrected')
        
    print('The final area is:', PeakAreaMinusBackground)
    BackgroundProportion = (Background/RawAreaOfPeak)*100
    print("The background proportion of raw area is:", BackgroundProportion)
    
    return PeakAreaMinusBackground

def BackgroundCalculationSeries(Wavelength,InterpolatedSpectrum):
    EmissionPeakArray = [696,706,714,727,738,794,750,763,772,383,360,480]
    EmissionPeakLimits = [(694.8,700.3),(704.8,709.4),(713.7,717.5),(725,730),(736.1,742),(791.1,797.29),(745.1,757.9),(758,769.5),(770.7,775),(382.17,383.4),(360,361.2),(478.45,483.07)]
    EmissionPeakPlotLimits = [(692,701),(704,710),(711,720),(724.5,731),(732,746),(785,810),(744,760),(756,770),(767,780),(381,384.5),(358,363),(475,485)]
    
    for EmissionPeak, PlotLimit, PeakLimit in zip(EmissionPeakArray,EmissionPeakPlotLimits,EmissionPeakLimits):
        LowerPeakLimit = PeakLimit[0]
        
        UpperPeakLimit = PeakLimit[1]
        print("")
        print("Peak:",EmissionPeak,"nm")
        print("")
        PeakAreaMinusBackground = Background_Calculation(LowerPeakLimit, UpperPeakLimit,InterpolatedSpectrum)
        
        EmissionPeakPlotName = str_to_class('ax'+str(EmissionPeak))
       
        EmissionPeakPlotName.plot(Wavelength,InterpolatedSpectrum(Wavelength),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
        EmissionPeakPlotName.set_xlim(PlotLimit)
        EmissionPeakPlotName.set_ylim(0,0.05)
        EmissionPeakPlotName.set_title(str(EmissionPeak))
        
        EmissionPeakPlotName.axvline(x=LowerPeakLimit,lw=0.2,color='black')
        EmissionPeakPlotName.axvline(x=UpperPeakLimit,lw=0.2,color='black')
            
        ix = LowerPeakLimit,UpperPeakLimit
        iy = InterpolatedSpectrum(LowerPeakLimit), InterpolatedSpectrum(UpperPeakLimit)
        
        verts = [(LowerPeakLimit,0), *zip(ix,iy), (UpperPeakLimit,0)]
        poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
        EmissionPeakPlotName.add_patch(poly)
        plt.show
        
    return PeakAreaMinusBackground
    

# EmissionPeakBackgroundLimitsArray = [(694.7,700.3),(704.8,709.4),(713.5,717.5),(725,730),(732.9,745),(791.1,797.29),(745,758),(758,769.5),(769.8,779.5),(382.17,383.92),(360,361.5),(478.45,483.07)]
# EmissionPeakArray = [696,706,714,727,738,794,750,763,772,383,360,480]
# EmissionPeakLimits = [(694.8,700.3),(704.8,709.4),(713.7,717.5),(725,730),(736.1,742),(791.1,797.29),(745.1,757.9),(758,769.5),(770.7,775),(382.17,383.4),(360,361.2),(478.45,483.07)]
# EmissionPeakPlotLimits = [(692,701),(704,710),(711,720),(724.5,731),(732,746),(785,810),(744,760),(756,770),(767,780),(381,384.5),(358,363),(475,485)]

PeakAreaMinusBackgroundArray = []
MetastablePeakFig, axes = plt.subplots(nrows=2,ncols=3,sharey=True)
ax696, ax706, ax714, ax727, ax738, ax794 = axes.flatten()

AllOtherPeaksFig, axes = plt.subplots(nrows=2,ncols=3,sharey=True)
ax750, ax763, ax772, ax383, ax360, ax480 = axes.flatten()

SpectraResult = FetchSpectra('RFPOWER0005')
Wavelength = SpectraResult[0]
InterpolatedSpectrum = SpectraResult[1]
BackgroundCalculationSeries(Wavelength, InterpolatedSpectrum)








