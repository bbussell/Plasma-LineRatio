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

def FetchSpectra(filename,path):
    os.chdir(path)
    NormalisedSpectraData = ReadAndNormaliseIntegrationTime(filename)
    NormalisedIrradiance = NormalisedSpectraData[0]
    Wavelength = NormalisedSpectraData[1]
    InterpolatedSpectrum = InterpolateNormalisedSpectrum(NormalisedIrradiance,Wavelength)
    
    return Wavelength, InterpolatedSpectrum

def Background_Calculation(LowerPeakLimit,UpperPeakLimit,InterpolatedSpectrum,LowerBackgroundPeakLimit,UpperBackgroundPeakLimit):
    RawAreaOfPeak = InterpolatedSpectrum.integral(LowerPeakLimit,UpperPeakLimit)
    PeakAreaMinusBackground = 0
    # if include772 == 0:
    #     PeakAreaMinusBackground = RawAreaOfPeak
    #     print("The Raw Area of Peak 772 is: ", PeakAreaMinusBackground)
    #     LowerIrradianceWithEdgeCorrection = 0
    #     UpperIrradianceWithEdgeCorrection = 0 
    # else:
    LowerIrradianceWithEdgeCorrection = 0.97*(InterpolatedSpectrum(LowerBackgroundPeakLimit))
    UpperIrradianceWithEdgeCorrection = 0.97*(InterpolatedSpectrum(UpperBackgroundPeakLimit))
    # print("The lower peak limit is: ",InterpolatedSpectrum(LowerPeakLimit))
    # print("The lower peak limit with edge correction is: ",LowerIrradianceWithEdgeCorrection)
    # print("")
    # print("The upper peak limit is: ",InterpolatedSpectrum(UpperPeakLimit))
    # print("The upper peak limit with edge correction is: ",UpperIrradianceWithEdgeCorrection)
    # print("")
    
    Background= 0.5*(LowerIrradianceWithEdgeCorrection+UpperIrradianceWithEdgeCorrection)*(UpperPeakLimit-LowerPeakLimit)
    print("The background is ", Background)
    print("The raw area is ", RawAreaOfPeak)
    BackgroundProportion = (Background/RawAreaOfPeak)*100
    # if Background > 0.7*RawAreaOfPeak:
    #     PeakAreaMinusBackground = RawAreaOfPeak
    #     print('The background is too large')
    
    if Background <= 0:
        PeakAreaMinusBackground = RawAreaOfPeak 
        print('The background is negative, no background substracted.')
    # elif BackgroundProportion > 50:
    #     PeakAreaMinusBackground = RawAreaOfPeak
    #     print("The BackgroundProportion was too large.")
    else:
        PeakAreaMinusBackground = RawAreaOfPeak - Background
        print('The background is positive and has been corrected')
        
    # if BackgroundProportion > 50:
    #     PeakAreaMinusBackground = RawAreaOfPeak
    
    if PeakAreaMinusBackground <= 0 :
        PeakAreaMinusBackground = RawAreaOfPeak
        print('The background was larger than the raw area. The background has')
        print('not been subtracted.')      
        
    print('The final area is:', PeakAreaMinusBackground)
    
    BackgroundProportion = (Background/RawAreaOfPeak)*100
    print("The background proportion of raw area is:", BackgroundProportion)
    
    return PeakAreaMinusBackground, LowerIrradianceWithEdgeCorrection, UpperIrradianceWithEdgeCorrection, RawAreaOfPeak

def SaveBackgroundCorrectedIrradianceToFile(EmissionPeakArray,array,File):
    #print("The final irradiance array sent to this function is: ", array)
    BackgroundCorrectedSpectra = pd.DataFrame(EmissionPeakArray,columns=["Wavelength"])
   # print("The background corrected wavelength is: ", BackgroundCorrectedSpectra)
    BackgroundCorrectedSpectra['Background Corrected Irradiance']=pd.Series(array)
    #print("The irradiance array to be printed is: ", array)
    #print("")
    print("The backround corrected spectra, including irradiance is: ", BackgroundCorrectedSpectra)
    
    with open(File+'_IntegratedIntensity.txt', 'w') as resultsfile:
        resultsfile.write('Datafile: '+File+'\n')   
    BackgroundCorrectedSpectra.to_csv(File+'_IntegratedIntensity.txt', index=False,mode="a")  
    return BackgroundCorrectedSpectra

def PlotPeakAndBackground(Wavelength, InterpolatedSpectrum,EmissionPeakArray,array):
    EmissionPeakLimits = [(694.8,700.3),(704.8,709.4),(713.7,716.7),(725.5,730),(736.1,740.9),(793.2,797.29),(745.1,757.9),(760.4,769.3),(770.7,775),(382.17,383.4),(360,361.2),(478.45,483.07)]
    EmissionPeakPlotLimits = [(691,701),(704,710),(711,720),(724.5,731),(732,746),(785,810),(744,760),(756,770),(767,780),(381,384.5),(358,363),(475,485)]
    InitialBackgroundPeakLimits = [(694.8,700.3),(704.8,709.4),(713.7,716.7),(725.5,730),(736.1,740.9),(793.2,797.29),(745.1,757.9),(760.4,769.3),(770.7,775),(382.17,383.4),(360,361.2),(478.45,483.07)]
    
    PeakPlots_x = [Wavelength,Wavelength,Wavelength,Wavelength,Wavelength,Wavelength,Wavelength,Wavelength,Wavelength,Wavelength,Wavelength,Wavelength]
    y = InterpolatedSpectrum(Wavelength)
    print(y)
    PeakPlots_y = [y,y,y,y,y,y,y,y,y,y,y,y]
    
    for i, plotlimit,PeakLimit,EmissionPeak,InitialPeakBackgroundLimit in zip(range(len(PeakPlots_x)),EmissionPeakPlotLimits,EmissionPeakLimits,EmissionPeakArray,InitialBackgroundPeakLimits):
        LowerPeakLimit = PeakLimit[0]
        UpperPeakLimit = PeakLimit[1]
        LowerBackgroundPeakLimit = InitialPeakBackgroundLimit[0]
        UpperBackgroundPeakLimit = InitialPeakBackgroundLimit[1]
        
        print("")
        print("This is the initial peak background result for Emission Peak: ",EmissionPeak)
        print("")
        InitialBackgroundCalculationResult = (Background_Calculation(LowerPeakLimit,UpperPeakLimit,InterpolatedSpectrum,LowerBackgroundPeakLimit,UpperBackgroundPeakLimit))
        PeakAreaMinusBackground = InitialBackgroundCalculationResult[0]
        LowerIrradianceWithEdgeCorrection = InitialBackgroundCalculationResult[1]
         #BackgroundCalculationResult[1]
        UpperIrradianceWithEdgeCorrection = InitialBackgroundCalculationResult[2]    
    
        fig = plt.figure(figsize=(10,8))
        axes = plt.subplot()
        axes.plot(PeakPlots_x[i],PeakPlots_y[i])
        axes.set_xlim(plotlimit[0],plotlimit[1])
        axes.set_ylim(0,0.02)
        axes.set_title(' Initial Background Fitting for '+ (str(EmissionPeak)))
        axes.set_xlabel("Wavelength (nm)", fontsize=10)
        axes.set_ylabel("Irradiance", fontsize=10)
       # axes.set_xticklabels(fontsize=8)
        
        for label in (axes.get_xticklabels() + axes.get_yticklabels()):
            label.set_fontsize(9)
        
        axes.axvline(x=LowerPeakLimit,lw=0.2,color='red')
        axes.axvline(x=UpperPeakLimit,lw=0.2,color='red')
       
        xx = LowerPeakLimit,UpperPeakLimit
        yy = (LowerIrradianceWithEdgeCorrection, UpperIrradianceWithEdgeCorrection)
        verticals = [(LowerPeakLimit,0), *zip(xx,yy), (UpperPeakLimit,0)]
        polys = Polygon(verticals,facecolor='0.9',edgecolor='0.5')
        axes.add_patch(polys)
        axes.grid(axis='x',linewidth=0.5)
    
        axes.xaxis.set_major_locator(plt.MaxNLocator(20))
        axes.xaxis.set_minor_locator(plt.MultipleLocator(5))
        plt.show()
        
        
        print("The initial background corrected area is: ", PeakAreaMinusBackground)
        print("")
        print("The limits used for the initial background calculationare: ")
        print("")
        print("Initial lower limit: ", LowerBackgroundPeakLimit)
        print("")
        print("Initial Upper Limit: ", UpperBackgroundPeakLimit)
        print("")
        UserDecisionOnBackground = input("Are you happy with the background area limits? Answer yes OR no and press Enter" )
        
        if UserDecisionOnBackground == 'yes':
            pass
        elif UserDecisionOnBackground == 'no':
            NewLowerBackgroundPeakLimit = float(input("Please enter your choice for the lower peak limit"))
            print("The lower peak limit you have chosen is: ", NewLowerBackgroundPeakLimit)
            NewUpperBackgroundPeakLimit = float(input("Please enter your choice for the upper peak limit")) 
            print("The upper peak limit you have chosen is: ", NewUpperBackgroundPeakLimit)
            print()
            print("Calculating new background with your chosen limits")
            print("")
            print("The second background calculation result is: ")
            BackgroundCalculationResult = (Background_Calculation(LowerPeakLimit,UpperPeakLimit,InterpolatedSpectrum,NewLowerBackgroundPeakLimit,NewUpperBackgroundPeakLimit))
            PeakAreaMinusBackground = BackgroundCalculationResult[0]
            LowerIrradianceWithEdgeCorrection = BackgroundCalculationResult[2]
         #BackgroundCalculationResult[1]
            UpperIrradianceWithEdgeCorrection = BackgroundCalculationResult[2]
            
            Newfigure = plt.figure(figsize=(10,8))
            Newaxes = plt.subplot()
            Newaxes.plot(PeakPlots_x[i],PeakPlots_y[i])
            Newaxes.set_xlim(plotlimit[0],plotlimit[1])
            Newaxes.set_ylim(0,0.02)
            Newaxes.set_title('New Background Choice for '+str(EmissionPeak))
            Newaxes.set_xlabel("Wavelength (nm)", fontsize=10)
            Newaxes.set_ylabel("Irradiance", fontsize=10)
       # axes.set_xticklabels(fontsize=8)
        
            for label in (Newaxes.get_xticklabels() + Newaxes.get_yticklabels()):
                label.set_fontsize(9)
            
            Newaxes.axvline(x=LowerPeakLimit,lw=0.2,color='red')
            Newaxes.axvline(x=UpperPeakLimit,lw=0.2,color='red')
       
            xx = LowerPeakLimit,UpperPeakLimit
            yy = (LowerIrradianceWithEdgeCorrection, UpperIrradianceWithEdgeCorrection)
            verticals = [(LowerPeakLimit,0), *zip(xx,yy), (UpperPeakLimit,0)]
            polys = Polygon(verticals,facecolor='0.9',edgecolor='0.5')
            Newaxes.add_patch(polys)
            Newaxes.grid(axis='x',linewidth=0.5)
    
            Newaxes.xaxis.set_major_locator(plt.MaxNLocator(20))
            plt.show()
            
        else:
            print("Your answer was not 'yes' or 'no', please amend your response")
            UserDecisionOnBackground = input("Are you happy with the background area limits? Answer yes OR no and press Enter" )
            print("")
            if UserDecisionOnBackground == 'yes':
                break
            elif UserDecisionOnBackground == 'no':
                NewLowerBackgroundPeakLimit = float(input("Please enter your choice for the lower peak limit"))
                print("The lower peak limit you have chosen is: ", NewLowerBackgroundPeakLimit)
                NewUpperBackgroundPeakLimit = float(input("Please enter your choice for the upper peak limit")) 
                print("The upper peak limit you have chosen is: ", NewUpperBackgroundPeakLimit)
                print()
                print("Calculating new background with your chosen limits")
                print("")
                print("The second background calculation result is: ")
                BackgroundCalculationResult = (Background_Calculation(LowerPeakLimit,UpperPeakLimit,InterpolatedSpectrum,NewLowerBackgroundPeakLimit,NewUpperBackgroundPeakLimit))
                PeakAreaMinusBackground = BackgroundCalculationResult[0]
                LowerIrradianceWithEdgeCorrection = BackgroundCalculationResult[2]
             #BackgroundCalculationResult[1]
                UpperIrradianceWithEdgeCorrection = BackgroundCalculationResult[2]
        array.append(PeakAreaMinusBackground)
        print("")
        #print("The Irradiance array for this result is: ",array)
        
        
        
    return PeakAreaMinusBackground, LowerIrradianceWithEdgeCorrection, UpperIrradianceWithEdgeCorrection, array

def BackgroundCalculationSeries(Wavelength,InterpolatedSpectrum,File):
    array = []
    EmissionPeakArray = [696,706,714,727,738,794,750,763,772,383,360,480]
    BackgroundCheckResult = PlotPeakAndBackground(Wavelength, InterpolatedSpectrum,EmissionPeakArray,array)
    PeakAreaMinusBackground = BackgroundCheckResult[0]
    FinalArray = BackgroundCheckResult[3]
    BackgroundCorrectedSpectra = SaveBackgroundCorrectedIrradianceToFile(EmissionPeakArray, FinalArray, File)

    return PeakAreaMinusBackground, BackgroundCorrectedSpectra

# MetastablePeakFig, axes = plt.subplots(nrows=2,ncols=3,sharey=True)
# ax696, ax706, ax714, ax727, ax738, ax794 = axes.flatten()

# AllOtherPeaksFig, axes = plt.subplots(nrows=2,ncols=3,sharey=True)
# ax750, ax763, ax772, ax383, ax360, ax480 = axes.flatten()

# filename='RFPOWER0001'
# path = r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\181220'
# SpectraResult = FetchSpectra(filename,path)
# Wavelength = SpectraResult[0]
# InterpolatedSpectrum = SpectraResult[1]
# BackgroundResult = BackgroundCalculationSeries(Wavelength, InterpolatedSpectrum,filename)
# PeakAreaMinusBackground = BackgroundResult[0]
# BackgroundCorrectedSpectra = BackgroundResult[1]


