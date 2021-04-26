# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 09:50:34 2021

@author: beaub
"""

import time
import os
import sys
start_time = time.time()
from math import exp,sqrt
import numpy as np
from background import FetchSpectra, BackgroundCalculationSeries
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from datetime import datetime
import pandas as pd
from astro_nist import Fetch_NIST

datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')

BoltzmannConstant = 1.38E-23*(1E4) #boltzman constant
#AtomicMassArgon = 9.109E-31
AtomicMassArgon = 39.948# 6.6335209E-23 #kg
#AtomicMassArgon = 6.6335209E-26 #atomic mass of Ar(kg)
GasConstant = 8.31446261815324*(1E3)
GasTemperatureInKelvin = 600 #757 = 15mtorr #gas temperature (K)
CharacteristicLength = 5.0 #charactersitic readsorption length (cm)

GasTemperatureInKelvin_str = str(GasTemperatureInKelvin)
CharacteristicLength_str = str(CharacteristicLength)

NIST_data = Fetch_NIST()

Aki = NIST_data['Aki']
#696 j=1s5
kij_696 = 1.43E-11
Aij_696 = Aki[0]
#727 j=1s4
kij_727 = 7.74E-12
Aij_727 = Aki[7]
#738 j=1s4
kij_738 = 6.25E-11
Aij_738 = Aki[9]
#706 j=1s5
kij_706 = 1.47E-11
Aij_706 = Aki[10]
#794 j=1s3
kij_794 = 3.07E-10
Aij_794 = Aki[12]
#714 j=1s5
kij_714 = 1.51E-12
Aij_714 = Aki[14]

def pause():
    programPause = input("Press the <ENTER> key to continue...")
    
def ExtractDensityInputs():
    
    with open('Density_inputs.txt', 'r') as DensityInputFile:
        ResonantDensityIncm3, MetastableDensityIncm3 = np.loadtxt("Density_inputs.txt", delimiter=';', 
                            unpack=True, usecols=(0, 1))
    return ResonantDensityIncm3, MetastableDensityIncm3

def CalculateReabsorptionCoefficient(k,n):
    return k*n*(GasTemperatureInKelvin**-0.5)

def CalculateEscapeFactor(k,n,CharacteristicLength):
    return (2-exp((-CalculateReabsorptionCoefficient(k,n)*CharacteristicLength)/1000))/(1+(CalculateReabsorptionCoefficient(k,n)*CharacteristicLength))

def CalculateModelBF_696(MetastableDensity,ResonantDensity,CharacteristicLength):
    
    return (Aij_696*CalculateEscapeFactor(kij_696,MetastableDensity,CharacteristicLength))/(Aij_727*CalculateEscapeFactor(kij_727,ResonantDensity,CharacteristicLength))

def CalculateModelBF_738(ResonantDensity,MetastableDensity,CharacteristicLength):
    
    return (Aij_738*CalculateEscapeFactor(kij_738,ResonantDensity,CharacteristicLength))/(Aij_706*CalculateEscapeFactor(kij_706,MetastableDensity,CharacteristicLength))

def CalculateModelBF_794(MetastableDensityn1s3,MetastableDensity,CharacteristicLength):
    return (Aij_794*CalculateEscapeFactor(kij_794,MetastableDensityn1s3,CharacteristicLength))/(Aij_714*CalculateEscapeFactor(kij_714,MetastableDensity,CharacteristicLength))


def CalculateExperimentalBranchingFraction(Iij,Iik):
    return Iij/Iik

def chi_squared(LRm,LRe,err):
    return ((LRe-LRm)/(err*LRe))**2

def raw(textfile):
    NormalisedIrradiance = np.loadtxt("RawInputData/" + textfile + ".txt", delimiter=' ', 
                  unpack=True, usecols=(1), skiprows=1) 
    return NormalisedIrradiance

def PrepareResultsFile(File,ExperimentalParameter):
    with open('mod_results.txt', 'w') as resultsfile:
        resultsfile.write('Metastable and Resonant Density Model Results\n')
        resultsfile.write('Datetime (Y-M-S) = ' + datestring + '\n')
        resultsfile.write('Filename:' + File + '\n')
        resultsfile.write('Experiment Name: ' + ExperimentalParameter + '\n')
        resultsfile.write('Gas Temperature (K) = ' + GasTemperatureInKelvin_str + '\n')
        resultsfile.write('Characteristic length (cm) = '+ CharacteristicLength_str + '\n')
        resultsfile.write('ResonantDensityIncm3 MetastableDensityIncm3 Chi-Squared \n')
        
def CalculateNeutralDensity():
    ProcessPressureInmbar = float(input("What is the process pressure in mbar?"))
    ProcessPressureinPa= ProcessPressureInmbar*100 #argon partial pressure in pascals
    ProcessPressureInmtorr = ProcessPressureInmbar/0.00133
    NeutralDensity = ((3.54E13)*ProcessPressureInmtorr*(273/GasTemperatureInKelvin)) #PP_mtorr/(BoltzmannConstant*GasTemperatureInKelvin)#6E13 
    print("the neutral density, calculated using PP, is: ", "%4.2e" % NeutralDensity)
    return NeutralDensity

def PrepareIntegratedIntensity(File,path):
    SpectraResult = FetchSpectra(File,path)
    Wavelength = SpectraResult[0]
    InterpolatedSpectrum = SpectraResult[1]
    BackgroundCalculationSeries(Wavelength, InterpolatedSpectrum, File)    
    WavelengthInNm, NormalisedIrradiance = np.loadtxt(File+'_IntegratedIntensity.txt', comments='#', delimiter=',', skiprows=2, unpack=True, 
                                        usecols=(0,1))
    os.chdir(r'..\..\..\Plasma Spec Code\RefactDevelopment\Plasma-LineRatio')
    return NormalisedIrradiance    

def CompareModelAndExperimentalLR(n_1s3,n_1s4,n_1s5,CharacteristicLength,File,NormalisedIrradiance):  
    I_696 = NormalisedIrradiance[0]
    I_706 = NormalisedIrradiance[1]
    I_714 = NormalisedIrradiance[2]
    I_727 = NormalisedIrradiance[3]
    I_738 = NormalisedIrradiance[4]
    I_794 = NormalisedIrradiance[8]
    AllModelLineRatios = [CalculateModelBF_696(n_1s5,n_1s4,CharacteristicLength),CalculateModelBF_738(n_1s4,n_1s5,CharacteristicLength),CalculateModelBF_794(n_1s3,n_1s5,CharacteristicLength)]
    AllExperimentalLineRatios = [CalculateExperimentalBranchingFraction(I_696,I_727),CalculateExperimentalBranchingFraction(I_738,I_706),CalculateExperimentalBranchingFraction(I_794,I_714)]
    
    return AllModelLineRatios, AllExperimentalLineRatios
    
def CalculateLoss(AllModelLineRatios,AllExperimentalLineRatios):
    Loss_696 = chi_squared(AllModelLineRatios[0],AllExperimentalLineRatios[0],0.05)
    Loss_738 = chi_squared(AllModelLineRatios[1],AllExperimentalLineRatios[1],0.05)
    Loss_794 = chi_squared(AllModelLineRatios[2], AllExperimentalLineRatios[2],0.05)
    TotalBFLoss = sum([Loss_696,Loss_738,Loss_794])
    return TotalBFLoss

def FindMinimumLossFromAllResults(ModelResonantDensity,ModelMetastableDensity,FullChiSquared,File,ExperimentalParameter):
    
    FullChiSquaredList = list(FullChiSquared)
    ResonantDensityList = list(ModelResonantDensity)
    MetastableDensityList = list(ModelMetastableDensity)
    MinimumChiSquaredIndex = FullChiSquaredList.index(min(FullChiSquaredList))

    print('Minimum chi-squared of ', min(FullChiSquared), 'occurs at position ', MinimumChiSquaredIndex, 
          'where n1s4 = ',"%7.1e" % ResonantDensityList[MinimumChiSquaredIndex], 'cm\u00b3', 'and MetastableDensityIncm3 = ',
          "%7.1e" % MetastableDensityList[MinimumChiSquaredIndex], 'cm\u00b3') 
    
    ResonantDensity_n1s2 = ResonantDensityList[MinimumChiSquaredIndex]
    ResonantDensity_n1s4 = ResonantDensityList[MinimumChiSquaredIndex]
    MetastableDensity_n1s3 = (MetastableDensityList[MinimumChiSquaredIndex])/6.5
    MetastableDensity_n1s5 = MetastableDensityList[MinimumChiSquaredIndex]
    FinalResonantDensity = ResonantDensity_n1s4 + ResonantDensity_n1s2 
    FinalMetastableDensity = MetastableDensity_n1s5 + MetastableDensity_n1s3 
    

    FinalCalculatedDensity = [ResonantDensity_n1s4,MetastableDensity_n1s5]
    FinalLoss = min(FullChiSquared)
    
    with open('MetaRes_Results.txt', 'w+') as finalresults:
        #for i in range(len(ResonantDensityIncm3)):
         #   mod_results.write("%d %d %d\n" % (ResonantDensityIncm3[i],MetastableDensityIncm3[i],chi_sum))
         finalresults.write('Metastable and Resonant Density Model Results\n')
         finalresults.write('Datetime (Y-M-S) = ' + datestring + '\n')
         finalresults.write('Filename:' + File + '\n')
         finalresults.write('Experiment Name: ' + ExperimentalParameter + '\n')
         finalresults.write('Gas Temperature (K) = ' + GasTemperatureInKelvin_str + '\n')
         finalresults.write('Characteristic length (cm) = '+ CharacteristicLength_str + '\n')
         finalresults.write('ResonantDensityIncm3_n1s2 ResonantDensityIncm3_n1s4 MetastableDensityIncm3_n1s3 MetastableDensity_n1s5 Chi-Squared \n')
         finalresults.write("%7.2e %7.2e %7.2e %7.2e %4.2f\n" % (ResonantDensity_n1s2,ResonantDensity_n1s4,MetastableDensity_n1s3,MetastableDensity_n1s5,FinalLoss))
    
    return FinalCalculatedDensity


def ModelMetastableAndResonantDensity(ResonantDensityIncm3,MetastableDensityIncm3,File,ExperimentalParameter,path):
    NormalisedIrradiance = PrepareIntegratedIntensity(File,path)
    for ResonantModelValue, MetastableModelValue in zip(ResonantDensityIncm3,MetastableDensityIncm3):
        
        n_1s4 = ResonantModelValue #resonant density (cm^-3)      
        n_1s5 = MetastableModelValue #metastable density (cm^-3)
        
        #-----------------------------------
        #Relation between total m and r 1sx levels
        n_1s3 = n_1s5/6.5
        n_1s2 = n_1s4
        #Setting model LR's
        
        #TESTING PURPOSES
        # print("For a resonant density of: ", n_1s4)
        # print("And a metastable density of: ", n_1s5)
        # print("")
               
        ModelAndExperimentalLineRatios = CompareModelAndExperimentalLR(n_1s3, n_1s4, n_1s5, CharacteristicLength,File,NormalisedIrradiance)
        AllModelLineRatios = ModelAndExperimentalLineRatios[0]
        AllExperimentalLineRatios = ModelAndExperimentalLineRatios[1]
        #print('The 696/727 model LR is:')
        #print(CalculateModelBranchingFraction(kij_696,kij_727,Aij_696,Aij_727,n_1s5,n_1s4,CharacteristicLength))    
       
        TotalBFLoss = CalculateLoss(AllModelLineRatios,AllExperimentalLineRatios)
        # print("The model line ratios are: ", AllModelLineRatios)
        # print("")
        # print("The experimental lin ratios are:", AllExperimentalLineRatios)
        # print("")
        # #print('Chi squared for n_r= ',n_1s4, 'and n_m= ',n_1s5, 'is: ',chi_sum)
        # print("Printing model result to file ")
        
        with open('mod_results.txt', 'a+') as mod_results:
        #for i in range(len(ResonantDensityIncm3)):
         #   mod_results.write("%d %d %d\n" % (ResonantDensityIncm3[i],MetastableDensityIncm3[i],chi_sum))
            mod_results.write("%7.2e %7.2e %6.2f\n" % (ResonantModelValue,MetastableModelValue,TotalBFLoss))
            
            #printing results to screen
        #plotting chi-sum as 3D plot
        #fig3d = plt.figure()
        
        #renaming the file with the current date and time
    with open('mod_results.txt', 'a+') as mod_results:
        ResonantDensity, MetastableDensity, BFLoss = np.loadtxt("mod_results.txt",unpack=True,usecols=(0, 1, 2),skiprows=7)
    MetastableResonantPlot3D = plt.axes(projection="3d")
    MetastableResonantPlot3D.plot3D(ResonantDensity,MetastableDensity,BFLoss,'black')
    plt.show()
    
    FinalCalculatedDensity = FindMinimumLossFromAllResults(ResonantDensity,MetastableDensity,BFLoss,File,ExperimentalParameter) 
    os.rename("mod_results.txt", time.strftime("MetaResResults/n1sDEN_"+File+ExperimentalParameter+"_%Y%m%d%H%M.txt"))  #calling function to print results 
    return FinalCalculatedDensity, TotalBFLoss, NormalisedIrradiance
    
# filelist = ['RFPOWER0005']

# NIST_data = Fetch_NIST()
# Aki = NIST_data['Aki']
# #696 j=1s5
# kij_696 = 1.43E-11
# Aij_696 = Aki[0]
# #727 j=1s4
# kij_727 = 7.74E-12
# Aij_727 = Aki[7]
# #738 j=1s4
# kij_738 = 6.25E-11
# Aij_738 = Aki[9]
# #706 j=1s5
# kij_706 = 1.47E-11
# Aij_706 = Aki[10]
# #794 j=1s3
# kij_794 = 3.07E-10
# Aij_794 = Aki[12]
# #714 j=1s5
# kij_714 = 1.51E-12
# Aij_714 = Aki[14]

# ExperimentalParameter = input("What system parameter was used during this experimental work? is this assessment for? E.g. 2kW or 0.0050 PP. Please respond and press Enter: ")

# print("Thank you. The parameter you are about to evaluate is: ",ExperimentalParameter)
# AllDensityInputValues = ExtractDensityInputs()
# ResonantDensityIncm3 = AllDensityInputValues[0]
# MetastableDensityIncm3 = AllDensityInputValues[1]

# for File in filelist:
#     os.chdir(r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\181220')
#     NormalisedIrradiance = PrepareIntegratedIntensity(File)
    
#     I_696 = NormalisedIrradiance[0]
#     I_706 = NormalisedIrradiance[1]
#     I_714 = NormalisedIrradiance[2]
#     I_727 = NormalisedIrradiance[3]
#     I_738 = NormalisedIrradiance[4]
#     I_794 = NormalisedIrradiance[8]

#     CalculateNeutralDensity()
#     PrepareResultsFile()
#     FinalCalculatedDensity = ModelMetastableAndResonantDensity(ResonantDensityIncm3,MetastableDensityIncm3)[0]
    
# current_time = time.time() - start_time
# print("The time taken to measure resonant and metastable density is:", current_time)
    