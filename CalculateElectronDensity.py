# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 19:05:15 2021

@author: beaub
"""

from math import exp
import os
import time
from datetime import datetime

def CalculateOpticalExcitationRateTemperatureDependence(ElectronTemperature):
    return 567 - (160*ElectronTemperature)+(16*(ElectronTemperature**2))

def CalculateElectronDensityAtom(Irradiance_2p,Irradiance_4p,TemperatureDependence,CriticalDensity_2p,CriticalDensity_4p):
    return (1-((Irradiance_2p/Irradiance_4p)/(TemperatureDependence)))/(((Irradiance_2p/Irradiance_4p)/(TemperatureDependence*CriticalDensity_2p))-(1/CriticalDensity_4p))

def CalculateElectronDensityIon(ElectronTemperature,Irradiance_Ion,Irradiance_Atom,NeutralDensity):
    # Irradiance_Ion = float(Irradiance_Ion)
    # Irradiance_Atom = float(Irradiance_Atom)
    # NeutralDensity = float(NeutralDensity)
    print("")
    print("For an electron temperature of ",ElectronTemperature)
    ExcitationRate_Atom = float((2.54E-10)*(ElectronTemperature**1.07)*(exp(-9.39/ElectronTemperature)))
    ExcitationRate_Ion = float(exp(-23.44-37.25*(exp(-ElectronTemperature/1.32))))
    return (ExcitationRate_Atom/ExcitationRate_Ion)*(Irradiance_Ion/Irradiance_Atom)*NeutralDensity  
    
def CalculateAllElectronDensity(NormalisedIrradiance,ElectronTemperature,NeutralDensity):
    nec_750 = 3E12
    nec_383 = 9E10
    nec_360 = 1E11
    
    I_750 = NormalisedIrradiance[6]
    I_383 = NormalisedIrradiance[9]
    I_360 = NormalisedIrradiance[10]
    I_480 = NormalisedIrradiance[11]
    
    ElectronDensity_383 = CalculateElectronDensityAtom(I_750, I_383, CalculateOpticalExcitationRateTemperatureDependence(ElectronTemperature), nec_750, nec_383)
    ElectronDensity_360 = CalculateElectronDensityAtom(I_750, I_360, CalculateOpticalExcitationRateTemperatureDependence(ElectronTemperature), nec_750, nec_360)
    ElectronDensity_480 = CalculateElectronDensityIon(ElectronTemperature,I_480,I_750,NeutralDensity)

    return ElectronDensity_383,ElectronDensity_360,ElectronDensity_480

def PrintElectronDensityResultsToScreen(AllElectronDensities,ElectronTemperature,File): 
    print("the electron density for using 383nm and T=",ElectronTemperature,"and file", File, " is: ", "%6.2e" % AllElectronDensities[0], "cm-3")
    print("the electron density for using 360nm and T=",ElectronTemperature,"and file", File, " is: ", "%6.2e" % AllElectronDensities[1], "cm-3")
    print("the electron density for the 480nm/750nm pair, with T=",ElectronTemperature, "and file", File, "is ", "%6.2e" % AllElectronDensities[2], "cm-3")
    
def PrintElectronDensityResultsToFile(AllElectronDensities,ElectronTemperature,ExperimentalParameter,File):
    ElectronDensity_383_str = str("%8.2e" %AllElectronDensities[0])
    ElectronDensity_360_str = str("%8.2e" %AllElectronDensities[1])
    ElectronDensity_480_str = str("%8.2e" %AllElectronDensities[2])
    ElectronTemperature_str = str(ElectronTemperature)
    
    with open(File+ElectronTemperature_str+ExperimentalParameter+'E_density.txt', 'w+') as resultsfile:
        resultsfile.write('Datafile: '+File+'\n')
        resultsfile.write('Experiment Name: '+ExperimentalParameter+'\n')
        resultsfile.write('Temperature value: '+ElectronTemperature_str+'\n')
        resultsfile.write('Emission Line (nm);Electron density (cm^-3) \n')
        resultsfile.write('383;'+ ElectronDensity_383_str+'\n')
        resultsfile.write('360;'+ ElectronDensity_360_str+'\n')
        resultsfile.write('480;' + ElectronDensity_480_str+'\n')

    os.rename(File+ElectronTemperature_str+ExperimentalParameter+'E_density.txt', time.strftime("EDensity_Results/"+File+ElectronTemperature_str+ExperimentalParameter+"K_%Y%m%d%H%M%S.txt")) 
     
    # N_e383 = abs(n_e383)
    # N_e360 = abs(n_e360)
    

# AllElectronDensities = CalculateAllElectronDensity(NormalisedIrradiance,FinalElectronTemperature[1],File)
# PrintElectronDensityResultsToScreen
# PrintElectronDensityResultsToFile