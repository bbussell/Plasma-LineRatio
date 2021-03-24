# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 20:13:41 2021

@author: beaub
"""
import time
start_time = time.time()
from datetime import datetime
import os
from MetaResonantDensityCalculation import ExtractDensityInputs, CalculateNeutralDensity, PrepareResultsFile, ModelMetastableAndResonantDensity
from ElectronTemperatureCalculation import FetchModelElectronTemperatureInputs,CalculateElectronTemperature,PlotLoss

def pause():
    programPause = input("Press the <ENTER> key to continue...")

datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')
BoltzmannConstant = 1.38E-23*(1E4) 
AtomicMassArgon = 39.948
GasConstant = 8.31446261815324*(1E3)
GasTemperatureInKelvin = 600
CharacteristicLength = 5

GasTemperatureInKelvin_str = str(GasTemperatureInKelvin)
CharacteristicLength_str = str(CharacteristicLength)

filelist = ['RFPOWER0005']

ExperimentalParameter = input("What system parameter was used during this experimental work? is this assessment for? E.g. 2kW or 0.0050 PP. Please respond and press Enter: ")

print("Thank you. The parameter you are about to evaluate is: ",ExperimentalParameter)
AllDensityInputValues = ExtractDensityInputs()
ResonantDensityIncm3 = AllDensityInputValues[0]
MetastableDensityIncm3 = AllDensityInputValues[1]

for File in filelist:
    os.chdir(r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\181220')
    NeutralDensity = CalculateNeutralDensity()
    PrepareResultsFile(File,ExperimentalParameter)
    print("Results file prepared")
    pause()    
    FullMetastableModelResults = ModelMetastableAndResonantDensity(ResonantDensityIncm3,MetastableDensityIncm3,File,ExperimentalParameter)
    FinalCalculatedDensity = FullMetastableModelResults[0]
    NormalisedIrradiance = FullMetastableModelResults[2]
    print("")
    print("Resonant and Metastable Modelling Complete.")
    print("The final density is: ", FinalCalculatedDensity)
    
    print("")
    print("Beginning Electron Temperature Modelling...")
    ElectronTemperatureInputs = FetchModelElectronTemperatureInputs()
    CalculateElectronTemperature(ElectronTemperatureInputs,FinalCalculatedDensity,NormalisedIrradiance,ExperimentalParameter,File,GasTemperatureInKelvin_str,datestring,CharacteristicLength_str,NeutralDensity)
    print("")   
    print("Finding minimum in X^2 and corresponding Electron Temperature...")
    PlotLoss(File,ExperimentalParameter,GasTemperatureInKelvin_str)
    print("Minimum found.")
    print("Electron Temperature Modelling Complete")

current_time = time.time() - start_time
print("The time taken to measure resonant and metastable density is:", current_time)
    