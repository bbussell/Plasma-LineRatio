# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 20:13:41 2021

@author: beaub
"""
import time
start_time = time.time()
from datetime import datetime
import numpy as np
import os
from MetaResonantDensityCalculation import ExtractDensityInputs, CalculateNeutralDensity, PrepareResultsFile, ModelMetastableAndResonantDensity
from ElectronTemperatureCalculation import FetchModelElectronTemperatureInputs,CalculateElectronTemperature,PlotLoss
from CalculateElectronDensity import CalculateAllElectronDensity, PrintElectronDensityResultsToScreen, PrintElectronDensityResultsToFile
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

filelist = [#'5050PLSCHAMBERINLET_CHAMBERPOS_RF0001',
#              '5050PLSCHAMBERINLET_CHAMBERPOS_RF0002',
#              '5050PLSCHAMBERINLET_CHAMBERPOS_RF0003',
#              '5050PLSCHAMBERINLET_CHAMBERPOS_RF0004',
#              '5050PLSCHAMBERINLET_CHAMBERPOS_RF0005',
              '5050PLSCHAMBERINLET_CHAMBERPOS_RF0005',
#              '5050PLSCHAMBERINLET_CHAMBERPOS_RF0007',
             # '5050PLSCHAMBERINLET_CHAMBERPOS_RF0008',
              '5050PLSCHAMBERINLET_CHAMBERPOS_RF0008']
             # '5050PLSCHAMBERINLET_CHAMBERPOS_RF0010']
            # 'PLS-RF0011',
            # 'PLS-RF0012']
            # 'FRONT_RF0002',
            # 'FRONT_RF0003',
            # 'FRONT_RF0004',
            # 'FRONT_RF0005',
            # 'FRONT_RF0006',
            # 'FRONT_RF0007',
            # 'FRONT_RF0008']
            # 'FRONT_RF0009',
            # 'FRONT_RF0010',
            # 'FRONT_RF0011']

AllDensityInputValues = ExtractDensityInputs()
ResonantDensityIncm3 = AllDensityInputValues[0]
MetastableDensityIncm3 = AllDensityInputValues[1]
ExperimentalParameter = input("What system parameter was used during this experimental work? is this assessment for? E.g. 2kW or 0.0050 PP. Please respond and press Enter: ")
print("Thank you. The parameter you are about to evaluate is: ",ExperimentalParameter)

for File in filelist:
    path = r'C:\Users\beau.bussell\Google Drive\EngD\Research Data\OES\Calibrated\130421'
    os.chdir(path)
    NeutralDensity = CalculateNeutralDensity()
    PrepareResultsFile(File,ExperimentalParameter)
    
    print("Results file prepared")
    pause()    
    FullMetastableModelResults = ModelMetastableAndResonantDensity(ResonantDensityIncm3,MetastableDensityIncm3,File,ExperimentalParameter,path)
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
    FinalElectronTemperature = PlotLoss(File,ExperimentalParameter,GasTemperatureInKelvin_str)
    print("Minimum found.")
    print("Electron Temperature Modelling Complete")
    
    print("")
    print("Calculating Electron Density...")
    FinalElectronTemperature = np.array(FinalElectronTemperature)
    AllElectronDensities = CalculateAllElectronDensity(NormalisedIrradiance,FinalElectronTemperature[1],NeutralDensity)
    PrintElectronDensityResultsToScreen(AllElectronDensities,FinalElectronTemperature[1],File)
    PrintElectronDensityResultsToFile(AllElectronDensities,FinalElectronTemperature[1],ExperimentalParameter,File)
    print("")
    print("Electron Density Calculated. Results Printed to file.")
    
current_time = time.time() - start_time
print("The time taken to measure everything is:", current_time)
    