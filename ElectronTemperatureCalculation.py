# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 17:36:18 2021

@author: beaub
"""

import time
import os
import sys
from math import exp
start_time = time.time()
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
#import RadTrap_No2 as RT # import rad trap code and run
from RadTrapAndEscapeFactorCalculation import radtrap, print_escape

#738 j=1s4 File=2p3
kG_738 = 2.20E-10
alphaG_738 = 0.59
EG_738 = 14.90

kM_738 = 6.09E-9
alphaM_738 = -0.11
EM_738 = 2.11

kR_738 = 8.23E-8
alphaR_738 = 0.01
ER_738 = 2.02

#750
kG_750 = 1.10E-9
alphaG_750 = 0.57
EG_750 = 14.96

kM_750 = 1.91E-8
alphaM_750 = -0.65
EM_750 = 2.60

kR_750 = 4.62E-8
alphaR_750 = 0.01
ER_750 = 2.02

#751.47
kG_751 = 4.49E-10
alphaG_751 = 0.58
EG_751 = 14.18

kM_751 = 7.87e-10
alphaM_751 = -0.59
EM_751 = 1.99

kR_751 = 6.29E-8
alphaR_751 = 0.01
ER_751 = 2.02

#763
kG_763 = 2.08E-10
alphaG_763 = 1.18
EG_763 = 12.74

kM_763 = 8.83E-8
alphaM_763 = 0.20
EM_763 = 1.88

kR_763 = 1.03E-7
alphaR_763 = 0.01
ER_763 = 2.02

#772.38
kG_772 = 1.04E-10
alphaG_772 = 0.69
EG_772 = 14.49

kM_772 = 1.03E-8
alphaM_772 = -0.08
EM_772 = 1.62

kR_772 = 3.85E-8
alphaR_772 = 0.01
ER_772 = 2.02

#772.42
kG_772_4 = 1.53E-10
alphaG_772_4 = 0.78
EG_772_4 = 14.61

kM_772_4 = 1.32E-8
alphaM_772_4 = 0.16
EM_772_4 = 1.91

kR_772_4 = 3.85E-8
alphaR_772_4 = 0.01
ER_772_4 = 2.02

#795
kG_795 = 2.39E-10
alphaG_795 = 0.73
EG_795 = 13.93

kM_795 = 3.93E-8
alphaM_795 = 0.08
EM_795 = 2.11

kR_795 = 5.15E-8
alphaR_795 = 0.01
ER_795 = 2.02

def str_to_class(str):
    return getattr(sys.modules[__name__], str)

def pause():
    programPause = input("Press the <ENTER> key to continue...")
    
def FetchCalculatedAtomDensities(FinalCalculatedDensity,File,ExperimentalParameter,GasTemperatureInKelvin_str):
    ResonantDensity_n1s4 = FinalCalculatedDensity[0]
    MetastableDensity_n1s5 = FinalCalculatedDensity[1]
    #LR's are: 738/750 ; 763/750 ; 772/750 ; 794/750
    n_ij, A, g_i, g_j = np.loadtxt("line_data2.txt", comments='#', delimiter=';', skiprows=2, unpack=True, 
                                        usecols=(0,1,2,3))
    Lower1sLevel = np.loadtxt("line_data2.txt", comments='#',dtype=str, delimiter=';', skiprows=2, unpack=True, 
                                        usecols=(4))
    #Taking n_m and n_r from calculation 
    MetastableDensity_n1s3 = MetastableDensity_n1s5/6.5
    ResonantDensity_n1s2 = ResonantDensity_n1s4
    for a, b, c, d, e in zip(n_ij,A,g_i,g_j,Lower1sLevel):
        if e=="1s2":
            N_j = ResonantDensity_n1s2
        elif e=="1s3":
            N_j = MetastableDensity_n1s3
        elif e=="1s4":
            N_j = ResonantDensity_n1s4
        elif e=="1s5":
            N_j = MetastableDensity_n1s5
        with open('live_density.txt', 'a+') as resultsfile:
            resultsfile.write("%4.2f %4.2f %3.1f %3.1f %4.2e \n" % (a,b,c,d,N_j))   
            
    n_ij, A, g_i, g_j, n_j = np.loadtxt("live_density.txt", comments='#', delimiter=' ', unpack=True, 
                                        usecols=(0,1,2,3,4))
    
def FetchModelElectronTemperatureInputs():
    ElectronTemperatureInputs = np.loadtxt("Te_Intervals.txt", unpack=True,
                              usecols=(0))
    return ElectronTemperatureInputs

def PrepareResultsFiles(datestring,File,ExperimentalParameter,GasTemperatureInKelvin_str,CharacteristicLength_str):
    with open('Te_results.txt', 'w') as resultsfile:
        resultsfile.write('Electron Temperature Results\n')
        resultsfile.write('Datetime (Y-M-S) = ' + datestring + '\n')
        resultsfile.write('Filename: ' + File + '\n')
        resultsfile.write('Experiment Name: ' + ExperimentalParameter + '\n')
        resultsfile.write('Gas Temperature (K) = ' + GasTemperatureInKelvin_str + '\n')
        resultsfile.write('Characteristic length (cm) = '+ CharacteristicLength_str + '\n')
        
    with open('LR_results.txt', 'w') as resultsfile:
        resultsfile.write('Electron Temperature Results - Raw Line Ratio Data\n')
        resultsfile.write('Datetime (Y-M-S) = ' + datestring + '\n')
        resultsfile.write('Filename: ' + File + '\n')
        resultsfile.write('Experiment Name: ' + ExperimentalParameter + '\n')
        resultsfile.write('Gas Temperature (K) = ' + GasTemperatureInKelvin_str + '\n')
        resultsfile.write('Characteristic length (cm) = '+ CharacteristicLength_str + '\n')    

def CalculateRadTrapCoefficients(ExperimentalParameter,CharacteristicLength,GasTemperatureInKelvin,File):
    AtomicMassArgon = 39.948# 6.6335209E-23 #kg
    GasConstant = 8.31446261815324*(1E3)
    GasTemperatureInKelvin_str = str(GasTemperatureInKelvin)
    GasTemperatureInKelvin = float(GasTemperatureInKelvin)
    
    n_ij, A, g_i, g_j, n_j = np.loadtxt("live_density.txt", comments='#', delimiter=' ', unpack=True,usecols=(0,1,2,3,4))
    os.rename("live_density.txt", ("live_density"+File+ExperimentalParameter+GasTemperatureInKelvin_str+'.txt'))
    radtrap(n_ij,A,g_i,g_j,n_j,CharacteristicLength,GasTemperatureInKelvin,AtomicMassArgon,GasConstant)
    n_ij, g_i, g_j, A, n_j, ko, k_ij = np.loadtxt("line_data_full.txt", comments='#', delimiter=' ', unpack=True, 
                                   usecols=(0,1,2,3,4,5,6))  

    print("")
    print("radtrap complete")
    print("")
    
    LineData = [n_ij,g_i,g_j,A,n_j,ko,k_ij]
    R_lam, Rad = np.genfromtxt("Rad_TrapCoeff.csv", delimiter=';', unpack=True, skip_header=1, usecols=(0,1))
    return Rad, LineData
    
def CalculateEscapeFactors(ExperimentalParameter,CharacteristicLength,GasTemperatureInKelvin,File):
    RadTrapResult = CalculateRadTrapCoefficients(ExperimentalParameter,CharacteristicLength,GasTemperatureInKelvin,File)
    
    Rad = RadTrapResult[0]
    LineData = RadTrapResult[1]
    
    n_ij = LineData[0]
    g_i = LineData[1]
    g_j = LineData[2]
    A = LineData[3]
    n_j = LineData[4]
    ko = LineData[5]
    k_ij = LineData[6]
    print_escape(n_ij,g_i,g_j,A,n_j,ko,k_ij,CharacteristicLength)
    print("")
    print("escape complete")
    print("")
    return Rad, LineData
    
def ex_rates(k,T,a,E):
    return k*((T)**a)*exp(-E/T)

def sum_nk(k,T,a,E,km,am,Em,kr,ar,Er,n_m,n_r,NeutralDensity):
    return NeutralDensity*ex_rates(k,T,a,E) + n_m*ex_rates(km,T,am,Em) +  n_r*ex_rates(kr,T,ar,Er)

def sum_750(ElectronTemperature,ResonantDensity,MetastableDensity,NeutralDensity,Rad):
    ResonantDensity_n1s2 = ResonantDensity
    MetastableDensity_n1s3 = MetastableDensity/6.5
    EffectiveResonantDensity = ResonantDensity + ResonantDensity_n1s2
    EffectiveMetastableDensity = MetastableDensity + MetastableDensity_n1s3
    sum_750=sum_nk(kG_750,ElectronTemperature,alphaG_750,EG_750,kM_750,alphaM_750,EM_750,kR_750,alphaR_750,ER_750,MetastableDensity,ResonantDensity,NeutralDensity)
    sum_751=sum_nk(kG_751,ElectronTemperature,alphaG_751,EG_751,kM_751,alphaM_751,EM_751,kR_751,alphaR_751,ER_751,MetastableDensity,ResonantDensity,NeutralDensity)
    R_750 = Rad[0]
    R_751 = Rad[4]
    BF_750_751 = (R_750*sum_750) + (R_751*sum_751)
    return BF_750_751

def LR_738(ElectronTemperature,ResonantDensity,MetastableDensity,Rad,NeutralDensity):
    ResonantDensity_n1s2 = ResonantDensity
    MetastableDensity_n1s3 = MetastableDensity/6.5
    EffectiveResonantDensity = ResonantDensity + ResonantDensity_n1s2
    EffectiveMetastableDensity = MetastableDensity + MetastableDensity_n1s3
    R_738 = Rad[2]
    BF_738 = R_738*(sum_nk(kG_738,ElectronTemperature,alphaG_738,EG_738,kM_738,alphaM_738,EM_738,kR_738,alphaR_738,ER_738,MetastableDensity,ResonantDensity,NeutralDensity))
    return ((BF_738/sum_750(ElectronTemperature,ResonantDensity,MetastableDensity,NeutralDensity,Rad)))
    
def LR_763(ElectronTemperature,ResonantDensity,MetastableDensity,Rad,NeutralDensity):
    ResonantDensity_n1s2 = ResonantDensity
    MetastableDensity_n1s3 = MetastableDensity/6.5
    EffectiveResonantDensity = ResonantDensity + ResonantDensity_n1s2
    EffectiveMetastableDensity = MetastableDensity + MetastableDensity_n1s3
    R_763 = Rad[5]
    BF_763 = R_763*(sum_nk(kG_763,ElectronTemperature,alphaG_763,EG_763,kM_763,alphaM_763,EM_763,kR_763,alphaR_763,ER_763,MetastableDensity,ResonantDensity,NeutralDensity))
    return BF_763/(sum_750(ElectronTemperature,ResonantDensity,MetastableDensity,NeutralDensity,Rad))

def LR_772(ElectronTemperature,ResonantDensity,MetastableDensity,Rad,NeutralDensity):
    ResonantDensity_n1s2 = ResonantDensity
    MetastableDensity_n1s3 = MetastableDensity/6.5
    EffectiveResonantDensity = ResonantDensity + ResonantDensity_n1s2
    EffectiveMetastableDensity = MetastableDensity + MetastableDensity_n1s3
    R_772_3 = Rad[6]
    R_772_4 = Rad[1]
    BF_772_3 = R_772_3*(sum_nk(kG_772,ElectronTemperature,alphaG_772,EG_772,kM_772,alphaM_772,EM_772,kR_772,alphaR_772,ER_772,MetastableDensity,ResonantDensity,NeutralDensity))
    BF_772_4 = R_772_4*(sum_nk(kG_772_4,ElectronTemperature,alphaG_772_4,EG_772_4,kM_772_4,alphaM_772_4,EM_772_4,kR_772_4,alphaR_772_4,ER_772_4,MetastableDensity,ResonantDensity,NeutralDensity))
    BF_772 = BF_772_3 + BF_772_4
    return BF_772/(sum_750(ElectronTemperature,ResonantDensity,MetastableDensity,NeutralDensity,Rad))

def LR_794(ElectronTemperature,ResonantDensity,MetastableDensity,Rad,NeutralDensity):
    ResonantDensity_n1s2 = ResonantDensity
    MetastableDensity_n1s3 = MetastableDensity/6.5
    EffectiveResonantDensity = ResonantDensity + ResonantDensity_n1s2
    EffectiveMetastableDensity = MetastableDensity + MetastableDensity_n1s3
    R_794 = Rad[3]
    BF_794 = R_794*(sum_nk(kG_795,ElectronTemperature,alphaG_795,EG_795,kM_795,alphaM_795,EM_795,kR_795,alphaR_795,ER_795,MetastableDensity,ResonantDensity,NeutralDensity))
    return BF_794/(sum_750(ElectronTemperature,ResonantDensity,MetastableDensity,NeutralDensity,Rad))

ElectronTemperatureLineRatios = [(738,750),(763,750),(772,750)]    

def chi_squared(LRm,LRe,err):
    return ((LRe-LRm)/(err*LRe))**2
    
def CalculateElectronTemperature(ElectronTemperatureInputs,FinalCalculatedDensity,NormalisedIrradiance,ExperimentalParameter,File,GasTemperatureInKelvin,datestring,CharacteristicLength,NeutralDensity):
    CharacteristicLength_str = str(CharacteristicLength)
    GasTemperatureInKelvin_str = str(GasTemperatureInKelvin)
    FetchCalculatedAtomDensities(FinalCalculatedDensity, File, ExperimentalParameter, GasTemperatureInKelvin_str)
    
    PrepareResultsFiles(datestring, File, ExperimentalParameter, GasTemperatureInKelvin_str, CharacteristicLength_str)
    ResonantDensity = FinalCalculatedDensity[0]
    MetastableDensity = FinalCalculatedDensity[1]
    Rad = CalculateEscapeFactors(ExperimentalParameter,CharacteristicLength,GasTemperatureInKelvin,File)[0]
    
    I_738 = NormalisedIrradiance[4]
    I_750 = NormalisedIrradiance[5]
    I_763 = NormalisedIrradiance[6]
    I_772 = NormalisedIrradiance[7]
    I_794 = NormalisedIrradiance[8]
    
    ExperimentalLineRatio_738 = I_738/I_750
    ExperimentalLineRatio_763 = I_763/I_750
    ExperimentalLineRatio_772 = I_772/I_750
    ExperimentalLineRatio_794 = I_794/I_750
    
    for ElectronTemperatureValue in ElectronTemperatureInputs:
        # for LineRatio in ElectronTemperatureLineRatios:
            
        #     UpperEmissionLine = LineRatio[0]
        #     LowerEmissionLine = LineRatio[1]
            
        ModelLineRatio_738 = LR_738(ElectronTemperatureValue,ResonantDensity,MetastableDensity,Rad,NeutralDensity)
        ModelLineRatio_763 = LR_763(ElectronTemperatureValue,ResonantDensity,MetastableDensity,Rad,NeutralDensity)
        ModelLineRatio_772 = LR_772(ElectronTemperatureValue,ResonantDensity,MetastableDensity,Rad,NeutralDensity)
        ModelLineRatio_794 = LR_794(ElectronTemperatureValue,ResonantDensity,MetastableDensity,Rad,NeutralDensity)
    
        ChiSqr_738 = chi_squared(ModelLineRatio_738,ExperimentalLineRatio_738,0.05)
        ChiSqr_763 = chi_squared(ModelLineRatio_763,ExperimentalLineRatio_763,0.05)
        ChiSqr_772 = chi_squared(ModelLineRatio_772,ExperimentalLineRatio_772,0.1)
        ChiSqr_794 = chi_squared(ModelLineRatio_794,ExperimentalLineRatio_794,0.05)
        ChiSqr_All = ChiSqr_738 + ChiSqr_763 + ChiSqr_772 + ChiSqr_794
        ChiSqr_Excluding772 = ChiSqr_738 + ChiSqr_763 + ChiSqr_794
        
        with open('Te_results.txt', 'a+') as te_results:
            te_results.write("%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n" % (ElectronTemperatureValue,ChiSqr_738, ChiSqr_763, ChiSqr_772, ChiSqr_794, ChiSqr_All, ChiSqr_Excluding772))
            
        with open("LR_results.txt", 'a+') as lr_results:
            lr_results.write("%4.2f %5.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f \n" % (ElectronTemperatureValue,ModelLineRatio_738,ExperimentalLineRatio_738,ModelLineRatio_763,ExperimentalLineRatio_763,ModelLineRatio_772,ExperimentalLineRatio_772,ModelLineRatio_794,ExperimentalLineRatio_794))
                
def FindMinimumLoss():
    
    ElectronTemperature, ChiSqr_738, ChiSqr_763, ChiSqr_772, ChiSqr_794, ChiSqr_All, ChiSqr_Excluding772 = np.genfromtxt("Te_results.txt", delimiter=" ", skip_header=6, unpack=True, usecols=(0,1,2,3,4,5,6))
    ElectronTemperatureLst = list(ElectronTemperature)
    ChiSqrList_738 = list(ChiSqr_738)
    ChiSqrList_763 = list(ChiSqr_763)
    ChiSqrList_772 = list(ChiSqr_772)
    ChiSqrList_794 = list(ChiSqr_794)    
    ChiSqrList_All = list(ChiSqr_All)
    ChiSqrList_Excluding772 = list(ChiSqr_Excluding772)
    
    SecondMinima_Te = ElectronTemperatureLst[70:]
    SecondMinimaLst_Excluding772 = ChiSqrList_Excluding772[70:]
    # print("second minima list:")
    # print(SecondMinimaLst_Excluding772)
    # print("")
    SecondMinimaArray_Excluding772 = np.array(SecondMinimaLst_Excluding772)
    
    print('------------------------------------------------------------------')
    print("The minimum in X^2 for all lines is: ", ChiSqrList_All[ChiSqrList_All.index(min(ChiSqrList_All))], "occurring when Te=", ElectronTemperature[ChiSqrList_All.index(min(ChiSqrList_All))])
    print("The minimum in X^2 for all excluding 772nm is: ", SecondMinimaLst_Excluding772[SecondMinimaLst_Excluding772.index(min(SecondMinimaLst_Excluding772))], "occurring when Te=", SecondMinima_Te[SecondMinimaLst_Excluding772.index(min(SecondMinimaLst_Excluding772))])
    print("The minimum in X^2 for 738nm is: ", ChiSqrList_738[ChiSqrList_738.index(min(ChiSqrList_738))], "occuring when Te=", ElectronTemperature[ChiSqrList_738.index(min(ChiSqrList_738))])
    print("The minimum in X^2 for 763nm is: ", ChiSqrList_763[ChiSqrList_763.index(min(ChiSqrList_763))], "occuring when Te=", ElectronTemperature[ChiSqrList_763.index(min(ChiSqrList_763))])
    print("The minimum in X^2 for 772nm is: ", ChiSqrList_772[ChiSqrList_772.index(min(ChiSqrList_772))], "occuring when Te=", ElectronTemperature[ChiSqrList_772.index(min(ChiSqrList_772))])
    print("The minimum in X^2 for 794nm is: ", ChiSqrList_794[ChiSqrList_794.index(min(ChiSqrList_794))], "occuring when Te=", ElectronTemperature[ChiSqrList_794.index(min(ChiSqrList_794))])
    
    ModelResultsInArrayFormat = [ElectronTemperature,ChiSqr_738,ChiSqr_763,ChiSqr_772, ChiSqr_794, ChiSqr_All,ChiSqr_Excluding772,SecondMinimaArray_Excluding772]
    ModelResultsInListFormat = [ElectronTemperature,ChiSqrList_738,ChiSqrList_763,ChiSqrList_772,ChiSqrList_794,ChiSqrList_All,ChiSqrList_Excluding772,SecondMinimaLst_Excluding772,SecondMinima_Te]
    
    return ModelResultsInArrayFormat, ModelResultsInListFormat

def PlotLoss(File,ExperimentalParameter,GasTemperatureInKelvin_str):
    AllElectronTemperatureModelResults = FindMinimumLoss()
    ModelResultsInArrayFormat = AllElectronTemperatureModelResults[0]
    ModelResultsInListFormat = AllElectronTemperatureModelResults[1]
    ElectronTemperature = ModelResultsInArrayFormat[0]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ElectronTemperature,ModelResultsInArrayFormat[1],c='b', label='X^2 for 738/750nm')
    ax.plot(ElectronTemperature,ModelResultsInArrayFormat[2],c='r', label='X^2 for 763/750nm')
    ax.plot(ElectronTemperature,ModelResultsInArrayFormat[3],c='g',label='X^2 for 772/750nm')
    ax.plot(ElectronTemperature,ModelResultsInArrayFormat[4],c='y',label='X^2 for 794/750nm')
    ax.plot(ElectronTemperature,ModelResultsInArrayFormat[5],c='magenta',label='X^2 for all lines')
    ax.plot(ElectronTemperature,ModelResultsInArrayFormat[6],c='orange',label='X^2 for all excl. 772nm')    
    ax.set_title('Effective Electron Temperature Chi-Sqruared '+ (str(File)))
    ax.set_xlabel("Electron Temperature (eV)", fontsize=10)
    ax.set_ylabel("Chi-Squared", fontsize=10)
    plt.legend(loc='upper right'); 
    ax.set_ylim([0,100])
    plt.show()
    
    os.rename("Te_results.txt", time.strftime("TeResults/te_results"+File+ExperimentalParameter+GasTemperatureInKelvin_str+"K_%Y%m%d%H%M%S.txt")) 
    os.rename("LR_results.txt", time.strftime("TeResults/LR_results"+File+ExperimentalParameter+GasTemperatureInKelvin_str+"_%Y%m%d%H%M%S.txt"))    
    SecondMinima_Te = ModelResultsInListFormat[8] 
    #print(SecondMinima_Te)
    FinalCalculatedElectronTemperature = [ElectronTemperature[ModelResultsInListFormat[5].index(min(ModelResultsInListFormat[5]))], SecondMinima_Te[ModelResultsInListFormat[7].index(min(ModelResultsInListFormat[7]))]]
    
    return FinalCalculatedElectronTemperature


# file='100CHAMBERINLET_RF0010'
# path = r'C:\Users\beau.bussell\Google Drive\EngD\Research Data\OES\Calibrated\130421'
# os.chdir(path)
# ExperimentalParameter = 'test'
# GasTemperatureInKelvin_str = str(600)

# PlotLoss(file,ExperimentalParameter,GasTemperatureInKelvin_str)
