# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 21:35:25 2020

@author: beaub
"""

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import auc

from scipy.interpolate import InterpolatedUnivariateSpline as IUS

#change directory path to one where files are located.

#os.chdir(r'..\..\OES\Calibrated\061020')

#wrapping the integral into a function so that any file
#can be used.

def integrate(fileN):
    #filename = input("What is the name of the file you want to analyse?")
    file = fileN+'.tit'
    
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
        
    f = IUS(lamda, irr, k=1)  # k=1 gives linear interpolation

    fig = plt.figure()
    est = fig.add_subplot(1,1,1)
    est.plot(lamda,f(lamda), color='tab:orange',lw=1.0,alpha=0.5)
    
    real=fig.add_subplot(1,1,1)
    
    real.plot(lamda,irr, color='tab:blue',lw=1.0,linestyle='dashed')
    est.set_xlim([350,450])
    real.set_xlim([350,450])
    est.set_ylim([0,0.4])

    #df2 = pd.DataFrame((lamda,Sp1,)
    
    #abc = IUS._call_(f,694.78,nu=0)
    
    RawArea_696 = f.integral(694.78, 698.15)
    bg_696 = 0.5*(f(694.78)+f(698.15))*(698.15-694.78)      
    print("The background for 696nm is:", bg_696)
    PeakA_696 = RawArea_696 - bg_696
    print('The original area is:',PeakA_696)
    print("Area minus background is:")
    print(PeakA_696)
    
    
    RawArea_705= f.integral(705.44,708.24)
    bg_705 = 0.5*(f(705.44)+f(708.24))*(708.24-704.44)      
    PeakA_705 = RawArea_705 - bg_705
    
    RawArea_714 = f.integral(713.84,716.64)
    bg_714 = 0.5*(f(713.84)+f(716.64))*(716.64-713.84)
    PeakA_714 = RawArea_714 - bg_714
    
    RawArea_727 = f.integral(725.58,728.94)
    bg_727 = 0.5*(f(725.58)+f(728.94))*(728.94-725.58)
    PeakA_727 = RawArea_727 - bg_727
    
    RawArea_738 = f.integral(736.19,740.65)
    bg_738 = 0.5*(f(736.19)+f(740.65))*(740.65-736.19)
    PeakA_738 = RawArea_738 - bg_738
    
    RawArea_750 = f.integral(747.9,754.02)
    bg_750 = 0.5*(f(747.9)+f(754.02))*(754.02-747.02)
    PeakA_750 = RawArea_750 - bg_750
    
    RawArea_763 = f.integral(761.25,766.25)
    bg_763 = 0.5*(f(761.25)+f(766.25))*(766.25-761.25)
    PeakA_763 = RawArea_763 - bg_763
    
    RawArea_772 = f.integral(770.14,775.13)
    bg_772 = 0.5*(f(770.14)+f(775.13))*(775.13-770.14)
    PeakA_772 = RawArea_772 - bg_772
    
    RawArea_794 = f.integral(792.86,797.29)
    bg_794 = 0.5*(f(792.86)+f(797.29))*(797.29-792.86)
    PeakA_794 = RawArea_794 - bg_794
    
    RawArea_451 = f.integral(449.5,451.88)
    bg_451 = 0.5*(f(449.5)+f(451.88))*(451.88-449.5)
    PeakA_451 = RawArea_451 - bg_451
    
    RawArea_419 = f.integral(418.3,421.1)
    bg_419 = 0.5*(f(418.3)+f(421.1))*(421.1-418.3)
    PeakA_419 = RawArea_419 - bg_419
    
    RawArea_367 = f.integral(366.42,368.75)
    bg_367 = 0.5*(f(366.42)+f(368.75))*(368.75-366.42)
    PeakA_367 = RawArea_367 - bg_367
    
    RawArea_415 = f.integral(414.77,417.09)
    bg_415 = 0.5*(f(414.77)+f(417.09))*(417.09-414.77)
    PeakA_415 = RawArea_415 - bg_415
    
    RawArea_383 = f.integral(382.17,383.92)
    bg_383 = 0.5*(f(382.17)+f(383.92))*(383.92-382.17)
    PeakA_383 = RawArea_383 - bg_383
    
    RawArea_365 = f.integral(364.66,366.42)
    bg_365 = 0.5*(f(364.66)+f(366.42))*(366.42-364.66)
    PeakA_365 = RawArea_365 - bg_365
    
    RawArea_360 = f.integral(359.99,361.74)
    bg_360 = 0.5*(f(359.99)+f(361.74))*(359.99-361.74)
    PeakA_360 = RawArea_360 - bg_360
    
    RawArea_488 = f.integral(487.1,488.83)
    bg_488 = 0.5*(f(487.1)+f(488.83))*(487.1-488.83)
    PeakA_488 = RawArea_488 - bg_488
    
    RawArea_425 = 0.73*(f.integral(424.6,428.7))
    bg_425 = 0.5*(f(424.6)+f(428.7))*(428.7-424.6)
    PeakA_425 = RawArea_425 - bg_425
    
    RawArea_480 = f.integral(478.45,483.07)
    bg_480 = 0.5*(f(478.45)+f(483.07))*(483.07-478.45)
    PeakA_480 = RawArea_480 - bg_480
    
    I = [PeakA_696,PeakA_705,PeakA_714, PeakA_727,PeakA_738,PeakA_750,PeakA_763,PeakA_772,PeakA_794,PeakA_451,PeakA_419,PeakA_367,PeakA_415,PeakA_383,PeakA_365,PeakA_360,PeakA_488,PeakA_425,PeakA_480]
    I_data = {'Wavelength (nm)':[696,706,714,727,738,750,763,772,794,451,419,367,415,383,365,360,488,425,480],'Integrated Intensity':I}
    #print("")
    #print("The result of using an integrate function is: ")
    #print("")
    df = pd.DataFrame(I_data,index=None,columns=["Wavelength (nm)","Integrated Intensity"])
    with open(fileN+'_IntegratedIntensity.txt', 'w') as resultsfile:
        resultsfile.write('Datafile: '+fileN+'\n')   
    df.to_csv(fileN+'_IntegratedIntensity.txt', index=False,mode="a")
    return fileN #for testing purposes swap fileN for df

#for testing purposes
# os.chdir(r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\181220')
# integral_result = integrate('RFPOWER0005')
# print(integral_result)










