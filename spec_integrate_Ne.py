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

    #Plotting the raw and interpolated data on the same figure
    fig, ax = plt.subplots(1)
    fig.suptitle('Raw Spectra and Linear Inerpolation')
    
    ax.plot(lamda,irr,color='blue', lw=1.0, label='Raw irradiance spectra')
    ax.plot(lamda,f(lamda),color='orange',lw=1.0, label='Interpolated spectra',linestyle='dashed')
    
    plt.legend(loc='upper left',frameon=False)
    
    ax.set_xlim(350,450)
    ax.set_ylim(0,0.4)

    
    #abc = IUS._call_(f,694.78,nu=0)
    
    #wave = ['696','705','714','727','738','750','763','772','794','451','419','367','415','383','365','360','488','425','480']
    
    #Backgorund calculation and subtraction
    RawArea_696 = f.integral(694.78, 698.15)
    bg_696 = 0.5*(f(694.78)+f(698.15))*(698.15-694.78)      
    print("The background for 696nm is:", bg_696)
    PeakA_696 = 0
    if bg_696 > 0:
        print('The background is positive')
        PeakA_696 = RawArea_696 - bg_696
        print('The original area is:',RawArea_696)
    else: 
        PeakA_696 = RawArea_696 
        print('The background is negative')
    PeakA_696 = PeakA_696
    print('The final area is:', PeakA_696)
    
    PeakA_705 = 0
    RawArea_705= f.integral(705.44,708.24)
    bg_705 = 0.5*(f(705.44)+f(708.24))*(708.24-704.44)   
    if bg_705 > 0:
        PeakA_705 = RawArea_705 - bg_705
    else: PeakA_705 = RawArea_705
    
    PeakA_714 = 0
    RawArea_714 = f.integral(713.84,716.64)
    bg_714 = 0.5*(f(713.84)+f(716.64))*(716.64-713.84)
    if bg_714 > 0:
        PeakA_714 = RawArea_714 - bg_714
    else: PeakA_714 = RawArea_714
    
    PeakA_727 = 0
    RawArea_727 = f.integral(725.58,728.94)
    bg_727 = 0.5*(f(725.58)+f(728.94))*(728.94-725.58)
    if bg_727 > 0:
        PeakA_727 = RawArea_727 - bg_727
    else: PeakA_727 = RawArea_727
    
    PeakA_738 = 0
    RawArea_738 = f.integral(736.19,740.65)
    bg_738 = 0.5*(f(736.19)+f(740.65))*(740.65-736.19)
    if bg_738 > 0:
        PeakA_738 = RawArea_738 - bg_738
    else: PeakA_738 = RawArea_738
    
    PeakA_750 = 0
    RawArea_750 = f.integral(747.9,754.02)
    bg_750 = 0.5*(f(747.9)+f(754.02))*(754.02-747.02)
    if bg_750 > 0:
        PeakA_750 = RawArea_750 - bg_750
    else: PeakA_750 = RawArea_750
    
    PeakA_763 = 0
    RawArea_763 = f.integral(761.25,766.25)
    bg_763 = 0.5*(f(761.25)+f(766.25))*(766.25-761.25)
    if bg_763 > 0:
        PeakA_763 = RawArea_763 - bg_763
    else: PeakA_763 = RawArea_763
    
    PeakA_772 = 0
    RawArea_772 = f.integral(770.14,775.13)
    bg_772 = 0.5*(f(770.14)+f(775.13))*(775.13-770.14)
    if bg_772 > 0:
        PeakA_772 = RawArea_772 - bg_772
    else: PeakA_772 = RawArea_772
    
    PeakA_794 = 0
    RawArea_794 = f.integral(792.86,797.29)
    bg_794 = 0.5*(f(792.86)+f(797.29))*(797.29-792.86)
    if bg_794 > 0:
        PeakA_794 = RawArea_794 - bg_794
    else: PeakA_794 = RawArea_794
    
    PeakA_451 = 0
    RawArea_451 = f.integral(449.5,451.88)
    bg_451 = 0.5*(f(449.5)+f(451.88))*(451.88-449.5)
    if bg_451 > 0:
        PeakA_451 = RawArea_451 - bg_451
    else: PeakA_451 = RawArea_451
    
    PeakA_419 = 0
    RawArea_419 = f.integral(418.3,421.1)
    bg_419 = 0.5*(f(418.3)+f(421.1))*(421.1-418.3)
    if bg_419 > 0:
        PeakA_419 = RawArea_419 - bg_419
    else: PeakA_419 = RawArea_419
    
    PeakA_367 = 0
    RawArea_367 = f.integral(366.42,368.75)
    bg_367 = 0.5*(f(366.42)+f(368.75))*(368.75-366.42)
    if bg_367 > 0:
        PeakA_367 = RawArea_367 - bg_367
    else: PeakA_367 = RawArea_367
    
    PeakA_415 = 0
    RawArea_415 = f.integral(414.77,417.09)
    bg_415 = 0.5*(f(414.77)+f(417.09))*(417.09-414.77)
    if bg_415 > 0:
        PeakA_415 = RawArea_415 - bg_415
    else: PeakA_415 = RawArea_415
    
    PeakA_383 = 0
    RawArea_383 = f.integral(382.17,383.92)
    bg_383 = 0.5*(f(382.17)+f(383.92))*(383.92-382.17)
    if bg_383 > 0:
        PeakA_383 = RawArea_383 - bg_383
    else: PeakA_383 = RawArea_383
    
    PeakA_365 = 0
    RawArea_365 = f.integral(364.66,366.42)
    bg_365 = 0.5*(f(364.66)+f(366.42))*(366.42-364.66)
    if bg_365 > 0:
        PeakA_365 = RawArea_365 - bg_365
    else: PeakA_365 = RawArea_365
    
    PeakA_360 = 0
    RawArea_360 = f.integral(359.99,361.74)
    bg_360 = 0.5*(f(359.99)+f(361.74))*(359.99-361.74)
    if bg_360 > 0:
        PeakA_360 = RawArea_360 - bg_360
    else: PeakA_360 = RawArea_360
    
    PeakA_488 = 0
    RawArea_488 = f.integral(487.1,488.83)
    bg_488 = 0.5*(f(487.1)+f(488.83))*(487.1-488.83)
    if bg_488 > 0:
        PeakA_488 = RawArea_488 - bg_488
    else: PeakA_488 = RawArea_488
    
    PeakA_425 = 0
    RawArea_425 = 0.73*(f.integral(424.6,428.7))
    bg_425 = 0.5*(f(424.6)+f(428.7))*(428.7-424.6)
    if bg_425 > 0:
        PeakA_425 = RawArea_425 - bg_425
    else: PeakA_425 = RawArea_425
    
    PeakA_480 = 0
    RawArea_480 = f.integral(478.45,483.07)
    bg_480 = 0.5*(f(478.45)+f(483.07))*(483.07-478.45)
    if bg_480 > 0:
        PeakA_480 = RawArea_480 - bg_480
    else: PeakA_480 = RawArea_480
    
    I = [PeakA_696,PeakA_705,PeakA_714, PeakA_727,PeakA_738,PeakA_750,PeakA_763,PeakA_772,PeakA_794,PeakA_451,PeakA_419,PeakA_367,PeakA_415,PeakA_383,PeakA_365,PeakA_360,PeakA_488,PeakA_425,PeakA_480]
    I_data = {'Wavelength (nm)':[696,706,714,727,738,750,763,772,794,451,419,367,415,383,365,360,488,425,480],'Integrated Intensity':I}
    #print("")
    #print("The result of using an integrate function is: ")
    #print("")
    df = pd.DataFrame(I_data,index=None,columns=["Wavelength (nm)","Integrated Intensity"])
    with open(fileN+'_IntegratedIntensity.txt', 'w') as resultsfile:
        resultsfile.write('Datafile: '+fileN+'\n')   
    df.to_csv(fileN+'_IntegratedIntensity.txt', index=False,mode="a")
    
    #Plotting the raw and interpolated data on the same figure
    #fig2, axes = plt.subplots(nrows=5,ncols=4)
    fig2, axes = plt.subplots(nrows=2,ncols=3,sharey=True)
    fig.suptitle('Background subtraction')
    
    #ax696,ax706,ax714,ax727,ax738,ax750,ax763,ax772,ax794,ax451,ax419,ax367,ax415,ax383,ax365,ax360,ax488,ax425,ax480,axblank = axes.flatten()
    ax696, ax706, ax714, ax727, ax738, ax794 = axes.flatten()
      
    ax696.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    # ax696.add_patch(
    #     patches.Rectangle(
    #         xy=(694.78,0),
    #         width=width,
    #         height=height,
    #         linewidth=1,
    #         color='red',
    #         fill=False
    #         )
    #     )
    ax696.set_xlim(694,701)
    ax696.set_ylim(0,0.06)
    ax696.set_title('696nm Background Subtraction')
    #--- poly plot
    ix = 694.78,700.5
    iy = f(700.5),f(700.5)
    verts = [(694.78,0), *zip(ix,iy), (700.5,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax696.add_patch(poly)
    
    ax706.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax706.set_xlim(704,710)
    ax706.set_ylim(0,0.06)
    ax706.set_title('706nm Background Subtraction')
    #--- poly plot
    ix = 704.8,709.4
    iy = f(704.8),f(709.4)
    verts = [(704.8,0), *zip(ix,iy), (709.4,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax706.add_patch(poly)
    
    ax714.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax714.set_xlim(713,718)
    ax714.set_ylim(0,0.06)
    ax714.set_title('706nm Background Subtraction')
    #--- poly plot
    ix = 713.8,717.5
    iy = f(717.9),f(717.9)
    verts = [(713.8,0), *zip(ix,iy), (717.5,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax714.add_patch(poly)
    
    ax727.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax727.set_xlim(724.5,730.5)
    ax727.set_ylim(0,0.06)
    ax727.set_title('706nm Background Subtraction')
    #--- poly plot
    ix = 725,730
    iy = f(725),f(730)
    verts = [(725,0), *zip(ix,iy), (730,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax727.add_patch(poly)
    
    ax738.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax738.set_xlim(733,746)
    ax738.set_ylim(0,0.06)
    ax738.set_title('706nm Background Subtraction')
    #--- poly plot
    ix = 736,742
    iy = f(733.5),f(733.5)
    verts = [(736,0), *zip(ix,iy), (742,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax738.add_patch(poly)
    
    ax794.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax794.set_xlim(790,800)
    ax794.set_ylim(0,0.06)
    ax794.set_title('706nm Background Subtraction')
    #--- poly plot
    ix = 791.2,797.5
    iy = f(791.2),f(791.2)
    verts = [(791.2,0), *zip(ix,iy), (797.5,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax794.add_patch(poly)
    
    plt.show
    
    fig2.tight_layout(h_pad=2)
    
    #Electron Temperature Background plots
    
    #Plotting the raw and interpolated data on the same figure
    #fig2, axes = plt.subplots(nrows=5,ncols=4)
    fig3, axes_3 = plt.subplots(nrows=2,ncols=2,sharey=True)
    fig3.suptitle('Te lines Background subtraction')
    
    #ax696,ax706,ax714,ax727,ax738,ax750,ax763,ax772,ax794,ax451,ax419,ax367,ax415,ax383,ax365,ax360,ax488,ax425,ax480,axblank = axes.flatten()
    ax750, ax763, ax772, axblank = axes_3.flatten()
    
    ax750.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax750.set_xlim(744,760)
    ax750.set_ylim(0,0.06)
    ax750.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 745,758
    iy = f(745),f(758)
    verts = [(745,0), *zip(ix,iy), (758,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax750.add_patch(poly)
    
    ax763.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax763.set_xlim(756,770)
    ax763.set_ylim(0,0.06)
    ax763.set_title('763nm Background Subtraction')
    #--- poly plot
    ix = 758,769.5
    iy = f(758),f(769.5)
    verts = [(758,0), *zip(ix,iy), (769.5,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax763.add_patch(poly)
    
    ax772.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax772.set_xlim(767,780)
    ax772.set_ylim(0,0.06)
    ax772.set_title('772nm Background Subtraction')
    #--- poly plot
    ix = 769.5,776
    iy = f(779.5),f(779.5)
    verts = [(769.5,0), *zip(ix,iy), (776,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax772.add_patch(poly)
        
    fig3.tight_layout()
    
    #Ne lines background plots
    
    #Plotting the raw and interpolated data on the same figure
    #fig2, axes = plt.subplots(nrows=5,ncols=4)
    fig4, axes_4 = plt.subplots(nrows=2,ncols=2,sharey=True)
    fig4.suptitle('Electron Density lines Background subtraction')
    
    #ax696,ax706,ax714,ax727,ax738,ax750,ax763,ax772,ax794,ax451,ax419,ax367,ax415,ax383,ax365,ax360,ax488,ax425,ax480,axblank = axes.flatten()
    ax451, ax419, ax415, ax383  = axes_4.flatten()
    #ax383, ax365, ax360, ax488, ax425, ax480
    
    ax451.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax451.set_xlim(449,455)
    ax451.set_ylim(0,0.01)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 449.6,453.1
    iy = f(449.6),f(453.1)
    verts = [(449.6,0), *zip(ix,iy), (453.1,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax451.add_patch(poly)
    
    ax419.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax419.set_xlim(416,422)
    ax419.set_ylim(0,0.06)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 417.1,421.2
    iy = f(417.1),f(421.2)
    verts = [(417.1,0), *zip(ix,iy), (421.2,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax419.add_patch(poly)
    
    # ax367.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    # ax367.set_xlim(365,369)
    # ax367.set_ylim(0,0.01)
    # #ax451.set_title('750nm Background Subtraction')
    # #--- poly plot
    # ix = 367,368.75
    # iy = f(367),f(368.75)
    # verts = [(367,0), *zip(ix,iy), (368.75,0)]
    # poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    # ax367.add_patch(poly)
        
    ax415.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax415.set_xlim(414,418)
    ax415.set_ylim(0,0.06)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 414.75,417.1
    iy = f(414.75),f(417.1)
    verts = [(414.75,0), *zip(ix,iy), (417.1,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax415.add_patch(poly)
    
    ax383.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax383.set_xlim(381,384.5)
    ax383.set_ylim(0,0.01)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 382.15,383.95
    iy = f(370.5),f(390.2)
    verts = [(382.15,0), *zip(ix,iy), (383.95,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax383.add_patch(poly)
    
    fig4.tight_layout()

    fig5, axes_5 = plt.subplots(nrows=2,ncols=2,sharey=True)
    fig4.suptitle('Electron Density lines Background subtraction')
    
    #ax696,ax706,ax714,ax727,ax738,ax750,ax763,ax772,ax794,ax451,ax419,ax367,ax415,ax383,ax365,ax360,ax488,ax425,ax480,axblank = axes.flatten()
    ax365, ax360, ax488, ax480  = axes_5.flatten()
    
    ax365.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax365.set_xlim(364,368)
    ax365.set_ylim(0,0.01)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 364.6,367
    iy = f(367),f(367)
    verts = [(364.6,0), *zip(ix,iy), (367,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax365.add_patch(poly)
    
    ax360.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax360.set_xlim(350,363)
    ax360.set_ylim(0,0.06)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 359,360.5
    iy = f(353),f(362)
    verts = [(359,0), *zip(ix,iy), (360.5,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax360.add_patch(poly)
    
    ax488.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax488.set_xlim(484,495)
    ax488.set_ylim(0,0.06)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 486,490
    iy = f(491.5),f(491.5)
    verts = [(486,0), *zip(ix,iy), (490,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax488.add_patch(poly)
    
    ax480.plot(lamda,f(lamda),color='orange',lw=0.75,label='Interpolated spectra',linestyle='dashed')
    ax480.set_xlim(475,485)
    ax480.set_ylim(0,0.01)
    #ax451.set_title('750nm Background Subtraction')
    #--- poly plot
    ix = 478.4,483.1
    iy = f(483.1),f(483.1)
    verts = [(478.4,0), *zip(ix,iy), (483.1,0)]
    poly = Polygon(verts,facecolor='0.9',edgecolor='0.5')
    ax480.add_patch(poly)
    
    return verts,fileN #for testing purposes swap fileN for df

#for testing purposes
# os.chdir(r'C:\Users\beau.bussell\Google Drive\EngD\Research Data\OES\Calibrated\181220')
# integral_result = integrate('RFPOWER0005')
# print(integral_result)











