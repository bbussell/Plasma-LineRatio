# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:06:52 2020

@author: beaub
"""

#Program to calculate the branching fractions of Ar emission lines such that the 
#metastable and resonant densities can be computed.

import time
import os

start_time = time.time()
from math import exp,sqrt
import numpy as np
from spec_integrate_Ne import integrate
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import csv
from datetime import datetime
#import RadTrap_No2 as RT # import rad trap code and run
from RadTrap_No2 import radtrap, print_escape
import pandas as pd


datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')

Kb = 1.38E-23*(1E4) #boltzman constant
#M = 9.109E-31
M = 39.948# 6.6335209E-23 #kg
#M = 6.6335209E-26 #atomic mass of Ar(kg)
R = 8.31446261815324*(1E4)
T_g = 600 #757 = 15mtorr #gas temperature (K)
p = 5 #charactersitic readsorption length (cm)

T_g_str = str(T_g)
p_str = str(p)
Tg = str(T_g)

p_str = str(p)

#!!!!!!!!!!!!!!!!!!!!!!!
#Before running, ensure all parameters are correct. Especially:
#Process Pressure, n_m, n_r

#List of datafiles to be used in electron density calculations.
filelist = ['LOCAL_RF_POS10006']
            # 'LOCAL_RF_POS10002',
            # 'LOCAL_RF_POS10003',
            # 'LOCAL_RF_POS10004',
            # 'LOCAL_RF_POS10005',
            # 'LOCAL_RF_POS10006',
            # 'LOCAL_RF_POS10007',
            # 'LOCAL_RF_POS10008',
            # 'LOCAL_RF_POS10009',
            # 'LOCAL_RF_POS10010',
            # 'LOCAL_RF_POS10011',
            # 'LOCAL_RF_POS10012']

#Transition Probabilties and reabsorption coeffcients for 6 Ar lines
#units for kij0 = cm^2K^0.5
#units for Aij = s^-1

#Fetching NIST data

#696 j=1s5
kij_696 = 1.43E-11
Aij_696 = 6.39E6

#727 j=1s4
kij_727 = 7.74E-12
Aij_727 = 1.83E6

#738 j=1s4
kij_738 = 6.25E-11
Aij_738 = 8.47E6

#706 j=1s5
kij_706 = 1.47E-11
Aij_706 = 3.8E6

#794 j=1s3
kij_794 = 3.07E-10
Aij_794 = 18.6E6

#714 j=1s5
kij_714 = 1.51E-12
Aij_714 = 0.63E6

#Transition Probabilties and reabsorption coeffcients for 6 Ar lines for Te calculation
#units for kij0 = cm^2K^0.5
#units for Aij = s^-1

#k = np.array([(2.20E-10),6.09E-9,8,23E-8],[1.10E-9,1.91E-8,4.62E-8],
           #  [2.08E-10,8.83E-8,1.03E-7],[1.04E-10,1.03E-8,3.85E-8],
          #   [2.39E-10,3.93E-8,5.15E-8])


#738 j=1s4 i=2p3
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

#772
kG_772 = 1.04E-10
alphaG_772 = 0.69
EG_772 = 14.49

kM_772 = 1.03E-8
alphaM_772 = -0.08
EM_772 = 1.62

kR_772 = 3.85E-8
alphaR_772 = 0.01
ER_772 = 2.02

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

#-------------------------------------------------------
with open('Density_inputs.txt', 'r') as inputs:

#Extract data and declare variables
    n1s4, n1s5 = np.loadtxt("Density_inputs.txt", delimiter=';', 
                            unpack=True, usecols=(0, 1))

#test
# n1s4 =[1,2,3,4,5]
# n1s5 = [10,20,30,40,50]

#exec(open('astro_nist.py').read())

param = input("What system parameter was used during this experimental work? is this assessment for? E.g. 2kW or 0.0050 PP. Please respond and press Enter: ")

print("Thank you. The parameter you are about to evaluate is: ",param)

def print_results(a,b,c):
    
    c_list = list(c)
    a_list = list(a)
    b_list = list(b)

    min_pos = c_list.index(min(c_list))

    print('Minimum chi-squared of ', min(c), 'occurs at position ', min_pos, 
          'where n1s4 = ',"%7.1e" % a_list[min_pos], 'cm\u00b3', 'and n1s5 = ',
          "%7.1e" % b_list[min_pos], 'cm\u00b3') 
    
    final_density = [a_list[min_pos],b_list[min_pos]]
    
    return final_density
    
#Function for caluclation of reabsorption coefficients 

def reabs_coeff(k,n):
    return k*n*(T_g**-0.5)

#Function for calculation of escape factors

def escape_factor(k,n,p):
    return (2-exp(-reabs_coeff(k,n)*p/1000))/(1+(reabs_coeff(k,n)*p))

#Function for calculation of Model Line Ratios
                
def LR_model(kij,kik,Aij,Aik,nij,nik,p):
    return (Aij*escape_factor(kij,nij,p))/(Aik*escape_factor(kik,nik,p))

#Function for calculation of experimental line ratios

def LR_exp(Iij,Iik):
    return Iij/Iik

#Calculation of chi-squared term without summing

def chi_squared(LRm,LRe,err):
    return ((LRe-LRm)/(err*LRe))**2

def raw(textfile):
    I = np.loadtxt("RawInputData/" + textfile + ".txt", delimiter=' ', 
                  unpack=True, usecols=(1), skiprows=1) 
    
    return I
#Te functions

#Function for caluclation of Excitation rates
#k=k_o, T=electron temperature, a=

def ex_rates(k,T,a,E):
    return k*((T)**a)*exp(-E/T)
#test of excitation rates

#def test_answer():
 #   assert ex_rates(kG_738,3.5,alphaG_738,EG_738) == 6.52E-12

def sum_nk(k,T,a,E,km,am,Em,kr,ar,Er,n_m,n_r):
    return n_g*ex_rates(k,T,a,E) + n_m*ex_rates(km,T,am,Em) +  n_r*ex_rates(kr,T,ar,Er)

def k_o(lamda,gi,gj,Aij):  
    term1 = ((lamda*(1E-7))**3)/(8*(np.pi**(1.5)))
    #print("term 1 =",term1)
    term2 = gi/gj
    #print("term 2 =", term2)    
    term3 = Aij*(sqrt(M))/((sqrt(2)*(sqrt(R))))
    #print("term 3 =", term3)    
    result = term1*term2*term3    
    #arr_k = np.array(result)    
    return result

def sum_750(T):
    sum_750=sum_nk(kG_750,T,alphaG_750,EG_750,kM_750,alphaM_750,EM_750,kR_750,alphaR_750,ER_750,n_m,n_r)
    return sum_750

def LR_mod(k,T,a,E,km,am,Em,kr,ar,Er,Rij,Rmn,n_m,n_r):
    LR = (Rij/Rmn)*((sum_nk(k,T,a,E,km,am,Em,kr,ar,Er,n_m,n_r)/sum_750(T)))
    return LR

def LR(n1s4,n1s5):
    
        for a, b in zip(n1s4,n1s5):
        
                #-----------------------------------
                #Experimental Results from J Boffard
                #1mtorr
                
                n_1s5 = b #metastable density (cm^-3)
                n_1s4 = a #resonant density (cm^-3)
                
                #-----------------------------------
                #Relation between total m and r 1sx levels
                n_1s3 = n_1s5/6.5
                n_1s2 = n_1s4
                    
                #Setting model LR's
                
                #print('The 696/727 model LR is:')
                #print(LR_model(kij_696,kij_727,Aij_696,Aij_727,n_1s5,n_1s4,p))
                
                LRm1=LR_model(kij_696,kij_727,Aij_696,Aij_727,n_1s5,n_1s4,p)
                
                #print('The 738/706 model LR is:')
                #print(LR_model(kij_738,kij_706,Aij_738,Aij_706,n_1s4,n_1s5,p))
                
                LRm2=LR_model(kij_738,kij_706,Aij_738,Aij_706,n_1s4,n_1s5,p)
                
                #print('The 794/714 model LR is:')
                #print(LR_model(kij_794,kij_714,Aij_794,Aij_714,n_1s3,n_1s5,p))
                
                LRm3=LR_model(kij_794,kij_714,Aij_794,Aij_714,n_1s3,n_1s5,p)
                
                #Setting experimental LR's
                
                #print('The 696/727 exp LR is:')
                #print(LR_exp(I_696,I_727))
                LRe1=LR_exp(I_696,I_727)
                
                #print('The 738/706 exp LR is:')
                #print(LR_exp(I_738,I_706))
                LRe2=LR_exp(I_738,I_706)
                
                #print('The 794/714 exp LR is:')
                #print(LR_exp(I_794,I_714))
                LRe3=LR_exp(I_794,I_714)
                
                #Summing all the chi-squared terms for each LR
                
                chi_sum=chi_squared(LRm1,LRe1,0.05)+chi_squared(LRm2,LRe2,0.05)+chi_squared(LRm3,LRe3,0.05)
                #print('Chi squared for n_r= ',n_1s4, 'and n_m= ',n_1s5, 'is: ',chi_sum)
                
                with open('mod_results.txt', 'a+') as mod_results:
                #for i in range(len(n1s4)):
                 #   mod_results.write("%d %d %d\n" % (n1s4[i],n1s5[i],chi_sum))
                    mod_results.write("%7.2e %7.2e %6.2f\n" % (a,b,chi_sum))
                    
                    #printing results to screen
        with open('mod_results.txt', 'a+') as mod_results:
            d, e, f = np.loadtxt("mod_results.txt",unpack=True,usecols=(0, 1, 2),skiprows=7)
        
            final_density = print_results(d,e,f) #calling function to print results
        
        #plotting chi-sum as 3D plot
        
        #fig3d = plt.figure()
        ax1 = plt.axes(projection="3d")
        
        ax1.plot3D(d,e,f,'black')
        plt.show()
        
        #renaming the file with the current date and time
        os.rename("mod_results.txt", time.strftime("MetaResResults/n1sDEN_"+i+param+"_%Y%m%d%H%M.txt")) 

        return final_density
def Te(T,n_m,n_r):
    #n_m = 3.0E10
    #n_r = 9.3E9
    #Ask user for n_m and n_r inputs
         
    I_738 = I[4]
    I_750 = I[5]
    I_763 = I[6]
    I_772 = I[7]
    I_794 = I[8]
    
    Exp_738 = I_738/I_750
    Exp_763 = I_763/I_750
    Exp_772 = I_772/I_750
    Exp_794 = I_794/I_750
    
    #LR's are:
    # 738/750 ; 763/750 ; 772/750 ; 794/750
    
    n_ij, A, g_i, g_j = np.loadtxt("line_data2.txt", comments='#', delimiter=';', skiprows=2, unpack=True, 
                                        usecols=(0,1,2,3))
    
    n_i = np.loadtxt("line_data2.txt", comments='#',dtype=str, delimiter=';', skiprows=2, unpack=True, 
                                        usecols=(4))

    #Taking n_m and n_r from calculation 
    n_1s3 = n_m/6.5
    n_1s2 = n_r
    
    for a, b in zip(n_ij,n_i):
        print("Wavelength: ", a)
        print("Lower level: ", b)
        if b=="1s2":
            N_j = n_1s2
        elif b=="1s3":
            N_j = n_1s3
        elif b=="1s4":
            N_j = n_r
        elif b=="1s5":
            N_j = n_m
        
        print("Density of ", b, "is ", N_j)
            
        with open('live_density.txt', 'a+') as resultsfile:
            resultsfile.write("%4.2f %s %4.2e \n" % (a,b,N_j))
    
    n_ij, n_j = np.loadtxt("live_density.txt", comments='#', delimiter=' ', unpack=True, 
                                        usecols=(0,2))  
    
    
    radtrap(n_ij,A,g_i,g_j,n_j)
    print("")
    print("radtrap complete")
    print("")
    
    n_ij, g_i, g_j, A, n_j, ko, k_ij = np.loadtxt("line_data_full.txt", comments='#', delimiter=' ', unpack=True, 
                                   usecols=(0,1,2,3,4,5,6))   
    print_escape(n_ij,g_i,g_j,A,n_j,ko,k_ij)
    print("")
    print("escape complete")
    print("")
    
    #print("The Te values being modelled are:")
    #print(T)
    
    with open('Te_results.txt', 'w') as resultsfile:
        resultsfile.write('Electron Temperature Results\n')
        resultsfile.write('Datetime (Y-M-S) = ' + datestring + '\n')
        resultsfile.write('Filename: ' + i + '\n')
        resultsfile.write('Experiment Name: ' + param + '\n')
        resultsfile.write('Gas Temperature (K) = ' + Tg + '\n')
        resultsfile.write('Characteristic length (cm) = '+ p_str + '\n')
        
    with open('LR_results.txt', 'w') as resultsfile:
        resultsfile.write('Electron Temperature Results - Raw Line Ratio Data\n')
        resultsfile.write('Datetime (Y-M-S) = ' + datestring + '\n')
        resultsfile.write('Filename: ' + i + '\n')
        resultsfile.write('Experiment Name: ' + param + '\n')
        resultsfile.write('Gas Temperature (K) = ' + Tg + '\n')
        resultsfile.write('Characteristic length (cm) = '+ p_str + '\n')
    
    for a in T:
                
        #-----------------------------------   
        #print('For an electron temperature of Te=',a)
        #print('')
        R_lam, Rad = np.genfromtxt("Rad_TrapCoeff.csv", delimiter=';', unpack=True, skip_header=1, usecols=(0,1))
        
        R_750 = Rad[0]
        R_772_4 = Rad[1]
        R_738 = Rad[2]
        R_794 = Rad[3]
        R_751 = Rad[4]
        R_763 = Rad[5]
        R_772_3 = Rad[6]
        
        #738/750
        LR_738 = LR_mod(kG_738,a,alphaG_738,EG_738,kM_738,alphaM_738,EM_738,kR_738,alphaR_738,ER_738,R_738,R_750,n_m,n_r)
        #print("The 738/750 line ratio is ", LR_738, "for Te =", a)
        #print("-------------------------")
        LR_763 = LR_mod(kG_763,a,alphaG_763,EG_763,kM_763,alphaM_763,EM_763,kR_763,alphaR_763,ER_763,R_763,R_750,n_m,n_r)
        #print("The 763/750 line ratio is, ", LR_763, "for Te =", a)
        
            #print("-------------------------")
    
        LR_772 = LR_mod(kG_772,a,alphaG_772,EG_772,kM_772,alphaM_772,EM_772,kR_772,alphaR_772,ER_772,R_772_3,R_750,n_m,n_r)
        #print("The 772/750 line ratio is, ", LR_772, "for Te =", a)
        
        #print("-------------------------")
        
        LR_794 = LR_mod(kG_795,a,alphaG_795,EG_795,kM_795,alphaM_795,EM_795,kR_795,alphaR_795,ER_795,R_794,R_750,n_m,n_r)
        #print("The 794/750 line ratio is, ", LR_794, "for Te =", a)
        
        chi_sum = chi_squared(LR_738,Exp_738,0.05) + chi_squared(LR_763,Exp_763,0.05) + chi_squared(LR_772,Exp_772,0.1) + chi_squared(LR_794,Exp_794,0.1)
        # print("")
        # print("chi_sum for T =", a,"is: ",chi_sum)
        # print("738 model LR is: ", LR_738, "and experimental 738 LR is: ", Exp_738)
        chi_738 = chi_squared(LR_738,Exp_738,0.05)
        
        #print("763 model LR is: ", LR_763, "and experimental 763 LR is: ", Exp_763)
        chi_763 = chi_squared(LR_763,Exp_763,0.05)
        
        #print("772 model LR is: ", LR_772, "and experimental 772 LR is: ", Exp_772)
        chi_772 = chi_squared(LR_772,Exp_772,0.1)
        
        #print("794 model LR is: ", LR_794, "and experimental 794 LR is: ", Exp_794)
        chi_794 = chi_squared(LR_794,Exp_794,0.1)
        
        #print("738 chi is: ", chi_738)
        # print("763 chi is: ", chi_763)
        # print("772 chi is: ", chi_772)
        # print("794 chi is: ", chi_794)
        
        #Adding calculation of chi-sum without including 763nm line
        
        chi_excl = chi_738 + chi_772 + chi_794
        
        chi_s_e_738 = chi_763 + chi_772 + chi_794
        
        chi_s_e_both = chi_772 + chi_794
        
        with open('Te_results.txt', 'a+') as te_results:
            #for i in range(len(n1s4)):
             #   mod_results.write("%d %d %d\n" % (n1s4[i],n1s5[i],chi_sum))
            te_results.write("%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n" % (a,chi_738, chi_763, chi_772, chi_794, chi_sum, chi_excl, chi_s_e_738, chi_s_e_both))
            
        with open("LR_results.txt", 'a+') as lr_results:
            lr_results.write("%4.2f %5.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n" % (a,LR_738,Exp_738,LR_763,Exp_763,LR_772,Exp_772,LR_794,Exp_794))
        
    Chi_df=pd.read_csv("Te_results.txt",sep=" ",names=['Electron Temperature (eV)','738 Chi','763 Chi','772 Chi','794 Chi','Chi Sum', 'Chi Sum excl. 763nm', 'Chi Sum excl. 738nm', 'Chi Sum excl 738 and 763nm'],comment="#",skiprows=6)
    Chi_df.to_csv('Te_results.csv',sep=';',index=False)
    
    LR_df=pd.read_csv("LR_results.txt",sep=" ",header=None,names=['Electron Temperature (eV)','738 Model LR','738 Exp LR','763 Model LR','763 Exp LR','772 Model LR','772 Exp LR','794 Model LR','794 Exp LR'],comment='#',skiprows=6)
    LR_df.to_csv('LR_results.csv',sep=";",index=False)
                
    Te, chi_738, chi_763, chi_772, chi_794, chi_s, chi_s_excl, chi_s_excl_738, chi_s_excl_both = np.genfromtxt("Te_results.csv", delimiter=";", skip_header=1, unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
    
    Te_l_lst = list(Te)
    chi_s_lst = list(chi_s)
    
    chi_s_e_lst = list(chi_s_excl)
    chi_s_e_738_lst = list(chi_s_excl_738)
    chi_s_e_both_lst = list(chi_s_excl_both)
    chi_763_lst = list(chi_763)
    
    min_pos = chi_s_lst.index(min(chi_s_lst))
    
    min_pos_excl = chi_s_e_lst.index(min(chi_s_e_lst))
    
    min_pos_738 = chi_s_e_738_lst.index(min(chi_s_e_738_lst))
    
    min_pos_both = chi_s_e_both_lst.index(min(chi_s_e_both_lst))
    
    min_763 = chi_763_lst.index(min(chi_763_lst))
    
    
    #plotting chi-sum 
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    ax1.plot(Te, chi_738, c='b', label='chi-738')
    ax2 = fig.add_subplot(111)
    
    ax2.plot(Te,chi_763, c='r', label='chi-763')
    
    ax3 = fig.add_subplot(111)
    ax3.plot(Te,chi_772, c='g',label='chi-772')
    
    ax4 = fig.add_subplot(111)
    ax4.plot(Te,chi_794, c='y',label='chi-794')
    
    ax5 = fig.add_subplot(111)
    ax5.plot(Te,chi_s, c='magenta',label='chi-sum')
    
    ax6 = fig.add_subplot(111)
    ax6.plot(Te,chi_s_excl, c='orange',label='chi-sum excl 763')
    
    ax7 = fig.add_subplot(111)
    ax7.plot(Te,chi_s_excl_738, c='purple',label='chi-sum excl 738')
    
    ax8 = fig.add_subplot(111)
    ax8.plot(Te,chi_s_excl_both, c='black', label='chi-sum excl both')
    
    plt.legend(loc='upper right');
    
    ax1.set_ylim([0,200])
    ax2.set_ylim([0,200])
    ax3.set_ylim([0,200])
    ax4.set_ylim([0,200])
    ax5.set_ylim([0,200])
    ax6.set_ylim([0,200])
    plt.show()
    
    plt.show()
    
    print('------------------------------------------------------------------')
    
    print("The minimum in chi-squared was", chi_s_lst[min_pos], "occurring when Te=", Te[min_pos])
    
    print("The minimum in chi-squared, excluding 763 is:", chi_s_e_lst[min_pos_excl], "occuring when Te=", Te[min_pos_excl])
    
    print("The minimum in chi-squared, excluding 738 is:", chi_s_e_738_lst[min_pos_738], "occuring when Te=", Te[min_pos_738])
    
    print("The minimum in chi-squared, excluding both 763 and 738 is:", chi_s_e_both_lst[min_pos_both], "occuring when Te=", Te[min_pos_both])
    
    print("The minimum in 763 chi is:", chi_763_lst[min_763], "occuring when Te=", Te[min_763])
    
    os.rename("Te_results.txt", time.strftime("TeResults/te_results"+i+param+Tg+"K_%Y%m%d%H%M%S.txt")) 
    os.rename("Te_results.csv", time.strftime("TeResults/te_results"+i+param+Tg+"_%Y%m%d%H%M%S.csv"))
    
    os.rename("LR_results.txt", time.strftime("TeResults/LR_results"+i+param+Tg+"_%Y%m%d%H%M%S.txt")) 
    os.rename("LR_results.csv", time.strftime("TeResults/LR_results"+i+param+Tg+"_%Y%m%d%H%M%S.csv"))
    
    #os.remove("te_results.csv")
    
    final_T = Te[min_pos_both]
    
    return final_T

def e_density(T,i):
    
    nec_750 = 3E12
    nec_451 = 4E11
    nec_419 = 3.2E11 #calculated from reciprocal of transition probability
    nec_367 = 1E11 #estimate based on Zhu 2007
    nec_415 = 5E11
    nec_383 = 1E11 #Zhu 2007
    nec_365 = 1E11
    nec_360 = 1E11
    
    K = 567 - (160*T)+(16*(T**2))
    
    #filename = input("What is the name of the file you want to analyse?")
    #file = filename+'_IntegratedIntensity.txt'
    #change directory 
    #os.chdir(r'..\..\OES\Calibrated\081020')
    
    #os.chdir(r'..\..\..\Plasma Spec Code\Plasma-LineRatio')
    
    I_750 = I[5]
    I_451 = I[9]
    I_419 = I[10]
    I_367 = I[11]
    I_415 = I[12]
    I_383 = I[13]
    I_365 = I[14]
    I_360 = I[15]
    
    # n_e = (1-((I_750/I_365)/(K)))/(((I_750/I_383)/(K*nec_750))-(1/nec_383))
    n_e = (1-((I_750/I_383)/K)) / (((I_750/I_383)/(K*nec_750))-(1/nec_383))
    N_e = abs(n_e)
    N_e_str = str(N_e)
    print("the electron density for", i, " is: ", "%6.2e" % N_e, "cm-3")
    
    with open(i+'383_E_density.txt', 'w+') as resultsfile:
        resultsfile.write('Datafile: '+i+'\n')
        resultsfile.write('Experiment Name: '+param+'\n')
        resultsfile.write('Electron density (cm^-3) \n')
        resultsfile.write(N_e_str)
    
    
    return N_e

#------------------------------------------------------
#Extract raw Intensity values for parameter chosen. 

#Integrating the spectral data to measure the intensity of each peak used

#filename = input("What is the name of the file you want to analyse?")
    
for i in filelist:
        # PP_mbar = float(input("What is the process pressure in mbar?"))
        # #PP_mbar = 0.0020 #argon partial pressure in mbar
        # PP_pa = PP_mbar*100 #argon partial pressure in pascals
        # n_g = PP_pa/(Kb*T_g)#6E13 
        
        #New n_g
        #Assuming that p=0.0050mbar, which is 3.75mtorr, then at 0.5kW, n_g = 6E13
        N_g = float(input("What is the ground state density, in x10^13 cm^-3"))
        n_g = N_g*(1E13)
        
        #Ask user for electon temperature value to be used in electron density calculation
        #calling integrate function to determine intensity of 451 and 750 emission lines
        os.chdir(r'C:\Users\beaub\Google Drive\EngD\Research Data\OES\Calibrated\061020')
        #integrating spectra in file i
        integrate(i)

        lamda, I = np.loadtxt(i+'_IntegratedIntensity.txt', comments='#', delimiter=',', skiprows=2, unpack=True, 
                                            usecols=(0,1))
        I_corr = I
        
        os.chdir(r'..\..\..\Plasma Spec Code\Plasma-LineRatio')
        
        #setting intensity values for met-res calculation
        #I = raw(filename)
        I_696 = I[0]
        I_706 = I[1]
        I_714 = I[2]
        I_727 = I[3]
        I_738 = I[4]
        I_794 = I[8]        
    
        #Creating headers in final results file
        with open('mod_results.txt', 'w') as resultsfile:
            resultsfile.write('Metastable and Resonant Density Model Results\n')
            resultsfile.write('Datetime (Y-M-S) = ' + datestring + '\n')
            resultsfile.write('Filename:' + i + '\n')
            resultsfile.write('Experiment Name: ' + param + '\n')
            resultsfile.write('Gas Temperature (K) = ' + T_g_str + '\n')
            resultsfile.write('Characteristic length (cm) = '+ p_str + '\n')
            resultsfile.write('n1s4 ns15 Chi-Squared \n')
        #Executing function to calculate n_m and n_r
        final_density = LR(n1s4,n1s5)
        
        print("--------------------------------------------------------------")
        print("")
        print("You are now calculating electron temperature")
        print("")
        
        #------------------
        #Te Calculation
               
        #-------------------------------------------------------
        #Loading Te model data 
        
        #Extract data and declare variables
        n_m = final_density[1]
        n_r = final_density[0]
        
        T = np.loadtxt("Te_Intervals.txt", unpack=True,
                              usecols=(0))
        #testing T values
        #T = [3.5,8]
        final_T = Te(T,n_m,n_r)
        
        print("--------------------------------------------------------------")
        print("")
        print("You are now calculating electron density")
        print("")
        
        e_density(final_T,i)


#----------------------------------------------
#END OF PROGRAM
#Calculating time program took to run
final_time = time.time() - start_time

print("This program took", "%5.3f" %  final_time,"s to run")




                
        
          
    
        
    
        
    
    
    
    
    
        
    
    
    
    
    



                
        
          
    
        
    
        
    
    
    
    
    
        
    
    
    
    
    
