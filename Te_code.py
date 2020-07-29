# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 17:26:05 2020

@author: beaub
"""
import numpy as np
from math import exp
from math import sqrt
import pandas as pd
import csv
import os
import time
import RadTrap_No2 as RT # import rad trap code and run

from datetime import datetime
start_time = time.time()
datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')

param = "1.5kW"


Kb = 1.38E-23*(1E4) #boltzman constant
T_g = 600 #757 = 15mtorr #gas temperature (K)
p = 10 #charactersitic readsorption length (cm)
#M = 9.109E-31
M = 39.948# 6.6335209E-23 #kg
#M = 6.6335209E-26 #atomic mass of Ar(kg)
R = 8.31446261815324*(1E4)

Tg = str(T_g)


PP_mbar = 0.0050 #argon partial pressure in mbar
PP_pa = PP_mbar*100 #argon partial pressure in pascals
n_g = 6E13 #PP_pa/(Kb*T_g)
#n_m = 3.0E10
#n_r = 9.3E9
n_m = 3.9E10 #n1s5
n_r = 1.2E10 #n1s4

#0.5kW RFPowerTopView004

#I_360 = 13.61 for n_e determination
#I_383 = 25.42 for n_e determination
#I_738 = 4051.31
#I_763 = 7209.19 
#I_750 = 8980.17  
#I_772 = 2532.46
#I_794 = 1979.53 #I_794 = 4595.85  

#3.0kW irradiance (integrated irradiance calculated from spectrum analyser)
#Note, checked these twice to ensure they were correct
# I_738 = 221.63
# I_763 = 622.25 
# I_750 = 599.88  
# I_772 = 277.20 #(772.38)
# I_794 = 310.63  

#1.5kW irradiance (integrated)
I_738 = 199.16
I_763 = 484.74 
I_750 = 505.74
I_772 = 240.35
I_794 = 249.30

#0.5kW irradiance (integrated)
# I_738 = 98.36
# I_763 = 234.89
# I_750 = 264.55
# I_772 = 111.54
# I_794 = 108.44

Exp_738 = I_738/I_750
Exp_763 = I_763/I_750
Exp_772 = I_772/I_750
Exp_794 = I_794/I_750

#LR's are:
# 738/750 ; 763/750 ; 772/750 ; 794/750

#Transition Probabilties and reabsorption coeffcients for 6 Ar lines
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
        
#Function for caluclation of Excitation rates
#k=k_o, T=electron temperature, a=

def ex_rates(k,T,a,E):
    return k*((T)**a)*exp(-E/T)
#test of excitation rates

#def test_answer():
 #   assert ex_rates(kG_738,3.5,alphaG_738,EG_738) == 6.52E-12

def sum_nk(k,T,a,E,km,am,Em,kr,ar,Er):
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

#Function for caluclation of reabsorption coefficients
#For use in radiation trapping term Rij

def reabs_coeff(k,n):
    return k*n*(T_g**-0.5)

#Function for calculation of escape factors

def escape_factor(k,n,p):
    return (2-exp(-reabs_coeff(k,n)*p/1000))/(1+(reabs_coeff(k,n)*p))
        
#Function for calculation of Model Line Ratios

n_ij, A, g_i, g_j, n_j = np.loadtxt("line_data2.txt", comments='#', delimiter=';', skiprows=2, unpack=True, 
                                        usecols=(0,1,2,3,4))

RT.radtrap(n_ij,A,g_i,g_j,n_j)

def sum_750(T):
    sum_750=sum_nk(kG_750,T,alphaG_750,EG_750,kM_750,alphaM_750,EM_750,kR_750,alphaR_750,ER_750)
    return sum_750

def LR_mod(k,T,a,E,km,am,Em,kr,ar,Er,Rij,Rmn):
    LR = (Rij/Rmn)*((sum_nk(k,T,a,E,km,am,Em,kr,ar,Er)/sum_750(T)))
    return LR
    
#-------------------------------------------------------
    
#-------------------------------------------------------

#Function for calculation of experimental line ratios

def LR_exp(Iij,Iik):
    return Iij/Iik

#Calculation of chi-squared term without summing

def chi_squared(LRm,LRe,err):
    return ((LRe-LRm)/(err*LRe))**2

#-------------------------------------------------------
#Loading Te model data 

#Extract data and declare variables
T = np.loadtxt("Te_Intervals.txt", unpack=True,
                      usecols=(0))
#testing T values
#T = [3.5,8]
print("The Te values being modelled are:")
print(T)

print("The line ratios are:")
#print(n1s4)
#print(n1s5)

#test
#n1s4 =[1,2,3,4,5]
#n1s5 = [10,20,30,40,50]

for a in T:
    #print("The ground, metastable and resonant excitation rate at 738nm, k, for T=",a," are:")
    #print(ex_rates(kG_738,a,alphaG_738,EG_738))
    #print(ex_rates(kM_738,a,alphaM_738,EM_738))
    #print(ex_rates(kR_738,a,alphaR_738,ER_738))
    
    #print("The sum of the excitation rates for 738nm is:")
    #print("%8.2e" % sum_nk(kG_738,a,alphaG_738,EG_738,kM_738,alphaM_738,EM_738,kR_738,alphaR_738,ER_738))
    
    #-----------------------------------
    
    #print("The ground, metastable and resonant excitation rates at 763nm, k, for T=",a," are:")
    #print(ex_rates(kG_763,a,alphaG_763,EG_763))
    #print(ex_rates(kM_763,a,alphaM_763,EM_763))
    #print(ex_rates(kR_763,a,alphaR_763,EM_763))
    
    #print("The sum of the excitation rates for 763nm is:")
    #print("%8.2e" % sum_nk(kG_763,a,alphaG_763,EG_763,kM_763,alphaM_763,EM_763,kR_763,alphaR_763,ER_763))
    
    #-------------------------------------
  
    #print("The ground, metastable and resonant excitation rate at 772nm, k, for T=",a," are:")
    #print(ex_rates(kG_772,a,alphaG_772,EG_772))
    #print(ex_rates(kM_772,a,alphaM_772,EM_772))
    #print(ex_rates(kR_772,a,alphaR_772,EM_772))
  
  
    #print("The sum of the excitation rates for 772nm is:")
    #print("%8.2e" % sum_nk(kG_772,a,alphaG_772,EG_772,kM_772,alphaM_772,EM_772,kR_772,alphaR_772,ER_772))
    
    #-------------------------------------
    
    #print("The ground, metastable and resonant excitation rate at 795nm, k, for T=",a," are:")
    #print(ex_rates(kG_795,a,alphaG_795,EG_795))
    #print(ex_rates(kM_795,a,alphaM_795,EM_795))
    #print(ex_rates(kR_795,a,alphaR_795,EM_795))
    
    #print("The sum of the excitation rates for 795nm is:")
    #print("%8.2e" % sum_nk(kG_795,a,alphaG_795,EG_795,kM_795,alphaM_795,EM_795,kR_795,alphaR_795,ER_795))
    
    #-----------------------------------   
    print('For an electron temperature of Te=',a)
    print('')
    R_lam, Rad = np.genfromtxt("Rad_TrapCoeff.csv", delimiter=';', unpack=True, skip_header=1, usecols=(0,1))
    
    R_750 = Rad[0]
    R_772_4 = Rad[1]
    R_738 = Rad[2]
    R_794 = Rad[3]
    R_751 = Rad[4]
    R_763 = Rad[5]
    R_772_3 = Rad[6]
    
    #738/750
    
    LR_738 = LR_mod(kG_738,a,alphaG_738,EG_738,kM_738,alphaM_738,EM_738,kR_738,alphaR_738,ER_738,R_738,R_750)
    #print("The 738/750 line ratio is ", LR_738, "for Te =", a)
    
    print("-------------------------")
    
    LR_763 = LR_mod(kG_763,a,alphaG_763,EG_763,kM_763,alphaM_763,EM_763,kR_763,alphaR_763,ER_763,R_763,R_750)
    #print("The 763/750 line ratio is, ", LR_763, "for Te =", a)
    
    
    #print("-------------------------")

    LR_772 = LR_mod(kG_772,a,alphaG_772,EG_772,kM_772,alphaM_772,EM_772,kR_772,alphaR_772,ER_772,R_772_3,R_750)
    #print("The 772/750 line ratio is, ", LR_772, "for Te =", a)
    
    #print("-------------------------")
    
    LR_794 = LR_mod(kG_795,a,alphaG_795,EG_795,kM_795,alphaM_795,EM_795,kR_795,alphaR_795,ER_795,R_794,R_750)
    #print("The 794/750 line ratio is, ", LR_794, "for Te =", a)
    
    chi_sum = chi_squared(LR_738,Exp_738,0.05) + chi_squared(LR_763,Exp_763,0.05) + chi_squared(LR_772,Exp_772,0.1) + chi_squared(LR_794,Exp_794,0.1)
    print("")
    print("chi_sum for T =", a,"is: ",chi_sum)
    print("738 model LR is: ", LR_738, "and experimental 738 LR is: ", Exp_738)
    chi_738 = chi_squared(LR_738,Exp_738,0.05)
    
    print("763 model LR is: ", LR_763, "and experimental 763 LR is: ", Exp_763)
    chi_763 = chi_squared(LR_763,Exp_763,0.05)
    
    print("772 model LR is: ", LR_772, "and experimental 772 LR is: ", Exp_772)
    chi_772 = chi_squared(LR_772,Exp_772,0.1)
    
    print("794 model LR is: ", LR_794, "and experimental 794 LR is: ", Exp_794)
    chi_794 = chi_squared(LR_794,Exp_794,0.1)
    
    print("738 chi is: ", chi_738)
    print("763 chi is: ", chi_763)
    print("772 chi is: ", chi_772)
    print("794 chi is: ", chi_794)
    
    with open('Te_results.txt', 'a+') as te_results:
        #for i in range(len(n1s4)):
         #   mod_results.write("%d %d %d\n" % (n1s4[i],n1s5[i],chi_sum))
        te_results.write("%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n" % (a,chi_738, chi_763, chi_772, chi_794, chi_sum))
        
    with open("LR_results.txt", 'a+') as lr_results:
        lr_results.write("%4.2f %5.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n" % (a,LR_738,Exp_738,LR_763,Exp_763,LR_772,Exp_772,LR_794,Exp_794))
            
Chi_df=pd.read_csv("te_results.txt",sep=" ",header=None,names=['Electron Temperature (eV)','738 Chi','763 Chi','772 Chi','794 Chi','Chi Sum'],comment="#")
Chi_df.to_csv('te_results.csv',sep=';',index=False)

LR_df=pd.read_csv("LR_results.txt",sep=" ",header=None,names=['Electron Temperature (eV)','738 Model LR','738 Exp LR','763 Model LR','763 Exp LR','772 Model LR','772 Exp LR','794 Model LR','794 Exp LR'],comment='#')
LR_df.to_csv('LR_results.csv',sep=";",index=False)
            
Te, chi_s = np.genfromtxt("te_results.csv", delimiter=";", skip_header=1, unpack=True, usecols=(0,5))

Te_l_lst = list(Te)
chi_s_lst = list(chi_s)

min_pos = chi_s_lst.index(min(chi_s_lst))

print("The minimum in chi-squared was", chi_s_lst[min_pos], "occurring when Te=", Te[min_pos])

#Calculating time program took to run
final_time = time.time() - start_time

print("This program took", "%5.3f" %  final_time,"s to run")

os.rename("te_results.txt", time.strftime("te_results"+param+Tg+"K_%Y%m%d%H%M%S.txt")) 
os.rename("te_results.csv", time.strftime("te_results"+param+Tg+"_%Y%m%d%H%M%S.csv"))

os.rename("LR_results.txt", time.strftime("LR_results"+param+Tg+"_%Y%m%d%H%M%S.txt")) 
os.rename("LR_results.csv", time.strftime("LR_results"+param+Tg+"_%Y%m%d%H%M%S.csv"))

#os.remove("te_results.csv")
                
        
          
    
        
    
        
    
    
    
    
    
        
    
    
    
    
    
