# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 17:26:05 2020

@author: beaub
"""
#Program to calculate the branching fractions of Ar emission lines such that the 
#metastable and resonant densities can be computed.

import time
import os

start_time = time.time()

from math import exp    
import numpy as np

from datetime import datetime
datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')

param = input("What system parameter was used during this experimental work? is this assessment for? E.g. 2kW or 0.0050 PP. Please respond and press Enter: ")

print("Thank you. The RF Power you are about to evaluate is: ",param)

def print_results(a,b,c):
    
    c_list = list(c)
    a_list = list(a)
    b_list = list(b)

    min_pos = c_list.index(min(c_list))

    print('Minimum chi-squared of ', min(c), 'occurs at position ', min_pos, 
          'where n1s4 = ',"%7.1e" % a_list[min_pos], 'cm\u00b3', 'and n1s5 = ',
          "%7.1e" % b_list[min_pos], 'cm\u00b3') 
    
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

def chi_squared(LRm,LRe):
    return ((LRe-LRm)/(0.05*LRe))**2


T_g = 600 #757 = 15mtorr #gas temperature (K)
p = 5 #charactersitic readsorption length (cm)

T_g_str = str(T_g)
p_str = str(p)


#param = "0.0050mbar"

#0.5kW RFPowerTopView004
# I_696 = 3955.84
# I_727 = 1146.58
# I_738 = 6848.47
# I_706 = 4028.31
# I_794 = 4595.85  
# I_714 = 634.94

#0.5kW irradiance 
# I_696 = 44.2
# I_727 = 15.01
# I_738 = 99
# I_706 = 46.26
# I_794 = 109.68  
# I_714 = 7.54

#0.75kW irradiance
# I_696 = 60.68
# I_727 = 20.12
# I_738 = 136.37
# I_706 = 63.50
# I_794 = 156.12  
# I_714 = 10.16

#1.0kW irradiance 
#I_696 = 61.07
#I_727 = 20.39
#I_738 = 135.51
#I_706 = 63.91
#I_794 = 156.79  
#I_714 = 10.16

#1.25kW irradiance
I_696 = 85.07
I_727 = 28.1
I_738 = 185.81
I_706 = 86.44
I_794 = 228.43  
I_714 = 13.52

#1.5kW irradiance
#I_696 = 94.65
#I_727 = 30.71
#I_738 = 204.38
#I_706 = 94.97
#I_794 = 259.22  
#I_714 = 15.17


#1.75kW irradiance
#I_696 = 97.76
#I_727 = 31.83
#I_738 = 212.52
#I_706 = 97.74
#I_794 = 273.19  
#I_714 = 15.78

#2kW irradiance
#I_696 = 100.34
#I_727 = 32.27
#I_738 = 214.45
#I_706 = 100.05
#I_794 = 280.27  
#I_714 = 16.05

#2.25kW irradiance
#I_696 = 105.95
#I_727 = 33.71
#I_738 = 225.37
#I_706 = 104.78
#I_794 = 301.28  
#I_714 = 16.40

#2.5kW irradiance
# I_696 = 109.80
# I_727 = 34.77
# I_738 = 226.72
# I_706 = 104.68
# I_794 = 308.39  
# I_714 = 16.15

#2.75kW irradiance
#I_696 = 109.22
#I_727 = 34.36
#I_738 = 220.40
#I_706 = 102.27
#I_794 = 310.10  
#I_714 = 16.63

# #3.0kW irradiance 
# I_696 = 107.5
# I_727 = 33.8
# I_738 = 223.2
# I_706 = 102.4
# I_794 = 309.39  
# I_714 = 16.28

# #0.0020mbar
# I_696 = 31.79
# I_727 = 9.02
# I_738 = 62
# I_706 = 26.73
# I_794 = 88.78  
# I_714 = 4.3

#0.0030mbar
# I_696 = 40.78
# I_727 = 11.6
# I_738 = 82.65
# I_706 = 36.57
# I_794 = 114.61  
# I_714 = 5.61

#0.0040mbar
# I_696 = 48.49
# I_727 = 13.63
# I_738 = 100.5
# I_706 = 45.11
# I_794 = 134.89  
# I_714 = 6.65

#0.0050mbar
# I_696 = 55.1
# I_727 = 15.74
# I_738 = 119.49
# I_706 = 55.33
# I_794 = 152.87  
# I_714 = 7.64

#1mtorr
#I_696 = 25.4
#I_727 = 10.1
#I_738 = 66.2
#I_706 = 18.7
#I_794 = 73
#I_714 = 4.81

#15mtorr
#I_696 = 267
#I_727 = 118
#I_738 = 711
#I_706 = 384
#I_794 = 581
#I_714 = 63.6

#Transition Probabilties and reabsorption coeffcients for 6 Ar lines
#units for kij0 = cm^2K^0.5
#units for Aij = s^-1

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

#-------------------------------------------------------
with open('Density_inputs.txt', 'r') as inputs:

#Extract data and declare variables
    n1s4, n1s5 = np.loadtxt("Density_inputs.txt", delimiter=';', 
                            unpack=True, usecols=(0, 1))
#print(n1s4)
#print(n1s5)
    
with open('mod_results.txt', 'w') as resultsfile:
    resultsfile.write('Metastable and Resonant Density Model Results\n')
    resultsfile.write('Datetime (Y-M-S) = ' + datestring + '\n')
    resultsfile.write('Experiment Name: ' + param + '\n')
    resultsfile.write('Gas Temperature (K) = ' + T_g_str + '\n')
    resultsfile.write('Characteristic length (cm) = '+ p_str + '\n')


#test
#n1s4 =[1,2,3,4,5]
#n1s5 = [10,20,30,40,50]

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
        
        chi_sum=chi_squared(LRm1,LRe1)+chi_squared(LRm2,LRe2)+chi_squared(LRm3,LRe3)
        #print('Chi squared for n_r= ',n_1s4, 'and n_m= ',n_1s5, 'is: ',chi_sum)
        
        with open('mod_results.txt', 'a+') as mod_results:
        #for i in range(len(n1s4)):
         #   mod_results.write("%d %d %d\n" % (n1s4[i],n1s5[i],chi_sum))
            mod_results.write("%7.2e %7.2e %6.2f\n" % (a,b,chi_sum))
        
       # mod_results = open("mod_results.txt",'a+')
        #mod_results.write(a,';',b,';',chi_sum)
        #mod_results.write()

#mod_results.close()

#printing results to screeen

with open('mod_results.txt', 'a+') as mod_results:
    
    d, e, f = np.loadtxt("mod_results.txt",unpack=True,usecols=(0, 1, 2),skiprows=5)
    
    #print(d)
    #print(e)
    #print(f)
    
    print_results(d,e,f) #calling function to print results
    
    #min_pos = f.index(min(f))
    
    #print('Minimum chi-squared of ', min(f), 'occurs at position ', min_pos, 
     #     'where n1s4 = ',"%7.1e" % d[63], 'cm\u00b3', 'and n1s5 = ',
      #    "%7.1e" % e[63], 'cm\u00b3') 

#renaming the file with the current date and time
os.rename("mod_results.txt", time.strftime("n1sDEN_"+param+"_%Y%m%d%H%M.txt")) 

#Calculating time program took to run
final_time = time.time() - start_time

print("This program took", "%5.3f" %  final_time,"s to run")


                
        
          
    
        
    
        
    
    
    
    
    
        
    
    
    
    
    
