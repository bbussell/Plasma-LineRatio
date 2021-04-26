# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 14:39:04 2020

@author: beau.bussell
"""
import numpy as np
from math import sqrt
import pandas as pd
import time

#Constants

# T_g = 600 #757 = 15mtorr #gas temperature (K)
# p = 5 #charactersitic readsorption length (cm)
# #M = 9.109E-31
# M = 39.948# 6.6335209E-23 #kg
# #M = 6.6335209E-26 #atomic mass of Ar(kg)
# R = 8.31446261815324*(1E4)


#Function to calculate the reabsorption coefficients k_o for each of the 22 lines
#that must be used to determine the contributuon by radiation trapping.

def ko_calc(lamda,gi,gj,Aij,M,R):
    term1 = ((lamda*(1E-7))**3)/(8*(np.pi**(1.5)))
    #print("term 1 =",term1)

    term2 = gi/gj
    #print("term 2 =", term2)
    
    term3 = Aij*(sqrt(M))/((sqrt(2)*(sqrt(R))))
    #print("term 3 =", term3)
    
    result = term1*term2*term3
    
    #arr_k = np.array(result)
    
    return result 

def reab_coeff_calc(n_ij,g_i,g_j,A,n_j,T_g,M,R):
        with open("line_data_full.txt","w+") as datafile:      
            for Wavelength, StatisticalWeightUpper,StatisticalWeightLower,Probability,LowerDensity in zip(n_ij,g_i,g_j,A,n_j):
                #print("ko for ", a, "nm is: ",ko_calc(a,b,c,d))
                ko = ko_calc(Wavelength,StatisticalWeightUpper,StatisticalWeightLower,Probability,M,R) 
                kij = ko*LowerDensity*((T_g)**(-0.5))
                #print( " and kij = ", kij)
                #print(" ")
                datafile.write("%6.2f %1.0f %1.0f %5.2e %7.2e %7.2e %7.5f \n" % (Wavelength,StatisticalWeightUpper,StatisticalWeightLower,Probability,LowerDensity,ko,kij))
    


        