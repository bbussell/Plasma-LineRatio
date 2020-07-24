# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:11:08 2020

@author: beaub
"""
import numpy as np
from math import exp
from math import sqrt
import pandas as pd
import csv
import os
import time
import Reab_Calc as RC

from datetime import datetime
datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')

T_g = 600 #757 = 15mtorr #gas temperature (K)
p = 10 #charactersitic readsorption length (cm)
#M = 9.109E-31
M = 39.948# 6.6335209E-23 #kg
#M = 6.6335209E-26 #atomic mass of Ar(kg)
R = 8.31446261815324*(1E4)

# n_1s3 = n_1s5/6.5
# n_1s2 = n_1s4

#Extract data from original file, calculate k_o values and then extract full data,
#including n_ij, A, g_i, g_j, n_j and k_o

#Function for caluclation of Excitation rates

def ex_rates(k,T,a,E):
    return k*((T)**a)*exp(-E/T)

def sum_nk(k,T,a,E,km,am,Em,kr,ar,Er):    
    return n_g*ex_rates(k,T,a,E) + n_m*ex_rates(km,T,am,Em) +  n_r*ex_rates(kr,T,ar,Er)

def reabs_coeff(k,n):
    return k*n*(T_g**-0.5)

#Function for calculation of escape factors, g

def escape_factor(k,n,p):
    return (2-exp(-reabs_coeff(k,n)*p/1000))/(1+(reabs_coeff(k,n)*p))

def rad_trap(k_ij,p,Ag,Am,Ar,kg,km,kr):
    return (escape_factor(k,n,p)*(Ag+Am+Ar))/((Ag*escape_factor(kg,n,p))+
            (Am*escape_factor(k,n,p))+(Ar*escape_factor(k,n,p)))
    
def LR_750(k,T,a,E,km,am,Em,kr,ar,Er):
    sum_750=sum_nk(kG_750,T,alphaG_750,EG_750,kM_750,alphaM_750,EM_750,kR_750,alphaR_750,ER_750) 
    
def LR_mod(k,T,a,E,km,am,Em,kr,ar,Er,Rij,Rmn):
    LR = Rij/Rmn(sum_nk(k,T,a,E,km,am,Em,kr,ar,Er)/LR_750(k,T,a,E,km,m,kr,ar,Er))

n_ij, A, g_i, g_j, n_j = np.loadtxt("line_data2.txt", comments='#', delimiter=';', skiprows=2, unpack=True, 
                                        usecols=(0,1,2,3,4))

#R738 = ((escape_factor(df.at[7,"k_o"],df.at[7,"n_j"],p))*(df.at[7,"A"] + 
#             df.at[6,"A"]+df.at[8,"A"]))/(((escape_factor(df.at[7,"k_o"],df.at[7,"n_j"],p))*df.at[7,"A"]) + 
#                  ((escape_factor(df.at[6,"k_o"],df.at[6,"n_j"],p))*df.at[6,"A"]) +
#                  (escape_factor(df.at[8,"k_o"],df.at[8,"n_j"],p)*df.at[8,"A"]))


def radtrap(n,A,gi,gj,nj):
    
    print("Calculating absorption coefficients...")
    print("Fetching ko data from Reab_Calc.py")
    
    RC.ko_data(n,gi,gj,A,nj)
    
    df=pd.read_csv("line_data_full.txt",sep=" ",header=None,names=['Wavelength (nm)','A','g_i','g_j','n_j','k_o','k_ij'],comment="#")
    df.to_csv('line_data_full.csv',sep=';',index=False)
    
    print("Absorption coefficients calculated.")
    
    #with open("k_oresults.txt","r") as f:
    
    #datafile = open('k_oresults.txt', 'a+')       
    
    #Calculation of radiation trapping coefficients for each upper level, taking into account all transitions to
    #the lower level of interest (with some minor exceptions)
    #2p1
    
    #with open("RadTrap_Coeff.csv","w+") as radfile:
    
    print("The escape factor for 738 is:",escape_factor(df.at[7,"k_o"],df.at[7,"n_j"],p))
    print("g*A for 738 is:",escape_factor(df.at[7,"k_o"],df.at[7,"n_j"],p)*(df.at[7,"A"] + 
             df.at[6,"A"]+df.at[8,"A"]))
        
    R_750 = (((escape_factor(df.at[0,"k_o"],df.at[0,"n_j"],p))*((df.at[0,"A"]) + 
             df.at[1,"A"])))/(((escape_factor(df.at[0,"k_o"],df.at[0,"k_o"],p))*df.at[0,"A"]) + 
                  (escape_factor(df.at[1,"k_o"],df.at[1,"n_j"],p))*df.at[1,"A"])
                              
    #print(R_750)
    
    #2p2
    R_772_4 = (((escape_factor(df.at[3,"k_o"],df.at[3,"n_j"],p))*(df.at[3,"A"] + 
             df.at[2,"A"] + df.at[4,"A"] + 
             df.at[5,"A"])))/(((escape_factor(df.at[3,"k_o"],df.at[3,"k_o"],p))*df.at[3,"A"]) + 
                  ((escape_factor(df.at[2,"k_o"],df.at[2,"n_j"],p))*df.at[2,"A"]) + 
                  ((escape_factor(df.at[4,"k_o"],df.at[4,"n_j"],p))*df.at[4,"A"]) + 
                  ((escape_factor(df.at[5,"k_o"],df.at[5,"n_j"],p))*df.at[5,"A"]))
    #print(R_772_4)
    
    #2p3
    
    R_738 = ((escape_factor(df.at[7,"k_o"],df.at[7,"n_j"],p))*(df.at[7,"A"] + 
             df.at[6,"A"]+df.at[8,"A"]))/(((escape_factor(df.at[7,"k_o"],df.at[7,"n_j"],p))*df.at[7,"A"]) + 
                  ((escape_factor(df.at[6,"k_o"],df.at[6,"n_j"],p))*df.at[6,"A"]) +
                  (escape_factor(df.at[8,"k_o"],df.at[8,"n_j"],p)*df.at[8,"A"]))
    #print(R_738)
    
    #2p4
    
    R_794 = (((escape_factor(df.at[10,"k_o"],df.at[10,"n_j"],p))*((df.at[10,"A"]) + 
             df.at[9,"A"] + df.at[11,"A"])))/(((escape_factor(df.at[10,"k_o"],df.at[10,"k_o"],p))*df.at[10,"A"]) + 
                  ((escape_factor(df.at[9,"k_o"],df.at[9,"n_j"],p))*df.at[9,"A"]) + 
                  ((escape_factor(df.at[11,"k_o"],df.at[11,"n_j"],p))*df.at[11,"A"]))
    #print(R_794)
    
    #2p5
    
    R_751 = ((escape_factor(df.at[12,"k_o"],df.at[12,"n_j"],p))*((df.at[12,"A"]))) / (((escape_factor(df.at[12,"k_o"],df.at[12,"k_o"],p))*df.at[12,"A"]))
    #print(R_751)
    
    #2p6
    
    R_763 = (((escape_factor(df.at[14,"k_o"],df.at[14,"n_j"],p))*(df.at[14,"A"] + 
             df.at[13,"A"])))/(((escape_factor(df.at[14,"k_o"],df.at[14,"k_o"],p))*df.at[14,"A"]) + 
                  (escape_factor(df.at[13,"k_o"],df.at[13,"n_j"],p))*df.at[13,"A"])
    #print(R_763)
    
    #2p7
    
    R_772_3 = ((escape_factor(df.at[17,"k_o"],df.at[17,"n_j"],p))*(df.at[17,"A"] + 
             df.at[15,"A"]+df.at[16,"A"]))/(((escape_factor(df.at[17,"k_o"],df.at[17,"n_j"],p))*df.at[17,"A"]) + 
                  ((escape_factor(df.at[15,"k_o"],df.at[15,"n_j"],p))*df.at[15,"A"]) +
                  (escape_factor(df.at[16,"k_o"],df.at[16,"n_j"],p)*df.at[16,"A"]))
    #print(R_772_3)
    
    #print("")
    Rad = [R_750,R_772_4,R_738,R_794,R_751,R_763,R_772_3]
    R_lam = [750.39,772.42,738.4,794.82,751.47,763.51,772.38]
    #for testing purposes
    
    Rad_arr = np.array(Rad)
    Rad_lam_arr = np.array(R_lam)
    
    df2 = pd.DataFrame(list(zip(R_lam,Rad)),columns=["Wavelength","Radation Trapping Coefficient"])
    df2.to_csv(r'Rad_TrapCoeff.csv', sep=';', index = False)
    print("Radiation trapping coefficients, along with their corresponding ")
    print("wavelength have been printed to Rad_TrapCoeff.csv")
    
    #with open("rad_trap.txt","w+") as datafile:
     #   for a, b in zip(Rad_lam_arr, Rad_arr):
      #      datafile.write("%3.1f %7.2e\n" % (a,b))
    
    #print(df)
    #print(df2)
    #print(EF_all)
    
    return Rad, df
    #return R_750, R_772_4, R_738, R_794, R_751, R_763, R_772_3
    
#radtrap(n_ij,A,g_i,g_j,n_j)

df=pd.read_csv("line_data_full.txt",sep=" ",header=None,names=['Wavelength (nm)','A','g_i','g_j','n_j','k_o'],comment="#")
                                                               
#print(df)
                                                               
#EF_all=[escape_factor(df.at[0,"k_ij"],df.at[0,"n_j"],p),escape_factor(df.at[3,"k_ij"],df.at[3,"n_j"],p),escape_factor(df.at[7,"k_ij"],df.at[7,"n_j"],p),escape_factor(df.at[10,"k_ij"],df.at[10,"n_j"],p),escape_factor(df.at[12,"k_ij"],df.at[12,"n_j"],p),escape_factor(df.at[14,"k_ij"],df.at[14,"n_j"],p),escape_factor(df.at[17,"k_ij"],df.at[17,"n_j"],p)]

#a = EF_all
#print(EF_all)

#print(R_list)
#print("")
#R = np.array(R_list)

#R_lam = np.array(R_lamda)

   # np.T([[R_arr]])
#np.T([[n_ij]])   


#for i in R_list:
 #   radfile.write("%7.2e\n" % (i))
#print("printed radiation trapping coefficients to RadTrap_Coeff.txt")
#print(R_arr)



  
  

