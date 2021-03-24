# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:11:08 2020

@author: beaub
"""
import numpy as np
from math import exp
from math import sqrt
import pandas as pd
from math import pi

from datetime import datetime
datestring = datetime.strftime(datetime.now(), '%Y-%m-%d-%H-%M-%S')


#Function for calculation of escape factors, g

#swapped reabs_coeff for k
def escape_factor(k,n,p):
    p = int(p)
    return (2-exp(-k*p/1000))/(1+(k*p))

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

def ko_explicit(lamda,gi,gj,Aij,M,R):
    T1 = (((lamda*1E-9)**3)/(8*(pi**(3/2))))*(gi/gj)*Aij*((sqrt(M)/sqrt(2)))
    T2 = (T1*((R)**(-0.5)))*(1E4)
    # print("T1: ", T1)
    # print("T2 = ", T2)
    
    return T2

def reab_coeff_calc(n_ij,g_i,g_j,A,n_j,T_g,M,R):
        with open("line_data_full.txt","a+") as datafile:      
            for a, b, c, d, e in zip(n_ij,g_i,g_j,A,n_j):
                #print("ko for ", a, "nm is: ",ko_calc(a,b,c,d))
                ko = ko_explicit(a,b,c,d,M,R) 
                kij = ko*e*((T_g)**(-0.5))
                #print( " and kij = ", kij)
                #print(" ")
                
                datafile.write("%6.2f %1.0f %1.0f %5.2e %7.2e %7.2e %7.5f \n" % (a,b,c,d,e,ko,kij))

def radtrap(n,A,gi,gj,nj,p,T_g,M,R):
    
    print("Calculating absorption coefficients...")
    print("Fetching ko data from Reab_Calc.py")
    
    reab_coeff_calc(n, gi, gj, A, nj,T_g,M,R)
    
    print("")
    print("ko and kij calculated and printed to file.")
    
    df=pd.read_csv("line_data_full.txt",sep=" ",header=None,names=['Wavelength (nm)','g_i','g_j','A','n_j','k_o','k_ij'],usecols=(0,1,2,3,4,5,6),comment="#")
    df.to_csv('line_data_full.csv',sep=';',index=True)
    
    print(" ")
    print("Absorption coefficients calculated.")
        
    R_750 = (((escape_factor(df.at[0,"k_ij"],df.at[0,"n_j"],p))*((df.at[0,"A"]) + 
             df.at[1,"A"])))/(((escape_factor(df.at[0,"k_ij"],df.at[0,"k_ij"],p))*df.at[0,"A"]) + 
                  (escape_factor(df.at[1,"k_ij"],df.at[1,"n_j"],p))*df.at[1,"A"])
    
    R_772_4 = (((escape_factor(df.at[3,"k_ij"],df.at[3,"n_j"],p))*(df.at[3,"A"] + 
             df.at[2,"A"] + df.at[4,"A"] + df.at[5,"A"])))/(((escape_factor(df.at[3,"k_ij"],df.at[3,"n_j"],p))*df.at[3,"A"]) + 
                  ((escape_factor(df.at[2,"k_ij"],df.at[2,"n_j"],p))*df.at[2,"A"]) + 
                  ((escape_factor(df.at[4,"k_ij"],df.at[4,"n_j"],p))*df.at[4,"A"]) + 
                  ((escape_factor(df.at[5,"k_ij"],df.at[5,"n_j"],p))*df.at[5,"A"]))
                                                            
    R_738 = ((escape_factor(df.at[7,"k_ij"],df.at[7,"n_j"],p))*(df.at[7,"A"] + 
             df.at[6,"A"]+df.at[8,"A"]))/(((escape_factor(df.at[7,"k_ij"],df.at[7,"n_j"],p))*df.at[7,"A"]) + 
                  ((escape_factor(df.at[6,"k_ij"],df.at[6,"n_j"],p))*df.at[6,"A"]) +
                  (escape_factor(df.at[8,"k_ij"],df.at[8,"n_j"],p)*df.at[8,"A"]))
    
    R_794 = (((escape_factor(df.at[10,"k_ij"],df.at[10,"n_j"],p))*((df.at[10,"A"]) + 
             df.at[9,"A"] + df.at[11,"A"] + df.at[12,"A"])))/(((escape_factor(df.at[10,"k_ij"],df.at[10,"n_j"],p))*df.at[10,"A"]) + 
                  ((escape_factor(df.at[9,"k_ij"],df.at[9,"n_j"],p))*df.at[9,"A"]) + 
                  ((escape_factor(df.at[11,"k_ij"],df.at[11,"n_j"],p))*df.at[11,"A"]) + escape_factor(df.at[12,"k_ij"],df.at[12,"n_j"],p))
                                                                 
    R_751 = (((escape_factor(df.at[14,"k_ij"],df.at[14,"n_j"],p))*(df.at[14,"A"] + 
             df.at[13,"A"])))/(((escape_factor(df.at[14,"k_ij"],df.at[14,"k_ij"],p))*df.at[14,"A"]) + 
                  (escape_factor(df.at[13,"k_ij"],df.at[13,"n_j"],p))*df.at[13,"A"])
  
    R_763 = (((escape_factor(df.at[17,"k_ij"],df.at[17,"n_j"],p))*(df.at[17,"A"] + 
             df.at[15,"A"] + df.at[16,"A"])))/(((escape_factor(df.at[17,"k_ij"],df.at[17,"k_ij"],p))*df.at[17,"A"]) + 
                  (escape_factor(df.at[15,"k_ij"],df.at[15,"n_j"],p))*df.at[15,"A"] + 
                  ((escape_factor(df.at[16,"k_ij"],df.at[16,"n_j"],p))*df.at[16,"A"]))
    
    R_772_3 = ((escape_factor(df.at[21,"k_ij"],df.at[21,"n_j"],p))*(df.at[21,"A"] + 
             df.at[18,"A"] + df.at[19,"A"] + df.at[20,"A"]))/(((escape_factor(df.at[21,"k_ij"],df.at[21,"n_j"],p))*df.at[21,"A"]) + 
                  ((escape_factor(df.at[18,"k_ij"],df.at[18,"n_j"],p))*df.at[18,"A"]) +
                  (escape_factor(df.at[19,"k_ij"],df.at[19,"n_j"],p)*df.at[19,"A"]) + 
                  (escape_factor(df.at[20,"k_ij"],df.at[20,"n_j"],p)*df.at[20,"A"]))
    
    #print("")
    Rad = [R_750,R_772_4,R_738,R_794,R_751,R_763,R_772_3]
    R_lam = [750.39,772.42,738.4,794.82,751.47,763.51,772.38]
    
    # print(" ")
    df2 = pd.DataFrame(list(zip(R_lam,Rad)),columns=["Wavelength","Radation Trapping Coefficient"])
    df2.to_csv(r'Rad_TrapCoeff.csv', sep=';', index = False)
    print("Radiation trapping coefficients, along with their corresponding ")
    print("wavelength have been printed to Rad_TrapCoeff.csv")
    
    return Rad, df
    
def print_escape(n_ij,g_i,g_j,A,n_j,ko,k_ij,p):
    
    with open("line_escape_factors.txt","w+") as datafile:
        for a, b, c, d, e, f, g in zip(n_ij,g_i,g_j,A,n_j,ko,k_ij):
            g_ij = escape_factor(g,e,p)
            gA = g_ij*d
            datafile.write("%6.2f %1.0f %1.0f %5.2e %7.2e %7.2e %7.5f %5.3f %3.5f \n" % (a,b,c,d,e,f,g,gA,g_ij))
         
    n_ij, g_i, g_j, A, n_j, ko, k_ij, gA, g_ij = np.loadtxt("line_escape_factors.txt", comments='#', delimiter=' ', unpack=True, 
                                             usecols=(0,1,2,3,4,5,6,7,8))
