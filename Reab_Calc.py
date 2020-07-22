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

T_g = 600 #757 = 15mtorr #gas temperature (K)
p = 10 #charactersitic readsorption length (cm)
#M = 9.109E-31
M = 39.948# 6.6335209E-23 #kg
#M = 6.6335209E-26 #atomic mass of Ar(kg)
R = 8.31446261815324*(1E4)


#Function to calculate the reabsorption coefficients k_o for each of the 22 lines
#that must be used to determine the contributuon by radiation trapping.

def ko_calc(lamda,gi,gj,Aij):
    term1 = ((lamda*(1E-7))**3)/(8*(np.pi**(1.5)))
    #print("term 1 =",term1)

    term2 = gi/gj
    #print("term 2 =", term2)
    
    term3 = Aij*(sqrt(M))/((sqrt(2)*(sqrt(R))))
    #print("term 3 =", term3)
    
    result = term1*term2*term3
    
    #arr_k = np.array(result)
    
    return result 

def ko_data(n,gi,gj,A,nj):

    with open("line_data_full.txt","w+") as datafile:
    #datafile = open('line_data_full.txt', 'w+')
    
        for a, b, c, d, e in zip(n,A,gi,gj,nj):
        
            #print("----------------------")
            #print("For spectral line n_ij =", a, "nm")
            #print("Where g_i =", b)
            #print("and g_j=", c)
            #print("n_j =", "%5.2e" % d)
            #print("And Aij =", "%5.2e" % d)
            #print("k_o = ", "%7.2e" % ko_calc(a,b,c,d))
            ko = ko_calc(a,b,c,d)
                #for i in range(len(n1s4)):
                 #   mod_results.write("%d %d %d\n" % (n1s4[i],n1s5[i],chi_sum))
            datafile.write("%6.2f %3.0f %3.1f %3.1f %7.2e %7.2e\n" % (a,b,c,d,e,ko))
                    #np.append(arr,arr_k,axis=1)
        
        print("")                
        print("Calculated reabsorption coeffients and printed to file:line_data_full.txt")
        
    #datafile.close()
    
    #n_ij, A, g_i, g_j, n_j, ko = np.loadtxt("line_data_full.txt", comments='#', delimiter=' ', unpack=True, 
     #                              usecols=(0,1,2,3,4,5))    

n_ij, A, g_i, g_j, n_j = np.loadtxt("line_data2.txt", comments='#', delimiter=';', skiprows=2, unpack=True, 
                                    usecols=(0,1,2,3,4))

print("g_i 738 =",g_i[7],"g_j_738= ",g_j[7])
print("loaded parameters from line_data2.txt")
    
#os.rename("line_data3.txt", "line_data_full.txt")                                                            
#os.rename("k_oresults.txt", time.strftime("Ko_%Y%m%d%H%M.txt")) 

ko_data(n_ij,g_i,g_j,A,n_j)
    
dfd=pd.read_csv("line_data_full.txt",sep=" ",header=None,names=['Wavelength (nm)','A','g_i','g_j','n_j','k_o'],comment="#")
dfd.to_csv('line_data_full.csv',sep=';',index=False)


        