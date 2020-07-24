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
            k_ij = ko*e*((T_g)**0.5)

            print("for lamda = ", a, "ko = ", ko,"k_ij = ", k_ij)
                #for i in range(len(n1s4)):
                 #   mod_results.write("%d %d %d\n" % (n1s4[i],n1s5[i],chi_sum))
            datafile.write("%6.2f %3.0f %3.1f %3.1f %7.2e %7.2e %7.5f\n" % (a,b,c,d,e,ko,k_ij))
                    #np.append(arr,arr_k,axis=1)       
        print("")                
        print("Calculated reabsorption coeffients an printed to file:line_data_full.txt")
        
    return ko
        
    #datafile.close()
    
    #n_ij, A, g_i, g_j, n_j, ko = np.loadtxt("line_data_full.txt", comments='#', delimiter=' ', unpack=True, 
     #                              usecols=(0,1,2,3,4,5))    

#n_ij, A, g_i, g_j, n_j, ko, k_ij = np.loadtxt("line_data_full.txt", comments='#', delimiter=' ', unpack=True, 
                                   #usecols=(0,1,2,3,4,5,6))   
n_ij, A, g_i, g_j, n_j = np.loadtxt("line_data2.txt", comments='#', delimiter=';', skiprows=2, unpack=True, 
                                    usecols=(0,1,2,3,4))

#n_ij = [750.39,667.73]
#A = [47200000,235000]
#g_i = [1.0,1.0]
#g_j = [3.0,3.0]
#nj = [1.20E10,1.20E10]

print("!!!!!!!!!!!!! TESTING !!!!!!!!!!!!!!")
print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

#print("ko data with 750 and 738")
#print(ko_data(n_ij,g_i,g_j,A,nj))
with open("line_data_full.txt","w+") as datafile:
    
    for a, b, c, d, e in zip(n_ij,g_i,g_j,A,n_j):
        print("ko for ", a, "nm is: ",ko_calc(a,b,c,d))
        ko = ko_calc(a,b,c,d) 
        kij = ko*e*((T_g)**(-0.5))
        print( " and kij = ", kij)
        print(" ")
        
        datafile.write("%6.2f %3.0f %3.1f %3.1f %7.2e %7.2e %7.5f \n" % (a,b,c,d,e,ko,kij))
    
print("")
print("ko and kij calculated and printed to file.")
    

# print("g_i 738 =",g_i[7],"g_j_738= ",g_j[7])
# print("A_i 738 =",A[7], "n_ij = ", n_ij[7])
# print("n_j 738 =", n_j[7])
# print("loaded parameters from line_data2.txt")
# print("ko738 = ")
# print(ko_calc(n_ij[7],g_i[7],g_j[7],A[7]))

#testing
#k738 = ko[7]*n_j[7]*(T_g**-0.5)
#print("k738 =", k738)
    
#os.rename("line_data3.txt", "line_data_full.txt")                                                            
#os.rename("k_oresults.txt", time.strftime("Ko_%Y%m%d%H%M.txt")) 

#ko_data(n_ij,g_i,g_j,A,n_j)
    
dfd=pd.read_csv("line_data_full.txt",sep=" ",header=None,names=['Wavelength (nm)','A','g_i','g_j','n_j','k_o', 'k_ij'],comment="#")
dfd.to_csv('line_data_full.csv',sep=';',index=False)


        