# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:56:10 2020

@author: beaub
"""
import numpy as np
import itertools as it
#import sys

#datafile = open("LR_data2.txt", "r")
with open('Density_Tg_inputs.txt', 'w') as density_tg_results:

    #Extract data and declare variables
    n1s4 = np.loadtxt("intervals.txt", delimiter=';', unpack=True,
                          usecols=(0))
    n1s5 = n1s4 
    
    tg = np.loadtxt("Tg_intervals.txt", delimiter=None, unpack=True, usecols=(0))
    #n1s5 = n1s4
    
    #Real data 
    #n1s4=[5E9,1.0E10,1.5E10,2E10,2.5E10,3.0E10,3.5E10,4.0E10]
    #n1s5=[1E10,2E10,3E10,4E10,5E10,6E10,7E10,8E10]
    
    #test data
    #n1s4=[1,2,3]
    #n1s5=[10,20,30]
    
    print('n1s4 =', n1s4)
    print('n1s5 =', n1s5)
    
    all_comb = []
    
    #In permutations(n1s4,r), r = the iterable, and the total number of 
    #input values, e.g. 3. 
    #or if the total number of values of each list is different then you can use 
    #len(n1s5)
    
    results = []
        
    for i in n1s4:
        for j in n1s5:
            for k in tg:
                results.append((i, j, k))
    
    print("Number of results:", len(results))
    
    #print(results)
    
    arr = np.asarray(results)
    #print(arr,sep=',')
    
    np.savetxt("Density_Tg_inputs.txt",arr,fmt=["%7.3e", "%7.3e", "%5.2f"],delimiter=";")
      
    #---------------------------------------------------------
    #Changing the list to a numpy array. Then reshaping this into a 2D array so 
    #it's easier to print to file.
    
    #print('\n'.join(map(str, all_comb)))     
    #print(*all_comb,sep="\n")
    #-----
    #arr = np.asarray(all_comb)
    
    #print(arr)
    #arrlst = arr.reshape(322560,2) #full No of combinations
    
    #---------------
    #arrlst = arr.reshape(18,2) #test No of combinations
    #print(arrlst)
    
    #np.savetxt("n_comb.txt",arrlst,fmt="%e",delimiter=",")
    
    #for row in arrlst:
        #print(row)
    
    #for j in list(all_comb):
        #print(j)    
    #print(*arr,sep='\n')
    
    #for row in arr:
    #    print(row)
    
    #np.savetxt(n_comb,arr)
    #print(all_comb[1])
    
    #print(all_comb)
    
    #theory_data = open('theory_data.txt', 'a+')
    #for i in all_comb:
        #theory_data.write("%f\n" % (str(all_comb)))
            #for i in range(len(n1s4)):
             #   mod_results.write("%d %d %d\n" % (n1s4[i],n1s5[i],chi_sum))
    #theory_data.write("%f\n" % (all_comb))
    
    #datafile.close()
    #theory_data.close()