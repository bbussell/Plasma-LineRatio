# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 14:51:18 2020

@author: beaub
"""


import pytest
from math import exp
from Te_code import ex_rates, sum_nk, k_o, reabs_coeff, escape_factor
from Te_code import p, kG_738, alphaG_738, EG_738, kM_738, alphaM_738,EM_738,kR_738,alphaR_738,ER_738,n_r,Rad
from RadTrap_No2 import radtrap
#import Te_code.py as T_code

lamda_738 = 738.4

A_738 = 8470000
g_i_738 = 5
g_j_738 = 3

lamda_738_ls = [738.4]
A_738_ls = [8470000]
g_i_738_ls = [5]
g_j_738_ls = [3]
n_r_ls = [12000000000.0]

#Testing functions in Reab_Calc
#!!fails because the function uses a list and 
#I'm only providing one value

def test_radtrap():
    assert radtrap(lamda_738_ls,A_738_ls,g_i_738_ls,g_j_738_ls,n_r_ls) == 1.277025826

#Testing functions in Te_code.py

#Testing excitation rates function
def test_answer():
    assert ex_rates(kG_738,3.5,alphaG_738,EG_738) == 6.52E-12
    
#Testing sumnk function
def test_sumnk():
    assert sum_nk(kG_738,3.5,alphaG_738,EG_738,kM_738,alphaM_738,EM_738,kR_738,alphaR_738,ER_738) == 1066.26

#Testing function that calculates reabsorption coefficients 
def test_ko():
    assert k_o(lamda_738,g_i_738,g_j_738,A_738) == 1.98E-09

#Testing function that calculates full reabsorption coefficients
def test_reabs_coeff():
    assert reabs_coeff(k_o(lamda_738,g_i_738,g_j_738,A_738),n_r) == 0.96878
 
#Testing function that calculates escape factors
def test_escape_factor():
    assert escape_factor(k_o(lamda_738,g_i_738,g_j_738,A_738),1.2E10,p) == 0.094466606
    
def test_rcode738():
    Rad738 = Rad[2]
    assert Rad738 == 1.277
    #Rad == [0.949336776,0.701471798,1.277025826,0.673757912,1,0.492865753,0.653786693]
    
def test_rcode772():
    Rad772 = Rad[1] 
    assert Rad772 == 0.701471798






