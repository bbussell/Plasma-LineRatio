# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 19:34:41 2020

@author: beaub
"""

from astroquery.nist import Nist 
from astropy.table import Table, join
import astropy.units as u

def Fetch_NIST():

    print("Fetching the latest NIST spectral data...")
    print("")
    
    #astropy NIST class, querying the NIST database
    #from 4000A to 7000A, Ar I spectrum, wavelengths in vacuum. 
    #BF emission lines
    sp_750 = Nist.query(750 * u.nm, 750.6 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_750 = sp_750.to_pandas(index=None)
    
    sp_667 = Nist.query(667 * u.nm, 668 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_667 = sp_667.to_pandas(index=None)
    
    sp_826 = Nist.query(826 * u.nm, 827 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_826 = sp_826.to_pandas(index=None)
    
    sp_772_4 = Nist.query(772.4 * u.nm, 773.45 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_772_4 = sp_772_4.to_pandas(index=None)
    
    sp_727 = Nist.query(727 * u.nm, 728 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_727 = sp_727.to_pandas(index=None)
    
    sp_696 = Nist.query(696.5 * u.nm, 697 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_696 = sp_696.to_pandas(index=None)
    
    # sp_840 = Nist.query(840 * u.nm, 841 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    # df_840 = sp_840.to_pandas(index=None)
    
    sp_738 = Nist.query(738 * u.nm, 739 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_738 = sp_738.to_pandas(index=None)
    
    sp_706 = Nist.query(706* u.nm, 707 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_706 = sp_706.to_pandas(index=None)
    
    sp_852 = Nist.query(852 * u.nm, 853 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_852 = sp_852.to_pandas(index=None)
    
    sp_794 = Nist.query(794 * u.nm, 796 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_794 = sp_794.to_pandas(index=None)
        
    sp_747 = Nist.query(747 * u.nm, 748 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_747 = sp_747.to_pandas(index=None)
        
    sp_714 = Nist.query(714 * u.nm, 715 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_714 =sp_714.to_pandas(index=None)
    
    # sp_857 = Nist.query(857 * u.nm, 858 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    # df_857 = sp_857.to_pandas(index=None)
    
    sp_751 = Nist.query(751 * u.nm, 752 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_751 = sp_751.to_pandas(index=None)
    
    sp_922 = Nist.query(922 * u.nm, 923 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_922 = sp_922.to_pandas(index=None)
    
    sp_800 = Nist.query(800 * u.nm, 801 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_800 = sp_800.to_pandas(index=None)
    
    sp_763 = Nist.query(763 * u.nm, 764 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_763 = sp_763.to_pandas(index=None)
    
    sp_935 = Nist.query(935 * u.nm, 936 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_935 = sp_935.to_pandas(index=None)
    
    # sp_866 = Nist.query(866 * u.nm, 867 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    # df_866 = sp_866.to_pandas(index=None)
    
    sp_810 = Nist.query(810 * u.nm, 811 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    df_810 = sp_810.to_pandas(index=None)
    
    # sp_772_3 = Nist.query(772.3 * u.nm, 772.4 * u.nm, linename="Ar I", wavelength_type="vacuum", output_order="wavelength", energy_level_unit="eV")
    # df_772_3 = sp_772_3.to_pandas(index=None)
    
    all_lines = df_696.append([df_750,df_667,df_826,df_772_4,df_727,df_696,df_738,df_706,df_852,df_794,df_747,df_714,df_751,df_922,df_800,df_763,df_935,df_810], ignore_index=True)

    final_lines = all_lines.drop(['Ritz','Acc.', 'TP', 'Line', 'Type'],axis=1)
    final_lines.to_csv('all_lines.txt',sep=";", index=False)

    return all_lines

# print("The spectral data is: ")
# print(" ")
# print(Fetch_NIST())

# fetched_lines = Fetch_NIST()

#total_table = join(table,table2)
#print(total_table)

