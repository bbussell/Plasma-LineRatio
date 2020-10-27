# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 21:24:44 2020

@author: beaub
"""

import PySimpleGUI as sg

sg.theme('Reddit')

layout =    [[sg.Text('Extended Corona Model for HiTUS'), sg.Text(size=(15,1))],
            [sg.Text('Neutral Gas Temperature (K)', size=(25, 1)), sg.InputText(key='-T_IN-')],
            [sg.Text('Reabsorption Length (cm)', size=(25, 1)), sg.InputText(key='-p_IN-')],
            [sg.Text('Process Pressure (mbar)', size=(25, 1)), sg.InputText(key='-PP_IN-')],
            [sg.Text('Neutral Gas Temperature (K)'), sg.Text(size=(25,1), key='-T_O-')],
            [sg.Text('Reabsorption Length (cm)'), sg.Text(size=(25,1), key='-p_O-')],
            [sg.Text('Process Pressure (mbar)'), sg.Text(size=(25,1), key='-PP_O-')],
            [sg.Button('Show'), sg.Button('Exit')]]

# input2 =    [sg.Text('Please enter your Name, Address, Phone')],
#             [sg.Text('Name', size=(15, 1)), sg.InputText()],
#             [sg.Text('Address', size=(15, 1)), sg.InputText()],
#             [sg.Text('Phone', size=(15, 1)), sg.InputText()],
#             [sg.Submit(), sg.Cancel()]

window = sg.Window('Extended Corona Model', layout,size=(800,400), auto_size_text=True, auto_size_buttons=True)

while True:  # Event Loop
    event, values = window.read()
    print(event, values)
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    if event == 'Show':
        # Update the "output" text element to be the value of "input" element
        window['-T_O-'].update(values['-T_IN-'])
        window['-p_O-'].update(values['-p_IN-'])
        window['-PP_O-'].update(values['-PP_IN-'])

window.close()