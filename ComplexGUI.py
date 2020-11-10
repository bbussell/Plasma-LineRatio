# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 15:03:12 2020

@author: beaub
"""


import PySimpleGUI as sg

superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

sg.theme('Reddit')
#sg.ChangeLookAndFeel('GreenTan')

#form = sg.FlexForm('Everything bagel', default_element_size=(40, 1))

# column1 = [[sg.Text('Column 1', background_color='#d3dfda', justification='center', size=(10,1))],
#            [sg.Spin(values=('Spin Box 1', '2', '3'), initial_value='Spin Box 1')],
#            [sg.Spin(values=('Spin Box 1', '2', '3'), initial_value='Spin Box 2')],
#            [sg.Spin(values=('Spin Box 1', '2', '3'), initial_value='Spin Box 3')]]
layout = [
    [sg.Text('Extended Corona Model', size=(30, 1), font=("Helvetica", 25))],
    [sg.Text('Experiment Name:'), sg.InputText(key='-Exp_Name-', justification='left')],
    [sg.Text('Choose an experimental file from your file directory:',auto_size_text=True)],
    [sg.Text('Chosen file:', size=(15, 1), auto_size_text=False, justification='right'),
     sg.InputText('Browse File Directory'), sg.FileBrowse()],
    [sg.HorizontalSeparator(color='black')],
    [sg.Text('Input System Parameters:',font=("Arial",12,"bold"))],
    [sg.Text('Neutral Gas Temperature (K): ', size=(25, 1)), sg.InputText(key='-T_IN-',justification='left')],
    [sg.Text('Reabsorption Length (cm): ', size=(25, 1)), sg.InputText(key='-p_IN-',justification='left')],
    [sg.Text('Process Pressure (mbar): ', size=(25, 1)), sg.InputText(key='-PP_IN-',justification='left')],
    [sg.Text('Neutral Gas Temperature (K): ', size=(25,1)), sg.InputText(key='-T_O-',justification='left')],
    [sg.Text('Reabsorption Length (cm): ', size=(25,1)), sg.InputText(key='-p_O-',justification='left')],
    [sg.Text('Process Pressure (mbar): ',size=(25,1)), sg.InputText(key='-PP_O-',justification='left')],
    [sg.HorizontalSeparator(color='black')],
    [sg.Text('Using a calibrated spectrum:'), sg.Radio('Yes     ', "RADIO1", default=True), sg.Radio('No', "RADIO1")],
    [sg.Text('Using a cosine corrector:'), sg.Radio('Yes     ', "RADIO2", default=True), sg.Radio('No', "RADIO2")],
    #[sg.Text('Here is some text.... and a place to enter text')],
    #[sg.Checkbox('My first checkbox!'), sg.Checkbox('My second checkbox!', default=True)],
    #[sg.Radio('My first Radio!     ', "RADIO1", default=True), sg.Radio('My second Radio!', "RADIO1")],
    #[sg.Multiline(default_text='This is the default Text should you decide not to type anything', size=(35, 3)),
     #sg.Multiline(default_text='A second multi-line', size=(35, 3))],
    #[sg.InputCombo(('Combobox 1', 'Combobox 2'), size=(20, 3)),
     #sg.Slider(range=(1, 100), orientation='h', size=(34, 20), default_value=85)],
#    [sg.Listbox(values=('Listbox 1', 'Listbox 2', 'Listbox 3'), size=(30, 3)),
#     sg.Slider(range=(1, 100), orientation='v', size=(5, 20), default_value=25),
#     sg.Slider(range=(1, 100), orientation='v', size=(5, 20), default_value=75),
#     sg.Slider(range=(1, 100), orientation='v', size=(5, 20), default_value=10)],
     #sliders for initial metastable/resonant range 
    [sg.Text('The metastable/resonant density starting point (x10' + '8 cm-3):'.translate(superscript))],
    [sg.Slider(range=(0, 1000), orientation='h', size=(20, 5), default_value=25, resolution=1)],
    [sg.Text('The metastable/resonant density end point (x10' + '8 cm-3):'.translate(superscript))], 
    [sg.Slider(range=(0, 1000), orientation='h', size=(20, 5), default_value=75, resolution=1)],
     #sg.Column(column1, background_color='#d3dfda')],
    [sg.Text('_'  * 80)],
    #[sg.Text('Choose A Folder', size=(35, 1))],
    #[sg.Text('Your Folder', size=(15, 1), auto_size_text=False, justification='right'),
    # sg.InputText('Default Folder'), sg.FolderBrowse()],
    [sg.Submit(), sg.Cancel()]
     ]

window = sg.Window('Extended Corona Model', layout, auto_size_text=True, auto_size_buttons=True)

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
#button, values = form.Layout(layout).Read()
##form.LayoutAndRead(layout)
#sg.Popup(button, values)