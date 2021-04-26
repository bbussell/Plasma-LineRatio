# Plasma-LineRatio
Programs that model emission spectra line ratio's against an extended corona model for argon plasmas. These programs allow the user to undertake plasma diagnostics and measure properties such as metastable density, resonant density, electron temperature and electron density. 

Fundamental files include:

CoronaModelExecution.py : Main mode execution file which calls all further files for specific functions. 

MetaResonantDensityCalculatio.py : First model execution which models the 600-800nm emission lines against an extended corona model that includes emission line branching fractions. As the first model, this also preprocesses the data by normalising the raw data, interpolating the spectra, integrating the relevant peaks and initiating a background calculation which requires user input for each emission.
Then it conducts the modelling against 2500 density trial values and outputs a predcited measuremnt for the metastable and resonant density.

ElectronTemperatureCalculation.py : Second model execution which models the 600-800nm emission lines against an extended corona model that includes emission line excitation rates and branching fractions that consider the effects of radiation trapping. This takes the preprocessed data from the 1st model as well as model 1 density results as inputs, before conducting the electron temperature modelling using approximately 1600-2000 electron temperature trial inputs.

CalculateElectronDensity.py : Third and final model execution which models the 300-400nm emission lines against an electron quenching model that includes the electron temperature dependence of excitation rates. This takes the preprocessed data from model 1 as well as the electron temperature results from model 2 as inputs.

Other important files:

IntTimeAndSpectaFit.py : Normalises and integrates the original raw emission spectra ready for background determination and subtraction. 

background.py : Takes in the preprocessed spectral data after normalisation and interpoloation. Uses the trapezium method to estimate the background. Plots the emission peak and the background for each peak and then asks for confirmation from the user as to whether they are happy to procede with the chosen background limits. If not, they can input the correct limits for the background calculation and continue with the processing.






