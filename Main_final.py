
"""
Main simluation

"""
#%% import libaries

import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from scipy import interpolate

plt.close('all')


#%% # --------------------------- import files ---------------------------

files = glob.glob('Data_set/*.tpd')


Data_set = []

# importing loop 
for file in files:
    
    Data = np.loadtxt(file,skiprows=3)   # Reads all numeric data
    
    with open(file) as f:
        [next(f) for i in range(2)]                     # skips first 2 lines
        titles = f.readline().strip('\n').split('\t')   # Reads tiltles

    
    # creates hashtable of spectrum called dictonary in python
    Spectrum = {}
    for column in range (Data.shape[1]): # iteraes over columns

        Spectrum.update({titles[column]: Data[:,column]})
        
        
    # appending spectrum to data set
    
    Data_set.append(Spectrum)
    

    
#%% #----------------------- initial calculations -------------------------
"""
This includes:
    - background removal
    - Spectrum integration
    - Energy estimations (Given Predefined prefactors)
                          
"""

from Pre_calculations import pre_calculations 

"""
Find full describtions in Pre_caculations.py
"""

Channel = 'UTI_17'

Experiment = pre_calculations(Data_set,Channel,titles)

"""
Sort spectra from low to high Coverages

"""

Experiment.sort_spectra()


"""
Define refence spectrum 

Will find the one with thre lowest area 
and use as reference if none is specified
"""
Experiment.find_ref_spectrum()
#Experiment.find_ref_spectrum(ref_spectrum)

"""
Correct background/ substract reference spectrum
"""
Experiment.Correct_background()

"""
Estimimate heating rate i.e Beta
"""

Experiment.Estimate_beta_value()

"""
Calculate Coverages as function Temperature/Time
"""
start_temp = 110
end_temp   = 500 

Experiment.Calculate_coverage(start_temp,end_temp)

#%% ------------- inspection of precalculated results --------------

"""
Inscpection of spectra 
"""
Fig = plt.figure()
ax = Fig.gca()

Experiment.plot_spectra(ax)

"""
Coverage inspection
"""

Fig = plt.figure()
ax = Fig.gca()

Experiment.plot_coverage(ax)


"""
Quick energy inspection
"""
Fig = plt.figure()
ax = Fig.gca()

Experiment.Estimate_energy(nu=10**14,ax=ax)



#%% ------------------------- Simulations  ---------------------



from Simulator import Simulation

Simulations = Simulation(Experiment,Channel)


nu_range =  10 ** np.arange(start = 8 , stop = 16,step = 2, dtype = np.float)

Coverages = []

for i in range (len(Experiment.Data_set)):
    Coverages.append ( Experiment.Data_set[i]['coverage'].max())

Coverages = np.array(Coverages)+.01

base_spectrum_number = 4 


start_coverages = Coverages  [base_spectrum_number +1 :]


u ,t , T, start,  end   = Simulations.Solve_polanyi_wigner_vary_prefactor(
    beta = 1,
    base_spectrum_number = base_spectrum_number ,
    nu_range    =  nu_range,
    start_coverages   = start_coverages,
    start_temp  = 110,
    end_temp    = 500,
    steps       = 250)


Signal = interpolate.interp1d(Experiment.Data_set[7].get('T'),Experiment.Data_set[7].get(Channel))(T)


    
Fig = plt.figure()
ax = Fig.gca()

m = 0 

colors = pl.cm.jet(np.linspace(0,1,len(nu_range))) 

for precfactor in u.keys():
    
    Sim_Data = u[precfactor]
    
    for i in range(Sim_Data.shape[0]):
        ax.plot(T,Sim_Data[i,:],color = colors [m])
        


    m = m +1
    
    
Experiment.plot_spectra(ax)







