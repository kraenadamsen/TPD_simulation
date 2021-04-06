# -*- coding: utf-8 -*-
"""
Simulation class

"""

from scipy import interpolate
import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt

class Simulation ():
    
    def __init__(self,experiment,Channel):
        
        self.Channel = Channel
        
        
        if experiment == None:
            pass
            
        else:          
            self.Data_set = experiment.Data_set
            
            
        global k_b  
        k_b = 8.617333262145*10**(-5)
        
   
    def Calculate_energy(self,nu,signal,Coverage,T):



        Energy =   - np.log(signal[:-1]/(nu*Coverage)) * k_b * T[:-1] 

        return Energy  
    
    
    def Calculate_coverage(self , start_temp , end_temp):
        
        Spectrum = self.ref_spectrum
        
        self.start_index = np.argmin(abs(Spectrum['T']-start_temp))
        self.end_index   = np.argmin(abs(Spectrum['T']-  end_temp)) 

        coverage = []
        for j in range(self.start_index , self.end_index-1):
            coverage.append(np.trapz(Spectrum.get(self.Channel)[j: self.end_index],Spectrum.get('sec')[j: self.end_index]))
       
        return  np.array([coverage])

    
    def Runge_kutta_calulator(self,f,u0,t0,tf,n):
        
        t = np.linspace(t0, tf, n+1)
        u = np.array((n+1)*[u0])
        h = t[1]-t[0]
        for i in range(n):
            k1 = h * f(u[i], t[i])    
            k2 = h * f(u[i] + 0.5 * k1, t[i] + 0.5*h)
            k3 = h * f(u[i] + 0.5 * k2, t[i] + 0.5*h)
            k4 = h * f(u[i] + k3, t[i] + h)
            u[i+1] = u[i] + (k1 + 2*(k2 + k3 ) + k4) / 6
            
        return u, t
    
    
    def Create_enegy_landscape(self,Coverage,Energy):
        
        E_des = interpolate.interp1d(Coverage,Energy,fill_value="extrapolate")
        
        return E_des
    
    def Solve_polanyi_wigner_vary_prefactor(self,beta, base_spectrum_number ,nu_range,start_coverages,start_temp , end_temp, steps ):
        
        
        
        self.ref_spectrum = self.Data_set[base_spectrum_number]        
        return_dict = {}
        
        T0 = start_temp
        
        Coverage = self.Calculate_coverage(start_temp , end_temp)
        
        t0 = 0
        tf = end_temp - start_temp       
        n  = steps
               
        for nu in nu_range:

            Energy = self.Calculate_energy(nu,
                                      self.ref_spectrum[self.Channel][self.start_index:self.end_index],
                                      Coverage,
                                      self.ref_spectrum['T'][self.start_index:self.end_index])
            
            E_des = self.Create_enegy_landscape(Coverage[0,:],Energy[0,:])
            
            dtheta_dt =  lambda  theta, t : - nu * theta * np.exp(-E_des(theta)/(k_b * (t * beta + T0)))    
            
            U = []
                   
            for start_coverage in start_coverages:
                
                u0 = start_coverage        
              
                u , t = self.Runge_kutta_calulator(dtheta_dt,u0,t0,tf,n) 
                
                U.append(u)           
                
            U = np.array(U)
            
            rate = -dtheta_dt(U,t)
            
            return_dict.update({f"{nu:.2E}": rate})
            
        T = (t * beta + T0)
        
        return [return_dict,t,T,self.start_index,self.end_index]
        
        
        
        
        
        
            