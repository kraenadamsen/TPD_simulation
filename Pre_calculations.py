
"""

Pre - calculations 

"""

import numpy as np



class pre_calculations():
    
    
    def __init__(self,Data_set,Channel,Channel_titles):
        
        self.Data_set = Data_set
        self.Channel = Channel        
        self.titles = Channel_titles
        
        
        self.cutting_spectra()
        
        pass
    
    
    # ----------- define reference spectrum  -----------
    def sort_spectra(self):
        """
        Sort spectra from low to high area

        """   

        Areas = [ np.trapz(self.Data_set[i].get(self.Channel)) for i in range(len(self.Data_set)) ]        
        

        self.Data_set =  [x for x, _ in sorted(zip(self.Data_set , Areas) , key = lambda x: x[1] , reverse = True)]
       
        
        
        
    
    # ----------- define reference spectrum  -----------
    def find_ref_spectrum(self,ref_spectrum = None):
        """
        Finds lowest area spectra an subtracted this from the rest

        """
        
        if ref_spectrum == None:
            self.ref_spectrum = self.Data_set[0].get(self.Channel)
            self.ref_area = np.trapz(self.ref_spectrum)
            
            for i in range(len(self.Data_set)):
                
                Spectrum = self.Data_set[i].get(self.Channel)
                
                self.area = np.trapz(Spectrum)
                
                if self.ref_area > self.area:
                    self.ref_spectrum = Spectrum
                    self.ref_area = self.area
                   
            
            
        else:
            
             self.ref_spectrum = ref_spectrum
             
             
    #-------- cutting spectra length --------
    
    def cutting_spectra(self):
    
        max_len = len(self.Data_set[0].get(self.Channel))
        
        for i in range(len(self.Data_set)):
            
            length = len(self.Data_set[i].get(self.Channel))
            
            if max_len > length:
                max_len = length
    
        for i in range(len(self.Data_set)):

            for column in range (len(self.titles)-1): # iteraes over columns
        
                 self.Data_set[i].update({self.titles[column]: self.Data_set[i].get(self.titles[column])[:max_len]})
    
    #-------- background removal ----------
    def Correct_background (self):
    
        for i in range(len(self.Data_set)):
        
            Spectrum = self.Data_set[i].get(self.Channel)  
            
            
            self.Data_set[i].update({self.Channel: Spectrum-self.ref_spectrum})
            
    #-------- Calculates the heating rate ----------       
    def Estimate_beta_value (self):
        
        for i in range(len(self.Data_set)):
            
            Spectrum = self.Data_set[i]
            
            fit =np.polyfit(Spectrum.get('sec'),Spectrum.get('T'), 1)

            self.Data_set[i].update({'beta': fit[0]})
 
    # ----------- caculates coverages in numerical integration -------------
    def Calculate_coverage(self, start_temp , end_temp):
        
        
        Spectrum  = self.Data_set[0]
        
        self.start_index = np.argmin(abs(Spectrum['T']-start_temp))
        self.end_index   = np.argmin(abs(Spectrum['T']-  end_temp))
         
        for i in range(len(self.Data_set)):
   
           Spectrum = self.Data_set[i]
           
           coverage = []
           
           for j in range(self.start_index , self.end_index-1):
               coverage.append(np.trapz(Spectrum.get(self.Channel)[j: self.end_index],Spectrum.get('sec')[j: self.end_index]))
       
           
           self.Data_set[i].update({'coverage': np.array(coverage)})
            
            
    def Estimate_energy(self,nu,ax):
        Energy = []
        k_b = 8.617333262145*10**(-5)
        
        for i in range (1,len(self.Data_set)):
           
                signal   = self.Data_set[i].get(self.Channel)[self.start_index:self.end_index-1]
                Coverage = self.Data_set[i].get('coverage')
                T        = self.Data_set[i].get('T')[self.start_index:self.end_index-1]
                energy =   - np.log(abs(signal/(nu*Coverage))) * k_b * T 

                
                ax.plot(Coverage,energy)
                Energy.append(energy)
                
        ax.set_xlabel('Coverage [Monolayers]')
        ax.set_ylabel('Energy [eV]')
                
        return Energy  
            
    # -----------plot functions ----------
    
    def plot_spectra(self,ax):
        for i in range(len(self.Data_set)):
            
            ax.plot(self.Data_set[i].get('T'),self.Data_set[i].get(self.Channel))
            ax.set_xlabel('Temperature [K]')
            ax.set_ylabel('Desorption Rate [Monolayers/s]')
            
    def plot_coverage(self,ax):
        for i in range(len(self.Data_set)):
            
            ax.plot(self.Data_set[i].get('T')[self.start_index:self.end_index-1],self.Data_set[i].get('coverage'))
            ax.set_xlabel('Temperature [K]')
            ax.set_ylabel('Coverage [Monolayers]')
        
        
        
        