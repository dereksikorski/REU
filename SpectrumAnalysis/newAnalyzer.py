## This file is analyzing CQ Spectra ##

## Imports
from re import X
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize, signal, ndimage
import os
import colorama as cl
import pandas as pd


class CQ:

    def __init__(self, param_path = None):
        """
        Create a "Cold Quasar" (CQ) object and analyze the object given the parameter file

        INPUT:
            - param_path (str)   :  None or the ABSOLUTE path to the parameter file. If no path given, uses "params.txt" in
            current directory
        """

        ## Store parameter path file:
        if param_path == None:      
            # Using the absolute path to params.txt
            self.parameters_path = os.getcwd() + '/params.txt' 
        else:                  
            self.parameters_path = param_path


        ## Initialize Colorama for Windows
        cl.init()

        ## Read parameters upon object creation
        self.readParameters()

        ## If User wants to just plot data, only plot data and end
        if self.plotSetting == "Data":
            self.plotData()


    def readParameters(self):
        """
        Read in the parameters for the Cold Quasar from the parameter file and store data
        """
        
        self.parameters_path = self.parameters_path.replace('\\', '/')      # Replace \ with / for windows users

        ## Load contents of the parameter file 
        if os.path.exists(self.parameters_path):    # Make sure file exists
            file_values_array = np.loadtxt(self.parameters_path, dtype=str, delimiter='=', usecols=(1))         
            print(cl.Fore.CYAN + "\nReading Parameters File Now! \n" +cl.Fore.WHITE)
        else:   # Raise error if file path does not exist
            raise Exception (cl.Fore.YELLOW + f"\nThe file path for the file containing the parameters:"
                                            + cl.Fore.RED + f"\n{self.parameters_path}\n"
                                            + cl.Fore.YELLOW + f"does not exist! Please enter a new path for the parameter file and try again!" + cl.Fore.WHITE)

        ## Read in parameters:
        self.data = pd.read_csv(file_values_array[0])     # Data from the .csv file
        
        # Read in spectrum data
        self.spectrum = [] 
        for index in list(self.data.index.values):
            self.spectrum.append([self.data["Wavelength"][index], self.data["Flux"][index]])
        self.spectrum = np.array(self.spectrum)

        # Redshift parameter
        self.redshift = float(file_values_array[1])
        
        # Number of prominences to fit
        self.prom_nums = int(file_values_array[2])

        # User plot setting
        self.plotSetting = file_values_array[3]


    def plotData(self, rs = True):
        """
        Plots the data with a redshift, unless specified

        INPUTS:
            - rs (bool)    : If True, applies redshift before plotting. If False, plots without redshift
        """
        ## Plot the data only
        if rs ==True:
            shifted_wl = self.DopplerShift(self.spectrum[:,0], self.redshift)
            self.spectrum = np.array(   [  [shifted_wl[i], self.spectrum[:,1][i]]  for i in range(len(shifted_wl))    ]   )

        plt.plot(self.spectrum[:,0], self.spectrum[:,1], 'b')
        plt.title("Flux vs Wavelength of Cold Quasar")
        plt.xlabel("Wavelength (Ang)")
        plt.ylabel("Flux")
        plt.show()


    def fitData(self):
        """
        Finds a Gaussian fit for the MgII, HB, and OIII lines as well as a fit for the continuum
        """
        ## Apply a redshift to the data:
        shifted_wl = self.DopplerShift(self.spectrum[:,0], self.redshift)

        self.spectrum = np.array(   [  [shifted_wl[i], self.spectrum[:,1][i]]  for i in range(len(shifted_wl))    ]   )

        ## Find peaks in data
        self.findPeaks()


        ## Fit Gaussians
        self.findGauss()

        ## Roughly fit Continuum
        self.roughContFit()

        self.fitGauss()
        ## Fit Continuum
        self.fitCont()



    def DopplerShift(self, wavelength_data, redshift_value):
        '''
        INPUTS:

            - wavelength_data (array)   :   Array of wavelengths
            - redshift_value (float)    :   Value of redshift to be applied to the wavelength array

        OUTPUT:

            - (array)  :  Wavelength array after the redshift has been applied

        This method helps "plotData" and "fitData" and is meant to shift the given data based on a given redshift value
        '''
        ## Apply a doppler shift to the array of wavelenghts based on the redshift value

        doppler_shift = lambda wl: wl / (1 + redshift_value)

        new_wavelength = list(map(doppler_shift, wavelength_data))

        return new_wavelength   # (array)  -->  The new array of wavelenghts after the doppler shift





    def findPeaks(self):
        """
        Finds the location of n-prominences in the spectrum, given a number n in the parameter file
        """
        ## Locate peaks

        peaks, feats = signal.find_peaks(self.spectrum[:,1], distance = 5, prominence=15, width=15)

        print(peaks)
        # Filter out any unwanted cosmic rays or false peaks that are too flat
        ind_to_delete = []
        for p_ind in peaks:
            y1, y2, y3 =  self.spectrum[:,1][p_ind-1],  self.spectrum[:,1][p_ind],  self.spectrum[:,1][p_ind+1]

            y_avgclose = (  (y2-y1) + (y2-y1)  )  /2     # Find average y-difference between the peak and surronding points

            if y_avgclose > y3/5:        # If the spike is far too sharp, cut it out
                ind_to_delete.append(np.where(peaks == p_ind))
        peaks = np.delete(peaks, np.array(ind_to_delete))

        # Filter out remaining peaks via signal to noise
        snr_list = []        # List to fill with signal to noise ratios
        for p_ind in peaks:
            signal_array = self.spectrum[:,1][p_ind-15:p_ind+15]
            snr = np.mean(signal_array) / np.std(signal_array)          # use mean/std for the signal-to-noise ratio
            snr_list.append(snr)

        ## Trim list until criteria is met:
        if len(snr_list) < self.prom_nums:  # Not enough peaks found
            print(cl.Fore.YELLOW + "Not enough prominences found! Try decreasing number of prominences to fit "+
                    "and check quality of data!" + cl.Fore.WHITE)
                   
        while len(snr_list) > self.prom_nums:
            max_ind = snr_list.index(max(snr_list))
            snr_list.remove(max(snr_list))
            peaks = np.delete(peaks, max_ind)
        
        ## Plot data with peak locations marked:
        if self.plotSetting == 'All':
            data = np.array([[self.spectrum[:,0][i], self.spectrum[:,1][i]] for i in peaks])
            plt.plot(self.spectrum[:,0], self.spectrum[:,1])
            plt.plot(data[:,0], data[:,1], 'ro')
            plt.xlabel("Wavelength (Ang)")
            plt.ylabel("Flux")
            plt.title("Peak Location")
            plt.show()

        self.peak_ind = peaks       # Define the indices of the peaks

        






    def findGauss(self):
        """
        Given the peak locations from 'self.findPeaks', identify data points used for Gaussian fitting
        """

        self.gaus_inds = []      # (list)  -->  Fill with lists of the gaussian indices 

        ## Loop through each peak in the emission spectrum
        for i, p_ind in enumerate(self.peak_ind):
            half_width = signal.peak_widths(self.spectrum[:,1], [p_ind], rel_height=0.5)[0][0]   /2     # (array) --> half width of peaks
            # MgII and H-Beta lines
            if i in (0,1):
                # append a smaller array of shape --> [ [index-1,....,index-n] ]  of the gaussian
                small_width = int(2.5*half_width)   # Gives the half width of the smaller gaussian in terms of indices
                large_width = int(7.5*half_width)   # Width of larger gaussian, in terms of indices (i.e. list of indices, not wavelengths!)
                small_list = list(range(p_ind-small_width, p_ind+small_width+1,1))
                large_list = list(range(p_ind-large_width,p_ind-small_width+1,1))
                large_list = list(range(p_ind-large_width, p_ind-small_width+1)) + list(range(p_ind+small_width, p_ind+large_width+1))
                self.gaus_inds.append( small_list )
                self.gaus_inds.append(large_list)

            # OIII line
            elif i == 2:
                small_width = int(2.5*half_width)   # Gives the half width of the smaller gaussian in terms of indices
                small_list = list(range(p_ind-small_width, p_ind+small_width+1,1))
                self.gaus_inds.append( small_list )
        print(self.gaus_inds)
            
        ## Plotting the points marked for Gaussian fitting:
        if self.plotSetting == "All":
            plt.plot(self.spectrum[:,0], self.spectrum[:,1])
            for c, idx_list in enumerate(self.gaus_inds):
                wl,f = [self.spectrum[:,0][i] for i in idx_list], [self.spectrum[:,1][i] for i in idx_list]
                if (c % 2)==0:
                    plt.plot(wl,f,'g.')
                else:
                    plt.plot(wl,f,'m.')
            plt.show()






    def roughContFit(self):
        """
        Provides a rough fit of the continuum after removing the Gaussians marked in 'self.findGauss'
        """
        ## Trim data exclude marked peaks:
        peak_pts = []
        for gauss in self.gaus_inds:
            for ind in gauss:
                peak_pts.append(ind)
        cont = np.delete(self.spectrum, peak_pts, axis = 0)

        ## Plot trimmed spectrum:
        if self.plotSetting == 'All':
            plt.plot(cont[:,0],cont[:,1],'b.')
            plt.show()
        

        ## Fit the curve with multiple power laws and store the error arrays:
        fit_errors = []
        params = []


        for ind, wl in enumerate(cont[:,0]):               # go from 50th --> -50th so that there are enough points to fit on either side of break
            if (ind<500) == False and (ind> (len(cont[:,0])-500) ) ==False:
                self.wlBreak = wl

                w1, w2 = cont[:,0][0:ind]  , cont[:,0][ind:len(cont[:,0])]                   # wl and flux to fit on left side of the wavelength break
                f1, f2 = cont[:,1][0:ind]  , cont[:,1][ind:len(cont[:,1])]  

                # Functions take flux/wl rather than flux:
                y1, y2 = f1/w1, f2/w2

                
                # wl to fit on the right side of the wavelength break

                p1_pars, p1_err = optimize.curve_fit(self.power1, w1, y1, p0=[1700, -0.6])

                self.temp_coeff, self.temp_power = p1_pars[0], p1_pars[1]

                p2_pars, p2_err = optimize.curve_fit(self.power2, w2, y2, p0=[-1])

                fit1 = w1*self.power1(w1, p1_pars[0], p1_pars[1])
                fit2 = w2*self.power2(w2, p2_pars[0])

                err = np.sum( (fit1-f1)**2 / fit1 ) + np.sum( (fit2-f2)**2 / fit2)

                params.append([ind, p1_pars,p2_pars, wl])
                fit_errors.append(err)

        ## Plot best fit
        min_ind = fit_errors.index(min(fit_errors))
        best_stats = params[min_ind]

        w1b, w2b = self.spectrum[:,0][0:best_stats[0]]  , self.spectrum[:,0][best_stats[0]:len(self.spectrum[:,0])] 
        self.wlBreak = best_stats[3]
        plt.plot(self.spectrum[:,0], self.spectrum[:,1])
        plt.plot(w1b, w1b*self.power1(w1b, best_stats[1][0], best_stats[1][1]), 'r.' )
        plt.plot(w2b, w2b*self.power2(w2b, best_stats[2][0]), 'r.')
        plt.show()
        

        

            
        




    def power1(self, w, a,b ):
        """
        Fit the continuum with a power low with a coefficient that is free to vary
        """
        return a*w**(-b-2)

    def power2(self,w, b):
        """
        Fit continuum, with a power law that does not have a coefficient that is free to vary
        """
        return  self.wlBreak**(b-self.temp_power)*(self.temp_coeff)*w**(-b-2)




    def fitGauss(self):
        """
        Given the peak locations from "findPeaks" function, fit Gaussians to each peak and store the gaussian data
        """

        ## Fit each prominence of the data:
        self.g_params , self.g_errors = [], []

        ## Loop through each peak in the emission spectrum
        for idx_list in self.gaus_inds:
           
            wl = [self.spectrum[:,0][i] for i in idx_list]
            flux = [self.spectrum[:,1][i] for i in idx_list]

            g_par, g_err = optimize.curve_fit(self.Gaussian, wl, flux, p0=[np.mean(flux), max(flux)-min(flux), np.mean(wl), np.std(flux)])


            self.g_params.append([g_par, wl])
            self.g_errors.append(g_err)


        if self.plotSetting=="All":
            plt.plot(self.spectrum[:,0], self.spectrum[:,1])
            for ind, parm in enumerate(self.g_params):
                color = ['r.-', 'g.-', 'm.-', 'b.-', 'r.-', 'g.-']
                plt.plot(parm[1], self.Gaussian(parm[1], parm[0][0], parm[0][1], parm[0][2], parm[0][3]), color[ind] )
            plt.xlabel("Wavelength (Ang)")
            plt.ylabel("Flux")
            plt.title("Gaussian fits on peaks")
            plt.show()


        
    def Gaussian(self, x, H, amp, mean, std):
        """
        INPUTS:
            
            - x     (array)   : Array of x-values
            - H     (float)   : The height of the Gaussian above 0
            - amp   (float)   : Amplitude of the gaussian (height)
            - mean  (float)   : Mean of the gaussian curve
            - std   (float)   : Standard Deviation of the Gaussian
        
        OUTPUTS:
            - (array)  --> A Gaussian distribution
        """
        return amp * np.exp(-((x-mean)**2) / (2*std**2) ) + H     # Gives a Gaussian distribution




    def fitCont(self):
        """
        Fits the continuum
        """
        f_copy = np.copy(self.spectrum[:,1])        # copy to avoid changing flux
        w_copy = np.copy(self.spectrum[:,0])

        ## Subtract the Gaussians from the continuum:
        all_ind = []

        for ind, gauss in enumerate(self.g_params):
        
            i_vals = [np.where(self.spectrum[:,0] == w)[0][0] for w in gauss[1]]    # trim spectrum to appropriate range of indices
            for i in i_vals:
                all_ind.append(i)
                
        f_copy = np.delete(f_copy, all_ind)
        w_copy = np.delete(w_copy, all_ind)


        exp_params, _ = optimize.curve_fit(self.exponential, w_copy, f_copy, p0=[0.005, 0.005, -30])
        power_params, _ = optimize.curve_fit(self.powerLaw, w_copy, f_copy)

        plt.plot(w_copy,f_copy)
        wspace = np.linspace(w_copy[0], w_copy[len(w_copy)-1], 5000)
        plt.plot(wspace, self.exponential(wspace, exp_params[0], exp_params[1], exp_params[2]), 'g')
        plt.plot(wspace, self.powerLaw(wspace, power_params[0], power_params[1]), 'r')
        plt.show()






    def exponential(self, x, a, b,c):
        """
        Fit the continuum with an exponential
        """
        return a*np.exp(-(b*x+c))

    def powerLaw(self,x,a,b):
        """
        Fit the continuum with a power law
        """
        return a*x**(-b)



if __name__ == "__main__":


    test = CQ('C:\\Users\sikor\OneDrive\Desktop\Research\KansasREU\REU\SpectrumAnalysis\params.txt')

    test.fitData()