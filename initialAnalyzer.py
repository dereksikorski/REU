## This file is analyzing CQ Spectra ##

## Imports
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

        ## Read in spectrum data:
        self.data = pd.read_csv(file_values_array[0])     # Data from the .csv file
        
        self.spectrum = [] 
        for index in list(self.data.index.values):
            self.spectrum.append([self.data["Wavelength"][index], self.data["Flux"][index]])
        self.spectrum = np.array(self.spectrum)

        self.redshift = float(file_values_array[1])
        
        self.prom_nums = int(file_values_array[2])

        self.plotSetting = file_values_array[3]


    def plotData(self, rs = True):
        """
        Plots the data with a redshift, unless specified

        INPUTS:
            - rs (bool)    : If True, applies redshift before plotting. If False, plots without redshift
        """
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
        self.fitGauss()

        ## Fit Continuum
        self.fitCont()




    def findPeaks(self):
        """
        Finds the location of n-prominences in the spectrum, given a number n in the parameter file
        """
        ## Locate peaks

        peaks, feats = signal.find_peaks(self.spectrum[:,1], distance = 5, prominence=15, width=15)

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

        # Trim list until criteria is met:
        if len(snr_list) < self.prom_nums:  # Not enough peaks found
            print(cl.Fore.YELLOW + "Not enough prominences found! Try decreasing number of prominences to fit "+
                    "and check quality of data!" + cl.Fore.WHITE)
        while len(snr_list) > self.prom_nums:
            max_ind = snr_list.index(max(snr_list))
            snr_list.remove(max(snr_list))
            peaks = np.delete(peaks, max_ind)
        
        ## Plot:
        if self.plotSetting == 'All':
            data = np.array([[self.spectrum[:,0][i], self.spectrum[:,1][i]] for i in peaks])
            plt.plot(self.spectrum[:,0], self.spectrum[:,1])
            plt.plot(data[:,0], data[:,1], 'ro')
            plt.xlabel("Wavelength (Ang)")
            plt.ylabel("Flux")
            plt.title("Peak Location")
            plt.show()

        self.peak_ind = peaks       # Define the indices of the peaks

        

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


    def fitGauss(self):
        """
        Given the peak locations from "findPeaks" function, fit Gaussians to each peak and store the gaussian data
        """

        ## Fit each prominence of the data:
        self.g_params , self.g_errors = [], []
        for i in self.peak_ind:
            half_width = signal.peak_widths(self.spectrum[:,1], [i], rel_height=0.5)[0][0]   /2
            
            wl1, f1 = self.spectrum[:,0][int(i- (2.5)*half_width) : int(i+ (2.5)*half_width)+1], self.spectrum[:,1][int(i- (2.5)*half_width) : int(i+ (2.5)*half_width)+1]
            wl2, f2 = self.spectrum[:,0][int(i- (1.0) * half_width) :int( i+(1.0)*half_width)+1], self.spectrum[:,1][int(i- (1.0)*half_width) : int(i+ (1.0)*half_width)+1]


            params1, g_errors1 = optimize.curve_fit(self.Gaussian, wl1, f1, p0 = [np.mean(f1), max(f1)-min(f1), np.mean(wl1), np.std(f1)])
            params1 = [params1, wl1]
            self.g_params.append(params1)
            self.g_errors.append(g_errors1)

            params2, g_errors2 = optimize.curve_fit(self.Gaussian, wl2, f2, p0 = [np.mean(f2), max(f1)-min(f1), np.mean(wl2), np.std(f2)])
            params2 = [params2, wl2]
            self.g_params.append(params2)
            self.g_errors.append(g_errors2)

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
            if ind % 2 != 0:            # Only work with even numbers as second curve is larger curve
                i_vals = [np.where(self.spectrum[:,0] == w)[0][0] for w in gauss[1]]    # trim spectrum to appropriate range of indices
                for i in i_vals:
                    all_ind.append(i)
                
        f_copy = np.delete(f_copy, all_ind)
        w_copy = np.delete(w_copy, all_ind)

        T, c = optimize.curve_fit(self.Planck, w_copy, f_copy, p0=500)

        if self.plotSetting in ['All', 'Result']:    
            plt.clf()
            plt.plot(self.spectrum[:,0], self.spectrum[:,1])
            ws = np.linspace(self.spectrum[:,0][0], self.spectrum[:,0][len(self.spectrum[:,0])-1], 5000)
            plt.plot(ws, self.Planck(ws,T), 'r')
            plt.plot(ws, self.Planck(ws, T+200), 'g')
            plt.plot(ws, self.Planck(ws, T-200), 'm')
            plt.legend(['Data', f'Best-fit Planck Funciton: {round(T[0], 2)}K', f'Planck Function: {round(T[0]+200,2)}K', f'Planck Function: {round(T[0]-200,2)}K'])
            plt.xlabel("Wavelength (Ang)")
            plt.ylabel("Flux")
            plt.title("Fitting the emission spectrum with Planck Curves")
            plt.show()


    def Planck(self, w, T):
        """
        Planck Function
        """
        c1= 2 * 6.626e-34 * 299792458
        c2 = 6.626e-34 * 299792458 /1.381e-23
        W = w * 10e-10
        B = (c1/W**5) * 1/((np.exp(c2/(W*T))) -1 )
        return B



if __name__ == "__main__":

    test = CQ()

    test.fitData()