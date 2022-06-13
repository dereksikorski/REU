## This file is analyzing CQ Spectra ##

## Imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
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





    def findPeaks(self):
        """
        Finds the location of n-prominences in the spectrum, given a number n in the parameter file
        """
        ## Locate peaks

        peaks, feats = signal.find_peaks(self.spectrum[:,1], distance = 5, prominence=15, width=15)

        # Filter out any unwanted cosmic rays
        ind_to_delete = []
        for p_ind in peaks:
            y1, y2, y3 = self.spectrum[:,1][p_ind-1],  self.spectrum[:,1][p_ind],  self.spectrum[:,1][p_ind+1]

            y_avg = (  (y2-y1) + (y2-y3)  )  /2     # Find average y-difference between the peak and surronding points

            if y_avg > y2/5:        # If the spike is far too sharp, cut it out
                ind_to_delete.append(np.where(peaks == p_ind))
        peaks = np.delete(peaks, np.array(ind_to_delete))


        data = np.array([[self.spectrum[:,0][i], self.spectrum[:,1][i]] for i in peaks])
        plt.plot(self.spectrum[:,0], self.spectrum[:,1])
        plt.plot(data[:,0], data[:,1], 'ro')
        plt.show()









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







if __name__ == "__main__":

    test = CQ()

    #test.plotData(True)

    test.fitData()
