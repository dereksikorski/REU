### This file is meant for importing in .fits files of 2D images and converting them into 1D spectra

### Imports:
import numpy as np
import os
import colorama as cl
from astropy.io import fits
import matplotlib.pyplot as plt


class spec:

    def __init__(self, param_path = None):
        """
        Create a "Data" object that holds information about the given science object

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
        
        # Path to directory containing data files
        self.specFolder = file_values_array[0]
        
        # Name of output file
        self.outName = file_values_array[1]

        # Whether or not to save plot
        self.plotOption = file_values_array[2]   


    def Combine(self):
        """
        Transform the data from a 2D Image to a 1D spectrum
        """
        # Read in the data from the input directory
        self.readData()

        # Combine data
        self.Combo()


    def readData(self):
        """
        Read and store data from the input spectra
        """
        ## Loop through the files in the directory and store data
        self.data_dict = {}
        for i, filename in enumerate(os.listdir(self.specFolder)):
            self.data_dict[str(i)] = np.loadtxt(self.specFolder+"\\" + filename, float, delimiter=",")
        print(self.data_dict)


    def Combo(self):
        """
        Combine the data using a weighted average
        """
        # Unpack the dictionary into a more usable form
        wls = []
        flux = []
        vars = []
        for key in self.data_dict.keys():
            data = self.data_dict[key]
            wls.append([d[0] for d in data])
            flux.append([d[1] for d in data])
            vars.append([d[2] for d in data])

        # Combine:
        


if __name__ == "__main__":



    test = spec("C:\\Users\sikor\OneDrive\Desktop\Research\KansasREU\REU\SpectrumCombination\comParams.txt")

    test.Combine()
