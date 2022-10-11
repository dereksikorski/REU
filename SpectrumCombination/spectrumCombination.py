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

        self.bin = float(file_values_array[3])  # Bin size in Angstroms


    def Combine(self):
        """
        Transform the data from a 2D Image to a 1D spectrum
        """
        # Read in the data from the input directory
        self.readData()

        # Do actual Combo
        self.Combo()

        # Plot the combined data
        self.plotCombo()

        self.writeFile()

    def readData(self):
        """
        Read and store data from the input spectra
        """
        ## Loop through the files in the directory and store data
        self.data_dict = {}
        for i, filename in enumerate(os.listdir(self.specFolder)):
            self.data_dict[str(i)] = np.loadtxt(self.specFolder+"\\" + filename, float, delimiter=",")


        # Unpack the dictionary into a more usable form
        self.wls = []
        self.flux = []
        self.vars = []
        for key in self.data_dict.keys():
            data = self.data_dict[key]
            for point in data:
                self.wls.append(point[0])
                self.flux.append(point[1])
                self.vars.append(point[2])




    def Combo(self):
        """
        Do an interperlation combo
        """
        # Find High and low wl
        low_wl = min(self.wls)
        high_wl = max(self.wls)

        # Create lists to fill
        self.comboWL = []
        self.comboFlux = []
        self.comboVar = []


        wl_bin = low_wl
        while wl_bin <= high_wl:
            wl_compare = np.array(self.wls) - wl_bin

            # Find values greater than 0:
            inds = np.where(wl_compare < 0)[0]
            
            wl_compare = np.delete(wl_compare, inds)
            wls = np.delete(self.wls, inds)
            fluxes = np.delete(self.flux, inds)
            vars = np.delete(self.vars, inds)

            # Find values less than bin size:
            inds = np.where(wl_compare > self.bin )[0]
            
            wls = np.delete(wls, inds)
            fluxes = np.delete(fluxes, inds)
            vars = np.delete(vars, inds)

            # Skip Chip gaps:
            if len(self.comboFlux) != 0:
                inds = np.where(fluxes < self.comboFlux[-1]/1.75)[0]
            else:
                inds = []

            wls = np.delete(wls, inds)
            fluxes = np.delete(fluxes, inds)
            vars = np.delete(vars, inds)
            
            try:
                self.comboWL.append(sum(wls)/len(wls))
                self.comboFlux.append(self.weightedAvg(fluxes, vars))
                self.comboVar.append(1/sum(1/vars))
            except ZeroDivisionError:
                break

            wl_bin += self.bin


    def weightedAvg(self,vals,vars):
        """
        Compute the weighted averages of arrays given the values and variances
        """
        num = sum(vals/vars)
        denom = sum(1/vars)

        return num/denom

            
    def plotCombo(self):
        """
        Plot the combined data
        """
        plt.plot(self.comboWL,self.comboFlux)

        plt.title(f"Combined Spectrum for {self.outName}")
        plt.xlabel(r"Wavelength ($\AA$)")
        plt.ylabel("Flux")
        # plt.annotate(text="MgII Emission", xy=(8000, 5e-17))
        # plt.arrow(7900, 5e-17, -450, 0, width = 1e-18, head_length=50)


        plt.savefig(self.specFolder + "\\" + self.outName + ".png")


    def writeFile(self):
        """
        Write the resulting files to a .txt files
        """

        # Data to write to the .txt file
        data = np.array([[self.comboWL[i], self.comboFlux[i], self.comboVar[i]] for i in range(len(self.comboWL))])

        # Create a filename based on the original input file
        newFile = self.specFolder + "\\" + self.outName  # Trim ending off of the input file
        newFile += ".txt"

        np.savetxt(newFile, data, delimiter=",")


if __name__ == "__main__":

    test = spec("C:\\Users\sikor\OneDrive\Desktop\Research\KansasREU\REU\SpectrumCombination\comParams.txt")

    test.Combine()
