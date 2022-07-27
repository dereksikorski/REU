### This file is meant for importing in .fits files of 2D images and converting them into 1D spectra

### Imports:
import numpy as np
import os
import colorama as cl
from astropy.io import fits
import matplotlib.pyplot as plt


class Data:

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
        
        # Path to the object .fits files
        self.objFilePath = file_values_array[0]
        
        # Path to the sensitivity file
        self.sensFilePath = file_values_array[1]

        # Center row of the image to use
        self.center = float(file_values_array[2]) - 1   # NOTE: This is the index of the python array and NOT the physical row.
                                                            # ex: If the center row of the image is 256 (starting at 1), then 
                                                            #     self.center = 255 as python indexing begins at 0, not 1.

        # Binning to use when finding flux:
        self.bin = float(file_values_array[3])



    def Transform(self):
        """
        Transform the data from a 2D Image to a 1D spectrum
        """

        ## Open the .fits files
        self.objData = fits.open(self.objFilePath)
        self.sensData = fits.open(self.sensFilePath)

        ## Load in data:
        print(cl.Fore.GREEN + "Object File Format table:" + cl.Fore.WHITE + "\n")
        self.objData.info()

        print()
        print(cl.Fore.GREEN + "Sensitivity Function File Format table:" + cl.Fore.WHITE + "\n")
        self.sensData.info()

        ## Pull and combine the needed flux values from the image:
        self.pullFlux()

        ## Preform wavelength transformation on pixel data:
        self.wlTransform()

        ## Plot final data:
        self.prePlot()

        ## Apply flux calibration:
        self.fluxCal()


        ## Trimmed final plot:
        self.trimmedPlot()

        ## Trimmed to MgII:
        self.trimmedtrimmedPlot()
        

        self.sensPlot()

        ## Close the .fits files
        self.objData.close()
        self.sensData.close()




    def pullFlux(self):
        """
        Pull the flux values from the 
        """

        ## Read in data as np.array objects
        sciData = self.objData["SCI"].data      # (array)  -->  Shape (#_columns, #_rows) containing flux values from 2D image
        varData = self.objData["VAR"].data      # (array)  -->  Shape (#_columns, #_rows) containing variance values from 2D image
        
        if self.bin == 0:
            usefulSci = [sciData[int(self.center)]]
            usefulVar = [varData[int(self.center)]]
        else:
            usefulSci = sciData[int(self.center - self.bin): int(self.center + self.bin )]   # Trim the array to only the useful rows
            usefulVar = varData[int(self.center - self.bin): int(self.center + self.bin )] 



        ## Loop through data and find a weighted value for the flux at each pixel value:

        self.fluxVals = []
        for px_num in range(len(usefulSci[0])):

            # Fill the lists with flux and var values for each pixel given the binned rows
            fluxes = np.array([ row[px_num] for row in usefulSci ])
            vars = np.array([  row[px_num] for row in usefulVar  ])
                       
            # Check if any of the variance values are 0
            inds = np.where(vars==0)[0] # List of indices where the variance is 0

            if len(inds) == 0:

                num = fluxes/vars
                denom = 1/vars

        

                self.fluxVals.append(sum(num)/sum(denom))
            
            elif len(inds) == len(fluxes):
                self.fluxVals.append(sum(fluxes)/len(fluxes))

            else:
                fluxes = np.delete(fluxes, inds)
                vars = np.delete(vars, inds)

                num = fluxes/vars
                denom = 1/vars


                self.fluxVals.append(sum(num)/sum(denom))

        

    
    def wlTransform(self):
        """
        Pull information from the image header and apply transformation to convert from pixels --> wl
        """

        ## Define a list of the pixel values:
        px_range  = np.array([p+1 for p in range(0,len(self.fluxVals))] ) # Need first pixel to be 1 not 0

        ## Pull transform info from header:
        p0 = self.objData["SCI"].header["CRPIX1"]  # Reference pixel number
        
        w0 = self.objData["SCI"].header["CRVAL1"]  # Wavelength of the initial pixel

        dw = self.objData["SCI"].header["CD1_1"]   # Change in wavelength between pixels


        ## Transfrom the pixel array to a wavelength array:

        self.wls = dw*(px_range - p0) + w0  # apply transformation -->   w = dw*(p - p0) + w0


            

    def fluxCal(self):
        """
        Calibrate the flux values based on the given sensitivity values
        """

        sens_vals = self.sensData[0].data

        self.fluxVals = self.fluxVals * sens_vals       # Multiply flux vals by sensitivity function to correct


    def prePlot(self):
        """
        Plot the final data
        """

        plt.plot(self.wls, self.fluxVals)

        plt.xlabel("Wavelength (Angstrom)")
        plt.ylabel("Flux")
        plt.xlim([7200,7800])
        plt.ylim([0,0.75e-16])

        plt.savefig("735CQ507.png")

    def trimmedPlot(self):
        """
        Plot the final trimmed data
        """
        plt.clf()
        plt.plot(self.wls, self.fluxVals)

        plt.xlabel("Wavelength (Angstrom)")
        plt.ylabel("Flux")
        plt.xlim([6000, 8300])
        plt.ylim([0,0.25e-14])

        plt.savefig("735CQ507_trimmed.png")
    
    def trimmedtrimmedPlot(self):
        """
        Trim to interesting region
        """
        plt.clf()
        plt.plot(self.wls, self.fluxVals)

        plt.xlabel("Wavelength (Angstrom)")
        plt.ylabel("Flux")
        plt.xlim([7200, 7800])
        plt.ylim([0,0.25e-14])

        plt.savefig("735CQ507MgII.png")

    def sensPlot(self):
        """
        Plot sens
        """
        plt.clf()
        sens_vals = self.sensData[0].data

        plt.plot(self.wls, sens_vals)
        plt.savefig("735sens.png")


if __name__ == "__main__":



    test = Data("C:\\Users\sikor\OneDrive\Desktop\Research\KansasREU\REU\SpectrumCreation\params.txt")

    test.Transform()
