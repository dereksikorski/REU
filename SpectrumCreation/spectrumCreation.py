### This file is meant for importing in .fits files of 2D images and converting them into 1D spectra

### Imports:
import numpy as np
import os
import colorama as cl
from astropy.io import fits
import matplotlib.pyplot as plt
import math


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
        
        # Central Wavelength
        self.centWL = file_values_array[1]

        # Center row of the image to use
        self.center = int(file_values_array[2]) - 1   # NOTE: This is the index of the python array and NOT the physical row.
                                                            # ex: If the center row of the image is 256 (starting at 1), then 
                                                            #     self.center = 255 as python indexing begins at 0, not 1.

        # Binning to use when finding flux:
        self.bin = float(file_values_array[3])

        # Plot setting:
        self.plot = file_values_array[4]

        # Ref Pixel:
        self.refPix = float(file_values_array[5])

        # Ref Wl:
        self.refWl = float(file_values_array[6])



    def Transform(self):
        """
        Transform the data from a 2D Image to a 1D spectrum
        """

        ## Open the .fits files
        self.objData = fits.open(self.objFilePath)

        ## Load in data:
        print(cl.Fore.GREEN + "Object File Format table:" + cl.Fore.WHITE + "\n")
        self.objData.info()


        ## Pull and combine the needed flux values from the image:
        self.pullFlux()


        ## Preform wavelength transformation on pixel data:
        self.wlTransform()


        ## Save file:
        self.writeFile()

        ## Plot:
        if self.plot in ("yes", "Yes", "y", "Y"):
            self.Plot()

        ## Close the .fits files
        self.objData.close()



    def pullFlux(self):
        """
        Pull the flux values from the 
        """

        ## Read in data as np.array objects
        sciData = self.objData["SCI"].data      # (array)  -->  Shape (#_columns, #_rows) containing flux values from 2D image
        varData = self.objData["VAR"].data      # (array)  -->  Shape (#_columns, #_rows) containing variance values from 2D image
        

        ## Pull the data from the desired rows from the data. This is based on the center aperture and binning amount specified in the parameters
        if self.bin == 0:
            # Only use the central row
            usefulSci = [sciData[int(self.center)]]
            usefulVar = [varData[int(self.center)]]
        else:
            # Use rows +/- the binning on either side of the central aperture
            usefulSci = sciData[int(self.center - self.bin): int(self.center + self.bin )]   # Trim the array to only the useful rows
            usefulVar = varData[int(self.center - self.bin): int(self.center + self.bin )] 


        ## Loop through data and find a weighted value for the flux at each pixel value:
        self.fluxVals = []      # List to fill with averaged flux values
        self.varVals = []       # List to fill with averaged variance values
        for px_num in range(len(usefulSci[0])):

            # Fill the lists with flux and var values for each pixel given the binned rows
            fluxes = np.array([ row[px_num] for row in usefulSci ])
            vars = np.array([  row[px_num] for row in usefulVar  ])
                       
            # Check if any of the variance values are 0
            inds = np.where(vars==0)[0]  # List of indices where the variance is 0


            if len(inds) != len(fluxes):
                # If some of the variances are 0, ignore those values and only use values with non-zero variance
                fluxes = np.delete(fluxes, inds)
                vars = np.delete(vars, inds)

                num = fluxes/vars
                denom = 1/vars


                self.fluxVals.append(sum(num)/sum(denom))
                self.varVals.append(1/sum(denom))
            

            else:
                # If all variance values are 0, then do simple unweighted average
                self.fluxVals.append(sum(fluxes)/len(fluxes))
                self.varVals.append(np.var(fluxes))

        

    
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


        ## Given wl of the reference sky emission line, apply a final shift on the entire spectrum to line them up

        dist = dw*(self.refPix - p0) + w0 - self.refWl      # This is how far the spectrum needs to be shifted

        self.wls -= dist


    def writeFile(self):
        """
        Write the calculated values to a .txt file
        """

        # Data to write to the .txt file
        data = np.array([[self.wls[i], self.fluxVals[i], self.varVals[i]] for i in range(len(self.wls))])

        # Create a filename based on the original input file
        newFile = self.objFilePath[:-5]  # Trim ending off of the input file
        newFile += ".txt"

        np.savetxt(newFile, data, delimiter=",")

    def Plot(self):
        """
        Save a plot of the spectra if plot setting specifies
        """
        # Plot full spectrum
        plt.clf()
        plt.plot(self.wls, self.fluxVals)

        plt.xlabel(r"Wavelength ($\AA$)")
        plt.ylabel("Flux")
        plt.title(r"CQ 4472 at $\lambda_{central}$=" f"{self.centWL}" +  r"$\mu$m")

        plt.savefig(self.objFilePath[:-5] + ".png")

        # Plot only MgII line
        plt.clf()
        plt.plot(self.wls, self.fluxVals)
        plt.xlim([8000,8600])
        plt.ylim([0,4e-18])

        plt.annotate("MgII Emission", (8200, 3.5e-18))
        plt.arrow(8270, 3.3e-18, -50, -0.5e-18, width = 1e-19, head_length=10)

        plt.xlabel(r"Wavelength ($\AA$)")
        plt.ylabel("Flux")
        plt.title(r"MgII Emission of CQ 4472 at $\lambda_{central}$=" f"{self.centWL}"+ r"$\mu$m")
        plt.savefig(self.objFilePath[:-5] + "_MgII.png")


    def truncate(self, num, n):
        '''
        Input:

            - num = Number to be truncated
            - n = Number of decimals to truncate to
        
        Output:

            - Truncated number
        '''
        stepper = 10.0 ** n
        return math.trunc(stepper * num) / stepper


if __name__ == "__main__":



    test = Data("C:\\Users\sikor\OneDrive\Desktop\Research\KansasREU\REU\SpectrumCreation\SCparams.txt")

    test.Transform()
