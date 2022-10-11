import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline



## First, import data in:
filepath = r"C:\Users\sikor\OneDrive\Desktop\Research\KansasREU\REU\NewAnalysis\CQ4472.txt"

data = np.loadtxt(filepath, delimiter=",")      # Load data into numpy array


## Define a gaussian fit:

def Gaussian(x, H, amp, mean, std):
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


plt.plot(data[:,0], data[:,1])
plt.show()

## Find all points in data that have flux below average flux:
avg_flux = np.average(data[:,1])
std_flux = np.std(data[:,1])



## Construct large gaus
large_gaus_data = []

for i, flux in enumerate(data[:,1]):
    if flux <= (avg_flux + 0.5*std_flux):

        large_gaus_data.append(data[i])

large_gaus_data = np.array(large_gaus_data)

plt.plot(data[:,0], data[:,1], "b.")
plt.plot(large_gaus_data[:,0], large_gaus_data[:,1], "g.")
plt.show()


lgaus, lgausError = curve_fit(Gaussian, large_gaus_data[:,0], large_gaus_data[:,1], p0 = [np.mean(large_gaus_data[:,1]),
                     np.max(large_gaus_data[:,1]) - np.min(large_gaus_data[:,1]), np.mean(large_gaus_data[:,0]), np.std(large_gaus_data[:,0])])

plt.plot(data[:,0], data[:,1], 'b.')
plt.plot(large_gaus_data[:,0], large_gaus_data[:,1], 'g.')
x = np.linspace(large_gaus_data[:,0][0], large_gaus_data[:,0][-1])
plt.plot(x, Gaussian(x, *lgaus), 'r')
plt.show()

## Subtract Gaus off of data:
sub_data = np.array([  [data[:,0][i], (data[:,1][i] - Gaussian(data[:,0][i], *lgaus))]     for i in range(len(data))])


plt.plot(sub_data[:,0], sub_data[:,1], 'b.')
plt.show()


## Fit narrow gaus

ngaus, ngausError = curve_fit(Gaussian, sub_data[:,0], sub_data[:,1], p0=[np.mean(sub_data[:,1]),
                        np.max(sub_data[:,1]) -np.min(sub_data[:,1]), np.mean(sub_data[:,0]), np.std(sub_data[:,0]) ])


plt.plot(sub_data[:,0], sub_data[:,1], 'b.')
x = np.linspace(sub_data[:,0][0], sub_data[:,0][-1], 10000)
plt.plot(x, Gaussian(x, *ngaus), 'r')
plt.show()



## Calculate residuals
diff = [sub_data[:,1][i] - Gaussian(sub_data[:,0][i], *ngaus) for i in range(len(sub_data))]

plt.scatter(data[:,0], diff, s=0.4)
plt.axhline(y=0, color="black")
plt.show()