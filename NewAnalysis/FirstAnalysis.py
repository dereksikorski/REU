import os
import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt


## First, import data in:
filepath = r"C:\Users\sikor\OneDrive\Desktop\Research\KansasREU\REU\SpectrumAnalysis\CQ507Trimmed.txt"

data = np.loadtxt(filepath, delimiter=",")


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


## Plot data to fit:
# plt.plot(data[:,0], data[:,1])
# plt.show()


## Defne index spots:
Li = 500
Ri = 675

## Plot narrow Gaussian

Firstwl = np.copy(data[:,0][Li:Ri])
Firstflux = np.copy(data[:,1][Li:Ri])


# Re-shift data:
FirstShift = np.mean(Firstwl)
ShiftedFirstwl = Firstwl - FirstShift
Firstflux *= 1E17

#FilteredFirst = gaussian_filter1d(Firstflux, sigma=5)

# Try to fit the narrow emission part with a gaussian
nGaus, nGaus_err = curve_fit(Gaussian, ShiftedFirstwl, Firstflux ,
                            p0 = [np.mean(Firstflux), np.max(Firstflux) - np.min(Firstflux), np.mean(Firstflux), np.std(Firstflux)])
    

print(f"Best fit params for narrow Gaussian: {nGaus}\n")

# Plot initial gauss fit:
# plt.plot(ShiftedFirstwl,Firstflux)
# plt.plot(ShiftedFirstwl, Gaussian(ShiftedFirstwl, nGaus[0], nGaus[1], nGaus[2], nGaus[3]), 'r')
# plt.show()



## Plot wide Gaussian

iToDelete = range(Li,Ri,1)    # indicies to delete (don't use in second fit)

Secondwl = np.delete( np.copy(data[:,0]), iToDelete)
SecondFlux = np.delete( np.copy(data[:,1]), iToDelete)


# Re-shift data:
SecondShift = np.mean(Secondwl)
ShiftedSecondwl = Secondwl- SecondShift
SecondFlux *= 1E17

# Pass through Gauss filter to smooth out noise

wGaus, wGaus_err = curve_fit(Gaussian, ShiftedSecondwl, SecondFlux,
                    p0 = [np.mean(SecondFlux), np.max(SecondFlux) - np.min(SecondFlux), np.mean(SecondFlux), np.std(SecondFlux)])

print(f"Best fit params for wide Gaussian: {wGaus}\n")


# Plot initial gauss fit:
# plt.plot(ShiftedSecondwl,SecondFlux)
# plt.plot(ShiftedSecondwl, Gaussian(ShiftedSecondwl, wGaus[0], wGaus[1], wGaus[2], wGaus[3]), 'r')
# plt.show()





## Create my total fit:



sect1 = [Gaussian(wl-SecondShift, wGaus[0], wGaus[1], wGaus[2], wGaus[3]) for wl in data[:,0][0:Li]]
sect2 = [Gaussian(wl-FirstShift, nGaus[0], nGaus[1], nGaus[2], nGaus[3]) for wl in data[:,0][Li:Ri]]
sect3 = [Gaussian(wl-SecondShift, wGaus[0], wGaus[1], wGaus[2], wGaus[3]) for wl in data[:,0][Ri:len(data[:,0])]]

fitFlux = np.concatenate((sect1,sect2,sect3))

# for i, wl in enumerate(data[:,0]):

#     if i <= Li:
#         nflux = Gaussian(wl-SecondShift, wGaus[0], wGaus[1], wGaus[2], wGaus[2])

#     elif i > Ri:
#         nflux = Gaussian(wl-SecondShift, wGaus[0], wGaus[1], wGaus[2], wGaus[2])
#     else:
#         nflux = Gaussian(wl-FirstShift, nGaus[0], nGaus[1], nGaus[2], nGaus[2])

#     fitFlux.append(nflux)


plt.plot(data[:,0], data[:,1]*1E17, 'b')
plt.plot(data[:,0], fitFlux, 'r')
plt.show()


## Find Residuals based on this Computation

diff = fitFlux - data[:,1]*1E17

plt.scatter(data[:,0], diff, s=0.4)
plt.axhline(y=0, color="black")
plt.show()