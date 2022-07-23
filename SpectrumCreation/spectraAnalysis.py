import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

### Inputs:

# Object Path:
objPath = "C:/Users/sikor/OneDrive/Desktop/Research/KansasREU/REU/cstObject.fits" 

# Transformation path:
transPath = "C:/Users/sikor/OneDrive/Desktop/Research/KansasREU/REU/sens.fits"















### Open the Object and Wavelength Transformation arrays:
objData = fits.open(objPath)
transData = fits.open(transPath)









### Close the files:
objData.close()
transData.close()









