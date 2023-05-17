###################
# Notes:
#       - This pipeline is meant to reduce Gemini data one waveband at a time
#       - Basic directory structure should be:
#           - input dir contains 
#                  a) A folder called "raw" containing all raw data for given observation night
#                     as well as a "Bad Mask Profile" (bpm) that matches the one in the bias 
#                     combination step below.
#                  b) A temp. folder will be created in the main directory that contains all temporary
#                     files created by gemini. This will automatically be deleted at the end of the pipeline
#                  c) The resulting image and log files will be stored in the main directory
#
#       - Before using, MUST fill out inputs below!!!
#
#       - This pipeline follows figure 2 on 
#         https://noirlab.edu/science/programs/csdc/usngo/gmos-cookbook/Processing/Concept_Overview.html
#         and each step in the process is numbered according to its position on the table.
###################
##
##
##
##
##
##
###############
# Define a task for the cosmic ray rejection
task lacos_spec = /N/u/dmsikors/Carbonate/Desktop/REUGemini/lacos_spec.cl

### Inputs:
       ## LINE 1: Enter the ABSOLUTE path of the desired night. The pasted
               ## directory should have a folder called "raw" that contains ALL raw data
       ## LINE 2: Enter the central wavelength of the images to reduce WITHOUT a "."
       ## LINE 3: Enter the IRAF name of the standard star being used, found at:
               ## https://noirlab.edu/science/programs/csdc/usngo/gmos-cookbook/Processing/Supplement.html#stdstar-list

cd /N/u/dmsikors/Carbonate/Desktop/REUGemini/Aug7/     # Enter path of data
string filterWl = "715"                      # Central wavelength of data
string stdName = "eg131"



################






###########   0   #############
### Create observation log and temp folder:
    ## Create an observation log of all data in given observing night and make
    ## a temp folder to store files away from the raw data.
    ## NOTE the temp folder will automatically be deleted at the end of the code,
    ## but this option can easily be removed at the end of the code.

    
mkdir temp      # Create a temp folder to store all temp data to be deleted

## Make a observation log in main directory
del "obsLog.txt"
cd raw

string fieldNames = "# File         ObsDate   ObsTime   Object  ObsType   ObsClass  CcdBin  ROI Filter1    Filter2    Disperser   SlitName  SlitType    Rotator CenWave T_exp   Airmass"
print (fieldNames, >"obsLog.txt")   # Create an observation log .txt file                     
string keyWords = "$I Date-Obs Time-Obs Object ObsType ObsClass Ccdsum DetRO1ys Filter1 Filter2 Grating MaskName MaskType PA GrWlen ExpTime Airmass"


hselect ("*.fits[1,inherit=yes]", keyWords, "obstype != 'acq'", >"obsLog.txt")
rename("obsLog.txt", "../")

###############################
###############################





##########  0 part II  #############
### Create lists of the relavant data:
    ## Make lists of different data types based on wavelength, and move to the
    ## temp directory. From here, all processes will be executed in temp rather
    ## than raw, so that the raw images will not be touched.
    

## Bias frames:
string biasCriteria = "obstype=='BIAS' && obsclass=='dayCal'"  # Criteria for hselect next line
hselect ("*.fits[1,inherit=yes]", "$I", biasCriteria, > "bias_tmp.txt") # Select names of bias images in "raw" and store in list
 
del "biasFiles.txt"     # Deletes file if repeated. Will give warning if file does not exist
gemextn ("@bias_tmp.txt", omit="index,kernel", outfile="biasFiles.txt")       # Copy the temporary file, removing uneeded info
del "bias_tmp.txt"
rename ("biasFiles.txt", "../temp/")     # Rename file so it is moved into temp directory

## Flats of given waveband:
string flatCriteria = "i_title=='GCALflat' && GrWlen=="   # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (flatCriteria // (filterWl//".")) , >> "flt_tmp.txt")   # Select flat field images based on wavelength

del "flatFiles.txt"  # Deletes file if repeated  
gemextn ("@flt_tmp.txt", omit="index,kernel", outfile= "flatFiles.txt" )  #Copy flats file, removing uneeded data
del "flt_tmp.txt" 
rename ("flatFiles.txt"  , "../temp/")     # Rename file to move into temp directory

## Arcs:
string arcCriteria = "obstype=='ARC' && GrWlen=="   # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (arcCriteria // (filterWl//".")) , > "arc_tmp.txt")     # Select arc files for given wl

del "arcFiles.txt"
gemextn ("@arc_tmp.txt", omit="index,kernel", outfile="arcFiles.txt")  # Trim the file names within the arc list
del "arc_tmp.txt"
rename ("arcFiles.txt", "../temp/")

## Science frames:
string objCriteria = "obstype=='OBJECT' && obsclass =='science' && GrWlen=="    # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (objCriteria // (filterWl//".")), >"obj_tmp.txt")

del "objFiles.txt"
gemextn ("@obj_tmp.txt", omit="index,kernel", outfile="objFiles.txt")
del "obj_tmp.txt"
rename ("objFiles.txt","../temp/")

## Standard Star:
string starCriteria = "obsclass=='partnerCal' && obstype == 'OBJECT'"
hselect ("*.fits[1,inherit=yes]", "$I", starCriteria, >"std_temp.txt")

del "stdFiles.txt"
gemextn ("@std_temp.txt", omit="index,kernel", outfile="stdFiles.txt")
del "std_temp.txt"
rename ("stdFiles.txt", "../temp/")

################################
################################


cd ..
cd temp
mkdir Logs      # Create a directory within temp to store relevant log files



#############  1  ##############
### Append MDF:
    ## Insert important metadata into file header and append appropriate MDF
    ## (Mask Definition File) as a table extension. The MDF relays info about 
    ## how the ccd is illuminated and how the fibers receive data given the
    ## imaging type in the header.
    
unlearn gprepare

gprepare.outpref = ""       # Rename the files after preparation to the original names
gprepare.rawpath = "../raw"
gprepare.logfile = "./Logs/gprepLog.txt"
    
## Bias:
gprepare ("@biasFiles.txt", outimages="@biasFiles.txt", fl_addm+)

## Flats:
gprepare ("@flatFiles.txt", outimages="@flatFiles.txt", fl_addm+)

## Arcs:
gprepare ("@arcFiles.txt", outimages="@arcFiles.txt", fl_addm+)

## Object:
gprepare ("@objFiles.txt", outimages="@objFiles.txt", fl_addm+)

## Std Star:
gprepare ("@stdFiles.txt", outimages="@stdFiles.txt", fl_addm+)

##############################
##############################





##############  2  ####################
### Create a bias Master Calibration image:
    ## Def:
        # Images have a "dark fixed-pattern noise" as a result of the manufacturing
        # of the CCD. This pattern is imprinted in every image the camera takes.
    ## Steps:
        # 1) Overscan correct and subtract overscan region of each image
        # 2) Combine each image to make a Master Bias frame
    ## Notes:
        # 1) Seems like this step does NOT need to be done interactively as the
        #    fit is only done for the overscan correction step.
        # 2) Applied order 11 to the fit after interactive examination as I felt it
        #    more accurately modeled the overscan regions

unlearn gbias

gbias.logfile = "./Logs/gbiasLog.txt"
gbias.order = 11 # After examination, this seems like a better fit than default.

gbias ("@biasFiles.txt", "MCbias.fits", fl_vardq+, fl_inte-)     # Combine and perform overscan corrections
    
    
######################################
######################################


############  3  #######################
### Create flat field Master Calibration
    ## Def:
        # CCDs have "dead pixels" and "hot pixels" where the sensitivity of the
        # pixel is higher/lower respectively. The Flat Filed image is an image of
        # a bright screen that is used to identify the pixel sensitivity.
    ## Steps:
        # 1) Apply overscan correction to each flat field frame
        # 2) Apply bias correction to each flat field frame using MasterCal Bias
        # 3) Combine the flat field frames together into one image
        # 4) Fit a response function to normalize the flat field frame for each
        #    extension separately
    ## Notes:
        # 1) Important to use fl_detec+, which means the extensions are fit
        #    individually rather than fit after being mosaicked. This then applies
        #    the fit to the un-mosaicked data images.
        # 2) Although it is recommended to fit the response functions interactively
        #    I'm not entirely sure how to do it.


unlearn gsflat
unlearn gemextn
gsflat.logfile = "./Logs/gsflatLog.txt"

gsflat ("@flatFiles.txt", "MCflat.fits", bias="MCbias.fits", \
   fl_vardq+, fl_detec+, fl_oversize-, order="13,11,28", fl_inte+)
   
######################################
######################################




###########  4  #####################
### Basic reductions:
    ## Def:
        # Apply initial reductions using the MasterCal calibrations to the Object,
        # Arc, and Std Star images
    ## Steps:
        # 1) Overscan correct all images
        # 2) Bias correct all images
        # 3) Flat-field correct ONLY the object images
        # 4) Apply the LA-Cosmic ray reject to the object images.
        # 5) Run gsappwave on reduce images to give initial wavelength calibration
    ## Notes:
        # 1) The cosmic ray rejection algorithm can be a bit finicky from
        #    time-to-time. As this is being used for bright quasars, the cutoff needs
        #    to be set high because the object is very bright. This means some 
        #    cosmic rays could slip through the algorithm.

unlearn gsreduce
unlearn gemextn

gsreduce.logfile= "./Logs/ReductionsLog.txt"
gsreduce.bias = "MCbias.fits"
gsreduce.cr_objlim = 1.    # Object has to be this times as bright as surronding objects to be flagged as a CR

## Object
gsreduce ("@objFiles.txt", outpref="o", flatim="MCflat.fits", \
          fl_crsp+, fl_fixpix-, fl_oversize-, fl_vardq+, fl_fulldq+, fl_inte+, ovs_fli-)
          
del "objFiles.txt"
hselect ("o*.fits", "$I", "obstype != 'acq'", >"objFiles.txt")

          
## Arc
gsreduce ("@arcFiles.txt", outpref="a", \
   fl_fixpix-, fl_flat-, fl_oversize-, fl_inte+, ovs_fli-)

del "arcFiles.txt"
hselect ("a*.fits", "$I", "obstype != 'acq'", >"arcFiles.txt")

## Std Star
gsreduce ("@stdFiles.txt", outpref="s", \
           fl_fixpix-, fl_flat-, fl_oversize-,fl_inte+, ovs_fli-)
          
del "stdFiles.txt"
hselect ("s*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")
################################
################################


############  5  ##############
### Combine same frames:
    ## Def:
        # Do a basic combination of images with the same data, grating, ROI, etc.
    ## Notes:
        # 1) It's important the images be taken at nearly the same time and at the
        #    same grating, central wavelength, ROI, ccdsum, etc. 
        # 2) It's recommended to check the images before and after a combination
        #    to make sure the data is sufficiently close.

#unlearn gemcombine
#gemcombine.logfile = "./Logs/Combination.txt"


#gemcombine ("@objFiles.txt", "Object.fits", fl_vardq+, fl_dqprop+)

rename ("o*.fits", "Object.fits")
#############################
############################



##########  6  ################
### Rectify Arc images:
    ## Def:
        # Determines the wavelength solution for arc images, starting with the 
        # initial values given by gsappwave (see step 5 in part 4)
    ## Steps:
        # 1) Run autoidentify on a defined section to identify spectral features
        #    based on the initial fit
        # 2) Run reidentify on full spectrum to establish wavelength calibration
        #    for all spatial points
    ## Notes:
        # 1) Although it is useful to check and make sure the correct lines have
        #    been identified, it is usually not necessary to interactively fit
        #    the spectral features
        # 2) It is best to use an arc image taken at approximately the same time
        #    as the image files for this step to avoid any inconsistancies based
        #    on how the CCDs flex or move as the telescope moves


unlearn gswavelength
unlearn gemextn

gswavelength.logfile="./Logs/gswaveLog.txt"
gswavelength.coordlist="gmos$data/CuAr_GMOS.dat"

gswavelength ("@arcFiles.txt", fwidth=6, order=5, nsum=50, fl_inte+)

##############################
#############################





##########  7  ###############
### Apply wavelength calibration:
    ## Def:
        # Defines a ploynomial to model the transformation from
        # pixel ---> wavelength scaling on the images
    ## Steps:
        # 1) Rectifies the bright sky emission lines in the science images. This
        #    means that the lines will be corrected so that they are vertically
        #    straight instead of being curved
        # 2) Applies the transformation found using gswavelength (above)
    ## Notes:
        # 1) It is worth looking at the images after this step to make sure the 
        #    sky emission lines have been correctly rectified so that they can 
        #    easily be removed in the next step


unlearn gstransform
gstransform.logfile="./Logs/gstransformLog.txt"

## Object:
gstransform ("Object.fits", wavtraname="@arcFiles.txt", fl_vardq+)

# Std Star
gstransform ("@stdFiles.txt", wavtraname="@arcFiles.txt")
          
del "stdFiles.txt"
hselect ("ts*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")

################################
################################



###########  8  ###############
### Sky sub:
    ## Def:
        # Remove the atmospheric emission lines from the image
    ## Steps:
        # 1) Define regions of the image to examine for atmospheric emission lines
        # 2) Model the emission lines with a function
        # 3) Apply corrections to the image to remove the emission lines
    ## Notes:
        # 1) This function essentially models the flux vs. y-value for each column
        #    of pixels in the image. It finds the continuum (highest flux) and
        #    marks an area significantly far away on either side of it. Outside the
        #    continuum, the sky emission should essentially be a flat shift off of
        #    0 flux. This function is then a simple shift of flux that removes the
        #    emission lines away from the continuum.
        # 2) It is worth checking the images before and after this step to see if 
        #    the emission lines have been removed correctly.
        # 3) MARKING RANGE --> Use form "r1:r2,r3:r4" to identify rows to sky subtract
        

unlearn gsskysub

gsskysub.logfile = "./Logs/gsskysub.txt"
gsskysub.long_sample = "10:200,300:500"     # Regions to define subtraction function
gsskysub.order = 5

## Object:
gsskysub ("tObject.fits", fl_answe-, fl_vardq+, fl_oversize-, fl_inte+)


## Std Star:
gsskysub ("@stdFiles.txt", fl_answe-,  fl_oversize-, fl_inte+)

del "stdFiles.txt"
hselect ("sts*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")

#############################
#############################



######### 9  ###############
### Extract Standard star spectrum:
    ## Def:
        # Extract a 1D spectrum from the standard star image
    ## Steps:
        # 1) Defines an "aperture", or region on the image over which the function
        #    looks to calculate the flux
        # 2) Using the wavelength calibration and the flux over the width of the 
        #    specified aperture, this function creates a 1D spectrum
    ## Notes:
        # 1) It is very important that the aperture is correctly defined so that
        #    the area immediately surronding and including the continuum is chosen.
        #    This may need to be done manually
        # 2) Weights should be added using the variance extension in order to
        #    eliminate noise on the flux values, raising the signal-to-noise ratio

unlearn gsextract

gsextract.logfile="./Logs/gsextractLog.txt"

gsextract ("@stdFiles.txt", tfunction="spline3", weights="variance", torder=9, \
   fl_inter-)

del "stdFiles.txt"
hselect ("ests*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")



###########  10  ##############
### Create Flux calibration:
    ## Def:
        # Creates a flux calibration based on the Std Star image to correct the
        # science image flux values
    ## Steps:
        # 1) Compares the flux value of the 1D Std star spectrum with a spectrum
        #    stored in the IRAF library
        # 2) Finds the flux ratio of the F_obs/F_lib and fits a low-order
        #    polynomial to the function
    ## Notes:
        # 1) Ideally, the polynomial should be a lower-order polynomial as the 
        #    flux calibration curve should be a somewhat smoother curve. This is
        #    then applied as a flux correction on the images later

unlearn gsstandard

gsstandard.logfile="./Logs/gsstandardLog.txt"
gsstandard.order = 3

gsstandard ("@stdFiles.txt", "std", "sens",  \
            starname = stdName, fl_inter-)


##################################
##################################




##########  11  ##############
### Flux calibration
    ## Def:
        # Applies a flux calibration and extinction correction to the science image
    ## Steps:
        # 1) Based on the flux calibration results in the last part, applies the
        #    the polynomial fit to the science data to correct the flux
        # 2) Applies an extinction fit to the data based on data from Manau Kea
   
unlearn gscalibrate

gscalibrate.logfile = "./Logs/gscalibrateLog.txt"
gscalibrate.extinction = "/N/u/dmsikors/Carbonate/Desktop/REUGemini/mk_extinct.dat"


gscalibrate ("@stdFiles.txt", sfunc="sens", fl_ext+, fl_scale-)

gscalibrate ("stObject.fits", sfunc="sens", fl_ext+, fl_scale-, fl_vardq+)



###############################
###############################


##########  12  ###############
### Extract spectrum
    ## Def:
        # Extract the spectrum
        
    
    
unlearn sarith


sarith ("cstObject.fits[SCI]", "copy", "", \
        "FinalSpectrum.fits", apertures = "240-270x5")













