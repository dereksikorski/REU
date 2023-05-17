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
### Inputs:
       ## LINE 1: Enter the ABSOLUTE path of the desired night. The pasted
       ## directory should have a folder called "raw" that contains ALL raw data
       ## LINE 2: Enter the central wavelength of the images to reduce WITHOUT a "."

cd /N/u/dmsikors/Carbonate/Desktop/REUGemini/Oct28/     # Enter path of data
string filterWl = "800"                      # Central wavelength of data

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
string objCriteria = "obstype=='OBJECT' && GrWlen=="    # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (objCriteria // (filterWl//".")), >"obj_tmp.txt")

del "objFiles.txt"
gemextn ("@obj_tmp.txt", omit="index,kernel", outfile="objFiles.txt")
del "obj_tmp.txt"
rename ("objFiles.txt","../temp/")

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

##############################
##############################





##############  2  ####################
### Create a bias Master Calibration image:
    ## Images have a "dark fixed-pattern noise" as a result of the manufacturing
    ## of the CCD. This pattern is imprinted in every image the camera takes.
    ## This step produces a Master Bias image by combining a number of bias images
    ## to approximate and remove the pattern from the other images.
    
unlearn gbias

gbias.logfile = "./Logs/gbiasLog.txt"

gbias ("@biasFiles.txt", "MCbias.fits")     # Combine and perform overscan corrections
    
    
######################################
######################################
    
    
    
    
###############  3  ####################
### Create a flat Master Calibration image:
    ## By overexposing the ccd with a bright light, the dead and hot pixels can
    ## be identified. This Master Cal is then applied to the other ccd images
    ## to remove different pixel sensitivities.

unlearn gsflat
    
    
gsflat.logfile = "./Logs/gsflat.txt"
gsflat.bias = "MCbias.fits"

gsflat ("@flatFiles.txt", "MCflat.fits", fl_vardq+, fl_fulldq+, \
        fl_fixpix-, fl_detec+, order="13,11,28")

#######################################
#######################################




##############  4  #####################
### Reduce the science flats:
    ## This step does a number of procedures:
    

unlearn gsreduce
unlearn gemextn


gsreduce.logfile = "./Logs/objReduction.txt"
gsreduce.outpref = "o"
gsreduce.bias = "MCbias.fits"
gsreduce.flatim = "MCflat.fits"

gsreduce ("@objFiles.txt",fl_crsp-, fl_gmos-, fl_fixp-, fl_gsap-, \
        fl_cut-, fl_titl-, fl_oversize-, fl_vardq+)

del "objFiles.txt"
hselect ("o*.fits[1,inherit=yes]", "$I", "obstype != 'acq'", >"temp.txt")
gemextn ( "@temp.txt", omit="index,kernel", outfile="objFiles.txt")
del "temp.txt"

#######################################
#######################################



###############  5  ##################
### Reduce the Arc images given the Master Cal Bias image:


unlearn gsreduce
unlearn gemextn


gsreduce.logfile = "./Logs/arcReduction.txt"
gsreduce.outpref = "a"
gsreduce.bias = "MCbias.fits"

gsreduce ("@arcFiles.txt", fl_flat-, fl_gmos-, fl_fixp-, fl_gsap-, \
            fl_cut-, fl_titl-, fl_vardq+)



del "arcFiles.txt"
hselect ("a*.fits[1,inherit=yes]", "$I", "obstype != 'acq'", >"temp.txt")
gemextn ( "@temp.txt", omit="index,kernel", outfile="arcFiles.txt")
del "temp.txt"




#################################
#################################


#############  6  ###############
### Mosaic images together:

unlearn gmosaic
gmosaic.logfile = "./Logs/gmosLog.txt"

## Science:

gmosaic ("@objFiles.txt", fl_vardq+, fl_full+)

del "objFiles.txt"
hselect ("mo*.fits", "$I", "obstype != 'acq'", >"objFiles.txt")


## Arcs:

gmosaic ("@arcFiles.txt", fl_fixp+, fl_clea+, fl_full+)

del "arcFiles.txt"
hselect ("ma*.fits", "$I", "obstype != 'acq'", >"objFiles.txt")













