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
       ## LINE 3: Enter the IRAF name of the standard star being used, found at:
               ## https: //noirlab.edu/science/programs/csdc/usngo/gmos-cookbook/Processing/Supplement.html#stdstar-list

cd /N/u/dmsikors/Carbonate/Desktop/REUGemini/Oct28/     # Enter path of data
string filterWl = "800"                      # Central wavelength of data
string stdName = "hilt600"

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

## Standard Star:
string starCriteria = "obsclass=='partnerCal' && obstype == 'OBJECT'"
hselect ("*.fits[1,inherit=yes]", "$I", starCriteria, >"std_temp.txt")

del "stdFiles.txt"
gemextn ("@std_tmp.txt", omit="index,kernel", outfile="stdFiles.txt")
del "std_tmp.txt"
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
    ## Images have a "dark fixed-pattern noise" as a result of the manufacturing
    ## of the CCD. This pattern is imprinted in every image the camera takes.
    ## This step produces a Master Bias image by combining a number of bias images
    ## to approximate and remove the pattern from the other images.
    
unlearn gbias

gbias.logfile = "./Logs/gbiasLog.txt"

gbias ("@biasFiles.txt", "MCbias.fits", fl_vardq+)     # Combine and perform overscan corrections
    
    
######################################
######################################
    

############  3  #######################
### Create flat field Master Calibration


unlearn gsflat
unlearn gemextn
gsflat.logfile = "./Logs/gsflatLog.txt"

gsflat ("@flatFiles.txt", "MCflat.fits", bias="MCbias.fits", \
   fl_vardq+, fl_detec+, fl_oversize-, fl_inter-, order="13,11,28")
   
######################################
######################################




###########  4  #####################
### Basic reductions:

unlearn gsreduce
unlearn gemextn

gsreduce.logfile= "./Logs/ReductionsLog.txt"
gsreduce.bias = "MCbias.fits"

## Object
gsreduce ("@objFiles.txt", outpref="o", flatim="MCflat.fits", \
          fl_crsp+, fl_fixpix-, fl_oversize-, fl_vardq+, fl_fulldq+)
          
del "objFiles.txt"
hselect ("o*.fits", "$I", "obstype != 'acq'", >"objFiles.txt")

          
## Arc
gsreduce ("@arcFiles.txt", outpref="a", \
   fl_fixpix-, fl_flat-, fl_oversize-)

del "arcFiles.txt"
hselect ("a*.fits", "$I", "obstype != 'acq'", >"arcFiles.txt")

## Std Star
gsreduce ("@stdFiles.txt", outpref="s", flatim="MCflat.fits", \
          fl_crsp+, fl_fixpix-, fl_oversize-, fl_vardq+, fl_fulldq+)
          
del "stdFiles.txt"
hselect ("s*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")
################################
################################


############  5  ##############
### Combine same frames:

unlearn gemcombine


gemcombine ("@objFiles.txt", "Object.fits", fl_vardq+, fl_dqprop+)


#############################
############################



##########  6  ################
### Rectify Arc images:


# Copy list for next steps
mkdir z
cp "arcFiles.txt" ./z
cd z
rename ("arcFiles.txt", "../wavtran.txt")
cd ..
!rm -r z
        ####
        
unlearn gswavelength
unlearn gemextn

gswavelength.logfile="./Logs/gswaveLog.txt"
gswavelength.coordlist="gmos$data/CuAr_GMOS.dat"

gswavelength ("@arcFiles.txt", fwidth=6, order=5, nsum=50)

##############################
#############################





##########  7  ###############
### Apply wavelength calibration


unlearn gstransform
gstransform.logfile="gstransformLog.txt"

## Object:
gstransform ("Object.fits", wavtraname="@wavtran.txt", fl_vardq+)

# Std Star
gstransform ("@stdFile.txt", wavtraname="@wavtran.txt", fl_vardq+)
          
del "stdFiles.txt"
hselect ("ts*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")

################################
################################



###########  8  ###############
### Sky sub
unlearn gsskysub

## Object:
gsskysub ("tObject.fits", fl_answe="no", fl_vardq+, fl_oversize-)


## Std Star:
gsskysub ("@stdFiles.txt", fl_answe="no", fl_vardq+, fl_oversize-)

del "stdFiles.txt"
hselect ("sts*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")

#############################
#############################



######### 9  ###############
### Extract Standard star spectrum:

unlearn gsextract

gsextract ("@stdFiles", tfunction="spline3", torder=9, \
   fl_inter-)

del "stdFiles.txt"
hselect ("ests*.fits", "$I", "obstype != 'acq'", >"stdFiles.txt")



###########  10  ##############
### Create Flux calibration:
unlearn gsstandard


gsstandard ("@stdFiles.fits", order="7", starname = stdName, fl_inter-)


##################################
##################################




##########  11  ##############
### Flux calibration
   
unlearn gscalibrate
gscalibrate.extinction = "../raw/mk_extinct.dat"


gscalibrate ("@stdFiles.fits", sfunc="sens", fl_ext+, fl_scale-)

gscalibrate ("stObject.fits", sfunc="sens", fl_ext+, fl_scale-, fl_vardq+)



###############################
###############################


##########  12  ###############
## Extract spectrum
unlearn sarith

sarith ("cstObject.fits[SCI]", "copy", "", "FinalSpectrum.fits")













