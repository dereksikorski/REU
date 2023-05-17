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
###################
##
##
##
##
##

###############
### Inputs:

cd /N/u/dmsikors/Carbonate/Desktop/REUGemini/Oct28/     # Enter path of data
string filterWl = "800"                      # Central wavelength of data

################




################
### Create temp folder in main directory and observation log for night

mkdir temp      # Create a temp folder to store all temp data to be deleted

## Make a observation log in main directory
del "obsLog.txt"
cd raw

string fieldNames = "# File         ObsDate   ObsTime   Object  ObsType   ObsClass  CcdBin  ROI Filter1    Filter2    Disperser   SlitName  SlitType    Rotator CenWave T_exp   Airmass"
print (fieldNames, >"obsLog.txt")   # Create an observation log .txt file                     
string keyWords = "$I Date-Obs Time-Obs Object ObsType ObsClass Ccdsum DetRO1ys Filter1 Filter2 Grating MaskName MaskType PA GrWlen ExpTime Airmass"


hselect ("*.fits[1,inherit=yes]", keyWords, "obstype != 'acq'", >"obsLog.txt")
rename("obsLog.txt", "../")
################




##################
### Make list of all frames in the raw folder:

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

del (filterWl//"flatFiles.txt"  )   # Deletes file if repeated  
gemextn ("@flt_tmp.txt", omit="index,kernel", outfile= filterWl // "flatFiles.txt" )  #Copy flats file, removing uneeded data
del "flt_tmp.txt" 
rename (filterWl//"flatFiles.txt"  , "../temp/")     # Rename file to move into temp directory

## Arcs:
string arcCriteria = "obstype=='ARC' && GrWlen=="   # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (arcCriteria // (filterWl//".")) , > "arc_tmp.txt")     # Select arc files for given wl

del (filterWl//"arcFiles.txt")
gemextn ("@arc_tmp.txt", omit="index,kernel", outfile=(filterWl//"arcFiles.txt"))  # Trim the file names within the arc list
del "arc_tmp.txt"
rename (filterWl//"arcFiles.txt", "../temp/")

## Science frames:
string objCriteria = "obstype=='OBJECT' && GrWlen=="    # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (objCriteria // (filterWl//".")), >"obj_tmp.txt")

del (filterWl//"objFiles.txt")
gemextn ("@obj_tmp.txt", omit="index,kernel", outfile=(filterWl//"objFiles.txt"))
del "obj_tmp.txt"
rename (filterWl//"objFiles.txt","../temp/")
#####################


cd ..
cd temp



####################
### Perform overscan corrections on all frames:

## Initialize parameters:
unlearn gsreduce
gsreduce.logfile="overscanLog.txt"
gsreduce.rawpath="../raw"
gsreduce.outpref=""

## Science frames:
gsreduce ("@"//(filterWl//"objFiles.txt"), outimag="@"//(filterWl//"objFiles.txt"), \
             fl_over+, fl_trim+, fl_bias-, fl_flat-, fl_gmos-, fl_fixp-, fl_cut-,fl_titl-, fl_vardq+)

## Arcs:
gsreduce ("@"//(filterWl//"arcFiles.txt"), outimag="@"//(filterWl//"arcFiles.txt"),  \
             fl_over+, fl_trim+, fl_bias-, fl_flat-, fl_gmos-, fl_fixp-, fl_cut-,fl_titl-, fl_vardq+)

## Flats:
gsreduce ("@"//(filterWl//"flatFiles.txt"), outimag="@"//(filterWl//"flatFiles.txt"),  \
             fl_over+, fl_trim+, fl_bias-, fl_flat-, fl_gmos-, fl_fixp-, fl_cut-,fl_titl-, fl_vardq+)

## Biases:
gsreduce ("@biasFiles.txt", outimag="@biasFiles.txt",  \
             fl_over+, fl_trim+, fl_bias-, fl_flat-, fl_gmos-, fl_fixp-, fl_cut-,fl_titl-, fl_vardq+)
###################





##################
### Combine Biases:

unlearn gbias

gbias.bpm="bpm_gmos-s_EEV_v1_2x2_img_MEF.fits"
gbias.logfile="gbiasLog.txt"

gbias ("@biasFiles.txt", "MCbias.fits", fl_over+, fl_trim+, fl_vardq+)

rename ("gbiasLog.txt", "../")

##################





#################
### Bias correct and combine flat fields:

unlearn gsflat
unlearn gemextn

gsflat.logfile = "gsflatLog.txt"
gsflat.bias = "MCbias.fits"

gsflat ("@"//(filterWl//"flatFiles.txt"), filterWl//"MCflat.fits", \
    fl_oversize-, fl_vardq+, fl_fulldq+, fl_detec+, order="13,11,28")
#################






################
### Reduce all of the science frames with the combined bias and flat field frames:

unlearn gsreduce
unlearn gemextn

gsreduce.outpref = "o"
gsreduce.bias = "MCbias.fits"
gsreduce.flatim = filterWl//("MCflat.fits")

gsreduce ("@"//(filterWl//"objFiles.txt"), fl_over-, fl_trim-, \
        fl_gmosaic-, fl_fixpix- , fl_gsappwave-, fl_cut- , fl_title- , fl_oversize- , fl_vardq+)

hselect ("o*.fits[1,inherit=yes]", "$I", "obstype != 'acq'", >"tempobj.txt")

gemextn ("@tempobj.txt", omit="index,kernel", outfile="obj.txt")  

del "tempobj.txt"

###############





###############
### Reduce the arc frames with bias corrections:

unlearn gsreduce
unlearn gemextn

gsreduce.outpref = "a"
gsreduce.bias = "MCbias.fits"

gsreduce ("@"//(filterWl//"arcFiles.txt"), fl_over-, fl_trim-,fl_flat-, \
            fl_gmosaic-, fl_fixpix-, fl_gsappwave-, fl_cut-, fl_title-, fl_oversize-, fl_vardq+)
            
            
            
hselect ("a*.fits[1,inherit=yes]", "$I", "obstype != 'acq'", >"temparc.txt")

gemextn ("@temparc.txt", omit="index,kernel", outfile="arc.txt")

del "temparc.txt"

##############



#############
### Mosaic science and arc frames:

unlearn gemextn
unlearn gmosaic

gmosaic.logfile = "gmosLog.txt"

## Science:

gmosaic ("@obj.txt", fl_paste-, fl_vardq+)

hselect ("mo*.fits", "$I", "obstype != 'acq'", >"tempobj.txt")

gemextn ("@tempobj.txt", omit="index,kernel", outfile="mobj.txt")
del "tempobj.txt"

## Arcs:

gmosaic ("@arc.txt", fl_paste-, fl_vardq+)

hselect ("ma*.fits", "$I", "obstype != 'acq'", >"temparc.txt")

gemextn ("@temparc.txt", omit="index,kernel", outfile="marc.txt")


#############








##############
### Apply wavelength calibrations with arcs
unlearn gswavelength

gswavelength ( "@marc.txt",  fwidth=6, order=7, nsum=50, step=2)






###############







