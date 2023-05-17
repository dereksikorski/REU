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


### Inputs:
cd /N/u/dmsikors/Carbonate/Desktop/REUGemini/Oct28/     # Enter path of data
string filterWl = "800"             # Central wavelength of data

###


### Create temp file in main directory and observation log for night

mkdir temp      # Create a temp folder to store all temp data to be deleted

    # Make a observation log in main directory
del "obsLog.txt"
cd raw

string fieldNames = "# File         ObsDate   ObsTime   Object  ObsType   ObsClass  CcdBin  ROI Filter1    Filter2    Disperser   SlitName  SlitType    Rotator CenWave T_exp   Airmass"
print (fieldNames, >"obsLog.txt")   # Create an observation log .txt file                     
string keyWords = "$I Date-Obs Time-Obs Object ObsType ObsClass Ccdsum DetRO1ys Filter1 Filter2 Grating MaskName MaskType PA GrWlen ExpTime Airmass"


hselect ("*.fits[1,inherit=yes]", keyWords, "obstype != 'acq'", >"obsLog.txt")
rename("obsLog.txt", "../")




### Make a list of the Bias frames in the directory:


    ## Identify header info and make list
string biasCriteria = "obstype=='BIAS' && obsclass=='dayCal'"   # Criteria for hselect next line
hselect ("*.fits[1,inherit=yes]", "$I", biasCriteria, > "bias_tmp.txt") # Select names of bias images in "raw" and store in list

    ## trim end of file names
string biasOmit = "index,kernel"
del "biasFiles.txt"     # Deletes file if repeated. Will give warning if file does not exist
gemextn ("@bias_tmp.txt", omit=biasOmit, outfile="biasFiles.txt")       # Copy the temporary file, removing uneeded info
del "bias_tmp.txt"
rename ("biasFiles.txt", "../temp/")     # Rename file so it is moved into temp directory







### Make list of the flat field frames in the directory:
    # NOTE: Still located in raw directory...

    ## Make list of all flats of given waveband
string flatCriteria = "i_title=='GCALflat' && GrWlen=="     # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (flatCriteria // (filterWl//".")) , >> "flt_tmp.txt")   # Select flat field images based on wavelength
    
    ## Trim file names in final list
string flatOmit = "index,kernel"        # info to omit
del (filterWl//"flatFiles.txt"  )   # Deletes file if repeated  
gemextn ("@flt_tmp.txt", omit=flatOmit, outfile= filterWl // "flatFiles.txt" )  #Copy flats file, removing uneeded data
del "flt_tmp.txt" 
rename (filterWl//"flatFiles.txt"  , "../temp/")     # Rename file to move into temp directory





### Make list of the arcs for the given waveband

    ## Make list of arc files within the "raw" folder
string arcCriteria = "obstype=='ARC' && GrWlen=="   # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (arcCriteria // (filterWl//".")) , > "arc_tmp.txt")     # Select arc files for given wl
    
    ## Trim the file names and move list
string arcOmit = "index,kernel" 
del (filterWl//"arcFiles.txt")
gemextn ("@arc_tmp.txt", omit=arcOmit, outfile=(filterWl//"arcFiles.txt"))  # Trim the file names within the arc list
del "arc_tmp.txt"
rename (filterWl//"arcFiles.txt", "../temp/")



### Make a list of object data for given waveband

    ## Make list of OBJECT files within "raw" folder
string objCriteria = "obstype=='OBJECT' && GrWlen=="    # hselect criteria
hselect ("*.fits[1,inherit=yes]", "$I", (objCriteria // (filterWl//".")), >"obj_tmp.txt")

    ## Trim file names and move list
string objOmit = "index,kernel"
del (filterWl//"objFiles.txt")
gemextn ("@obj_tmp.txt", omit=objOmit, outfile=(filterWl//"objFiles.txt"))
del "obj_tmp.txt"
rename (filterWl//"objFiles.txt","../temp/")


cd ..




### Combine the bias frames

cd temp     # change to temp dir to create files.

del "MCbias.fits"      # Delete MasterCal bias image in case name is taken
unlearn gbias          # Unlearn all params for gbias

gbias.bpm="bpm_gmos-s_EEV_v1_2x2_img_MEF.fits"      # Set the static BPM (Bad Pixel Mask)
gbias.logfile="gbiasLog.txt"            # Set a logfile in the data directory
gbias.rawpath = "../raw"             # Set path to the raw files

gbias ("@biasFiles.txt", "MCbias.fits", fl_vardq+, fl_over+)      # Run gbias making variance and data quality extensions
rename ("gbiasLog.txt", "../")


fxheader MCbias.fits       # List the file contents, display the image, etc., for QA







### Bias correct and combine Flats

del (filterWl//"MCflat.fits")       # Delete repeated file

unlearn gsflat
unlearn gemextn

gsflat.logfile = "gsflatLog.txt"        # Set logfile location
gsflat.rawpath = "../raw"               # Set location of raw images
gsflat.bias="MCbias.fits"               # Call Bias image (threw error when calling with gsflat, so did it before)


gsflat ("@"//(filterWl//"flatFiles.txt"), filterWl//("MCflat.fits"), fl_vardq+, fl_dete+, fl_over+, fl_inte-, order="13,11,28")

rename ("gsflatLog.txt", "../")





### Bias correct the arc images:

unlearn gsreduce
gsreduce.logfile="arc_gsreduceLog.txt"
gsreduce.rawpath="../raw"
gsreduce.bias="MCbias.fits"

gsreduce ("@"//(filterWl//"arcFiles.txt"), fl_fixp-, fl_over+,fl_trim+, fl_vard+, fl_full+,fl_flat-, fl_inte-)

rename ("arc_gsreduceLog.txt", "../")




### Process object files with bias and flat-fields

unlearn gsreduce

gsreduce.logfile="obj_gsreduceLog.txt"
gsreduce.rawpath="../raw"
gsreduce.bias="MCbias.fits"
gsreduce.flatim= (filterWl//("MCflat.fits"))

gsreduce ("@"//(filterWl//"objFiles.txt"), fl_fixp-, fl_over+,fl_trim+, fl_vard+, fl_full+, fl_inte-)

rename ("obj_gsreduceLog.txt", "../")




















