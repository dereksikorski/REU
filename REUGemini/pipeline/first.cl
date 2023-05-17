### Inputs:
cd /N/u/dmsikors/Carbonate/Desktop/REU/Aug7     # Enter path of data
string filterWl = "715."             # Central wavelength of data


filterName = "w"//filterWl


### Make a list of the Bias frames in the directory:

    # Identify header info and make list
string biasCriteria = "obstype=='BIAS' && obsclass=='dayCal'"
hselect ("*.fits[1,inherit=yes]", "$I", biasCriteria, > "bias_tmp.txt")
    # trim end of file names
string omit = "index,kernel"
del "biasFiles.txt"     # Deletes file if repeated
gemextn ("@bias_tmp.txt", omit=omit, outfile="biasFiles.txt")
del "bias_tmp.txt"







### Make list of the flat field frames in the directory:

    # Make list of all flats of given waveband
string flatCriteria = "i_title=='GCALflat' && GrWlen=="
hselect ("*.fits[1,inherit=yes]", "$I", (flatCriteria // filterWl) , > "flt_tmp.txt")
    # Trim file names in final list
string flatOmit = "index,kernel"
fileEnd = filterName // "txt"
flatFileName = "flatFiles"//fileEnd
del fileName      # Deletes file if repeated
gemextn ("@flt_tmp.txt", omit=omit, outfile=flatFileName)
del "flt_tmp.txt"





### Make list of the arcs for the given waveband

string arcCriteria = "obstype=='ARC' && GrWlen=="
hselect ("*.fits[1,inherit=yes]", "$I", (arcCriteria // filterWl) , > "arc_tmp.txt")
string arcOmit = "index,kernel"
arcFileName = "arcFiles"//fileEnd
del fileName
gemextn ("@arc_tmp.txt", omit=omit, outfile=arcFileName)
del "arc_tmp.txt"






### Combine the bias frames

del "MCbias.fits"
unlearn gbias
gbias.bpm="bpm_gmos-s_EEV_v1_2x2_img_MEF.fits"
gbias.logfile="gbiasLog.txt"

gbias ("@biasFiles.txt", "MCbias.fits", fl_vardq+)

# List the file contents, display the image, etc., for QA
fxheader MCbias.fits







### Bias correct and combine Flats

ender = filterWl//"fits"
flatFile = "MCflat"//ender

unlearn gsflat
unlearn gemextn

gsflat.logfile = "gsflatLog.txt"


gsflat ("@"//flatFileName, flatFile, bias="MCbias.fits", \
   fl_vard+, fl_dete+, fl_over-, fl_inte+, order="13,11,28")










