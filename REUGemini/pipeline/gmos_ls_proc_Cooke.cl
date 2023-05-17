# GMOS Data Reduction Cookbook companion script to the chapter:
#   "Reduction of Longslit Spectra with IRAF" 
#
# IRAF CL command script to:

#Originally from
# Select GMOS calibration images for CQ507, in program GN-2021B-FT-201.
#    2021-Aug-07  kcc7952@rit.edu
#
# The names for the relevant header keywords and their expected values are described 
# in the DRC chapter entitled "Supplementary Material"
#
# Perform the following starting in the parent work directory:
#cd /Volumes/Cooke/Gemini_2021_FT/CQ507/combined
#cd /path/to/work_directory
# 
# To run this script in an IRAF session, load the gemini, gemtools & gmos packages, 
# then enter the following at the IRAF prompt: 
#   cl < gmos_img_proc.cl


# Set up the log file for this reduction
gmos.logfile="GN-2021B-FT-201.log"

print ("### Begin Processing GMOS/MOS Spectra ###")
print ("###")

print ("===Exposure Selection===")




print (" --Selecting Bias Exposures--")
#
## Select bias exposures. 
# Bias exposures with a common observation class, RoI, and CCD binning:
s1 = "obstype?='BIAS' && obsclass?='dayCal' && detrO1ys==512 && ccdsum?='1 2'"
# Select bias exposures within ~2 months of the target observations:
s2 = " && @'date-obs' > '2021-07-07' && @'date-obs' < '2021-09-07'"
string biasSelect
biasSelect = s1 // s2
# Select exposures using information from both PHDU and HDU-1
hselect ("*.fits[1,inherit=yes]", "$I", biasSelect, > "bias_tmp.txt")
# Now remove the IRAF index and kernel parameters from the file names
gemextn ("@bias_tmp.txt", omit="index,kernel", outfile="biasFull.txt")

# Also prepare CenterSpec RoI for the standard star.
#biasSelect = "obstype?='BIAS' && obsclass?='dayCal' && detrO1ys==512 && ccdsum?='1 2'"
#hselect ("*.fits[1,inherit=yes]", "$I", biasSelect, > "bias_tmp.txt")
#gemextn ("@bias_tmp.txt", omit="index,kernel", outfile="biasCenSp.txt")





### Flat-fields
#
print (" --Selecting Flat-field Expsoures--")
#
## Flat-fields must match the observation type, RoI and CCD binning.
s1 = "obstype?='FLAT' && obsclass?='partnerCal' && detro1ys==512 && ccdsum?='1 2'"

# Must also match the grating, aperture, and central wavelength 725:
s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=725.0"
string flatSelect
flatSelect = s1 // s2
# Generate lists of calibration filenames using above criteria.
hselect ("N*.fits[1,inherit=yes]", "$I", flatSelect, > "flt_tmp725.txt")
gemextn ("@flt_tmp725.txt", omit="index,kernel", outfile="flatFull725.txt")

# Must also match the grating, aperture, and central wavelength 715:
s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=715.0"
flatSelect = s1 // s2
# Generate lists of calibration filenames using above criteria.
hselect ("N*.fits[1,inherit=yes]", "$I", flatSelect, > "flt_tmp715.txt")
gemextn ("@flt_tmp715.txt", omit="index,kernel", outfile="flatFull715.txt")

# Must also match the grating, aperture, and central wavelength 735:
s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=735.0"
flatSelect = s1 // s2
# Generate lists of calibration filenames using above criteria.
hselect ("N*.fits[1,inherit=yes]", "$I", flatSelect, > "flt_tmp735.txt")
gemextn ("@flt_tmp735.txt", omit="index,kernel", outfile="flatFull735.txt")

### Science, Arcs, Standards
#
print (" --Selecting Comparison Arc Expsoures--")
#
# Select comparison arc exposures via observation type, CCD RoI & binning:
s1 = "obstype?='ARC' && obsclass?='dayCal' && detro1ys==512 && ccdsum?='1 2'"
s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=725"
string arcSelect
arcSelect = s1 // s2
hselect ("N*.fits[1,inherit=yes]", "$I", arcSelect, > "arc_tmp725.txt")
gemextn ("@arc_tmp725.txt", omit="index,kernel", outfile="arcFull725.txt")

s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=715.0"
arcSelect = s1 // s2
hselect ("N*.fits[1,inherit=yes]", "$I", arcSelect, > "arc_tmp715.txt")
gemextn ("@arc_tmp715.txt", omit="index,kernel", outfile="arcFull715.txt")

s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=735.0"
arcSelect = s1 // s2
hselect ("N*.fits[1,inherit=yes]", "$I", arcSelect, > "arc_tmp735.txt")
gemextn ("@arc_tmp735.txt", omit="index,kernel", outfile="arcFull735.txt")

print (" --Selecting Science Expsoures--")
#
# Select science exposures via observation type, RoI and CCD binning.
# This will also select data from object CQ 507. Delete these filename if you want.
s1 = "obstype?='OBJECT' && obsclass?='science' && detro1ys==512 && ccdsum?='1 2'"
s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=725.0"
string sciSelect
sciSelect = s1 // s2
hselect ("N*.fits[1,inherit=yes]", "$I", sciSelect, > "sci_tmp725.txt")
gemextn ("@sci_tmp725.txt", omit="index,kernel", outfile="sciFiles725.txt")

s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=715.0"
sciSelect = s1 // s2
hselect ("N*.fits[1,inherit=yes]", "$I", sciSelect, > "sci_tmp715.txt")
gemextn ("@sci_tmp715.txt", omit="index,kernel", outfile="sciFiles715.txt")

s2 = "&& grating?='R831+_G5302' && maskname?='1.0arcsec' && grwlen=735.0"
sciSelect = s1 // s2
hselect ("N*.fits[1,inherit=yes]", "$I", sciSelect, > "sci_tmp735.txt")
gemextn ("@sci_tmp735.txt", omit="index,kernel", outfile="sciFiles735.txt")



#print (" --Selecting Standard Star Expsoures--")
#
#s1 = "obstype?='OBJECT' && obsclass?='partnerCal' && detro1ys==512 && ccdsum?='1 2'"
#string stdSelect
#stdSelect = s1 // s2
#hselect ("N*.fits[1,inherit=yes]", "$I", stdSelect, > "std_tmp.txt"
#gemextn ("@std_tmp.txt", omit="index,kernel", outfile="stdFiles.txt")

print ("===Creating MasterCals===")

print (" --Creating Bias Residual MasterCals--")

## Create Bias Residual, including VAR and DQ arrays.
# Use primarily the default task parameters.
unlearn gbias
unlearn gemextn    # Disarm a bug in gbias
gbias.log = "gbiasLog.txt"
#gbias.rawpath="./"
gbias.verbose=no

#gbias ("@biasCenSp.txt", "MCbiasCenSp", fl_vardq+)

#gbias is broken so the following line will fail, so just use the pre-prepared bias
#gbias ("@biasFull.txt", "MCbiasFull", fl_vardq+)

#use N20210726S0225_bias.fits for bias

print (" --Creating Flat-field MasterCals--")
#
## Process flat-field, including VAR and DQ arrays.
# Use primarily the default task parameters.
unlearn gireduce
unlearn gsflat
gsflat.logfile="gsflatLog.txt"
#gsflat.rawpath="./"
gsflat.verbose=no

# Normalize the spectral flats per CCD.
# Interactive curve fitting is better, but we won't interupt the flow here.
print ("--Flat-field normalization, non-interactive--")
#
gsflat ("@flatFull725.txt", "MCflatFull725.fits",  bias="N20210726S0225_bias.fits", fl_vardq+, fl_fulldq+, fl_detec+, fl_oversize-, fl_inter-, order="13,11,28")

gsflat ("@flatFull715.txt", "MCflatFull715.fits",  bias="N20210726S0225_bias.fits", fl_vardq+, fl_fulldq+, fl_detec+, fl_oversize-, fl_inter-, order="13,11,28")

gsflat ("@flatFull735.txt", "MCflatFull735.fits",  bias="N20210726S0225_bias.fits", fl_vardq+, fl_fulldq+, fl_detec+, fl_oversize-, fl_inter-, order="13,11,28")

#gsflat ("S20070623S0108", "MCflatCenSp.fits", bias="MCbiasCenSp.fits", \
   fl_vardq+, fl_fulldq+, fl_detec+, fl_oversize-, fl_inter-, order="13,11,28")

### Perform basic processing
print ("===Processing Science Files===")
print (" --Performing Basic Processing--")

# Use primarily the default task parameters.
unlearn gsreduce
gsreduce.logfile="gsreduceLog.txt"
#gsreduce.rawpath="./"
gsreduce.verbose=no

# Perform basic reductions on all exposures for science targets.
gsreduce ("@sciFiles725.txt", bias="N20210726S0225_bias.fits", flatim="MCflatFull725", fl_fixpix-, fl_oversize-, fl_vardq+, fl_fulldq+)
gsreduce ("@sciFiles715.txt", bias="N20210726S0225_bias.fits", flatim="MCflatFull715", fl_fixpix-, fl_oversize-, fl_vardq+, fl_fulldq+)
gsreduce ("@sciFiles735.txt", bias="N20210726S0225_bias.fits", flatim="MCflatFull735",
    fl_fixpix-, fl_oversize-, fl_vardq+, fl_fulldq+)

# Perform basic reductions on the Arcs and standard star.
gsreduce ("@arcFull725.txt", bias="N20210726S0225_bias.fits", fl_fixpix-, fl_flat-, fl_oversize-)
gsreduce ("@arcFull715.txt", bias="N20210726S0225_bias.fits", fl_fixpix-, fl_flat-, fl_oversize-)
gsreduce ("@arcFull735.txt", bias="N20210726S0225_bias.fits", fl_fixpix-, fl_flat-, fl_oversize-)
         
         
         
#starndard star prepare files
#gsreduce ("S20070623S0109", bias="MCbiasCenSp", \
   fl_fixpix-, fl_flat-, fl_oversize-)

#gsreduce ("@stdFiles.txt", bias="MCbiasCenSp", flatim="MCflatCenSp", \
   fl_fixpix+, fl_oversize-, fl_vardq-)

### End of basic processing. Continue with advanced processing.
## Cosmic ray rejection.
print (" --Perform multi-frame cosmic ray rejection--")
#
# Create a list of common slit alignments for the target & standard star.
#sections gs//@stdFiles.txt > gsStdFiles.txt

s1 = "obstype?='OBJECT' && obsclass?='science' && grwlen=725.0"
hselect ("gsN*.fits[0]", "$I", (s1 // "&& i_title?='CQ 507'"), > "CQ_507_725_tmp.txt")
gemextn ("@CQ_507_725_tmp.txt", omit="index", outfile="CQ_507_725.txt")
dele ("CQ_507_725_tmp.txt")
                                
s1 = "obstype?='OBJECT' && obsclass?='science' && grwlen=715.0"
hselect ("gsN*.fits[0]", "$I", (s1 // "&& i_title?='CQ 507'"), > "CQ_507_715_tmp.txt")
gemextn ("@CQ_507_715_tmp.txt", omit="index", outfile="CQ_507_715.txt")
dele ("CQ_507_715_tmp.txt")
                                                                
s1 = "obstype?='OBJECT' && obsclass?='science' && grwlen=735.0"
hselect ("gsN*.fits[0]", "$I", (s1 // "&& i_title?='CQ 507'"), > "CQ_507_735_tmp.txt")
gemextn ("@CQ_507_735_tmp.txt", omit="index", outfile="CQ_507_735.txt")
dele ("CQ_507_735_tmp.txt")

print (" --Combine science exposures--")

# Use primarily the default task parameters.
unlearn gemcombine
gemcombine.logfile="gemcombineLog.txt"
gemcombine.reject="ccdclip"
gemcombine.verbose=no

# Combine the exposures with outlier rejection for each orientation.
                                
#We only have one exposure at each central wavelength for this target, so we can comment this out for now, we'll have to combine the final wavelength calibrated exposures later for proper alignment
#gemcombine ("@gsStdFiles.txt", "LTT9239.fits", fl_vardq-, fl_dqprop+)
#gemcombine ("CQ_507.txt", "CQ507", fl_vardq+, fl_dqprop+)


# Clean up
#imdele ("gS20*.fits")

print (" --Perform wavelength calibration--")
#
# Use primarily the default task parameters.
unlearn gswavelength
gswavelength.logfile="gswaveLog.txt"

# In this case, the default medium-resolution line list will work well.
gswavelength.coordlist="gmos$data/CuAr_GMOS.dat"

# The fit to the dispersion relation should be performed interactively. 
# Here we will us a previously determined result.
#725
gswavelength ("gs//@arcFull725.txt", fl_inter-, fwidth=6, order=5, nsum=50)
#715
gswavelength ("gs//@arcFull715.txt", fl_inter-, fwidth=6, order=5, nsum=50)
#735
gswavelength ("gs//@arcFull735.txt", fl_inter-, fwidth=6, order=5, nsum=50)

print (" --Apply wavelength calibration--")
#
# Use primarily the default task parameters.
unlearn gstransform
gstransform.logfile="gstransformLog.txt"

#gstransform ("LTT9239", wavtraname="gsS20070623S0109", fl_vardq+)
gstransform ("@CQ_507_725.txt", wavtraname="gsN20210809S0319", fl_vardq+)

gstransform ("@CQ_507_715.txt", wavtraname="gsN20210809S0321", fl_vardq+)

gstransform ("@CQ_507_735.txt", wavtraname="gsN20210809S0320", fl_vardq+)


#Clean up.
#imdele ("gsS2007*.fits")

print (" --Perform sky subtraction--")
#
## Sky subtraction. This will require summing the spectra along columns, e.g.: 
#pcols tAM2306b.fits[SCI] 1100 2040 wy1=40 wy2=320

# Subtract sky spectrum using selected regions. 
# The regions should be selected with care, using e.g. prows/pcols. 
unlearn gsskysub
gsskysub.logfile="gsskysubLog.txt"

#gsskysub ("tLTT9239", fl_oversize-, fl_vardq-, long_sample="20:70,190:230")
#735
gsskysub ("tgsN20210807S0178", fl_oversize-, fl_vardq+)
#715
gsskysub ("tgsN20210807S0181", fl_oversize-, fl_vardq+)
#725
gsskysub ("tgsN20210807S0175", fl_oversize-, fl_vardq+)


#print (" --Extract Std spectrum--")
#
### Flux calibration with LTT9239.
# Extract the std spectrum using a large aperture.
# It's important to do this interactively. 
unlearn gsextract
#gsextract.logfile="gsextractLog.txt"
#gsextract ("stLTT9239", fl_inter-, apwidth=3., tfunction="spline3", torder=9)

#print (" --Perform Flux calibration--")
#
# Derive the sensitivity function.
# Be sure to download the custom Mauna Kea atmospheric exteinction function. 
#unlearn gsstandard
#gsstandard.logfile="gsstdLog.txt"
#gsstandard.caldir="onedstds$ctionewcal/"

#gsstandard ("estLTT9239", sfile="std", sfunction="sens", fl_inter-, \
   order=7, starname="l9239", extinction="./mk_extinct.txt")

## Apply the sensitivity function.
#unlearn gscalibrate
#gscalibrate.logfile="gscalibrateLog.txt"
#gscalibrate.extinction="./mk_extinct.txt"

#gscalibrate ("stAM2306*", sfunc="sens", fl_ext+, fl_scale-, fl_vardq+)
#gscalibrate ("estLTT9239", sfunc="sens", fl_ext+, fl_scale-)

print (" --Extract Target Spectra--")
#No idea what this method is, but for reference we'll keep it
#nsum=4
#sarith ("cstAM2306b.fits[SCI]", "copy", "", "ecstAM2306b.ms", apertures="222-346x4")

                                
#YOU REALLY SHOULD EXTRACT MANUALLY, the automatic target identification will often get it wrong
#715 good
gsextract stgsN20210807S0181 apwidth=3 fl_inter+
#725 iffy
gsextract stgsN20210807S0175 apwidth=2 fl_inter+
#735 good
gsextract stgsN20210807S0178 apwidth=3 fl_inter+
                                
                                
                                
                                
#Now combine the good exposures
gemscombine inimages=estgsN20210807S0178,estgsN20210807S0181 outimage=CQ507stack.fits combine=median
print ("===Finished Calibration Processing===")
