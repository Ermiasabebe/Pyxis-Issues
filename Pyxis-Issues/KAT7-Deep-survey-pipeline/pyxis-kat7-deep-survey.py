import Pyxis
import ms
import mqt
import im
import lsm
import stefcal
import os
import sys
import math
import tempfile
from Pyxis.ModSupport import *
import pyrap.tables
import Tigger
from Tigger.Coordinates import angular_dist_pos_angle
import im.argo as argo

from pyrap.tables import table

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pyfits as pf
from astLib.astWCS import WCS

## 2. Default variable assignments
# These are variables that control the behaviour of modules.
# Here we just set up some reasonable defaults; config files 
# will almost certainly override them.
# Set up destination directory for all plots, images, etc.
DESTDIR_Template = 'plots-${MS:BASE}${-stage<STAGE}'
#  up base filename for these files
OUTFILE_Template = '${DESTDIR>/}${MS:BASE}${_s<STEP}${_<LABEL}'

# Extract unpolarized sky models by default
lsm.PYBDSM_POLARIZED = False     
# 0 means we make a single image from all channels. Use 1 to
# make a per-channel cube.
# default place to look for MSs
MS_TARBALL_DIR = os.path.expanduser("~/data")


# Default initial sky model
LSM0 = 'Deep-de-initial.lsm.html'

LMS1 = 'Deep-de-cal1.lsm.html'

LSM2 = "Deep-de-cal2.lsm.html"

LSM3 = "Deep-de-cal3.lsm.html"

LSM4 = "Deep-de-cal4.lsm.html"

LSM5 = "Deep-de-cal5.lsm.html"

LSM_INTRINSIC = "INTRINSIC-flux.lsm.html"

# this is where we keep track of generated image names
# the use of a Safelist is encouraged, because otherwise parallel jobs might
# trample over the same file
#IMAGE_LIST = Safelist("$OUTDIR/image.list");

DESTDIR_Template = "${OUTDIR>/}plots-${MS:BASE}"



FITS_L_AXIS, FITS_M_AXIS = "-L", "M"

## 2. Procedures
# Procedures are just python functions
#ms.addbitflagcol(MS)

PI = numpy.pi

# Default dE solution or smoothing intervals
DE_TIME_INTERVAL = 50
DE_FREQ_INTERVAL = 30

# Default dE smoothing solution intervals


DEB_TIME_INTERVAL = 5
DEB_FREQ_INTERVAL = 3
# If True, default i smoothed solution, else make stepwise solution
DE_SMOOTHING     = True

# G + B solution intervals
DI_TIME_INTERVAL = 30
DI_FREQ_INTERVAL = 15


# If True, default is smoothed solution, else make stepwise solution
DI_SMOOTHING     = True

def reset_ms():
  """reset_ms: make a clean start. 
  If FULL_RESET is set (or MS does not exist): does a complete refresh. Removes the measurement set given by the MS variable, 
  untars a fresh copy from the data dir, and does various WSRT-specific initialization.
  If FULL_RESET is not set, then simply clears out flagsets raised during calibration."""
  # init MODEL_DATA/CORRECTED_DATA/etc. columns
  #pyrap.tables.addImagingColumns(MS)
  print "resetting"
  # add a bitflag column (used by MeqTrees)
  #cal.addbitflagcol(MS)
  ms.addbitflagcol(MS)
  # convert UVWs to J2000
  # change the WEIGHT column of redundant baselines, this reduces grating lobes
  # in images
  # use the flag-ms.py tool to copy the FLAG column to the "legacy" bitflag
  ms.flagms(MS,"-Y +L -f legacy -c")

# This steps will make a model image to extract the sources for the calibration 
rad = lambda d: d * np.pi/180.0 # degs to radians
deg = lambda r: r * 180.0/np.pi # radians to degrees


def correlation_factor(src,psf,img,pos_sky,step=None):

    image = pf.open(img)
    psf_image = pf.open(psf)

    hdr = image[0].header
    psf_hdr = psf_image[0].header


    ndim = len(image[0].shape)
    if ndim == 4: 
        image_data = image[0].data[0,0,:,:]
        psf_data_ = psf_image[0].data[0,0,:,:]
    if ndim == 3:
        image_data = image[0].data[0,:,:]
        psf_data_ = psf_image[0].data[0,:,:]
    if ndim == 2:
        image_data = image[0].data[:,:]
        psf_data_ = psf_image[0].data[:,:]
    elif ndim < 2:
        raise ValueError('The FITS file needs at least two dimensions')
    elif ndim > 4:
       raise ValueError('FITS file has more than 4 axes. Aborting')

    # image padding
    n_pix = hdr['NAXIS1']

  
    padding = np.zeros([n_pix+2.0*step,n_pix+2.0*step])        
    padding[step:-step,step:-step] = image_data    
    hdr['CRPIX1'],hdr['CRPIX2'] = (n_pix+2.0*step)/2.0,(n_pix+2.0*step)/2.0
    # taking the positions

    wcs = WCS(hdr,mode='pyfits')
        
    pos_sky = [pos_sky]
    pos = [wcs.wcs2pix(*position) for position in pos_sky]
    step = [step,step]
       
    wcs_psf = WCS(psf_hdr,mode='pyfits')

    #pf.writeto('expore.fits', padding, header=hdr, clobber=True)   
    center_psf_ra0,center_psf_dec0 = wcs_psf.wcs2pix(psf_hdr['CRVAL1'],psf_hdr['CRVAL2'])
    psf_region = psf_data_[center_psf_dec0-step[0]:center_psf_dec0+step[0],center_psf_ra0-step[1]:center_psf_ra0+step[1]]
     
    psf_data = psf_region.flatten()

    for (ra,dec) in pos:
        data = padding[dec-step[0]:dec+step[0],ra-step[1]:ra+step[1]]     
    data_region = data.flatten()
          
    #if len(data_region) and len(psf_data) > 1:
    if data_region.shape == psf_data.shape:
        # computing th correlation matrix and cf is the correlation factor
        cmatrix = np.corrcoef((data_region,psf_data))
        if len(cmatrix) > 0: 
        
            cf =  (np.diag((np.rot90(cmatrix))**2).sum())**0.5/2**0.5
    else:
          cf = 0.0000001
    src.setTag('cf',cf) 
    
    return cf




def make_clean_model(image="${im.RESTORED_IMAGE}",psf_image="${im.PSF_IMAGE}",lsm0="$LSM0",threshold=7):
	
    image, psf_image, lsm0 = interpolate_locals("image psf_image lsm0")
    lsm.pybdsm_search(image,output=lsm0,threshold=threshold)
    catalog = Tigger.load(lsm0)
    src = catalog.sources
    cs =[]
    for i in range(len(src)):
    	position = [deg(src[i].pos.ra),deg(src[i].pos.dec)]
    	c = correlation_factor(src[i],psf=psf_image,img=image,pos_sky=position,step=120)
    	cs.append(c)
    for i in range(len(src)):
    	if cs[i] > 0.6* max(cs):
    		src[i].setTag('dE',True)
    
    catalog.save(lsm0)
  	
	
# The direction independent calibration 
 
def calibrate_DI(msname='$MS', lsmname='$LSM', tdlcon='$TDLCONF', tdlsec='${stefcal.STEFCAL_SECTION}',timeint=None,freqint=None,smooth=None,
              column='$COLUMN', do_dE=False, options={},args=[],**kw):
    """ Calibrate MS """
    
    msname, lsmname, column, tdlcon, tdlsec = \
interpolate_locals('msname lsmname column tdlcon tdlsec')
    
    v.MS = msname
    v.LSM = lsmname 
    args = ["${ms.MS_TDL} ${lsm.LSM_TDL}"] + list(args)
    
    # Direction independent Calibration (Gain and Bandpass calibration)
    if timeint is None:
    	timeint = DI_TIME_INTERVAL;
    if freqint is None:
    	freqint = DI_FREQ_INTERVAL;
    if smooth is None:
    	smooth = DI_SMOOTHING;
    if smooth:
    	kw['stefcal_gain.freqint'] = kw['stefcal_gain.timeint'] = 0;
    	kw['stefcal_gain.freqsmooth'] = freqint;
    	kw['stefcal_gain.timesmooth'] = timeint;
    else:
    	kw['stefcal_gain.freqsmooth'] = kw['stefcal_gain.timesmooth'] = 0;
    	kw['stefcal_gain.freqint'] = freqint;
    	kw['stefcal_gain.timeint'] = timeint;
    	

    
    
    stefcal.stefcal(msname,section=tdlsec,diffgains=False, output="CORR_RES",restore=dict(niter=1000),
           restore_lsm=False,apply_only=False, flag_threshold=(1,.5), options=options, args=args,**kw)
           
# Change apparent flux to intrinsic flux for beam calibration    

def APP_INT(msname="$MS",lsmname="$LSM0",\
                     beams="$BEAM_PATTERN",output="$LSM_INTRINSIC"):
                     
	
                   
    msname,lsmname,beams,output = interpolate_locals('msname lsmname beams output')

    phs_dir = table(msname+"/FIELD").getcol("PHASE_DIR")[0,0,:]
    direction = "J2000,%.3fdeg,%.3fdeg" %((np.rad2deg(phs_dir[0])),(np.rad2deg(phs_dir[1])))
    
    
    
    
    
    
    
   
    
    x.sh("tigger-convert --app-to-int --beam-freq  ${ms.SPW_CENTRE_MHZ} \
         --beam-clip 0.0005  --center $direction --primary-beam='./kat_beams/KAT7_$$(XY)-chan-304-allfreq-$$(realimag).fits' \
          --fits-l-axis=$FITS_L_AXIS \
         --fits-m-axis=$FITS_M_AXIS -f --pa-from-ms $MS:${ms.FIELD} --pa-range -45.92,142.19 $lsmname $output")
         
# Calibate ucing the beam

def calibrate_beam(msname='$MS', lsmname='$LSM0', tdlcon='$TDLCONF', timeint=None,freqint=None,smooth=None,tdlsec='${stefcal.STEFCAL_SECTION3}',
              column='$COLUMN', do_dE=True,options={}, args=[],**kw):
    
    msname, lsmname, column, tdlcon, tdlsec = \
interpolate_locals('msname lsmname column tdlcon tdlsec')
    
    v.MS = msname 
    v.LSM = lsmname
    args = ["${ms.MS_TDL} ${lsm.LSM_TDL}"] + list(args)
    options['pybeams_fits.filename_pattern'] = './kat_beams/KAT7_$$(XY)-chan-304-allfreq-$$(realimag).fits'
    options['critical_flag_threshold']= 20.0
    options['me.e_enable']= 1.0
    options['pybeams_fits.l_axis']=FITS_L_AXIS
    options['pybeams_fits.m_axis']=FITS_M_AXIS
    
    options['pybeams_fits.ampl_interpolation']= 0.0
    options['pybeams_fits.l_beam_offset']= 0.0
    options['pybeams_fits.m_beam_offset'] = 0.0
    options['pybeams_fits.missing_is_null'] = 1
    options['pybeams_fits.normalize_gains'] = 0.0
    options['pybeams_fits.verbose_level']= None
    options['pybeams_fits.sky_rotation']=1
    options['pybeams_fits.spline_order'] = 3
   
    
    
    if timeint is None:
    	timeint = DEB_TIME_INTERVAL;
    if freqint is None:
    	freqint = DEB_FREQ_INTERVAL;
    if smooth is None:
    	smooth = DE_SMOOTHING;
    if smooth:
    	kw['stefcal_diffgain.freqint'] = kw['stefcal_diffgain.timeint'] = 5;
    	kw['stefcal_diffgain.freqsmooth'] = freqint;
    	kw['stefcal_diffgain.timesmooth'] = timeint;
    else:
    	kw['stefcal_diffgain.freqsmooth'] = kw['stefcal_diffgain.timesmooth'] = 3;
    	kw['stefcal_diffgain.freqint'] = freqint;
    	kw['stefcal_diffgain.timeint'] = timeint;
    	
           
    stefcal.stefcal(section=tdlsec,diffgains=False, output="CORR_RES",restore=dict(niter=3000),
           restore_lsm=True,apply_only=False, flag_threshold=(1,0.5),options=options, args=args,**kw) 

    #stefcal.stefcal(section=tdlsec,diffgains=True, output="CORR_RES",restore=dict(niter=3000),
           #restore_lsm=True,apply_only=False, flag_threshold=(1,0.5),options=options, args=args,**kw)

# Calibrate the direction dependent (differential gain calibration )          
           
def calibrate_DD(msname='$MS', lsmname='$LSM', tdlcon='$TDLCONF', timeint=None,freqint=None,smooth=None,tdlsec='${stefcal.STEFCAL_SECTION}',
              column='$COLUMN', do_dE=True,options={},args=[],**kw):
    
    msname, lsmname, column, tdlcon, tdlsec = \
interpolate_locals('msname lsmname column tdlcon tdlsec')
    
    v.MS = msname 
    v.LSM = lsmname
    args = ["${ms.MS_TDL} ${lsm.LSM_TDL}"] + list(args)
    
    if timeint is None:
    	timeint = DE_TIME_INTERVAL;
    if freqint is None:
    	freqint = DE_FREQ_INTERVAL;
    if smooth is None:
    	smooth = DE_SMOOTHING;
    if smooth:
    	kw['stefcal_diffgain.freqint'] = kw['stefcal_diffgain.timeint'] = 40;
    	kw['stefcal_diffgain.freqsmooth'] = freqint;
    	kw['stefcal_diffgain.timesmooth'] = timeint;
    else:
    	kw['stefcal_diffgain.freqsmooth'] = kw['stefcal_diffgain.timesmooth'] = 20;
    	kw['stefcal_diffgain.freqint'] = freqint;
    	kw['stefcal_diffgain.timeint'] = timeint;
    	
           
    stefcal.stefcal(msname,section=tdlsec,diffgains=True, output="CORR_RES",restore=dict(niter=3000),
           restore_lsm=True,apply_only=False, flag_threshold=(1,.5),options=options, args=args,**kw)            
           
           
                
    
def cal_DI(msname="$MS", lsm0='$LSM0', timeint=None,freqint=None,smooth=None,start=1, stop=5,**kw): 
    """ Run pipeline on a single MS"""
   
   
    msname, lsm0 = interpolate_locals("msname lsm0")
   

    # Use the initial model to do the direction independent calibration.
    
    
    calibrate_DI(msname,lsm0)     
    im.make_image(restore=True, psf=True, restore_lsm=False)

    lsm.pybdsm_search(thresh_pix=5 , thresh_isl=3)
    v.LSM = lsm.PYBDSM_OUTPUT
    lsm1= v.LSM
    
    x.sh("tigger-convert --append $LSM $lsm0 $LSMFINAL -f")

    x.sh("tigger-restore ${im.RESTORED_IMAGE} $LSMFINAL ${im.FULLREST_IMAGE} -f ")
    
    
    
    
def cal_BM(msname="$MS",lsm0="$LSM0", timeint=None,freqint=None,start=1, stop=5): 
    """ Run pipeline on a single MS"""
   
   
    msname, lsm0 = interpolate_locals("msname lsm0")
   
    calibrate_beam(msname,lsm0,tdlsec='${stefcal.STEFCAL_SECTION3}')       
    #im.make_image(restore=True, psf=True, restore_lsm=False) 
    #image=im.RESTORED_IMAGE

    #lsm.pybdsm_search(image,thresh_pix= 5, thresh_isl=3)
    #v.LSM = lsm.PYBDSM_OUTPUT
    #lsm1= v.LSM
    #x.sh("tigger-convert --append $LSM $lsm1 $LSM0 $LSMFINAL -f")
    #x.sh("tigger-convert --append $LSM $LSM0 $LSMFINAL -f")
    #x.sh("tigger-restore ${im.RESIDUAL_IMAGE} $LSMFINAL ${im.FULLREST_IMAGE} -f ")
    


def cal_DD(msname="$MS", lsm0='$LSM0', timeint=None,freqint=None,start=1, stop=5): # you can assign start and stop to None
    """ Run pipeline on a single MS"""
   
   
    msname, lsm0 = interpolate_locals("msname lsm0")
   

    # Direction dependent calibration stage and produce the sky model.
    
    calibrate_DD(msname,lsm0,tdlsec='${stefcal.STEFCAL_SECTION}') # these two are global variables      
    #im.make_image(restore=True, psf=True, restore_lsm=False)

    #lsm.pybdsm_search(thresh_pix= 5, thresh_isl=3)
    #v.LSM = lsm.PYBDSM_OUTPUT
    #lsm1= v.LSM
    
    #x.sh("tigger-convert --append $LSM $LSM0 $LSMFINAL -f")
    #x.sh("cp -r Deep_field_2.MS Deep_field-cal-10s2ch.MS")
    

#The statstical analysis at each extracted source position. It calculates the local variance at each source position.

  
def local_variance(lsm="$LSMFINAL",image="${im.RESIDUAL_IMAGE}"):
    lsm, image = interpolate_locals("lsm image")
    x.sh("python stats.py -fit gaussian -s -nb 200 -cat $lsm $image -z 0.9")


# The final call to run all together. Write "pyxis OUTDIR=GIVEAPPFILENAME run_all" to execute the script  
def run_all():
	
	info("############## Make image from the data column, Model the sky and dE tag ######################")
	im.make_image(column="DATA", restore=True, psf=True,restore_lsm=False)
	
	make_clean_model()
	info("######################## Calibrating G+B only ####################")
	
	
	cal_DI(lsm0=LSM0)
	
	im.make_image(column="CORRECTED_DATA", restore=True, psf=True,restore_lsm=False)
	image=im.RESTORED_IMAGE

	lsm.pybdsm_search(image,thresh_pix= 5, thresh_isl=3)
	v.LSM = lsm.PYBDSM_OUTPUT
	lsm1= v.LSM
	
	make_clean_model(image=im.RESTORED_IMAGE,lsm0=LSM,threshold=3)
	
	
	
	info("################# CALIBRATING DIFFENTIAL GAIN (dE): the flyswatter (O.M.Smirnov) ###################")
	cal_DD(lsm0=LSMFINAL) 
	im.make_image(column="CORRECTED_DATA", restore=True, psf=True,restore_lsm=False)
	image=im.RESTORED_IMAGE 
	lsm.pybdsm_search(image,thresh_pix= 5, thresh_isl=3)
	v.LSM = lsm.PYBDSM_OUTPUT
	lsm1= v.LSM
	make_clean_model(image=im.RESTORED_IMAGE,lsm0=LSM,threshold=3)
	
	info("############## Change the apparent flux to intrinsic flux for beam calibration ######################")
	
	APP_INT(lsmname=LSM0,output=LSM_INTRINSIC)
	
	info("############## Do only G+B with beam and with beam ######################")
	cal_BM(lsm0=LSM_INTRINSIC)
	
	
	
	x.sh("tigger-restore ${im.RESTORED_IMAGE} $LSM0 ${im.FULLREST_IMAGE} -f ")
	
	
	
   
 
     



