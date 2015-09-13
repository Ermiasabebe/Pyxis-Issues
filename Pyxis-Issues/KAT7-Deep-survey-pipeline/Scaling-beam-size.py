import pyfits
import numpy as np
from scipy.ndimage.interpolation import zoom
import glob

import os

import sys

fitfilenames = glob.glob("*.fits") 



def get_new_beam(filename):

	
	sourcefile = pyfits.open(filename)

# The input information about the simulated beam.

	header = sourcefile[0].header
	sourcedata =  sourcefile[0].data
	
	init_freq0  = header['CRVAL3'] # The initial frequency 
	d_initfreq  = header['CDELT3'] # Cahnnel width 
	nchan = header['NAXIS3'] # channel numbers. It is 304 in our case.
	newdata = []
	
# The input informatnion for the beam that will be scaled.	
	new_freq0 = 1707351562.5 # Intial frequency 
	d_newfreq = 781250.0 # Channel width 
	new_freqf = 1944070312.5 + d_newfreq # The final frequency 
	
	new_chan = np.arange(new_freq0,new_freqf,d_newfreq)
	
	z,x,y = sourcedata.shape
	
	header['CRVAL3'] = new_chan[0]
	header['CDELT3'] = 781250.0
	
	
	
	for ichan in range(len(new_chan)):

		tmpdata = np.zeros((x,y))
		
		factor = init_freq0/new_chan[ichan] # The factor to sacle the beam. If the factor is less than one, it will scale down.
		
		temp = zoom(sourcedata[ichan,:,:], factor)
		
		diff = np.array(sourcedata[ichan,:,:].shape) - np.array(temp.shape)
		
		tmpdata[diff[0]/2:-diff[0]/2,diff[0]/2:-diff[0]/2] = temp
		
		newdata.append(tmpdata)
		
	
	pyfits.writeto(filename.split('fits')[0]+'SCALED.fits', newdata, header=header, clobber=True)

for i in range(len(fitfilenames)):
    
    get_new_beam(fitfilenames[i])
    
    

	


	
	