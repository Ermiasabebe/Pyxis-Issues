import pyfits
import numpy as np
from scipy.ndimage.interpolation import zoom
import glob

import os

import sys

fitfilenames = glob.glob("*.fits") # The lists of fits files 



def Get_new_beam(filename):

	
	sourcefile = pyfits.open(filename) # Open the simulated fits files 

# The input information about the simulated beam.

	header = sourcefile[0].header
	sourcedata =  sourcefile[0].data
	
	init_freq0  = header['CRVAL3'] # The initial frequency 
	d_initfreq  = header['CDELT3'] # Channel width 
	nchan = header['NAXIS3'] # channel numbers. It is 304 in our case.
	newdata = []
	
# The input informatnion for the beam that will be scaled.	
	new_freq0 = 1707351562.5 # Intial frequency 
	d_newfreq = 781250.0 # Channel width 
	new_freqf = 1944070312.5 + d_newfreq # The final frequency 
	
	new_chan = np.arange(new_freq0,new_freqf,d_newfreq)
	
	z,x,y = sourcedata.shape
	
	header['CRVAL3'] = new_chan[0] # The newly formed fits file initial frequency 
	header['CDELT3'] = 781250.0 # The newly formed fits channel width
	
	
	
	for ichan in range(len(new_chan)):

		tmpdata = np.zeros((x,y))
	# The factor to scale the beam. If the factor is less than one, it will scale down.	
		factor = init_freq0/new_chan[ichan] 
		
		temp = zoom(sourcedata[ichan,:,:], factor)
		
		diff = np.array(sourcedata[ichan,:,:].shape) - np.array(temp.shape)
		
		tmpdata[diff[0]/2:-diff[0]/2,diff[0]/2:-diff[0]/2] = temp
		
		newdata.append(tmpdata)
		
	
	pyfits.writeto(filename.split('fits')[0]+'SCALED.fits', newdata, header=header, clobber=True)

for i in range(len(fitfilenames)):
    
    Get_new_beam(fitfilenames[i])
    
    

	


	
	
