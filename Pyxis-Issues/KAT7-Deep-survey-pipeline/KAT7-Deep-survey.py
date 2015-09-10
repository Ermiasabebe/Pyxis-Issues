#!/bin/env python
#
# Reduction script for KAT7 Deep survey data 

#
# These must be done before invoking the script/pipeline. For example

#
#
#
import time
import os.path, math
from sys import exit
import timeit

scriptmode = True


# Defining some parameters to be used later


c=2.99792458E8
D=12 # diameter of KAT-7 dish

# Boltzmann's Constant, Converting to Jy
k = 1.38e03;

# Efficiency of the telescope

eff = 0.7

tsys = 25

degrad=180/pi
arcsecond_deg=3600

# attempt at auto script
#
print "=============================================================="
print "KAT-7 Reduction Pipeline"
print "=============================================================="
print Version
print "=============================================================="
print "Warning: This script was written for CirX1 and may require some "
print "tweaking for your specificexperiment"
print ""
print "NOTE: that most input are expecting a string, so no need to use \'\' "
print "Avoid pressing Ctrl+C during a run...wait when for the return key"
print "=============================================================="

msfile_name = "1355337337_20121212_1_1822.hh_vv"

prefix = "Reduce" 

ref_ant='ant5'

subchans ="0:250~350"


try: msfile_name
except NameError:
    raise NameError('msfile variable not specified')
 
if not os.path.exists(msfile_name+'.ms'):
    raise IOError('msfile file does not exist.')
else:
    msfile = msfile_name+'.ms'
 


########################################################################
# definitions of standard amplitude calibrators
#
acal_std_name = {}
acal_std_name['3C286'] = frozenset(['3C286','3C 286',
                                    '1328+30','1328+307','B1328+307','B1328+30',
                                    '1331+305','J1331+305','J1331+3030'])
acal_std_name['3C138'] = frozenset(['3C138','3C 138',
                                    '0518+16','0518+165','0518+1635','B0518+16','B0518+165','B0518+1635',
                                    '0521+166','0521+1638','J0521+166','J0521+1638'])
acal_std_name['3C147'] = frozenset(['3C147','3C 147',
                                    '0538+49','0538+498','0538+4949','B0538+49','B0538+498','B0538+4949',
                                    '0542+498','0542+4951','J0542+498','J0542+4951'])
acal_std_name['3C48'] = frozenset(['3C48','3C 48',
                                   '0134+32','0134+329','0134+3254','B0134+32','B0134+329','B0134+3254',
                                   '0137+331','0137+3309','J0137+331','J0137+3309'])
acal_std_name['1934-638'] = frozenset(['1934-638','PKSJ1939-6342','PKS 1934-638'])


fluxcal = '0'

        
phasecal = '1'

srcfield = '2'

FOIfields = srcfield+','+phasecal+','+fluxcal #fields of interest

print "=============================================================="
print "The fields of interest are "+FOIfields
print "=============================================================="

       
     
########################################################################
#
# Initial stuff.

os.system('rm -rf '+prefix+'.*')

print "=============================================================="
print "STARTING REDUCTION OF "+msfile+" AT "
time.asctime
print ""
print '------------Cleaning up old files ('+prefix+'*)--------------'
print "=============================================================="
# clean up old files

os.system('rm -rf '+prefix+'_*')

#
#===================================================================
# listobs to find out obs srcs ...
#
print "=============================================================="
print "Listobs"
#print "=============================================================="

default('listobs')
listobs(vis = msfile)

#print "==============================================================" 
print "\033[1;48mUse listobs (outputted in the logger) to find the field names\033[1;m"
print "=============================================================="

#if scriptmode: 
    #user_check=raw_input('Return to continue script <==|\n')

print "=============================================================="
print " "
print "Extract various useful quantities ...."
print "These are not fully used right now but shall be in the future..."
print "Taken from Lazio's script and tweaked for KAT-7"
print " "
print "=============================================================="

tb.open(msfile)
uvw=tb.getcol('UVW')
interval=tb.getcol('INTERVAL')
tb.close()
 
source_table=msfile + '/SOURCE'
tb.open(source_table)
calcode=tb.getcol("CODE")
source_id=tb.getcol("SOURCE_ID")
source_name=tb.getcol("NAME")
tb.close()
 
spw_table=msfile + '/SPECTRAL_WINDOW'
tb.open(spw_table)
freq_list=tb.getcol("REF_FREQUENCY")
channel_width=tb.getcol("CHAN_WIDTH")
num_chan=tb.getcol("NUM_CHAN")
tb.close()


#plotants(vis=msfile,figfile=prefix+'_antenna_layout.png')
#plotuv(vis=msfile,figfile=prefix+'_uvcov.png',field='2')


freq=freq_list[0]
if (1E9 < freq) and (freq < 2E9):
    band = 'L'
elif (2E9 < freq) and (freq < 4E9):
    band = 'S'
elif (4E9 < freq) and (freq < 8E9):
    band='C'
elif (8E9 < freq) and (freq < 12E9):
    band='X'
elif (12E9 < freq) and (freq < 40E9):
    band = 'K'
elif freq > 40E9:
    band = 'Q'

print "=============================================================="
print "Observations are determined to be in the ", band, " band."
print "==============================================================" 



wavelength=c/freq
primary_beam=(1.02*wavelength/D)*degrad   # HPBW, from Napier (1999)
FoV=(primary_beam/math.sqrt(2))*arcsecond_deg
 
b=[]
for u, v, w in zip(uvw[0], uvw[1], uvw[2]):
    b.append(math.sqrt(u*u+v*v)/wavelength)
bmax=max(b)
                                          # HPBW, synthesized beam,
                                          # from Bridle & Schwab (1999)
synthesized_beam=(1.2/bmax)*degrad*arcsecond_deg

print "=============================================================="
print "Maximum baseline is ",round(bmax,3), "lambda"
print "The synthesized beam is ",round((synthesized_beam/60),3),"arcmin"
print "==============================================================" 





setjy(vis=msfile,field=fluxcal,spw='0',scalebychan=True,standard='Perley-Butler 2010', fluxdensity=-1)


gain_table0=prefix+'.gaincal0'

print 'Gain table is - '+gain_table0

default('gaincal')
vis=msfile
caltable=gain_table0
field=fluxcal
solint='int'	
spw=subchans
refant=ref_ant
minblperant=4
minsnr=3.0
calmode='p'
gaintable=['']

gaincal()

plotcal(caltable=gain_table0,xaxis='time',yaxis='phase',field=fluxcal,spw='0',plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-Gaincal0-'+fluxcal+'-phase-ants.png',iteration='antenna',subplot=334)
clearpanel


delay_table0=prefix+'.delaycal0'

print '\n---> Delay corrections - '+delay_table0

default('gaincal')

vis=msfile
caltable=delay_table0
field=fluxcal
refant=ref_ant
spw='0'
gaintype='K'
solint='inf'
combine='scan'
minblperant=4
minsnr=3.0
gaintable=[gain_table0]

gaincal()
        
plotcal(caltable=delay_table0,xaxis='antenna',yaxis='delay',figfile=prefix+'_plotcal-delay0-'+fluxcal+'.png')

clearpanel


print "=============================================================="
print"# Making the bandpass now"
print"=============================================================="
	
clearstat()
	
band_pass_table0=prefix+'.bandpasscal0'
	
print '\n---> Bandpass - '+band_pass_table0

default('bandpass')
vis=msfile
caltable=band_pass_table0
field=fluxcal
solint='inf'
combine='scan'
refant=ref_ant
minblperant=4
minsnr=3.0
bandtype = 'B'
interp =['nearest']
gaintable=[gain_table0,delay_table0]
solnorm=True

bandpass()

#plotbandpass(caltable = band_pass_table0, xaxis = 'freq',yaxis = 'both', subplot=42, figfile=prefix+'_plotcal-BPOLY-'+fluxcal+'-amp.png')

#clearpanel

#plotcal(caltable=band_pass_table0,xaxis='chan',yaxis='amp',field=fluxcal,iteration='antenna',figfile=prefix+'_plotcal-B0-'+fluxcal+'-amp.png',subplot=324)

clearpanel

print '\n---> Applycal: the bandpass... '+band_pass_table0
	
applycal(vis=msfile, gaintable=[band_pass_table0,delay_table0],calwt=False)

#clearcal()
print "=========================================================="
print "Gain calibration (GAINCAL) - starting ...."
print "=========================================================="

gain_table1=prefix+'.gaincal1'
	
print '\n---> Gain - '+gain_table1+' - '+fluxcal


clearstat()
	
default('gaincal')

vis=msfile
caltable=gain_table1
field=fluxcal
spw='0'
solint='int'
combine='scan'
refant=ref_ant
gaintype='G'
calmode='ap'
solnorm=F
gaintable=[band_pass_table0,delay_table0]

gaincal()


#plotcal(caltable=gain_table1,xaxis='time',yaxis='phase',field=fluxcal,spw='0',plotrange=[-1,-1,-180,180],
              #figfile=prefix+'_plotcal-G1'+fluxcal+'-phase-ants.png',iteration='antenna', subplot=324)

clearpanel



default('gaincal')

vis=msfile
caltable=gain_table1
field=phasecal
spw='0'
solint='int'
combine='scan'
refant=ref_ant
gaintype='G'
calmode='ap'
solnorm=F
append=T
minsnr=1.0
gaintable=[band_pass_table0,delay_table0]

gaincal()


plotcal(caltable=gain_table1,xaxis='time',yaxis='phase',field=phasecal,spw='0',plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-Gaincal1-'+phasecal+'-phase.png',iteration='antenna',subplot=334)
                
                
clearpanel                


flux_table1=prefix+'.fluxscalecal1'
print '---> Fluxscale - '+flux_table1

default('fluxscale')
vis=msfile
caltable=gain_table1
fluxtable=flux_table1
reference=[fluxcal]
transfer=[phasecal]

fluxscale()


# Applying calibrations solutions (APPLYCAL)
#
print '---> Applycal: '+fluxcal
applycal(vis=msfile,gaintable=[flux_table1,band_pass_table0,delay_table0],spw='0',parang=False,field=fluxcal,gainfield=[fluxcal,''],interp=['nearest','',''],calwt=F)

print '---> Applycal: '+phasecal
applycal(vis=msfile,gaintable=[flux_table1,band_pass_table0,delay_table0],spw='0',parang=False,field=phasecal,gainfield=[phasecal,''],interp=['nearest','',''],calwt=F)

print '---> Applycal: '+srcfield
applycal(vis=msfile,gaintable=[flux_table1,band_pass_table0,delay_table0],spw='0',parang=False,field=srcfield,gainfield=[fluxcal+','+phasecal,''],interp=['linear','',''],calwt=F)

print '\n**************************'
print '    END of Calibration'
print '**************************'



srcfield = srcfield

target_vis = "Deep_field" +'_'+srcfield +".MS"

os.system('rm -rf Deep_field*')


print '\n----Splitting the source data - for time bin of 10s and channel width of 2' +target_vis

# split for time bin of 10 s and every two channel

split(vis=msfile,outputvis=target_vis,datacolumn='corrected',field=srcfield, spw='0',timebin='10s', width='2')


plotms(vis = target_vis, xaxis='channel')





































