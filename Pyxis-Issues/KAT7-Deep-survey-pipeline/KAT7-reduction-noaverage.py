#!/bin/env python
#
# Reduction script for KAT7 data 
#
# These must be done before invoking the script/pipeline. For example

#
#
#
import time
import os.path, math
from sys import exit

scriptmode = True

#
Version='v. 0.0 NO  2012-08-22'
#
# Defining some parameters to be used later
#

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

msfile_name = raw_input('\nEnter the measurement set without dot ms: ') # '1341077121.full_pol.ms'
prefix = raw_input('\nEnter the prefix for this reduction (e.g Cirx1): ') # 'CirX1_school'

clip_sigma = raw_input('\nEnter sigma level above which you will consider RFIs: ')

ref_ant = raw_input('\nEnter your reference antenna (e.g. ant5): ')

#
# User specified parameters
#
try: clip_sigma
except NameError:
    clip_sigma=10.0
#
# Chwcking the existence of the msfile
# 
try: msfile_name
except NameError:
    raise NameError('msfile variable not specified')
 
if not os.path.exists(msfile_name+'.ms'):
    raise IOError('msfile file does not exist.')
else:
    msfile = msfile_name+'.ms'
 
try: ref_ant
except NameError:
    raise NameError('Reference antenna not specified in variable ref_ant.')

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

#
########################################################################
#
# Initial stuff.
#
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

if scriptmode: 
    user_check=raw_input('Return to continue script <==|\n')

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

#
# Selecting calibrators and sources
#

print "=============================================================="
print "Enter the fields of interest:"
print "=============================================================="

fluxcal = raw_input('\nEnter the src id of Primary calibrator: ') 
print 'Your calibrator is', source_name[int(fluxcal)]
user_check = 'y'

flatlist = []
while (user_check.lower() == 'y'):

    for a in acal_std_name.values(): flatlist += list(a)   # I am cheating here
    if source_name[int(fluxcal)] in flatlist:
        print "Your calibrator is in standard format"
    else:
        raise RuntimeError("Your calibrator does not seem to be a standard name, use browsetable to correct")
        exit()
    
    user_check = raw_input('Is there more than one primary calibrator? (Y/N): ')
    if user_check.lower() == 'y':
        newcal = raw_input(' What is it\'s code: ---> ')
        fluxcal = fluxcal+','+newcal
        
phasecal = raw_input('\nEnter the Gain calibrator: ')

user_check = raw_input('Is there more than one gain calibrator? (Y/N): ')
if user_check.lower() == 'y':
	newcal = raw_input(' What is it\'s code: ---> ')
	phasecal = phasecal+','+newcal

srcfield = raw_input('\nEnter the science Target field: ')

FOIfields = srcfield+','+phasecal+','+fluxcal #fields of interest

print "=============================================================="
print "The fields of interest are "+FOIfields
print "=============================================================="
print "  "
print "=============================================================="
print "Plotting the antenna layout"
print "=============================================================="

plotants(vis=msfile,figfile=prefix+'_antenna_layout.png')

if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
 
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
  
# 
# How much time-average smearing can be tolerated?
# Assume no more than time_smearing_loss peak
# intensity loss over half of primary beam.
# Follow Section 2 of Bridle & Schwab (1999).
# Note undocumented option that the amount of time-average
# smearing loss can be specified by user.
#
try: time_smearing_loss
except NameError:
    time_smearing_loss=0.01
 
tau=math.sqrt(time_smearing_loss*1E9)*(synthesized_beam/FoV)
 
arbitrary_maximum=30
if tau > arbitrary_maximum:
    tau=arbitrary_maximum
 
dt=min(interval)
#
# tau is allowed value, dt is actual (minimum)
# (could be an issue if baseline-dependent
#  correlator accumulation used)
# make sure that tau is an integer multiple of dt
# 
tau=dt*math.floor(tau/dt)
if tau < dt:
    tau=dt
print "=============================================================="  
print "Data will be averaged in time."
print "Original time sampling [s]  : ",round(dt,3)
print "Averaging time [s]          : ",round(tau,3)
#tau = tau*0.5
print "Averaging time to be used[s]: ",round(tau,3)
print "=============================================================="
 
#
# How much bandwidth smearing can be tolerated?
# Assume no more than band_smearing_loss peak
# intensity loss over half of primary beam.
# Follow Section 1 of Bridle & Schwab (1999).
# Assume square bandpass, no taper, expand resulting
# sine integral to lowest order
# Note undocumented option that the amount of time-average
# smearing loss can be specified by user.
#
try: band_smearing_loss
except NameError:
    band_smearing_loss=0.01
 
eta_band=3.79
delta_nu=freq*(2/eta_band)*(synthesized_beam/FoV)*math.sqrt(18*band_smearing_loss)
#
# delta_nu is allowed value,
# figure out even divisibles of the actual
# value, stored in num_chan that is smaller than
# delta_nu
# potential bug if uneven channel widths
# used in different spectral windows
# 
dnu=channel_width[0][0]
nchan_log2=math.log(num_chan[0],2)
for i in range(int(nchan_log2-1), 0, -1):
    if (dnu*math.pow(2,i)) < delta_nu:
        nchav=math.pow(2,i)
        break
 
if nchav < 1:
    nchav=1
 
print "=============================================================="
print "Data will be averaged in frequency."
print "Original channel width [kHz]: ",round(min(channel_width[0])/1E3,0)
print "Averaged channel width [kHz]: ",round((nchav*min(channel_width[0]))/1E3,0)
print "Number of channels: ", nchav
print "=============================================================="

if scriptmode: 
    user_check=raw_input('Return to continue script <==|\n')

#
# Compress data for faster processing.
# Throw away edge channels.
#
 
bchan=num_chan[0]*0.05
echan=num_chan[0]*0.95
chans='*:%d~%d'%(bchan, echan)
 
#Av_msfile = prefix +'_avsplit.ms'
tave='%.fs'%tau

scan_ave = raw_input('\nIf you want to average in scan, enter it now e.g. 100~200:')

print "=============================================================="
print "Splitting for time averaged data over..."
print '\033[1;41mDo NOT press Ctrl+C...wait for this to be over\033[1;m'
print "=============================================================="



default('flagdata')
flagdata(vis=msfile, mode='clip', clipzeros=True, flagbackup=False)

if scriptmode: 
    user_check=raw_input('Return to continue script <==|\n')


#------Look at the data to find bad channels - these are generally the edge channels where the signal drops off, and ch0 is always RFI

clearstat()

print "=============================================================="
print "Look at the data to find bad channels (keep this open)..."
print "Using correlation XY to check RFI. "
print "=============================================================="


plotms(vis=msfile,field=fluxcal,spw='0',xaxis='channel',yaxis='amp',correlation ='XY,YX')

if scriptmode: 
    user_check=raw_input('Return to continue script <==|\n')

#------Remove known bad channels (huge rfi or fall-off)
clearstat()

badchannels = raw_input('Bad channel range e.g. 0:0~10;40~50 OR a charcater for none: ')

print "=============================================================="
print "Removing bad channels, eg fall-off etc..."
print "=============================================================="

while (badchannels.isdigit()):

    default('flagdata')
    
    flagdata(vis=msfile,mode='manual',field=FOIfields,spw=badchannels,display='data',flagbackup=True)

    if scriptmode: 
        user_check=raw_input('Return to continue script <==|\n')
 
user_check = raw_input('From PLOTMS; Are there obvious baseline- or antenna-specific problems? (Y/N): ')

#possible while loop for flagging

while (user_check.lower() == 'y'):
	antORbase = raw_input('Faulty antenna or baseline: ')
	spwORchan = raw_input('Faulty spw or channel(s) for above: ')
	clearstat()
	print " =============================================================="
        print "-- Flagging '+antORbase+' on spw '+spwORchan+'..."
        print "=============================================================="
	flagdata(vis=msfile,mode='manual',field=FOIfields,spw=spwORchan,antenna=antORbase,display='data',flagbackup=True)
	
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')

        clearstat()
	user_check = raw_input('From PLOTMS; More baseline- or antenna-specific problems? (Y/N): ')
        
user_check = raw_input('Have you closed PLOTMS? (Y/N): ')

while (user_check.lower() == 'n'):
	user_check = raw_input('Have you closed PLOTMS? (Y/N): ')

clearstat()


#------Automatic RFI detection
print "=============================================================="
print "Will try Automatic RFI flagging at later stage..."
print "=============================================================="

# tflagdata(vis=Av_msfile,mode='tfcrop',field=FOIfields,spw='0',datacolumn='data',freqcutoff=2.0,timecutoff=2.0,display='both',flagbackup=True)

clearstat()

user_check = raw_input('Have you closed PLOTMS? (Y/N): ')
while (user_check.lower() == 'n'):
	user_check = raw_input('Have you closed PLOTMS? (Y/N): ')

clearstat()
default('plotms')
print '\n-- Use PLOTMS to choose a reference antenna...'

vis=msfile
xaxis='channel'
yaxis='amp'
averagedata=F
transform=F
extendflag=F
iteraxis='antenna'
selectdata=T
field=fluxcal
spw='0'
avgbaseline=T

if scriptmode: 
    user_check=raw_input('Return to continue script <==|\n')

plotms()

subchans = raw_input('Choose a subset of channels for bandpass calibration (eg 0:50~55): ') 

user_check = raw_input('Have you closed PLOTMS? (Y/N): ')
while (user_check.lower() == 'n'):
	user_check = raw_input('Have you closed PLOTMS? (Y/N): ')

clearstat()

#
# Use setjy to find the path to model images
#
#setjy(vis=msfile, listmodels=True)

print "=============================================================="
print "=============START CALIBRATION LOOP==========================="
print "=============================================================="

doCalibration = True
calstep = 0

while doCalibration:
	calstep = calstep + 1
	print '\n**************************'
	print '*** Calibration Run ',calstep
	print '**************************'
	clearstat()
        print "=============================================================="
        print "Ensure MODEL_DATA = unity and CORRECTED_DATA = DATA, before"
        print "attempting calibration"
        print "=============================================================="
	print "----Clearing calibrations----"
        print "=============================================================="
	clearcal(vis=msfile,field=FOIfields,spw='0')

        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
	print "=============================================================="
        print "Removing old tables"
        print "=============================================================="
	os.system('rm -rf *cal*')
	user_check=raw_input('If unable to remove, go do it manually. Press Enter to continue')
	moreflagging = False
        print "=============================================================="
	print "Setting core flux level"
        print "=============================================================="
#	setjy(vis=Av_msfile,field=fluxcal,spw='0', modimage=fluxcal+'.im',scalebychan=True,standard='Perley-Butler 2010', fluxdensity=-1)
	setjy(vis=msfile,field=fluxcal,spw='0',scalebychan=True,standard='Perley-Butler 2010', fluxdensity=-1)

        print "==============================================================" 
        print "Please do check your calibration results before proceeding"
        print "=============================================================="
        
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')

        print "=============================================================="    
        print "Bandpass calibration (GAINCAL)"
        print "=============================================================="
	gain_table0=prefix+'.gcal0'
	print 'gaincal table is - '+gain_table0

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

        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')

        gaincal()

	plotcal(caltable=gain_table0,xaxis='time',yaxis='phase',field=fluxcal,spw='0',plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-G0-'+fluxcal+'-phase-ants.png',iteration='antenna',subplot=334)
	
        user_check=raw_input('Press Enter to continue ')

	plotcal(caltable=gain_table0,xaxis='time',yaxis='phase',field=fluxcal,spw='0',plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-G0-'+fluxcal+'-phase.png',iteration='antenna',subplot=334)
	
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
 
	user_check=raw_input('Did you notice data that needs to be flagged? (Y/N): ')

	clearstat()
	if user_check.lower() == 'y':
            print"Look for 'dodgy' data.... take note of time slots and flag in plotms"
            clearstat()
            default('plotms')
            plotms(vis=msfile,field=fluxcal,xaxis='time',yaxis='amp')
            userinp = raw_input('Flag out the dodgy data. Press Enter to continue.')
            continue

        print "=============================================================="
        print "# Delay Calibration - should get rid of the slopes"
        print "=============================================================="

	delay_table0=prefix+'.kcal0'
	default('gaincal')
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

        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')

        gaincal()
        
	plotcal(caltable=delay_table0,xaxis='antenna',yaxis='delay',figfile=prefix+'_plotcal-K0-'+fluxcal+'.png',iteration='antenna',subplot=334)

	user_check=raw_input('Press Enter to continue')

        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
        
        print "=============================================================="
        print"# Making the bandpass now"
        print"=============================================================="
	
        clearstat()
	band_pass_table0=prefix+'.bcal0'
	default('bandpass')
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


        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')

        #bandpass()
	#plotbandpass(caltable = band_pass_table0, xaxis = 'freq',yaxis = 'both', subplot=42, figfile=prefix+'_plotcal-BPOLY-'+fluxcal+'-amp.png')

	plotcal(caltable=band_pass_table0,xaxis='chan',yaxis='amp',field=fluxcal,figfile=prefix+'_plotcal-B0-'+fluxcal+'-amp.png',iteration='antenna',subplot=334)
	
	
	user_check=raw_input('Press Enter to continue')

	plotcal(caltable=band_pass_table0,xaxis='chan',yaxis='phase',field=fluxcal,plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-B0-'+fluxcal+'-phase.png',iteration='antenna',subplot=334)
	user_check=raw_input('Press Enter to continue')

	plotcal(caltable=band_pass_table0,xaxis='chan',yaxis='phase',field=fluxcal,plotrange=[-1,-1,-180,180],
                iteration='antenna',figfile=prefix+'_plotcal-B0-'+fluxcal+'-phase-ants.png',subplot=334)
	
        user_check=raw_input('Press Enter to continue')

	user_check=raw_input('Did you notice data that needs to be flagged? (Y/N): ')

	clearstat()
        
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
            
	if user_check.lower() == 'y':
	# Look for 'dodgy' data.... take note of time slots and flag in plotms
		clearstat()
		default('plotms')
		plotms(vis=msfile,field=fluxcal,xaxis='time',yaxis='phase')
		userinp = raw_input('Flag out the dodgy data. Press Enter to continue.')
		continue

	print '\n---> Applycal: the bandpass... '+band_pass_table0
	applycal(vis=msfile, gaintable=[band_pass_table0,delay_table0],calwt=False)



        print "=========================================================="
        print "Gain calibration (GAINCAL) - starting ...."
        print "=========================================================="
	gain_table1=prefix+'.gcal1'
	print '\n---> Gain - '+gain_table1+' - '+fluxcal

#
# Determine the complex gains for the primary calibrator
#
	clearstat()
	default('gaincal')

	gaincal(vis=msfile,caltable=gain_table1,field=fluxcal,spw=subchans,solint='int',combine='scan',
                refant=ref_ant,gaintype='G',calmode='ap',solnorm=F,gaintable=[band_pass_table0,delay_table0])
        
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
            
	plotcal(caltable=gain_table1,xaxis='time',yaxis='phase',field=fluxcal,spw='0',plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-G1-'+fluxcal+'-phase.png',iteration='antenna',subplot=334)
	
	
	user_check=raw_input('Press Enter to continue')

	#plotcal(caltable=gain_table1,xaxis='time',yaxis='phase',field=fluxcal,spw='0',plotrange=[-1,-1,-180,180],
                #figfile=prefix+'_plotcal-G1'+fluxcal+'-phase-ants.png',iteration='antenna',subplot=334)
	
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
            
	user_check=raw_input('Did you notice data that needs to be flagged? (Y/N): ')

	clearstat()
	if user_check.lower() == 'y':
	# Look for 'dodgy' data.... take note of time slots and flag in plotms
		clearstat()
		default('plotms')
		plotms(vis=msfile,field=fluxcal,xaxis='time',yaxis='phase')
		userinp = raw_input('Flag out the dodgy data. Press Enter to continue.')
		continue

        print "==============================================================" 
        print "Determine the complex gains for the phase calibrator."
        print "==============================================================" 

	clearstat()
	print '\n---> Gain - '+gain_table1+' - '+phasecal

	gaincal(vis=msfile,caltable=gain_table1,field=phasecal,spw='0',solint='int',combine='scan',refant=ref_ant,
                gaintype='G',calmode='ap',solnorm=F,append=T,minsnr=1.0,gaintable=[band_pass_table0,delay_table0])
	
	plotcal(caltable=gain_table1,xaxis='time',yaxis='phase',field=phasecal,spw='0',plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-G1-'+phasecal+'-phase.png',iteration='antenna',subplot=334)
	
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')
          
	plotcal(caltable=gain_table1,xaxis='time',yaxis='phase',field=phasecal,spw='0',plotrange=[-1,-1,-180,180],
                figfile=prefix+'_plotcal-G1-'+phasecal+'-phase-ants.png',iteration='antenna',subplot=334)
	
        
        if scriptmode: 
            user_check=raw_input('Return to continue script <==|\n')

	user_check=raw_input('Did you notice data that needs to be flagged? (Y/N): ')

	clearstat()
	if user_check.lower() == 'y':
	# Look for 'dodgy' data.... take note of time slots and flag in plotms
		clearstat()
		default('plotms')
		plotms(vis=msfile,field=phasecal,xaxis='time',yaxis='phase')
		userinp = raw_input('Flag out the dodgy data. Press Enter to continue.')
		continue

	doCalibration = False
	

print "================================================================"
print "=============END OF CALIBRATION LOOP============================"
print "================================================================"


if scriptmode: 
    user_check=raw_input('Return to continue script <==|\n')
          
#========================================================================
# Using primary calibrator to know system response (FLUXSCALE)
#
flux_table1=prefix+'.fcal1'
print '---> Fluxscale - '+flux_table1

fluxscale(vis=msfile,caltable=gain_table1,fluxtable=flux_table1, reference=[fluxcal],transfer=[phasecal])

print "================================================================"
print "Check the logger that you are getting the right expected fluxes"
print "================================================================"

if scriptmode: 
    user_check=raw_input('Return to continue script <==|\n')


#========================================================================
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

print 'Check the corrected data for the calibrators - should have amp centred on the flux given in setjy and fluxscale, and phases centred on zero... If not, something has gone wrong!'
clearstat()
default('plotms')
plotms(vis=msfile,xaxis='phase',yaxis='amp',xdatacolumn='corrected',ydatacolumn='corrected',field=fluxcal+','+phasecal)
user_check=raw_input('\nPress Enter to continue...')
clearstat()
###################################################################################################################

# split the file based on the source ID
srcfield = raw_input('\nEnter the field number of the source to be splitted: ')
target_vis = 'Deep_survey' +'_'+srcfield +'.MS'

print '\n----Splitting the source data - '+target_vis
split(vis=msfile,outputvis=target_vis,datacolumn='corrected',field=srcfield, spw='0',timebin='20s', width='2')





clearstat()

plotms(vis=target_vis,field='',spw='',xaxis='uvwave',yaxis='amp')

#                                                                                                                                               
# Clip any egregious data.                                                                                                                      
#                                                                                                                                               
Vstat=visstat(vis=target_vis,selectdata=F)
Vclip=float(clip_sigma)*Vstat['DATA']['rms']

flagdata(vis=target_vis,mode='clip',correlation='ABS XX, ABS YY',clipminmax=[0.0,Vclip])

#
# to put the theoretical rms here: Borrowed from Brad the part below
#

msfile = target_vis

# Field where to calculate the expected rms 

field_id = srcfield

# Get Number of Antennas N:
tb.open(msfile+'/ANTENNA');
N = len(tb.getcol('NAME'));

 #Get Dish Diameter
D = tb.getcol('DISH_DIAMETER')[0]; 
A = pl.pi*(D/2.)**2;
tb.close();

# Get Total Bandwidth dv
tb.open(msfile+'/SPECTRAL_WINDOW');
dv = tb.getcol('TOTAL_BANDWIDTH')[0];
Nchan = tb.getcol('NUM_CHAN')[0];
dv_per_chan = tb.getcol('CHAN_WIDTH')[0];
tb.close();

field_id = 0

# Total Integration Time
tb.open(msfile);
Nints = len(pl.where(tb.getcol('FIELD_ID')==field_id)[0]);
dt = tb.getcol('EXPOSURE')[0];
T = Nints * dt; 

pss_N = pl.sqrt(2)*tsys*k;
pss_D = eff * A * pl.sqrt(N*(N-1)*dv*T);
pss = 1.0E3* pss_N / pss_D;

print "\n"
print "Point Source Sensitivity Parameters for Source "+str(field_id);
print "-------------------------------------------------"
print "Total Integration Time [s] = "+str(T);
print "Total Bandwidth [Hz] = "+str(dv);
print "Number of Channels = "+str(Nchan);
print "Number of Antennas = "+str(N);
print "Approximate Dish Area [m^2] = "+str(A);
print "User Input Tsys [K] = "+str(tsys);
print "\n"
print "Your estimated Point Source Sensitivity [mJy] = "+str(pss);
print "\n"
fluxthresh = raw_input('*** What flux threshold do you want to go down to (eg 0.15mJy for 15min)? ---> ')

imname =prefix+'_'+srcfield+'_I'

#
# First make a Hogbom CLEAN to pik up the mask
#

cells='%.3farcmin'%(synthesized_beam/(5*60))

print ""
print 'Making an image with '+str(cells)
print 'Image size will be about 3deg '
print 'Image will be saved to',imname
print ""

#clean(vis=target_vis,imagename=imname,field='',spw='',niter=2500,gain=0.1,threshold=fluxthresh,psfmode='hogbom',interactive=True,imsize=[512,512], cell=cells,stokes='I',weighting='natural',robust=0.5)

#viewer(imname+'.image')


