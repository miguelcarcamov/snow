######################################################################### 
#
#  plotbandpass3.py
#
#  This is meant to be a generic task to display CASA bandpass solutions
#  with options to overlay them in various combinations.  It is meant to
#  work the 'new' (casa 3.4) calibration tables, and (unlike plotbandpass
#  and plotbandpass2) to allow mixed spws (e.g. TDM and FDM for ALMA).  
#  This function uses the msmd tool when run in casa 4.1.0, but for older
#  versions, it still relies on
#  the ValueMapping class written by Stuartt Corder in analysisUtils.py.
#  It is designed to be called from inside the analysisUtils namespace:
#      import analysisUtils as au
#      au.plotbandpass()
#  or directly on its own:
#      import plotbandpass3 as p
#      p.plotbandpass()
#
#  (not to be confused with plotbpass.py, which was written for a specific
#   purpose of analyzing ALMA bandpass stability)
#
#  Todd R. Hunter  September 2012
#
# for regressions, see plotbandpass_regression.py
#
from __future__ import print_function  # prevents adding old-style print statements
PLOTBANDPASS_REVISION_STRING = "$Id: plotbandpass3.py,v 1.202 2020/03/23 20:18:30 thunter Exp $" 
import pylab as pb
import math, os, sys, re, inspect
import time as timeUtilities
import numpy as np
import analysisUtils as au
try:
    import casalith
    casaVersion = casalith.version_string()
except:
    import casadef     # necessary to read the casa version strings
    if casadef.casa_version >= '5.0.0':
        import casa as mycasa
        if 'cutool' in dir(mycasa):
            cu = mycasa.cutool()
            casaVersion = '.'.join([str(i) for i in cu.version()[:-1]]) + '-' + str(cu.version()[-1])
        else:
            casaVersion = mycasa.casa['build']['version'].split()[0]
    else:
        casaVersion = casadef.casa_version

try:
    # will work if casaVersion < '5.9.9' or there is a local file taskinit.py
    from taskinit import *  # needed for things like tb.open
#    print("plotbandpass3: imported casatools using taskinit *")
except:
    from casatasks import casalog
    from casatasks import gaincal
    from casatasks import tclean
    from casatools import measures as metool
    from casatools import table as tbtool
    from casatools import atmosphere as attool
    from casatools import msmetadata as msmdtool
    from casatools import image as iatool
    from casatools import ms as mstool
    from casatools import quanta as qatool
    from casaplotms import plotms
#    print("plotbandpass3: imported casatasks and casatools individually")
    
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter
import matplotlib.transforms
import inspect
try:
    from six.moves import input as raw_input
except:
    pass  # raw_input works in older versions where six.moves is unavailable

TOP_MARGIN  = 0.25   # Used if showatm=T or showtksy=T
BOTTOM_MARGIN = 0.25 # Used if showfdm=T
MAX_ATM_CALC_CHANNELS = 512

# This is a color sequence found online which has distinguishable colors
markeredgewidth=0.0
overlayColors = [
      [0.00,  0.00,  0.00], # black
      [0.00,  0.00,  1.00], # blue
      [0.00,  0.50,  0.00], # green
      [1.00,  0.00,  0.00], # red
      [0.00,  0.75,  0.75], # cyan
#      [0.75,  0.00,  0.75], # magenta, avoid because it's the same as atmcolor
      [0.25,  0.25,  0.25], # gray
      [0.75,  0.25,  0.25], # brick
      [0.95,  0.95,  0.00], # yellow
      [0.25,  0.25,  0.75], # light blue
#      [0.75,  0.75,  0.75],  this color is invisible for some reason
      [0.00,  1.00,  0.00], # bright green
      [0.76,  0.57,  0.17], #  olive
      [0.54,  0.63,  0.22],
      [0.34,  0.57,  0.92],
      [1.00,  0.10,  0.60],
#      [0.88,  0.75,  0.73],  invisible
      [0.10,  0.49,  0.47],
      [0.66,  0.34,  0.65],
      [0.99,  0.41,  0.23]]
overlayColors += overlayColors + overlayColors  # 17*3 = 51 total colors
overlayColors += overlayColors + overlayColors # try to support antenna,time
overlayColors += overlayColors + overlayColors # try to support antenna,time
overlayColors += overlayColors + overlayColors # try to support antenna,time
overlayColors += overlayColors + overlayColors # try to support antenna,time

# Enumeration to keep track of plot pages
PAGE_ANT = 0
PAGE_SPW = 1
PAGE_TIME = 2
PAGE_AP = 3

# Used to parse command line arguments
myValidCharacterList = ['~', ',', ' ', '*',] + [str(m) for m in range(10)]
myValidCharacterListWithBang = ['~', ',', ' ', '*', '!',] + [str(m) for m in range(10)]
LARGE_POSITIVE = +1e20
LARGE_NEGATIVE = -1e20
maxAntennaNamesAcrossTheTop = 17
maxTimesAcrossTheTop = 13 # 17 for HH:MM, reduced by 1 below for subplot=11
antennaVerticalSpacing = 0.018 # 0.016
antennaHorizontalSpacing = 0.05
xstartTitle = 0.07
ystartTitle = 0.955
xstartPolLabel = 0.05
ystartOverlayLegend = 0.931
opaqueSky = 270. # Kelvin, used for scaling TebbSky

developerEmail = "thunter@nrao.edu"

#class Polarization:
    # taken from Stokes.h in casa
#    (Undefined, I,Q,U,V,RR,RL,LR,LL,XX,XY,YX,YY) = range(13)  

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno
        
def println(linenumber, mystring):
    global terminal
    mystring = "v%s_%04d: " % (version(False).split()[2], linenumber) + mystring
    print(mystring)
    if (len(mystring) > 0):
        mystring += '\n'
    try:  # check if global is defined
        terminal.write(mystring)
        terminal.flush()
    except:
        return
    
def version(showfile=True):
    """
    Returns the CVS revision number.
    """
    myversion = "$Id: plotbandpass3.py,v 1.202 2020/03/23 20:18:30 thunter Exp $" 
    if (showfile):
        print("Loaded from %s" % (__file__))
    return myversion

def refTypeToString(rtype):
    rtypes = ['REST','LSRK','LSRD','BARY','GEO','TOPO','GALACTO','LGROUP','CMB']
    return(rtypes[rtype])

def corrTypeToString(ptype):
    ptypes = ['Undefined','I','Q','U','V','RR','RL','LR','LL','XX','XY','YX','YY']
    mystring = ptypes[ptype]
#    print "mystring = %s" % (mystring)
    return(mystring)
    
def buildAntString(antID,msFound,msAnt):
      if (msFound):
          antstring = msAnt[antID]
      else:
          antstring = '%02d' % (antID)
      return(antstring)
      
def makeplot(figfile,msFound,msAnt,overlayAntennas,pages,pagectr,density,
             interactive,antennasToPlot,spwsToPlot,overlayTimes,overlayBasebands,
             locationCalledFrom, xant, ispw, subplot, resample='1', debug=False, 
             figfileSequential=False, figfileNumber=0):
  if (type(figfile) == str):
      if (figfile.find('/')>=0):
          directories = figfile.split('/')
          directory = ''
          for d in range(len(directories)-1):
              directory += directories[d] + '/'
          if (os.path.exists(directory)==False):
              print("Making directory = ", directory)
              os.system("mkdir -p %s" % directory)
  if (debug):
      print("makeplot(%d): pagectr=%d, len(pages)=%d, len(spwsToPlot)=%d, len(figfile)=%d, figfileNumber=%d, xant=%d, msAnt=%s, antennasToPlot=%s, pages(ANT,SPW,TIME,AP)=" % (locationCalledFrom,
                                                            pagectr, len(pages),len(spwsToPlot), len(figfile), figfileNumber, xant, str(msAnt), str(antennasToPlot)), pages)
  if (pages[pagectr][PAGE_SPW] >= len(spwsToPlot)):
      # necessary for test86: overlay='spw' of spectral scan dataset.  to avoid indexing beyond the
      # end of the array in the the case that the final frame is of a baseband with n spw, and
      # earlier frames had >n spws   2014-04-08
      ispw = spwsToPlot[-1]
      if debug:
          print("setting ispw to final (because %d >= %d)" % (pages[pagectr][PAGE_SPW],len(spwsToPlot)))
  else:
      # CAS-8285: Added 'if' to be sure to use ispw passed in for single-panel plots, but 
      # use the original behavior for multi-panel plots simply to preserve the pngfile
      # naming scheme (i.e. including the spw name of lower right panel) to match old 
      # regressions.  Should probably remove this whole 'else' block someday, if I don't 
      # mind if future multi-panel filenames contain spw name of upper left panel.
      if (subplot != 11 or overlayBasebands):  # Add only this line for CAS-8285.
          ispw = spwsToPlot[pages[pagectr][PAGE_SPW]]
          if debug:
              print("setting ispw to spwsToPlot[pages[pagectr=%d]=%d[PAGE_SPW]] = %d" % (pagectr,pages[pagectr][PAGE_SPW],ispw))
  t = pages[pagectr][PAGE_TIME] #  + 1
  if (subplot == 11):
      antstring = buildAntString(xant, msFound, msAnt) # fix for CAS-8285
  else:
      # this causes png file to be named for the antenna in the upper left corner, rather than lower right corner
      antstring = buildAntString(antennasToPlot[pages[pagectr][PAGE_ANT]], msFound, msAnt)  # original behavior
  figfile = figfile.split('.png')[0]  # this will be added back later
  if (figfileSequential):
      plotfilename = figfile + '.%03d' % (figfileNumber)
  else:
      if (msFound):
          if (overlayAntennas and overlayTimes):
              plotfilename = figfile+'.spw%02d'%(ispw)
          elif (overlayAntennas):
              plotfilename = figfile+'.spw%02d'%(ispw)+'.t%02d'%(t)
          elif (overlayTimes):
              plotfilename = figfile+'.'+antstring+'.spw%02d'%(ispw)
          else:
              plotfilename = figfile+'.'+antstring+'.spw%02d'%(ispw)+'.t%02d'%(t)
      else:
          if (overlayAntennas and overlayTimes):
              plotfilename = figfile+'.spw%02d'%(ispw)
          elif (overlayAntennas):
              plotfilename = figfile+'.spw%02d'%(ispw)+'.t%02d'%(t)
          elif (overlayTimes):
              plotfilename = figfile+'.ant'+antstring+'.spw%02d'%(ispw)
          else:
              plotfilename = figfile+'.ant'+antstring+'.spw%02d'%(ispw)+'.t%02d'%(t)

  if (int(resample) > 1):
      plotfilename += '.resample%d.png' % (int(resample))
  else:
      plotfilename += '.png'
  if (interactive == False or 1==1):
      print("Building %s" % (plotfilename))
  pb.savefig(plotfilename, format='png', dpi=density)
  return(plotfilename)

def utdatestring(mjdsec):
    (mjd, dateTimeString) = au.mjdSecondsToMJDandUT(mjdsec)
    tokens = dateTimeString.split()
    return(tokens[0])

def mjdsecArrayToUTString(timerangeListTimes):
    """
    accepts [4866334935, 4866335281] etc.
    returns '08:04:10, 09:03:00' etc.
    """
    timerangeListTimesString = ''  
    for t in timerangeListTimes:
        timerangeListTimesString += utstring(t,3) + ' '
    return(timerangeListTimesString)

def utstring(mjdsec, xframeStart=110):
    (mjd, dateTimeString) = au.mjdSecondsToMJDandUT(mjdsec)
    tokens = dateTimeString.split()
    hoursMinutes = tokens[1][0:len(tokens[1])-3]
    hoursMinutesSeconds = tokens[1][0:len(tokens[1])]
    if (xframeStart == 110):  # 2011-01-01 UT 00:00
        return(tokens[0]+' '+tokens[2]+' '+hoursMinutes)
    elif (xframeStart == 3):
        return(hoursMinutesSeconds)
    else:  # 00:00
        return(hoursMinutes)
    
def openBpolyFile(caltable, debug=False):
    mytb = au.createCasaTool(tbtool)
    mytb.open(caltable)
    desc = mytb.getdesc()
    if ('POLY_MODE' in desc):
        polyMode = mytb.getcol('POLY_MODE')
        print("This is a BPOLY solution = %s" % (polyMode[0]))
        polyType = mytb.getcol('POLY_TYPE')
        scaleFactor = mytb.getcol('SCALE_FACTOR')
        antenna1 = mytb.getcol('ANTENNA1')
        times = mytb.getcol('TIME')
        cal_desc_id = mytb.getcol('CAL_DESC_ID')
        nRows = len(polyType)
        for pType in polyType:
            if (pType != 'CHEBYSHEV'):
                print("I do not recognized polynomial type = %s" % (pType))
                return
        # Here we assume that all spws have been solved with the same mode
        uniqueTimesBP = np.unique(mytb.getcol('TIME'))
        nUniqueTimesBP = len(uniqueTimesBP)
        print("There are %d unique times in the BPOLY solution:" % (nUniqueTimesBP))
        if (nUniqueTimesBP == 2):
            print("differing by %g seconds" % (uniqueTimesBP[1]-uniqueTimesBP[0]))
        mystring = ''
        for u in uniqueTimesBP:
            mystring += '%.3f, ' % (u)
        print(mystring)
        nPolyAmp = mytb.getcol('N_POLY_AMP')
        nPolyPhase = mytb.getcol('N_POLY_PHASE')
        frequencyLimits = mytb.getcol('VALID_DOMAIN')
        increments = 0.001*(frequencyLimits[1,:]-frequencyLimits[0,:])
        frequenciesGHz = []
        for i in range(len(increments)):
           freqs = (1e-9)*np.arange(frequencyLimits[0,i],frequencyLimits[1,i],increments[i])
           frequenciesGHz.append(freqs)
        polynomialAmplitude = []
        polynomialPhase = []
        for i in range(len(polyMode)):
            polynomialAmplitude.append([1])
            polynomialPhase.append([0])
            if (polyMode[i] == 'A&P' or polyMode[i] == 'A'):
                polynomialAmplitude[i]  = mytb.getcell('POLY_COEFF_AMP',i)[0][0][0]
            if (polyMode[i] == 'A&P' or polyMode[i] == 'P'):
                polynomialPhase[i] = mytb.getcell('POLY_COEFF_PHASE',i)[0][0][0]
  
        mytb.close()
        mytb.open(caltable+'/CAL_DESC')
        nSpws = len(mytb.getcol('NUM_SPW'))
        spws = mytb.getcol('SPECTRAL_WINDOW_ID')
        spwBP = []
        for c in cal_desc_id:
            spwBP.append(spws[0][c])
        mytb.close()
        nPolarizations = len(polynomialAmplitude[0]) / nPolyAmp[0]
        if (debug):
            print("(3)Set nPolarizations = %s" % (str(nPolarizations)))
        
        # This value is overridden by the new function doPolarizations in ValueMapping.
        # print "Inferring %d polarizations from size of polynomial array" % (nPolarizations)
        return([polyMode, polyType, nPolyAmp, nPolyPhase, scaleFactor, nRows, nSpws, nUniqueTimesBP,
                uniqueTimesBP, nPolarizations, frequencyLimits, increments, frequenciesGHz,
                polynomialPhase, polynomialAmplitude, times, antenna1, cal_desc_id, spwBP])
    else:
        mytb.close()
        return([])
   # end of openBpolyFile()

def displayTimesArray(uniqueTimesPerFieldPerSpw):
  """
  Displays an array of MJD second timestamps as UT timestamps
  """
  legendString = ''
  for s in uniqueTimesPerFieldPerSpw:
      legendString += "["
      for f in s:
          legendString += "["
          for t in f:
              legendString += "%s" % utstring(t,3)
              if (t != f[-1]):
                  legendString += ", "
          legendString += "]"
          if (f != s[-1]):
              legendString += ', '
      legendString += "], "
      if (s != uniqueTimesPerFieldPerSpw[-1]):
          legendString += ', '
  print(legendString)          

def checkPolsToPlot(polsToPlot, corr_type_string):
  firstFailure = 0
  for pol in polsToPlot:
      if ((pol in corr_type_string) == False):
          print("Polarization product %s is not in the ms" % (pol))
          firstFailure += 1
          if (pol in ['XX','YY']):
              polsToPlot = ['RR','LL']
          else:
              polsToPlot = ['XX','YY']
          break
  if (firstFailure>0):
     print("Looking for instead: ", polsToPlot)
     for pol in polsToPlot:
        if ((pol in corr_type_string) == False):
            print("Polarization product %s is not in the ms" % (pol))
            firstFailure += 1
            if (pol in ['XX']):
                polsToPlot = ['YY']
            elif (pol in ['YY']):
                polsToPlot = ['XX']
            elif (pol in ['RR']):
                polsToPlot = ['LL']
            elif (pol in ['LL']):
                polsToPlot = ['RR']
            break
  if (firstFailure > 1):
      print("Looking for instead: ", polsToPlot)
      for pol in polsToPlot:
          if ((pol in corr_type_string) == False):
              print("Polarization product %s is not in the ms" % (pol))
              return([])
  return(polsToPlot)

def getCorrType(msName, spwsToPlot, debug):
    """
    Open the DATA_DESCRIPTION_ID table.  Find the polarization_id of the first
    spw in the list of spwsToPlot, then read the CORR_TYPE from POLARIZATION
    table.
    """
    mytb = au.createCasaTool(tbtool)
    mytb.open(msName+'/DATA_DESCRIPTION')
    spws = mytb.getcol('SPECTRAL_WINDOW_ID')
    polarization_id = mytb.getcol('POLARIZATION_ID')
    mytb.close()
    pol_id = 0
    telescopeName = au.getObservatoryName(msName)
    mytb.open(msName+'/POLARIZATION')
    for myspw in spwsToPlot:
        if (debug):
            print("looking for %d in %s" % (myspw, str(spws)))
        row = list(spws).index(myspw)
        if (row >= 0):
            pol_id = polarization_id[row]
            corr_type = mytb.getcell('CORR_TYPE',pol_id)
            if (corr_type[0] >= 5 or (telescopeName.find('ALMA')<0 and telescopeName.find('VLA')<0)):
                # Undefined, I, Q, U, V, which ALMA and VLA never use
                # Need to allow non-VLA, non-ALMA to stop here
                break
#    num_corr = mytb.getcol('NUM_CORR')
    mytb.close()
    corr_type_string = []
    if (len(corr_type) == 4):
        print("This is a 4-polarization dataset.")
        if (corr_type[0] in [5,6,7,8]):
            corr_type = [5,8]
        elif (corr_type[0] in [9,10,11,12]):
            corr_type = [9,12]
        else:
            print("Unsupported polarization types = ", corr_type)
            return(corr_type, corr_type_string)
    # This overrides the len(gain_table) because it can have length=2 even when only 1 pol present
    nPolarizations = len(corr_type)
    if (debug):
        print("getCorrType():  (2)Set nPolarizations = %d" % nPolarizations)
    for ct in corr_type:
        corr_type_string.append(corrTypeToString(ct))
    print("corr_types = ", corr_type,  " = ", corr_type_string)
    return(corr_type, corr_type_string, nPolarizations)

def writeArgument(f,name,arg):
    if (type(arg) == str):
        s = "%-18s = '%s'" % (name,arg)
        t = "%s='%s'" % (name,arg)
    else:
        s = "%-18s = %s" % (name,str(arg))
        t = "%s=%s" % (name,arg)
    f.write(s+'\n')
    return(t)

def resampleSolution(x,y,resample=1):
    """
    Takes a solution (quantity y vs. channel x) and expands the number
    of points via linear interpolation
    """
    newx = []
    for i in range(len(x)):
        newx.append(x[i])
        if (i < len(x)-1):
            for j in range(1,resample):
                newx.append((x[i]*(resample-j) + x[i+1]*j)/(1.0*resample))
    newx = np.array(newx)
    if (False):
        f = scipy.interpolate.interp1d(x,y,kind='linear')
        newy = f(newx)
    else:
        newy = np.interp(newx, x, y)
    return(newx,newy)
        
def channelDifferences(y, x, resample=1):
    """
    Takes a vector, and computes the channel-to-channel derivative.
    Optionally, it will also resample the data and compute the
    derivative.
    - Todd Hunter
    """
    x = np.array(x)
    y = np.array(y)
    if (len(x) > 1):
        channelWidth = x[1]-x[0]
        d = (np.diff(y)/np.diff(x))
        newy = d*channelWidth 
        newx = (x[1:]+x[:-1])/2.  # midpoints of input x-axis
    else:
        newx = x
        newy = y
    if (resample > 1):
        x,y = resampleSolution(x,y,resample)
        if (len(x) > 1):
            channelWidth = x[1]-x[0]
            d = (np.diff(y)/np.diff(x))
            resy = d*channelWidth 
            resx = (x[1:]+x[:-1])/2.  # midpoints of input x-axis
        else:
            resx = x
            resy = y
    else:
        resy = newy
        resx = newx
    return(newy, newx, resy, resx)

def drawOverlayTimeLegends(xframe,firstFrame,xstartTitle,ystartTitle,caltable,titlesize,
                           fieldIndicesToPlot,ispwInCalTable,uniqueTimesPerFieldPerSpw,
                           timerangeListTimes, solutionTimeThresholdSeconds,debugSloppyMatch,
                           ystartOverlayLegend,debug,mysize, fieldsToPlot,myUniqueColor,
                           timeHorizontalSpacing, fieldIndex,overlayColors,
                           antennaVerticalSpacing, overlayAntennas, timerangeList, caltableTitle,
                           mytime, scansToPlot, scansForUniqueTimes):
    """
    Draws the legend at the top of the page, if it is the correct time to do so,
    including the overlayTimes, the 'UT' label, and the caltable name.
    """
    if (xframe == firstFrame):
        # draw title including caltable name
        pb.text(xstartTitle, ystartTitle, caltableTitle, size=titlesize,
                color='k', transform=pb.gcf().transFigure)
        # support multi-fields with overlay='time'
        uTPFPS = []  # stands for uniqueTimesPerFieldPerSpw
        uTPFPStimerange = []
        # Find all timerange integers for all fields, not just the ones that were plotted
        allTimeranges = []
        for f in range(len(uniqueTimesPerFieldPerSpw[ispwInCalTable])):
            for t in uniqueTimesPerFieldPerSpw[ispwInCalTable][f]:
                if (t in timerangeListTimes):
                    allTimeranges.append(list(timerangeListTimes).index(t))
        for f in fieldIndicesToPlot:
#            print "uniqueTimesPerFieldPerSpw[spw=%d][field=%d] = " % (ispwInCalTable,f), uniqueTimesPerFieldPerSpw[ispwInCalTable][f]
            for t in uniqueTimesPerFieldPerSpw[ispwInCalTable][f]:
                matched, mymatch = sloppyMatch(t, timerangeListTimes, solutionTimeThresholdSeconds,
                                               myprint=debugSloppyMatch, whichone=True)
                if (matched):
                    uTPFPS.append(t)
                    uTPFPStimerange.append(mymatch)

        allTimeranges = list(np.sort(np.unique(allTimeranges)))
        idx = np.argsort(uTPFPS)
        uTPFPStimerange = np.array(uTPFPStimerange)[idx]
        uTPFPS = np.sort(uTPFPS)
        if (debug):
            print("uTPFPS = ", np.array(uTPFPS, dtype='int'))
            print("uTPFPStimerange = ", np.array(uTPFPStimerange, dtype='int'))
            print("scansForUniqueTimes = ", scansForUniqueTimes)
            print("timerangeList = ", timerangeList)
            print("allTimeranges = ", allTimeranges)
            print("fieldsToPlot = ", fieldsToPlot)
            print("fieldIndex = ", fieldIndex)
            print("fieldsIndicesToPlot = ", fieldIndicesToPlot)
        timeFormat = 3  # HH:MM:SS
        maxTimesAcross = maxTimesAcrossTheTop
        if (firstFrame == 111):
            maxTimesAcross -= 2

        for a in range(len(uTPFPS)):
            legendString = utstring(uTPFPS[a],timeFormat)
            if (debug): print("----> %d: Defined legendString: %s" % (a,legendString))
            if (a==0):
                pb.text(xstartTitle-0.03, ystartOverlayLegend, 'UT',color='k',fontsize=mysize,
                        transform=pb.gcf().transFigure)
            if (a < maxTimesAcross):
                x0 = xstartTitle + (a*timeHorizontalSpacing)
                y0 = ystartOverlayLegend
            else:
                # start going down the righthand side
                x0 = xstartTitle + (maxTimesAcross*timeHorizontalSpacing)
                y0 = ystartOverlayLegend-(a-maxTimesAcross)*antennaVerticalSpacing
            if (True):
                if (debug):
                    print("3)checking time %d" % (a))
                if (sloppyMatch(uTPFPS[a],timerangeListTimes,solutionTimeThresholdSeconds,
                                mytime, scansToPlot, scansForUniqueTimes,
                                debugSloppyMatch)):
                    myUniqueTime = uTPFPS[a]
                    if (debug):
                        print("3)setting myUniqueTime to %d" % (myUniqueTime))
            if (debug): print("----> Drawing legendString: %s" % (legendString))
            if ((len(fieldsToPlot) > 1 or len(timerangeList) > 1) and overlayAntennas==False):
                # having overlayAntennas==False here will force all time labels to be black (as desired)
                if (debug):
                    print("allTimeranges.index(%d) = %d" % (a,allTimeranges.index(uTPFPStimerange[a])))
# old method
#                    print "len(uTPFPS)=%d, a=%d, len(myUniqueColor)=%d, overlayColors[%d]=%s" % (len(uTPFPS),a,len(myUniqueColor),timerangeList[uTPFPStimerange[a]],str(overlayColors[timerangeList[uTPFPStimerange[a]]]))
#                pb.text(x0, y0, legendString,color=overlayColors[timerangeList[uTPFPStimerange[a]]],fontsize=mysize,
#                        transform=pb.gcf().transFigure)
                    print("len(uTPFPS)=%d, a=%d, len(myUniqueColor)=%d, overlayColors[%d]=%s" % (len(uTPFPS),a,len(myUniqueColor),timerangeList[allTimeranges.index(uTPFPStimerange[a])],str(overlayColors[timerangeList[allTimeranges.index(uTPFPStimerange[a])]])))
# old method
#                pb.text(x0, y0, legendString,color=overlayColors[timerangeList[a]],fontsize=mysize,
#                        transform=pb.gcf().transFigure)
                pb.text(x0, y0, legendString,color=overlayColors[timerangeList[allTimeranges.index(uTPFPStimerange[a])]],fontsize=mysize,
                        transform=pb.gcf().transFigure)
            else:
                pb.text(x0, y0, legendString,fontsize=mysize, transform=pb.gcf().transFigure)

def lineNumber():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno
        
def drawAtmosphereAndFDM(showatm, showtsky, atmString, subplotRows, mysize, TebbSky,
                         TebbSkyImage,plotrange, xaxis, atmchan, atmfreq, transmission,
                         subplotCols, showatmPoints,xframe, channels,LO1,atmchanImage,
                         atmfreqImage,transmissionImage, firstFrame,showfdm,nChannels,
                         tableFormat,
                         originalSpw_casa33, chanFreqGHz_casa33,originalSpw,chanFreqGHz,
                         overlayTimes, overlayAntennas, xant, antennasToPlot, overlaySpws,
                         baseband, showBasebandNumber, basebandDict, overlayBasebands,
                         drewAtmosphere, showtsys=False, Trx=None):
    """
    If requested by the user at the command line, draw the atmospheric curve
    and the FDM window locations.
    """
    mylineno = lineNumber()
    ylim = pb.ylim()  # CAS-8655
    if ((showatm or showtsky) and len(atmString) > 0): 
        ylim = DrawAtmosphere(showatm, showtsky, subplotRows, atmString,
                       mysize, TebbSky, plotrange, xaxis, atmchan,
                       atmfreq, transmission, subplotCols,
                       showatmPoints=showatmPoints, xframe=xframe, channels=channels,
                       mylineno=mylineno,xant=xant, overlaySpws=overlaySpws,
                       overlayBasebands=overlayBasebands, drewAtmosphere=drewAtmosphere,
                       loc=201, showtsys=showtsys, Trx=Trx)
        if (LO1 is not None):
            # Now draw the image band
            ylim = DrawAtmosphere(showatm,showtsky, subplotRows, atmString,
                                  mysize, TebbSkyImage, plotrange, xaxis,
                                  atmchanImage, atmfreqImage, transmissionImage,
                                  subplotCols, LO1, xframe, firstFrame, showatmPoints, channels=channels,
                                  mylineno=mylineno,xant=xant, overlaySpws=overlaySpws,
                                  overlayBasebands=overlayBasebands,drewAtmosphere=drewAtmosphere,
                                  loc=202, showtsys=showtsys, Trx=Trx)
    if (xaxis.find('freq')>=0 and showfdm and nChannels <= 256):
        if (tableFormat == 33):
            showFDM(originalSpw_casa33, chanFreqGHz_casa33, baseband,
                    showBasebandNumber, basebandDict)
        else:
            showFDM(originalSpw, chanFreqGHz, baseband, showBasebandNumber,
                    basebandDict)
        ylim = pb.ylim()  # CAS-11062 need to pass the new wider limits back up to calling function
    return ylim  # CAS-8655

def DrawPolarizationLabelsForOverlayTime(xstartPolLabel,ystartPolLabel,corr_type,polsToPlot,
                                         channeldiff,ystartMadLabel,subplotRows,gamp_mad,
                                         mysize,
                                         ampmarkstyle,markersize,ampmarkstyle2, gamp_std):
    """
    Currently this is only called for amp vs. X plots. The corresponding code for phase
    vs. X plots is still inside plotbandpass().  But this is okay because overlay='time'
    is mainly intended for Tsys plots.
    """
#    print "DrawPolarizationLabelsForOverlayTime"
    x0 = xstartPolLabel
    y0 = ystartPolLabel
    if (corrTypeToString(corr_type[0]) in polsToPlot):
        if (channeldiff > 0):
            pb.text(x0, ystartMadLabel-0.03*subplotRows*0,
                    corrTypeToString(corr_type[0])+' MAD = %.4f, St.Dev = %.4f'%(gamp_mad[0]['mad'], gamp_std[0]['std']),
                    color='k',size=mysize, transform=pb.gca().transAxes)
        if (ampmarkstyle.find('-')>=0):
            pb.text(x0, y0, corrTypeToString(corr_type[0])+' solid', color='k',
                    size=mysize, transform=pb.gca().transAxes)
        else:
            pb.text(x0+0.02, y0, corrTypeToString(corr_type[0]), color='k',
                    size=mysize, transform=pb.gca().transAxes)
            pdesc = pb.plot([x0-0.1], [y0], '%sk'%ampmarkstyle, markersize=markersize,
                            scalex=False,scaley=False, transform=pb.gca().transAxes,markeredgewidth=markeredgewidth)
    if (len(corr_type) > 1):
        if (corrTypeToString(corr_type[1]) in polsToPlot):
            if (channeldiff > 0):
                pb.text(x0, ystartMadLabel-0.03*subplotRows*1,
                        corrTypeToString(corr_type[1])+' MAD = %.4f, St.Dev = %.4f'%(gamp_mad[1]['mad'], gamp_std[1]['std']),
                        color='k',size=mysize, transform=pb.gca().transAxes)
            if (ampmarkstyle2.find('--')>=0):
                pb.text(x0, y0-0.03*subplotRows, corrTypeToString(corr_type[1])+' dashed',
                        color='k', size=mysize, transform=pb.gca().transAxes)
            else:
                pb.text(x0, y0-0.03*subplotRows, corrTypeToString(corr_type[1]), # removed +0.02*xrange on 11-Mar-2014
                        color='k', size=mysize, transform=pb.gca().transAxes)
                pdesc = pb.plot([x0-0.1], [y0-0.03*subplotRows], '%sk'%ampmarkstyle2,
                                markersize=markersize, scalex=False,scaley=False, 
                                transform=pb.gca().transAxes,markeredgewidth=markeredgewidth)

def countDigitsInXTickLabels(adesc):
    """
    This was supposed to detect when x-axis labels overlap, and reduce the
    number of major ticks accordingly.
    I have not gotten this to work yet, it only reads back '' as each label.
    """
    xticklabels = adesc.get_xticklabels()  # both of these give blank labels?
#    xticklabels = pb.gca().xaxis.get_majorticklabels() # both of these give blank labels?
    print("xticklabels = %s" % (str(xticklabels)))
    for l in range(len(xticklabels)):
        print("l = ", xticklabels[l])
    digits = 0
    for l in range(len(xticklabels)):
        digits += len('%g' % (float(xticklabels[l].get_text())))
    print("%d digits in %d labels" % (digits,len(xticklabels)))

def GetFieldIdsForFieldName(token, vm, mymsmd, msFields):
    if (casaVersion >= '4.1.0'):
        if (mymsmd != '' and mymsmd is not None):
            myid = mymsmd.fieldsforname(token)[0]
        else:
            myid = list(msFields).index(token)
    else:
        myid = vm.getFieldIdsForFieldName(token)
    return(myid)

def GetFieldNamesForFieldId(u, vm, mymsmd, msFields):
    if (casaVersion >= '4.1.0'):
        if (mymsmd != '' and mymsmd is not None):
            myFieldName = mymsmd.namesforfields(u)[0]
        else:
            myFieldName = msFields[u]
    else:
        myFieldName = vm.getFieldNamesForFieldId(u)
    return(myFieldName)

def computeHighestSpwIndexInSpwsToPlotThatHasCurrentScan(spwsToPlot, scansToPlotPerSpw, scan):
    highestSpwIndex = -1
    for i,spw in enumerate(spwsToPlot):
        if (scan in scansToPlotPerSpw[spw]):
            highestSpwIndex = i
    return(highestSpwIndex)

def madOfDiff(solution):
    """
    This function is used to decide which of two curves has more scatter, and hence should
    be plotted first (i.e. shown in the background) when overlaying two solutions.
    Added as part of CAS-9474 to do a better job of the selection
    """
    if (len(solution) < 4):
        return au.MAD(np.diff(solution))
    else:
        start = len(solution)/4
        stop = len(solution)*3/4
        ydata = np.array(solution[start:stop+1])
        return au.MAD(np.diff(ydata))
        
DEFAULT_PLATFORMING_THRESHOLD = 10.0 # unused if platformingSigma != 0
def plotbandpass3(caltable='', antenna='', field='', spw='', yaxis='amp',
                  xaxis='chan', figfile='', plotrange=[0,0,0,0], help=False,
                  caltable2='', overlay='', showflagged=False, timeranges='',
                  buildpdf=False, caltable3='', markersize=3, density=108,
                  interactive=True, showpoints='auto', showlines='auto',   # 20 args
                  subplot='22', zoom='', poln='', showatm=False, pwv='auto',
                  gs='gs', convert='convert', chanrange='',
                  solutionTimeThresholdSeconds=30.0, debug=False, vm='',
                  phase='',  ms='', showtsky=False, showfdm=False,showatmfield='',
                  lo1=None, showimage=False, showatmPoints=False, parentms='', # 40 args
                  pdftk='pdftk', channeldiff=False, edge=8, resample=1, vis='',
                  platformingThreshold=DEFAULT_PLATFORMING_THRESHOLD,
                  platformingSigma=5.0, basebands=None, showBasebandNumber=False,
                  scans='', figfileSequential=False, groupByBaseband=False,
                  cleanup=False, caltable2amplitudeOffset=0, xcolor='b', 
                  ycolor='g', chanrangeSetXrange=False, 
                  overlaySpwDistinguish='', asciiFile=False, 
                  maxAtmCalcChannels=MAX_ATM_CALC_CHANNELS,maxAltitude=60,
                  firstPlot=0, Feff=0.99, SBGain=0.99, Trx='auto', showtsys=False): # 66 args
  """
  This is a tool to plot bandpass solutions faster than plotcal.  It is 
  designed to work on both the old style and new style cal tables.  The 
  source code is in plotbandpass3.py.  For more
  detailed help, run au.plotbandpass3(help=True) or see examples at:
  http://casaguides.nrao.edu/index.php?title=Plotbandpass
  -- Todd Hunter
  """
  print("%s" % (PLOTBANDPASS_REVISION_STRING))
  DEBUG = debug
  if (help):
      print("Usage: plotbandpass(caltable='', antenna='', field='', spw='', yaxis='amp',")
      print("   xaxis='chan', figfile='', plotrange=[0,0,0,0], help=False, caltable2='',") 
      print("   overlay='', showflagged=False, timeranges='', buildpdf=False, caltable3='',")
      print("   markersize=3, density=108, interactive=True, showpoints='auto',")
      print("   showlines='auto', subplot='22', zoom='', poln='', showatm=False, pwv='auto',")
      print("   gs='gs', convert='convert', chanrange='', debug=False, vm='',")
      print("   solutionTimeThresholdSeconds=30.0, phase='', ms='', showtsky=False,")
      print("   showfdm=False, showatmfield='', lo1=None, showimage=False,")
      print("   showatmPoints=False, parentms='', pdftk='pdftk', channeldiff=False,")
      print("   edge=8, resample=1, vis='',platformingThreshold=%f," % (DEFAULT_PLATFORMING_THRESHOLD))
      print("   platformingSigma=%.1f, basebands=None, showBasebandNumber=False," % (5.0))
      print("   scans='', figfileSequential=False, groupByBaseband=False,")
      print("   cleanup=False, caltable2amplitudeOffset=0, xcolor='b',")
      print("   ycolor='g', chanrangeSetXrange=False)")
      print(" antenna: must be ID (int or string or list), or a single antenna name or list")
      print(" asciiFile: if True, then also dump an ascii file of the amplitude spectrum")
      print(" basebands: show only spws from the specified baseband or list of basebands (default:None=all)")
      print(" buildpdf: True/False, if True and figfile is set, assemble pngs into a pdf")
      print(" caltable: a bandpass table, of type B or BPOLY")
      print(" caltable2: a second cal table, of type BPOLY or B, to overlay on a B table")
      print(" caltable2amplitudeOffset: constant value to add to caltable2")
      print(" caltable3: a third cal table, of type BPOLY, to overlay on the first two")
      print(" chanrange: set xrange ('5~100' or '80%') over which to autoscale y-axis for xaxis='freq'")
      print(" chanrangeSetXrange: if True, then chanrange also sets the xrange to display")
      print(" channeldiff: set to value > 0 (sigma) to plot derivatives of amplitude")
      print(" cleanup: remove pngs after making pdf when buildpdf=True")
      print(" convert: full path for convert command (in case it's not found)")
      print(" density: dpi to use in creating PNGs and PDFs (default=108)")
      print(" edge: the number of edge channels to ignore in finding outliers (for channeldiff>0)")
      print(" field: must be an ID, source name, or list thereof; can use trailing *: 'J*'")
      print(" figfile: the base_name of the png files to save: base_name.antX.spwY.png")
      print(" figfileSequential: naming scheme, False: name by spw/antenna (default)")
      print("                    True: figfile.1.png, figfile.2.png, etc.")
      print(" firstPlot: relevant for overlaying caltables; 1 -> plot first solution 1st, etc.")
      print(" groupByBaseband: group spws for display by baseband") 
      print(" gs: full path for ghostscript command (in case it's not found)")
      print(" help: print this message")
      print(" interactive: if False, then figfile will run to completion automatically")
      print(" lo1: specify the LO1 setting (in GHz) for the observation (used if showimage=T)")
      print(" overlay: 'antenna','time','antenna,time','spw', or 'baseband'")
      print("         makes 1 plot with different items in different colors")
      print(" overlaySpwDistinguish: '','color','width','color,width'; use 'width2' for wider")
      print(" showflagged:  show the values of data, even if flagged")
      print(" markersize: size of points (default=3)")
      print(" ms: name of the ms for this table, in case it does not match the string in the caltable")
      print(" parentms: name of the parent ms, in case the ms has been previously split")
      print(" pdftk: full path for pdftk command (in case it's not found)")
      print(" phase: the y-axis limits to use for phase plots when yaxis='both'")
      print(" platformingSigma: declare platforming if the amplitude derivative exceeds this many times the MAD")
      print(" platformingThreshold: if platformingSigma=0, then declare platforming if the amplitude")
      print("                       derivative exceeds this percentage of the median")
      print(" plotrange: define axis limits: [x0,x1,y0,y1] where 0,0 means auto")
      print(" poln: polarizations to plot (e.g. 'XX','YY','RR','LL' or '' for both)")
      print(" pwv: define the pwv to use for the showatm option: 'auto' or value in mm")
      print(" resample: channel expansion factor to use when computing MAD of derivative (for channeldiff>0)")
      print(" scans: show only solutions for the specified scans (int, list, or string)")
      print(" showatm: compute and overlay the atmospheric transmission curve")
      print(" showatmfield: for overlay='time', use first observation of this fieldID or name")
      print(" showatmPoints: draw atmospheric curve with points instead of a line")
      print(" showBasebandNumber: put the BBC_NO in the title of each plot")
      print(" showfdm: when showing TDM spws with xaxis='freq', draw locations of FDM spws.")
      print("     If showBasebandNumber=True, then show all FDM spws regardless of baseband")
      print(" showimage: also show the atmospheric curve for the image sideband (in black)")
      print(" showtsky: compute and overlay the sky temperature curve instead of transmission")
      print(" showlines: draw lines connecting the data (default=T for amp, F for phase)")
      print(" showpoints: draw points for the data (default=F for amp, T for phase)")
      print(" solutionTimeThresholdSeconds: consider 2 solutions simultaneous if within this interval (default=30)")
      print(" spw: must be single ID or list or range (e.g. 0~4, not the original ID)")
      print(" subplot: 11..81,22,32 or 42 for RowsxColumns (default=22), any 3rd digit is ignored")
      print(" timeranges: show only these timeranges, the first timerange being 0")
      print(" vm: the result from ValueMapping('my.ms'), or as returned from a previous call to plotbandpass")
      print(" xaxis: 'chan' or 'freq'")
      print(" xcolor: color for XX polarization points (default = blue)")
      print(" yaxis: 'amp', 'tsys', 'phase', or 'both' amp&phase == 'ap'; append 'db' for dB")
      print(" ycolor: color for YY polarization points (default = green)")
      print(" zoom: 'intersect' will zoom to overlap region of caltable with caltable2")
      return(vm)
  mytimestamp = timeUtilities.time()
  if type(poln) == list:
      # fix for CASA6 which forces string to be a list of strings of length one
      if len(poln) == 1:
          poln = poln[0]
  debugSloppyMatch = debug
  doneOverlayTime = False  # changed from True on 08-nov-2012
  missingCalWVRErrorPrinted = False
  if (ms == '' and vis != ''):
      # make au.plotbandpass compatible with casapy's argument list (which uses vis instead of ms)
      ms = vis  # this is the only use of 'vis' from the command line
  # initialize the arguments to DrawAtmosphereAndFDM() 
  TebbSky = None      
  TebbSkyImage = None 
  Tsys = None      
  TsysImage = None 
  atmchan = None      
  atmfreq = None      
  transmission = None 
  atmchanImage = None      
  atmfreqImage = None      
  transmissionImage = None 
  originalSpw_casa33 = None
  originalSpw = None       
  chanFreqGHz_casa33 = None
  chanFreqGHz = None
  # initialize arguments to DrawPolarizationLabelsForOverlayTime()
  gamp_mad = None
  gamp_std = None
  figfileNumber = 0  # only used if figfileSequential == True
  
  # Write a .last file
  cmd = 'plotbandpass'
  if (os.access(os.getcwd(),os.W_OK)):
    if (os.path.exists('plotbandpass.last') == False or os.access('plotbandpass.last',os.W_OK)):
      lastfile = open('%s.last'%cmd, 'w')
      lastfile.write('taskname           = "%s"\n'%cmd)
      cmd += '(' + writeArgument(lastfile, "caltable", caltable)
      cmd += ',' + writeArgument(lastfile, "antenna" , antenna)
      cmd += ',' + writeArgument(lastfile, "field" , field)
      cmd += ',' + writeArgument(lastfile, "spw" , spw)
      cmd += ',' + writeArgument(lastfile, "yaxis", yaxis)
      cmd += ',' + writeArgument(lastfile, "xaxis", xaxis)
      cmd += ',' + writeArgument(lastfile, "figfile", figfile)
      cmd += ',' + writeArgument(lastfile, "plotrange" , plotrange)
      cmd += ',' + writeArgument(lastfile, "help", help)
      cmd += ',' + writeArgument(lastfile, "caltable2", caltable2)
      cmd += ',' + writeArgument(lastfile, "overlay", overlay)
      cmd += ',' + writeArgument(lastfile, "showflagged", showflagged)
      cmd += ',' + writeArgument(lastfile, "timeranges", timeranges)
      cmd += ',' + writeArgument(lastfile, "buildpdf", buildpdf)
      cmd += ',' + writeArgument(lastfile, "caltable3", caltable3)
      cmd += ',' + writeArgument(lastfile, "markersize", markersize)
      cmd += ',' + writeArgument(lastfile, "density", density)
      cmd += ',' + writeArgument(lastfile, "interactive", interactive)
      cmd += ',' + writeArgument(lastfile, "showpoints", showpoints)
      cmd += ',' + writeArgument(lastfile, "showlines", showlines)
      cmd += ',' + writeArgument(lastfile, "subplot", subplot)
      cmd += ',' + writeArgument(lastfile, "zoom", zoom)
      cmd += ',' + writeArgument(lastfile, "poln", poln)
      cmd += ',' + writeArgument(lastfile, "showatm", showatm)
      cmd += ',' + writeArgument(lastfile, "showatmfield", showatmfield)
      cmd += ',' + writeArgument(lastfile, "pwv", pwv)
      cmd += ',' + writeArgument(lastfile, "gs", gs)
      cmd += ',' + writeArgument(lastfile, "convert", convert)
      cmd += ',' + writeArgument(lastfile, "chanrange", chanrange)
      cmd += ',' + writeArgument(lastfile, "chanrangeSetXrange", chanrangeSetXrange)
      cmd += ',' + writeArgument(lastfile, "solutionTimeThresholdSeconds", solutionTimeThresholdSeconds)
      cmd += ',' + writeArgument(lastfile, "debug", debug)
      cmd += ',' + writeArgument(lastfile, "vm", vm)
      cmd += ',' + writeArgument(lastfile, "phase", phase)
      cmd += ',' + writeArgument(lastfile, "ms", ms)
      cmd += ',' + writeArgument(lastfile, "parentms", parentms)
      cmd += ',' + writeArgument(lastfile, "lo1", lo1)
      cmd += ',' + writeArgument(lastfile, "showimage", showimage)
      cmd += ',' + writeArgument(lastfile, "showtsky", showtsky)
      cmd += ',' + writeArgument(lastfile, "showatmPoints", showatmPoints)
      cmd += ',' + writeArgument(lastfile, "showfdm", showfdm)
      cmd += ',' + writeArgument(lastfile, "pdftk", pdftk)
      cmd += ',' + writeArgument(lastfile, "channeldiff", channeldiff)
      cmd += ',' + writeArgument(lastfile, "edge", edge)
      cmd += ',' + writeArgument(lastfile, "resample", resample)
      cmd += ',' + writeArgument(lastfile, "vis", vis)
      cmd += ',' + writeArgument(lastfile, "platformingThreshold", platformingThreshold)
      cmd += ',' + writeArgument(lastfile, "platformingSigma", platformingSigma)
      cmd += ',' + writeArgument(lastfile, "basebands", basebands)
      cmd += ',' + writeArgument(lastfile, "showBasebandNumber", showBasebandNumber)
      cmd += ',' + writeArgument(lastfile, "scans", scans)
      cmd += ',' + writeArgument(lastfile, "figfileSequential", figfileSequential)
      cmd += ',' + writeArgument(lastfile, "groupByBaseband", groupByBaseband)
      cmd += ',' + writeArgument(lastfile, "cleanup", cleanup)
      cmd += ',' + writeArgument(lastfile, "caltable2amplitudeOffset", caltable2amplitudeOffset)
      cmd += ',' + writeArgument(lastfile, "xcolor", xcolor)
      cmd += ',' + writeArgument(lastfile, "ycolor", ycolor) + ')'
      lastfile.write('#%s\n'%(cmd))
      lastfile.close()

#  if (platformingThreshold != DEFAULT_PLATFORMING_THRESHOLD and channeldiff==False):
#      channeldiff = 10
  LO1 = None  # Fix for SCOPS-4877
  lo1s = None # Fix for SCOPS-4877
  if (showimage == False):
      LO1 = lo1 = None
  elif (lo1 is not None):
      if (lo1 > 1e6):
          # convert from Hz to GHz
          lo1 *= 1e-9
  if (showatm and showtsky):
      print("You have selected both showatm and showtsky!  Defaulting to showatm=True only.")
      showtsky = False
  if (showatm==False and showtsky==False and showatmfield!=''):
      print("Defaulting to showatm=True because showatmfield was specified.")
      showatm = True
  if (showatm==False and showtsky==False and showimage==True):
      print("Defaulting to showatm=True because showimage was True.")
      showatm = True
  if showtsys:
      showtsky = True

  if (overlay.find('time') < 0 and showatmfield != ''):
      print("The showatmfield only has meaning for overlay='time'.")
      return(vm)
  
  if (plotrange=='' or plotrange==[]):
      plotrange = [0,0,0,0]
  if (type(plotrange) != list):
      print("plotrange must be an array:  [0,1,-180,180]")
      return(vm)
  if (len(plotrange) < 4):
      print("plotrange must be an array:  [0,1,-180,180]")
      return(vm)
  if (phase != ''):
      if (type(phase) != list):
          print("phase must be either '' or 2 values: [x,y]")
          return(vm)
      if (len(phase) != 2):
          print("phase must be either '' or 2 values: [x,y]")
          return(vm)

  if (edge < 0):
      print("edge must be >= 0")
      return(vm)
  
  if (resample < 1):
      print("resample must be an integer >= 1")
      return(vm)
  resample = int(resample)
  solutionTimeThresholdSeconds = float(solutionTimeThresholdSeconds)  # so it also accepts a string
  
  if (buildpdf and figfile==''):
      print("With buildPDF=True, you must specify figfile='yourFileName' (.png will be appended if necessary).")
      return(vm)

  if (interactive==False and figfile=='' and channeldiff == False):
      print("With interactive=False and channeldiff=False, you must specify figfile='yourFileName' (.png will be appended if necessary).")
      return(vm)

  pxl = 0 # polarization number to use for setting xlimits if plotrange=[0,0...]
  chanrangePercent = None
  if (type(chanrange) != str):
      if (type(chanrange) != list):
          print("Chanrange must be a string or list:  '8~120' or [8,120]")
          return(vm)
      elif (len(chanrange) != 2):
          print("Chanrange must be a string or list:  '8~120' or [8,120]")
          return(vm)
      elif ((type(chanrange[0]) != int) or (type(chanrange[1]) != int)):
          print("Chanrange list members must be integers, not ", type(chanrange[0]), type(chanrange[1]))
          return
  elif (len(chanrange) < 1):
      chanrange = [0,0]
  else:
      if (chanrange.find('%')>0):
          chanrangePercent = float(chanrange.split('%')[0])
          if (debug): print("******************* Set chanrangePercent to %f, chanrangeSetXrange=" % (chanrangePercent), chanrangeSetXrange)
          if (chanrangePercent >= 100 or chanrangePercent <= 0):
              chanrangePercent = None
          chanrange = [0,0]
      elif (chanrange.find('~')>=0):
          tokens = chanrange.split('~')
          if (len(tokens) < 2):
              print("Invalid chanrange, too few tokens")
              return(vm)
          try:
              chanrange = [int(tokens[0]),int(tokens[1])]
              if (DEBUG):
                  print("Using chanrange = ", chanrange)
          except:
              print("Invalid chanrange, not integers")
              return(vm)
      else:
          print("Invalid chanrange, no tilde or percent sign found")
          return(vm)
      if (xaxis.find('chan')>=0):
              print("The chanrange parameter is only valid for xaxis='freq', and only if the plotrange is [0,0,0,0].")
              return(vm)
  if (chanrange[0] < 0):
      print("Invalid chanrange, cannot be negative")
      return(vm)
  if ((chanrange[0] != 0 or chanrange[1] != 0 or chanrangePercent is not None) and
      (plotrange[0] != 0 or plotrange[1] != 0 or plotrange[2] != 0 or plotrange[3] != 0)):
      print("If chanrange is specified, then plotrange must be all zeros.")
      return(vm)
      
  if (pwv==''):
      pwv = 1.0
  if (type(poln) != list):
        poln = poln.upper()
  if (poln == 'X'):
        poln = 'XX'
  if (poln == 'Y'):
        poln = 'YY'
  if (poln == 'X,Y' or poln=='Y,X'):
        poln = 'XX,YY'
  if (poln == 'R'):
        poln = 'RR'
  if (poln == 'L'):
        poln = 'LL'
  if (poln == 'R,L' or poln=='L,R'):
        poln = 'RR,LL'

  # Parse the polarizations to plot from the command line
  # Prior to opening the .ms (later), we cannot tell which products are actually present
  useAllPols = False
  if (poln == ''):
        useAllPols = True
        polsToPlot = ['XX','YY']  # assume ALMA initially
  elif (type(poln) == list):
        polsToPlot = poln.split(',')
  else:
        if ((poln in ['','RR','RL','LR','LL','XX','XY','YX','YY','RR,LL','XX,YY']) == False):
              print("Unrecognized polarization option = ", poln)
              return(vm)
        if (poln.find(',')>0):
            polsToPlot = poln.split(',')
        else:
            polsToPlot = [poln]
      
        
  if ((overlay in ['antenna', 'spw', 'time', 'baseband', '', 'antenna,time', 'time,antenna']) == False): # try to support antenna,time
     print("Unrecognized option for overlay: only 'antenna', 'spw', 'baseband', 'time' and 'antenna,time' are supported.")
     return(vm)
     
  allowedFrames = [11,21,31,41,51,61,71,81,22,32,42] # [11,22,32,42]
  if (int(subplot) > 100):
      # This will accept 111, 221, 321, 421, etc.
      subplot = int(subplot/10)
  if ((int(subplot) in allowedFrames)==False):
    print("Subplot choice (rows x columns) must be one of ", allowedFrames)
    print("(with an optional trailing digit that is ignored).")
    return(vm)

  if ((int(subplot) % 2) == 1):
      timeHorizontalSpacing = 0.06*1.3  # *1.3 is for HH:MM:SS
  else:
      timeHorizontalSpacing = 0.05*1.3  # *1.3 is for HH:MM:SS

  if (yaxis.find('both')<0 and yaxis.find('ap')<0 and yaxis.find('tsys')<0 and
      yaxis.find('amp')<0 and yaxis.find('phase')<0):
      print("Invalid yaxis.  Must be 'amp', 'tsys', 'phase' or 'both'.")
      return(vm)

  if (yaxis.find('tsys')>=0):
      yaxis = 'amp'

  if (xaxis.find('chan')<0 and xaxis.find('freq')<0):
      print("Invalid xaxis.  Must be 'chan' or 'freq'.")
      return(vm)

  if (showatm and showtsky):
      print("showatm=True and showtsky=True are mutually exclusive options")
      return(vm)

  if (showfdm and xaxis.find('freq')<0):
      print("The option showfdm=True requires xaxis='freq'.")
      return(vm)

  # Plotting settings
  minPhaseRange = 0.2
  plotfiles = []
  if (int(subplot) % 2 == 1):
    mysize = '10'
    titlesize = 10
  elif (int(subplot) == 22 or int(subplot) == 32):
    mysize = '8'
    titlesize = 8
  else:
    mysize = '7'
    titlesize = 8
  maxCharsBeforeReducingTitleFontSize = 72
  
  if (type(subplot) == str):
      subplot = int(subplot)
  if (subplot in allowedFrames == False):
      print("Invalid subplot = %d.  Valid options are: " % (subplot), allowedFrames)
      return(vm)
  xframeStart = int(subplot)*10  # i.e. 110 or 220 or 420
  firstFrame = xframeStart + 1
  lastFrame = xframeStart + int(subplot/10)*(subplot%10)
  bottomRowFrames = [111,212,313,414,515,616,717,818,223,224,325,326,427,428]  # try to make this more general
  leftColumnFrames = [111,211,212,311,312,313,411,412,413,414,511,512,513,514,515,611,612,613,614,615,616,
                      711,712,713,714,715,716,717,811,812,813,814,815,816,817,818,221,223,321,323,325,421,423,425,427]
  rightColumnFrames = [111,211,212,311,312,313,411,412,413,414,511,512,513,514,515,611,612,613,614,615,616,
                       711,712,713,714,715,716,717,811,812,813,814,815,816,817,818,222,224,322,324,326,422,424,426,428]
  subplotCols = subplot % 10
  subplotRows = int(subplot/10)
  ystartPolLabel = 1.0-0.04*subplotRows
  ystartMadLabel = 0.04*subplotRows
  if (subplotCols == 1):
      fstringLimit = 40 # character length of multi-field overlay title string
  elif (subplotCols == 2):
      fstringLimit = 12 # character length of multi-field overlay title string

  if (debug):
      print("0)setting xframe to %d" % (xframeStart))
  xframe = xframeStart
  previousSubplot = xframe
#  print "Using markersize = ", markersize
  pcolor = [xcolor,ycolor]
  x2color = 'k'
  y2color = 'c'
  p2color = ['k','c']
  x3color = 'm'
  y3color = 'r'
  p3color = ['m','r']
  if (showpoints == 'auto'):
      if (showlines == 'auto'):
          ampmarkstyle = '-'
          phasemarkstyle = '.'
          if (len(polsToPlot) == 1):
                ampmarkstyle2 = '-'
          else:
                ampmarkstyle2 = '--'
          phasemarkstyle2 = 'o'
      elif (showlines == False):
          ampmarkstyle = '.'
          ampmarkstyle2 = 'o'
          phasemarkstyle = '.'
          phasemarkstyle2 = 'o'
      else:
          ampmarkstyle = '-'
          phasemarkstyle = '-'
          if (len(polsToPlot) == 1):
                ampmarkstyle2 = '-'
                phasemarkstyle2 = '-'
          else:
                ampmarkstyle2 = '--'
                phasemarkstyle2 = '--'
  elif (showpoints == True):
      if (showlines == 'auto'):
          ampmarkstyle = '.-'
          phasemarkstyle = '.'
          if (len(polsToPlot) == 1):
                ampmarkstyle2 = 'o-'
          else:
                ampmarkstyle2 = 'o--'
          phasemarkstyle2 = 'o'
      elif (showlines == False):
          ampmarkstyle = '.'
          ampmarkstyle2 = 'o'
          phasemarkstyle = '.'
          phasemarkstyle2 = 'o'
      else:
          ampmarkstyle = '.-'
          phasemarkstyle = '.-'
          if (len(polsToPlot) == 1):
                ampmarkstyle2 = 'o-'
                phasemarkstyle2 = 'o-'
          else:
                ampmarkstyle2 = 'o--'
                phasemarkstyle2 = 'o--'
  else:  # showpoints == False
      if (showlines == False):
          print('You must have either showpoints or showlines set True or auto, assuming showlines=T')
      ampmarkstyle = '-'
      phasemarkstyle = '-'
      if (len(polsToPlot) == 1):
            ampmarkstyle2 = '-'
            phasemarkstyle2 = '-'
      else:
            ampmarkstyle2 = '--'
            phasemarkstyle2 = '--'

  ampmarkstyles = [ampmarkstyle,ampmarkstyle2]
  phasemarkstyles = [phasemarkstyle,phasemarkstyle2]
  # bpoly solutions should always be shown as lines, not dots or dots+lines
  bpolymarkstyle = '-'

  amplitudeWithPhase = (yaxis.find('both')>=0 or yaxis.find('ap')>=0)
  if (amplitudeWithPhase):
      myhspace = 0.30
      if (overlay.find('antenna')>=0 or overlay.find('time')>=0  or overlay.find('spw')>=0):
          print("Option overlay='antenna' or 'time' is incompatible with yaxis='both'.  Pick either amp or phase.")
          return(vm)
      if (subplotRows % 2 == 1):
          print("Option yaxis=='both' is incompatible with odd numbers of rows, change subplot")
          return(vm)
  else:
      myhspace = 0.30
  if (subplot/10 > 2):
      myhspace = 0.4
  if (subplot/10 > 3):
      myhspace = 0.6
  mywspace = 0.25
  
  # Now open the Bandpass solution table
  if (len(caltable) < 1):
      print("You need to specify a caltable.")
      return(vm)
  if (caltable[-1] == '/'):
      print("Stripping off the trailing '/' from the caltable name.")
      caltable = caltable[:-1]
  if not os.path.exists(caltable):
      print("Caltable does not exist.")
      return
  mytb = au.createCasaTool(tbtool)
  try:
      mytb.open(caltable)
  except:
      print("Could not open the caltable = %s" % (caltable))
      return(vm)
  if (caltable[0] != '/'):
      # print this so when someone sends me a bug report I can find their data!
      try:
          print("caltable = %s:%s/%s" % (os.uname()[1], os.getcwd(), caltable))
      except:
          print("caltable = localhost:%s/%s" % (os.getcwd(), caltable))
  else:
      try:
          print("caltable = %s:%s" % (os.uname()[1], caltable))
      except:
          print("caltable = localhost:%s" % (caltable))
  
  if (len(caltable) > 90):
      caltableTitle = '...' + caltable[-90:]
  else:
      caltableTitle = caltable
  names = mytb.colnames()
  if ('UVW' in names):
      print("This appears to be a measurement set, not a caltable. Aborting.")
      return(vm)
  casalog.post(cmd)
  ant = mytb.getcol('ANTENNA1')
  fields = mytb.getcol('FIELD_ID')
#  if (DEBUG):
#      print "FIELD_ID column = ", fields
  validFields = False
  for f in fields:
      if (f != -1):
          validFields = True
  if (validFields == False):
      print("The field_id is -1 (invalid) for all rows of this caltable.")
      print("Did you remember to run interpTsys.assignFieldAndScanToSolution()?")
      return(vm)
  try:
      flags = {}
      complete = True
      for f in range(len(fields)):
          if mytb.iscelldefined('FLAG',f):
              flags[f] = mytb.getcell('FLAG',f)
          else:
              complete = False
              print("Warning: Missing data in the FLAG column, table may not be complete.")
              break
  except:
      pass
  if not complete:
      print("Missing data in the FLAG column.") 
      print("If it is a solution file, does it contain solutions for both TDM and FDM spws?")
      if 0 not in flags:
          print("Are you sure this is a bandpass solution file, or is it the .ms?")
          return(vm)
#  if (debug): print "%d: flags.keys() = %s" % (len(flags.keys()), str(flags.keys()))

  times = mytb.getcol('TIME')
  intervals = mytb.getcol('INTERVAL')
  if ('SPECTRAL_WINDOW_ID' not in names):
      # This is an old-style CASA cal table.
      tableFormat = 33  
      msAnt = []
      cal_desc_id = mytb.getcol('CAL_DESC_ID')
      VisCal = (mytb.info())['subType']
      if (VisCal == "BPOLY"):
          print("This appears to be a BPOLY cal table written in the casa 3.3/3.4 style.")
      else:
          print("This appears to be an old-format cal table from casa 3.3 or earlier.")
      if (debug): print("VisCal = ", VisCal)
      mytb.close()
      ParType = "unknown"  # i.e. not Complex
      calDesc = mytb.open(caltable+'/CAL_DESC')
      originalSpws = mytb.getcol('SPECTRAL_WINDOW_ID')  # [[0,1,2,3]]
      if debug: print("originalSpws = ", originalSpws)
      originalSpw = originalSpws[0]                   # [0,1,2,3]
      if debug: print("originalSpw = ", originalSpw)
      msName = mytb.getcol('MS_NAME')[0]
      if debug: print("msName in table = ", msName)
      if (ms != ''):
          msName = ms
      # This appears to be the channel range extracted from the original spw, but is
      # only present in B solutions.  
      if (VisCal == "BPOLY"):
          originalChannelStart = np.zeros(len(originalSpw))
      else:
          originalChannelRange = mytb.getcol('CHAN_RANGE')
          originalChannelStart = originalChannelRange[0][0][:][0]
      mytb.close()
      try:
          mytb.open(msName+'/SPECTRAL_WINDOW')
          refFreq = mytb.getcol('REF_FREQUENCY')    
          net_sideband = mytb.getcol('NET_SIDEBAND')
          measFreqRef = mytb.getcol('MEAS_FREQ_REF')
          originalSpw_casa33 = range(len(measFreqRef))
          chanFreqGHz_casa33 = []     # used by showFDM
          for i in originalSpw_casa33:
              # The array shapes can vary.
              chanFreqGHz_casa33.append(1e-9 * mytb.getcell('CHAN_FREQ',i))
          mytb.close()
      except:
          print("2) Could not open the associated measurement set tables (%s). Will not translate antenna names." % (msName))
#          print "I will assume ALMA data: XX, YY, and refFreq=first channel."
#          corr_type_string = ['XX','YY']
#          corr_type = [9,12]
  else:  # 3.4
      tableFormat = 34
      cal_desc_id = mytb.getcol('SPECTRAL_WINDOW_ID')
      cal_scans = mytb.getcol('SCAN_NUMBER')
      unique_cal_scans = np.unique(cal_scans)
      cal_scans_per_spw = {}
      for myspw in np.unique(cal_desc_id):
          cal_scans_per_spw[myspw] = np.unique(cal_scans[np.where(myspw == cal_desc_id)[0]])
          if (debug):
              print("spw %d: scans %s" % (myspw,str(cal_scans_per_spw[myspw])))
      ParType = mytb.getkeyword('ParType')    # string = 'Complex'
      msName = mytb.getkeyword('MSName')      
      VisCal = mytb.getkeyword('VisCal')      # string = 'B TSYS'
      PolBasis = mytb.getkeyword('PolBasis')  # string = 'LINEAR'
      spectralWindowTable = mytb.getkeyword('SPECTRAL_WINDOW').split()[1]
      antennaTable = mytb.getkeyword('ANTENNA').split()[1]
      fieldTable = mytb.getkeyword('FIELD').split()[1]
      mytb.close()
      mytb.open(spectralWindowTable)
      chanFreqGHz = []
      originalSpws = range(len(mytb.getcol('MEAS_FREQ_REF')))
      originalSpw = originalSpws  # may need to do a global replace of this <------------------------
      originalSpwNames = mytb.getcol('NAME')
      for i in originalSpws:
          # The array shapes can vary.
          chanFreqGHz.append(1e-9 * mytb.getcell('CHAN_FREQ',i))
      mytb.close()
      #      CAS-6801 changes
      mytb.open(antennaTable)
      msAnt = mytb.getcol('NAME')
      mytb.close()
      mytb.open(fieldTable)
      msFields = mytb.getcol('NAME')
      mytb.close()
      
  if (VisCal == 'K Jones'):
      delay = True
      showpoints = True
      ampmarkstyle = '.'
      ampmarkstyle2 = 'o'
      if (markersize < 8): markersize = 8
  else:
      delay = False
  # Now open the associated ms tables via ValueMapping
  mymsmd = ''
#  msAnt = []  # comment this out when CAS-6801 changes are in place
  observatoryName = ''
  if (vm == ''):
    if (debug): print("msName = %s." % (msName))
    if (os.path.exists(msName) or os.path.exists(os.path.dirname(caltable)+'/'+msName)):
      if (os.path.exists(msName) == False):
          msName = os.path.dirname(caltable)+'/'+msName
          if (debug): print("found msName = %s." % (msName))
      if (casaVersion >= '4.1.0'):
          mymsmd = au.createCasaTool(msmdtool)
          if (debug): print("Running mymsmd on %s..." % (msName))
          mymsmd.open(msName)
          donetime = timeUtilities.time()
          if (debug): print("%.1f sec elapsed" % (donetime-mytimestamp))
          mytimestamp = timeUtilities.time()
          if (debug): print("time = %s" % (str(mytimestamp)))
          msAnt = mymsmd.antennanames(range(mymsmd.nantennas()))
          if (debug): print("msAnt = %s" % (str(msAnt)))
#          msFields = mymsmd.namesforfields(range(mymsmd.nfields())) # bombs if split has been run on subset of fields
          msFields = mymsmd.fieldnames()  # used to be namesforfields(), but that is broken in early CASA 6
          observatoryName = mymsmd.observatorynames()[0]
          print("Available antennas = %s" % (str(msAnt)))
      else:  # old CASA
        try:
          print("Running ValueMapping on %s..." % (msName))
          print("(This can take awhile, try the vm option next time: vm=au.plotbandpass(..)")
          print("                         then, on subsequent calls:    au.plotbandpass(..,vm=vm)")
          vm = au.ValueMapping(msName)
          donetime = timeUtilities.time()
          if (debug): print("%.1f sec elapsed" % (donetime-mytimestamp))
          mytimestamp = timeUtilities.time()
          msAnt = vm.antennaNamesForAntennaIds
          msFields = vm.fieldNamesForFieldIds
          print("Available antennas = ", msAnt)
          observatoryName = au.getObservatoryName(msName)
        except:
          print("1)Could not open the associated measurement set tables (%s). Will not translate antenna names or frequencies." % (msName))
    else: # ms does not exist
      if (ms=='' and tableFormat < 34):
          print("Could not find the associated measurement set (%s). Will not translate antenna names or frequencies." % (msName))
      elif (ms != ''):
          # Use the ms name passed in from the command line
          msName = ms
#          print "************* 2) Set msName to ", msName
          if (casaVersion >= '4.1.0'):
              mymsmd = au.createCasaTool(msmdtool)
              if (debug): print("Running mymsmd on %s..." % (msName))
              mymsmd.open(msName)
              donetime = timeUtilities.time()
              if (debug): print("%.1f sec elapsed" % (donetime-mytimestamp))
              mytimestamp = timeUtilities.time()
              if (debug): print("time = %s" % (str(mytimestamp)))
              msAnt = mymsmd.antennanames(range(mymsmd.nantennas()))
              if (debug): print("msAnt = %s" % (str(msAnt)))
          #          msFields = mymsmd.namesforfields(range(mymsmd.nfields())) # bombs if split has been run on subset of fields
              msFields = mymsmd.namesforfields()
              observatoryName = mymsmd.observatorynames()[0]
              print("Available antennas = %s" % (str(msAnt)))
          else:
            try:
              print("Running ValueMapping on %s..." % (msName))
              print("(This can take a minute, try using the vm option next time)")
              vm = au.ValueMapping(msName)
              donetime = timeUtilities.time()
              if (debug): print("%.1f sec elapsed" % (donetime-mytimestamp))
              mytimestamp = timeUtilities.time()
              msAnt = vm.antennaNamesForAntennaIds
              msFields = vm.fieldNamesForFieldIds
              print("Available antennas = ", msAnt)
              observatoryName = au.getObservatoryName(msName)
            except:
              print("1b) Could not open the associated measurement set tables (%s). Will not translate antenna names or channels to frequencies." % (msName))
  else:
      # vm was specified on the command line
      if (msName.find(vm.getInputMs()) < 0):
          if (msName.find(vm.getInputMs().split('/')[-1]) < 0):
              print("WARNING:  There is a mismatch between the ms name in the ValueMapping")
              print("structure provided and the ms name in the bandpass table:")
              print("   %s vs. %s" % (vm.getInputMs(), msName))
              print("Using the name passed in via the ms parameter.")
          else:
              msName = vm.getInputMs()
              print("WARNING:  There is a mismatch in the directory path for the ms between")
              print("the ValueMapping structure and the bandpass table.")
              print("Using the path from the ValueMapping result provided via the vm parameter.")
      msAnt = vm.antennaNamesForAntennaIds
      msFields = vm.fieldNamesForFieldIds
      if (casaVersion >= '4.1.0'):
          # This is necessary to get the baseband labelling correct.
          mymsmd = au.createCasaTool(msmdtool)
          if (debug): print("Running mymsmd on %s..." % (msName))
          mymsmd.open(msName)

  msFound =  False
  if (len(msAnt) > 0):
      msFields = list(msFields) # necessary to avoid having to index with extra 0: msFields[fieldIndex][0]
      msFieldsUnique = np.unique(msFields)
      msFound = True  # will also be true if necessary information is found in the caltable subtables
      print("Fields in ms  = ", msFields)
  else:
      msFields = []
      

  if (tableFormat == 33 and msFound):  # casa 3.3
      # Now open the associated ms tables via ValueMapping to figure out channel freqs
      chanFreqGHz = []
      for ictr in range(len(originalSpw)):
        if debug: print("ictr = ", ictr)
        if (casaVersion >= '4.1.0'):
            if debug: print("nspw = %d, np.max(originalSpw) = %d" % (mymsmd.nspw(False),np.max(originalSpw)))
            if (mymsmd.nspw(False) < np.max(originalSpw)): # waiting on CAS-4285
                # Then there was an extra split
                i = ictr
            else:
                i = originalSpw[ictr]
            nchan = mymsmd.nchan(i)
            if (nchan > 1):
                missingFrequencyWidth = originalChannelStart[ictr]*(mymsmd.chanfreqs(i)[-1]-mymsmd.chanfreqs(i)[0])/(nchan-1)
            else:
                missingFrequencyWidth = 0
            if (missingFrequencyWidth > 0):
                if (DEBUG):
                    print("Correcting for channels flagged prior to running bandpass by %f GHz" % (missingFrequencyWidth*1e-9))
            newfreqs = 1e-9*(mymsmd.chanfreqs(i)) + missingFrequencyWidth*1e-9
        else:
            if debug: print("len(vm.spwInfo) = %d, np.max(originalSpw) = %d" % (len(vm.spwInfo),np.max(originalSpw)))
            if (len(vm.spwInfo) < np.max(originalSpw)):
                # Then there was an extra split
                i = ictr
            else:
                i = originalSpw[ictr]
            nchan = len(vm.spwInfo[i]["chanFreqs"])
            if (nchan > 1):
                missingFrequencyWidth = originalChannelStart[ictr]*(vm.spwInfo[i]["chanFreqs"][-1]-vm.spwInfo[i]["chanFreqs"][0])/(nchan-1)
            else:
                missingFrequencyWidth = 0
            if (missingFrequencyWidth > 0):
                if (DEBUG):
                    print("Correcting for channels flagged prior to running bandpass by %f GHz" % (missingFrequencyWidth*1e-9))
            newfreqs = 1e-9*(vm.spwInfo[i]["chanFreqs"]) + missingFrequencyWidth*1e-9
#        if (debug): print "Appending onto chanFreqGHz: ", newfreqs
        chanFreqGHz.append(newfreqs)

  uniqueSpwsInCalTable = np.unique(cal_desc_id)

  # initial calculation for final message if not all spws appear with overlay='antenna'
  uniqueTimes = sloppyUnique(np.unique(times), 1.0)
  nUniqueTimes = len(uniqueTimes)
  if (nUniqueTimes == 1):
      solutionTimeSpread = 0
  else:
      solutionTimeSpread = np.max(uniqueTimes)-np.min(uniqueTimes)
  print("Found solutions with %d unique times (within a threshold of 1.0 second)." % (nUniqueTimes))


  uniqueTimes = sloppyUnique(np.unique(times), solutionTimeThresholdSeconds)
  nUniqueTimes = len(uniqueTimes)
  if (nUniqueTimes == 1):
      print("Found solutions with %d unique time (within a threshold of %d seconds)." % (nUniqueTimes,solutionTimeThresholdSeconds))
  else:
      print("Found solutions with %d unique times (within a threshold of %d seconds)." % (nUniqueTimes,solutionTimeThresholdSeconds))
  displayTimesArray([[uniqueTimes]])
  scansForUniqueTimes = []
  if (tableFormat >= 34):
      if (len(unique_cal_scans) == 1):
          print("Found solutions with %d unique scan number %s" % (len(unique_cal_scans), str(unique_cal_scans)))
      else:
          print("Found solutions with %d unique scan numbers %s" % (len(unique_cal_scans), str(unique_cal_scans)))

      scansForUniqueTimes, nUniqueTimes = computeScansForUniqueTimes(uniqueTimes, cal_scans, times, unique_cal_scans)
  elif (scans != ''):
      print("Selection by scan is not support for old-style tables that do not have the scan number filled.")
      return
  uniqueTimesCopy = uniqueTimes[:]

  mystring = ''
  if (debug):
     for u in uniqueTimes:
         mystring += '%.6f, ' % (u)
     print(mystring)
  uniqueAntennaIds = np.unique(ant)
  uniqueFields = np.unique(fields)
  nFields = len(uniqueFields)
  spwlist = []
  uniqueTimesPerFieldPerSpw = []
  for s in uniqueSpwsInCalTable:
        uniqueTimesPerField = []
        for f in uniqueFields:
            timelist = []
            for row in range(len(fields)):
                if (fields[row] == f and cal_desc_id[row] == s):
                    if (sloppyMatch(times[row], timelist, solutionTimeThresholdSeconds) == False):
                        timelist.append(times[row])
                        spwlist.append(cal_desc_id)
            uniqueTimesPerField.append(timelist)
        uniqueTimesPerFieldPerSpw.append(uniqueTimesPerField)
#  print "uniqueTimesPerFieldPerSpw (len=%d) = %s" % (len(uniqueTimesPerFieldPerSpw), str(uniqueTimesPerFieldPerSpw))

  # Parse the spws to plot from the command line
  if (spw==''):
     spwsToPlot = list(uniqueSpwsInCalTable)
  else:
     if (type(spw) == str):
           if (spw.find('!')>=0):
                 print("The ! modifier is not (yet) supported")
                 return(vm)
           tokens = spw.split(',')
           spwsToPlot = []
           for token in tokens:
                 if (len(token) > 0):
                       if (token.find('*')>=0):
                             spwsToPlot = list(uniqueSpwsInCalTable)
                             break
                       elif (token.find('~')>0):
                             (start,finish) = token.split('~')
                             spwsToPlot +=  range(int(start),int(finish)+1)
                       else:
                             spwsToPlot.append(int(token))
     elif (type(spw) == list):
         spwsToPlot = np.sort(spw)
     else:
         spwsToPlot = [spw]

  if (len(uniqueSpwsInCalTable) > 1):
      print("%d spws in the solution = " % (len(uniqueSpwsInCalTable)), uniqueSpwsInCalTable)
  else:
      print("%d spw in the solution = " % (len(uniqueSpwsInCalTable)), uniqueSpwsInCalTable)
  keepSpwsToPlot = spwsToPlot[:]
  for myspw in spwsToPlot:
      if (myspw not in uniqueSpwsInCalTable):
          print("WARNING: spw %d is not in the solution. Removing it from the list to plot." % (myspw))
          print("Available spws = ", uniqueSpwsInCalTable)
          keepSpwsToPlot.remove(myspw)
          if (casaVersion >= '4.1.0' and mymsmd != ''):
#              nonwvrspws = list(set(range(mymsmd.nspw())).difference(set(mymsmd.wvrspws())))
              if (myspw not in range(mymsmd.nspw())):
                  print("FATAL: spw %d is not even in the ms.  There might be a bug in your script." % (myspw))
                  return
              elif (myspw in mymsmd.wvrspws()):
                  print("WARNING: spw %d is a WVR spw." % (myspw))
                  return
  spwsToPlot = keepSpwsToPlot[:]
  if (spwsToPlot == []):
      print("FATAL: no spws to plot")
      return
  originalSpwsToPlot = computeOriginalSpwsToPlot(spwsToPlot, originalSpw, tableFormat, debug)
         
  # Now generate the list of minimal basebands that contain the spws to be plotted
  if (casaVersion >= '4.1.0' and msFound):
      allBasebands = []
      if (mymsmd != ''):
          try:
              for spw in originalSpwsToPlot:
                  mybaseband = mymsmd.baseband(spw)
                  if (debug): print("appending: spw=%d -> bb=%d" % (spw,mybaseband))
                  allBasebands.append(mybaseband)
              allBasebands = np.unique(allBasebands)
              basebandDict = getBasebandDict(msName,caltable=caltable,mymsmd=mymsmd)  # needed later by showFDM()
          except:
              basebandDict = {}
              print("This dataset (%s) does not have a BBC_NO column in the SPECTRAL_WINDOW_TABLE." % (msName))
      else:
          basebandDict = {}
          telescopeName = getTelescopeNameFromCaltable(caltable)
          print("Measurement set not found.")
      if (basebandDict == {}):
          if (overlay.find('spw') >= 0):
              print("As such, since the ms cannot be found, overlay='spw' is not supported, but overlay='baseband' should work.")
              return
          if (showfdm):
              print("As such, since the ms cannot be found, showfdm=True is not supported.")
              showfdm = False
          if (showBasebandNumber):
              print("As such, since the ms cannot be found, showBasebandNumber=True is not supported.")
              showBasebandNumber = False

  elif (msFound==False):
      allBasebands = [1,2,3,4]
  else:
      basebandDict = getBasebandDict(msName,caltable=caltable,mymsmd=mymsmd)  # needed later by showFDM()
      allBasebands = []
      for spw in originalSpwsToPlot:
          mybaseband = [key for key in basebandDict if spw in basebandDict[key]]
          if (len(mybaseband)>0): allBasebands.append(mybaseband[0])
      allBasebands = np.unique(allBasebands)
      if (allBasebands == []):
          allBasebands = [1,2,3,4]
  if (debug):
      print("================ allBasebands = ", allBasebands)
      
  if (basebands is None or basebands == ''):
      basebands = allBasebands
  elif (type(basebands) == str):
      basebands = [int(s) for s in basebands.split(',')]
  elif (type(basebands) != list):
      # it is a single integer
      basebands = [basebands]
  for baseband in basebands:
      if (baseband not in allBasebands):
          print("Baseband %d is not in the dataset (only %s)" % (baseband,str(allBasebands)))
          return
  if (msFound):
      msFieldsList = str(np.array(msFields)[uniqueFields])
  else:
      msFieldsList = 'unknown'
  print("%d field(s) in the solution = %s = %s" % (len(uniqueFields), uniqueFields,msFieldsList))
  
  # Figure out which kind of Bandpass solution this is.
  bOverlay = False  # Am I trying to overlay a second B-type solution?
  if (os.path.exists(caltable) == False):
        print("Caltable does not exist = %s" % (caltable))
        return(vm)
  try:
      ([polyMode, polyType, nPolyAmp, nPolyPhase, scaleFactor, nRows, nSpws, nUniqueTimesBP, uniqueTimesBP,
        nPolarizations, frequencyLimits, increments, frequenciesGHz, polynomialPhase,
        polynomialAmplitude, timesBP, antennasBP, cal_desc_idBP, spwBP]) = openBpolyFile(caltable, debug)
      bpoly = True
      bpolyOverlay = bpolyOverlay2 = False
      if (xaxis.find('chan') >= 0):
          print("Sorry, but BPOLY solutions cannot be plotted with xaxis='chan'. Proceeding with xaxis='freq'.")
          xaxis = 'freq'
      if (chanrange[0] != 0 or chanrange[1] != 0 or chanrangePercent is not None):
          print("The chanrange parameter only applies if the first caltable is a B solution, not a BPOLY.")
          return(vm)
      if (len(caltable2) > 0):
          try:
              # figure out if the next file is a BPOLY or another B solution to pick the proper error message.
              ([polyMode, polyType, nPolyAmp, nPolyPhase, scaleFactor, nRows, nSpws, nUniqueTimesBP, uniqueTimesBP,
                nPolarizations, frequencyLimits, increments, frequenciesGHz, polynomialPhase,
                polynomialAmplitude, timesBP, antennasBP, cal_desc_idBP, spwBP]) = openBpolyFile(caltable2, debug)
              print("Sorry, but you cannot overlay two BPOLY solutions (unless caltable is a B solution and caltable2 and 3 are BPOLYs).")
          except:
              print("Sorry, but for overlays, caltable must be a B solution, while caltable2 and 3 can be either type.")
          return(vm)
  except:
      print("caltable: This is a %s solution." % (VisCal))
      bpoly = bpolyOverlay = bpolyOverlay2 = False

      # Now check if there is a second file to overlay
      if (len(caltable2) > 0):
        if (os.path.exists(caltable2) == False):
              print("Caltable2 does not exist = %s" % (caltable2))
              return(vm)
        try:
          # figure out if the next file is a BPOLY or another B solution
          if (debug): print("Calling openBpolyFile('%s')" % (caltable2))
          ([polyMode, polyType, nPolyAmp, nPolyPhase, scaleFactor, nRows, nSpws, nUniqueTimesBP, uniqueTimesBP,
            nPolarizations, frequencyLimits, increments, frequenciesGHz, polynomialPhase,
            polynomialAmplitude, timesBP, antennasBP, cal_desc_idBP, spwBP]) = openBpolyFile(caltable2, debug)
          if (debug): print("Done")
          bpolyOverlay = True
          if (debug): print("Overlay the BPOLY solution")
          if (xaxis.find('chan')>=0):
              print("Sorry, but overlap of BPOLY is currently possible only with xaxis='freq'")
              return(vm)
          if (len(caltable3) > 0):
             if (os.path.exists(caltable3) == False):
                   print("Caltable3 does not exist = %s" % (caltable3))
                   return(vm)
             bpolyOverlay2 = True
             if (debug): print("Overlay the second BPOLY solution")
             ([polyMode2, polyType2, nPolyAmp2, nPolyPhase2, scaleFactor2, nRows2, nSpws2,
               nUniqueTimesBP2, uniqueTimesBP2,
               nPolarizations2, frequencyLimits2, increments2, frequenciesGHz2, polynomialPhase2,
               polynomialAmplitude2, timesBP2, antennasBP2, cal_desc_idBP2, spwBP2]) = openBpolyFile(caltable3, debug)
        except:
            # this is another B solution
            print("Overlay another %s solution" % (VisCal))
            bOverlay = True
            if (xaxis.find('freq')<0):
                  print("Currently, you must use xaxis='freq' to overlay two B solutions.")
                  return(vm)
            if (len(caltable3) > 0):
                  print("You cannot overlay caltable3 because caltable2 is a B solution.")
                  return(vm)
      elif (len(caltable3) > 0):
          print("You cannot have a caltable3 argument without a caltable2 argument.")
          return(vm)
          
  if (overlay.find('antenna')>=0):
      overlayAntennas = True
      if (bpoly == True):
            print("The overlay of times or antennas is not supported with BPOLY solutions")
            return(vm)
      if (len(caltable2)>0):
            print("The overlay of times or antennas not supported when overlaying a B or BPOLY solution")
            return(vm)
      print("Will overlay solutions from different antennas")
  else:
      overlayAntennas = False

  if (overlay.find('time')>=0):
      overlayTimes = True
      if (bpoly == True):
            print("The overlay of times or antennas is not supported with BPOLY solutions")
            return(vm)
      if (len(caltable2)>0):
            print("The overlay of times or antennas not supported when overlaying a B or BPOLY solution")
            return(vm)
      print("Will overlay solutions from different times")
  else:
      overlayTimes = False
      
  if (overlay.find('spw')>=0):
      if (tableFormat < 34):
          print("Overlay spw may not work reliably for old cal tables")
      overlaySpws = True
      if (bpoly == True):
            print("The overlay of times, antennas, or spws is not supported with BPOLY solutions")
            return(vm)
      if (len(caltable2)>0):
            print("The overlay of times, antennas, or spws not supported when overlaying a B or BPOLY solution")
            return(vm)
      print("Will overlay solutions from different spws within a baseband")
  else:
      overlaySpws = False
      
  if (overlay.find('baseband')>=0):
      if (tableFormat < 34):
          print("Overlay baseband may not work reliably for old cal tables")
      overlayBasebands = True
      if (bpoly == True):
            print("The overlay of times, antennas, spws, or basebands is not supported with BPOLY solutions")
            return(vm)
      if (len(caltable2)>0):
            print("The overlay of times, antennas, spws, or basebands not supported when overlaying a B or BPOLY solution")
            return(vm)
      print("Will overlay solutions from all spws regardless of baseband")
  else:
      overlayBasebands = False

  if (bOverlay):        
        # Now open the Bandpass solution table
        try:
              mytb.open(caltable2)
        except:
              print("Could not open the second caltable = %s" % (caltable2))
              return(vm)
        names = mytb.colnames()
        ant2 = mytb.getcol('ANTENNA1')
        fields2 = mytb.getcol('FIELD_ID')
        times2 = mytb.getcol('TIME')
        if ('SPECTRAL_WINDOW_ID' not in names):
            if ('SNR' not in names):
                print("This does not appear to be a cal table.")
                return(vm)
            else:
                tableFormat2 = 33
                print("This appears to be an old-format cal table from casa 3.3 or earlier.")
                cal_desc_id2 = mytb.getcol('CAL_DESC_ID')
                VisCal2 = (mytb.info())['subType']
                mytb.close()
                ParType = "unknown"  # i.e. not Complex
                calDesc2 = mytb.open(caltable2+'/CAL_DESC')
                originalSpws2 = mytb.getcol('SPECTRAL_WINDOW_ID')  # [[0,1,2,3]]
                originalSpw2 = originalSpws2[0]                   # [0,1,2,3]
                msName2 = mytb.getcol('MS_NAME')[0]
                mytb.close()
                # Now open the associated ms tables via ValueMapping to figure out channel freqs
                chanFreqGHz2 = []
                for ictr in range(len(originalSpw2)):
                    if debug: print("ictr = ", ictr)
                    if (casaVersion >= '4.1.0'):
                        if debug: print("mymsmd.nspw() = %d, np.max(originalSpw) = %d" % (mymsmd.nspw(),np.max(originalSpw2)))
                        if (mymsmd.nspw() < np.max(originalSpw2)):
                            # Then there was an extra split
                            i = ictr
                        else:
                            i = originalSpw2[ictr]
                        nchan = len(mymsmd.chanfreqs(i))
                        if (nchan > 1):
                            missingFrequencyWidth = originalChannelStart[ictr]*(mymsmd.chanfreqs(i)[-1]-mymsmd.chanfreqs(i)[0])/(mymsmd.nchan(i))
                        else:
                            missingFrequencyWidth = 0
                        if (missingFrequencyWidth > 0):
                            if (DEBUG):
                                print("Correcting for channels flagged prior to running bandpass by %f GHz" % (missingFrequencyWidth*1e-9))
                        newfreqs = 1e-9*(mymsmd.chanfreqs(i)) + missingFrequencyWidth*1e-9
#                        if debug: print "Appending onto chanFreqGHz: ", newfreqs
                        chanFreqGHz2.append(newfreqs)
                    else:
                        if debug: print("len(vm.spwInfo) = %d, np.max(originalSpw) = %d" % (len(vm.spwInfo),np.max(originalSpw2)))
                        if (len(vm.spwInfo) < np.max(originalSpw2)):
                            # Then there was an extra split
                            i = ictr
                        else:
                            i = originalSpw2[ictr]
                        nchan = len(vm.spwInfo[i]["chanFreqs"])
                        if (nchan > 1):
                            missingFrequencyWidth = originalChannelStart[ictr]*(vm.spwInfo[i]["chanFreqs"][-1]-vm.spwInfo[i]["chanFreqs"][0])/(nchan-1)
                        else:
                            missingFrequencyWidth = 0
                        if (missingFrequencyWidth > 0):
                            if (DEBUG):
                                print("Correcting for channels flagged prior to running bandpass by %f GHz" % (missingFrequencyWidth*1e-9))
                        newfreqs = 1e-9*(vm.spwInfo[i]["chanFreqs"]) + missingFrequencyWidth*1e-9
#                        if debug: print "Appending onto chanFreqGHz: ", newfreqs
                        chanFreqGHz2.append(newfreqs)
                    
        else:
            tableFormat2 = 34
            cal_desc_id2 = mytb.getcol('SPECTRAL_WINDOW_ID')
            msName2 = mytb.getkeyword('MSName')      
            ParType2 = mytb.getkeyword('ParType')    # string = 'Complex'
            VisCal2 = mytb.getkeyword('VisCal')      # string = 'B TSYS'
            PolBasis2 = mytb.getkeyword('PolBasis')  # string = 'LINEAR'
            spectralWindowTable2 = mytb.getkeyword('SPECTRAL_WINDOW').split()[1]
            mytb.close()
            mytb.open(spectralWindowTable2)
            chanFreqGHz2 = []
            originalSpws2 = range(len(mytb.getcol('MEAS_FREQ_REF')))
            for i in originalSpws2:
                # The array shapes can vary.
                chanFreqGHz2.append(1e-9 * mytb.getcell('CHAN_FREQ',i))
            originalSpws2 = range(len(mytb.getcol('MEAS_FREQ_REF')))
            originalSpw2 = originalSpws2  # may want to do a global replace of this <----------------------------------

        uniqueSpwsInCalTable2 = np.unique(cal_desc_id2)
        mytb.open(caltable2)
        if 'FLAG' in mytb.colnames():
            flags2 = {}
            for f in range(len(fields2)):
                if mytb.iscelldefined('FLAG',f):
                    flags2[f] = mytb.getcell('FLAG',f)
                else:
                    flags2[f] = np.zeros(len(chanFreqGHz[f]))
        else:
            print("bOverlay: No Flag column found. Are you sure this is a bandpass solution file, or is it the .ms?")
            print("If it is a solution file, does it contain solutions for both TDM and FDM spws?")
            return(vm)
        uniqueTimes2 = sloppyUnique(np.unique(times2), solutionTimeThresholdSeconds)
        nUniqueTimes2 = len(uniqueTimes2)
#        print "Found %d solutions in time: MJD seconds = " % (nUniqueTimes2), uniqueTimes2
        spacing = ''
        for i in range(1,nUniqueTimes2):
            spacing += '%.0f, ' % (np.abs(uniqueTimes2[i]-uniqueTimes2[i-1]))
        print("Found %d solutions in time, spaced by seconds: " % (nUniqueTimes2), spacing)
        displayTimesArray([[uniqueTimes2]])
        uniqueAntennaIds2 = np.unique(ant2)
        uniqueFields2 = np.unique(fields2)
        nFields2 = len(uniqueFields2)

        if (debug): print("(boverlay) original unique spws in the second dataset = ", np.unique(originalSpw2))

        uniqueTimesPerFieldPerSpw2 = []
        for s in uniqueSpwsInCalTable2:
            uniqueTimesPerField2 = []
            for f in uniqueFields2:
                timelist2 = []
                for row in range(len(fields2)):
                      if (fields2[row] == f and cal_desc_id2[row] == s):
                          if (sloppyMatch(times2[row], timelist2, solutionTimeThresholdSeconds,
                                          myprint=False) == False):
                              timelist2.append(times2[row])
                uniqueTimesPerField2.append(timelist2)
            uniqueTimesPerFieldPerSpw2.append(uniqueTimesPerField2)
        if debug: 
            print("uniqueTimesPerFieldPerSpw2 = ") #, uniqueTimesPerFieldPerSpw2
        displayTimesArray(uniqueTimesPerFieldPerSpw2)
          
        if debug: 
            print("%d spw(s) in the second solution = " % (len(uniqueSpwsInCalTable2)), uniqueSpwsInCalTable2)
        if (msFound):
            msFieldsList = str(np.array(msFields)[uniqueFields2])
        else:
            msFieldsList = 'unknown'
        print("%d field(s) in the solution = %s = %s" % (len(uniqueFields2), uniqueFields2, msFieldsList))

  # Parse the timeranges field from the command line
  if timeranges != '':
      timerangesWasSpecified = True
  else:
      timerangesWasSpecified = False
  if (type(timeranges) == str):
         # a list of timeranges was given
         tokens = timeranges.split(',')
         timerangeList = []
         removeTime = []
         for token in tokens:
             if (len(token) > 0):
                 if (token.find('!')==0):
                     timerangeList = range(len(uniqueTimes))
                     removeTime.append(int(token[1:]))
                 elif (token.find('~')>0):
                     (start,finish) = token.split('~')
                     timerangeList +=  range(int(start),int(finish)+1)
                 else:
                     timerangeList.append(int(token))
         timerangeList = np.array(timerangeList)
         for rt in removeTime:
             timerangeList = timerangeList[np.where(timerangeList != rt)[0]]
         timerangeList = list(timerangeList)
         if (len(timerangeList) < 1):
            if (len(removeTime) > 0):
                print("Too many negated timeranges -- there are none left to plot.")
                return
            else:
                # then a blank list was specified
                timerangeList = range(len(uniqueTimes))
  elif (type(timeranges) == list):
      # it's already a list of integers
      timerangeList = timeranges
  else:
      # It's a single, integer entry
      timerangeList = [timeranges]

  if (timerangesWasSpecified and scans != ''):  # CAS-8489
      if (type(scans) == list or type(scans) == np.ndarray):
          myscan = int(scans[0])
      else:
          myscan = int(str(scans).split(',')[0])
      if (myscan not in scansForUniqueTimes):
          print("No rows for scan %d, only " % (myscan), np.unique(scansForUniqueTimes))
          return
      timerangeOffset = scansForUniqueTimes.index(myscan)
      timerangeList = np.array(timerangeList) + timerangeOffset
      print("Since both timeranges and scans was specified, generated new effective timerangeList: ", timerangeList)
  if (max(timerangeList) >= len(uniqueTimes)):
      print("Invalid timerange.  Solution has %d times (%d~%d)" % (len(uniqueTimes),0,len(uniqueTimes)-1))
      return
  timerangeListTimes = np.array(uniqueTimes)[timerangeList]
  if (tableFormat == 33 or scansForUniqueTimes == []):
      # SMA data with scan numbers of -1 has empty list for scansForUniqueTimes
      scansToPlot = []
      if (scans != ''):
          print("Selection by scan is not possible for this dataset.")
          return
  else:
      if (debug): print("scansForUniqueTimes = %s" % (str(scansForUniqueTimes)))
      scansToPlot = np.array(scansForUniqueTimes)[timerangeList]
      if (np.unique(scansToPlot)[0] == -1):
          # scan numbers are not correct in this new-style cal table
          scansToPlot = []
          if (scans != ''):
              print("Selection by scan number is not possible with this dataset.")
              return
      if (scans != '' and scans != []):
          if (type(scans) == list):
              scansToPlot = scans
          elif (type(scans) == str):
              scansToPlot = [int(a) for a in scans.split(',')]
          else:
              scansToPlot = [scans]
          for scan in scansToPlot:
              if (scan not in scansForUniqueTimes):
                  print("Scan %d is not in any solution" % (scan))
                  return
      print("scans to plot: %s" % (str(scansToPlot)))
  scansToPlotPerSpw = {}
  for myspw in np.unique(cal_desc_id):
      scansToPlotPerSpw[myspw] = []
  for scan in scansToPlot:
      for myspw in np.unique(cal_desc_id):
          if (scan in cal_scans_per_spw[myspw]):
              scansToPlotPerSpw[myspw].append(scan)

  # remove spws that do not have any scans to be plotted
  # but only for tables that have a scan number column, and not filled with all -1
  if (tableFormat > 33 and scansForUniqueTimes != []):
      for myspw in np.unique(cal_desc_id):
          if (debug):
              print("scans to plot for spw %d: %s" % (myspw, scansToPlotPerSpw[myspw]))
          if (scansToPlotPerSpw[myspw] == []):
              indexDelete = np.where(np.array(spwsToPlot)==myspw)[0]
#              print "indexDelete = ", indexDelete
              if (len(indexDelete) > 0):
                  spwsToPlot = list(np.delete(spwsToPlot, indexDelete[0]))
      print("spwsToPlot = ", spwsToPlot)
  timerangeListTimesString = mjdsecArrayToUTString(timerangeListTimes)
  print("UT times to plot:  ", timerangeListTimesString)
  print("%d MJDSecond to plot: " % (len(timerangeListTimes)), str([int(b) for b in timerangeListTimes]))
  if (len(timerangeListTimes) > len(np.unique(scansToPlot))):
      # fix for CAS-9474
      uniqueScansToPlot, idx = np.unique(scansToPlot, return_index=True)
      if (len(uniqueScansToPlot) < len(scansToPlot)):
          # If the solution time for one spw differs by more than solutionTimeThresholdSeconds from
          # another spw, then we will get 2 identical entries for the same scan, and thus duplicate 
          # plots.  So, remove one.
          print("Engaging fix for CAS-9474")
          scansToPlot = uniqueScansToPlot
          timerangeListTimes = list(np.array(timerangeListTimes)[idx])
          timerangeList = list(np.array(timerangeList)[idx])
          timerangeListTimesString = mjdsecArrayToUTString(timerangeListTimes)
          print("Revised scans to plot: %s" % (str(scansToPlot)))
          print("Revised UT times to plot:  ", timerangeListTimesString)
          print("%d MJDSecond to plot: " % (len(timerangeListTimes)), str([int(b) for b in timerangeListTimes]))
          print("Corresponding time IDs (0-based):", timerangeList) # OneBased 

  # Check for mismatch
  if (bpolyOverlay):
      if (len(timerangeListTimes) > nUniqueTimesBP):
          print("There are more timeranges (%d) to plot from %s than exist in the caltable2=%s (%d)" % (len(timerangeListTimes), caltable,caltable2, nUniqueTimesBP))
          for i in timerangeList:
              if (sloppyMatch(timerangeListTimes[i],uniqueTimesBP[0],
                              solutionTimeThresholdSeconds, mytime,
                              scansToPlot, scansForUniqueTimes, False)):
                  print("Try adding 'timeranges=%d'" % (i+1))
          return(vm)
      if (bpolyOverlay2):
          if (len(timerangeListTimes) > nUniqueTimesBP2):
              print("There are more timeranges to plot (%d) from %s than exist in the caltable3=%s (%d)" % (len(timerangeListTimes), caltable, caltable3, nUniqueTimesBP2))
              return(vm)
          
  # Parse the antenna string to emulate plotms
  if (type(antenna) == str):
     if (len(antenna) == sum([m in myValidCharacterListWithBang for m in antenna])):
         # a simple list of antenna numbers was given 
         tokens = antenna.split(',')
         antlist = []
         removeAntenna = []
         for token in tokens:
             if (len(token) > 0):
                 if (token.find('*')==0 and len(token)==1):
                     antlist = uniqueAntennaIds
                     break
                 elif (token.find('!')==0):
                     antlist = uniqueAntennaIds
                     removeAntenna.append(int(token[1:]))
                 elif (token.find('~')>0):
                     (start,finish) = token.split('~')
                     antlist +=  range(int(start),int(finish)+1)
                 else:
                     antlist.append(int(token))
         antlist = np.array(antlist)
         for rm in removeAntenna:
             antlist = antlist[np.where(antlist != rm)[0]]
         antlist = list(antlist)
         if (len(antlist) < 1 and len(removeAntenna)>0):
             print("Too many negated antennas -- there are no antennas left to plot.")
             return(vm)
     else:
         # The antenna name (or list of names) was specified
         tokens = antenna.split(',')
         if (msFound):
             antlist = []
             removeAntenna = []
             for token in tokens:
               if (casaVersion >= '4.1.0'):
                 if (token in msAnt): 
                     antlist = list(antlist)  # needed in case preceding antenna had ! modifier
#                     antlist.append(mymsmd.antennaids(token)[0]) # This crashes if ms was not actually found
                     antlist.append(list(msAnt).index(token)) # This alternative should work in all cases
                 elif (token[0] == '!'):
                     if (token[1:] in msAnt):
                         antlist = uniqueAntennaIds
#                         removeAntenna.append(mymsmd.antennaids(token[1:])) # This crashes if ms was not actually found
                         removeAntenna.append(list(msAnt).index(token[1:])) # This alternative should work in all cases
                     else:
                         print("Antenna %s is not in the ms. It contains: " % (token), msAnt)
                         return(vm)
                 else:
                     print("Antenna %s is not in the ms. It contains: " % (token), msAnt)
                     return(vm)
               else:
                 if (token in vm.uniqueAntennas):
                     antlist = list(antlist)  # needed in case preceding antenna had ! modifier
                     antlist.append(vm.getAntennaIdsForAntennaName(token))
                 elif (token[0] == '!'):
                     if (token[1:] in vm.uniqueAntennas):
                         antlist = uniqueAntennaIds
                         removeAntenna.append(vm.getAntennaIdsForAntennaName(token[1:]))
                     else:
                         print("Antenna %s is not in the ms. It contains: " % (token), vm.uniqueAntennas)
                         return(vm)
                 else:
                     print("Antenna %s is not in the ms. It contains: " % (token), vm.uniqueAntennas)
                     return(vm)
             antlist = np.array(antlist)
             for rm in removeAntenna:
                 antlist = antlist[np.where(antlist != rm)[0]]
             antlist = list(antlist)
             if (len(antlist) < 1 and len(removeAntenna)>0):
                 print("Too many negated antennas -- there are no antennas left to plot.")
                 return(vm)
         else:
             print("Antennas cannot be specified my name if the ms is not found.")
             return(vm)
  elif (type(antenna) == list):
      # it's a list of integers
      antlist = antenna
  else:
      # It's a single, integer entry
      antlist = [antenna]

  if (len(antlist) > 0):
      antennasToPlot = np.intersect1d(uniqueAntennaIds,antlist)
  else:
      antennasToPlot = uniqueAntennaIds
  if (len(antennasToPlot) < 2 and overlayAntennas):
      print("More than 1 antenna is required for overlay='antenna'.")
      return
  
  # Parse the field string to emulate plotms
  removeField = []
  if (type(field) == str):
     if (len(field) == sum([m in myValidCharacterListWithBang for m in field])):
         # a list of field numbers was given
         tokens = field.split(',')
         fieldlist = []
         for token in tokens:
             if (token.find('*')>=0):
                 fieldlist = uniqueFields
                 break
             elif (token.find('!')==0):
                 fieldlist = uniqueFields
                 removeField.append(int(token[1:]))
             elif (len(token) > 0):
                 if (token.find('~')>0):
                     (start,finish) = token.split('~')
                     fieldlist +=  range(int(start),int(finish)+1)
                 else:
                     fieldlist.append(int(token))
         fieldlist = np.array(fieldlist)
         for rm in removeField:
             fieldlist = fieldlist[np.where(fieldlist != rm)[0]]
         fieldlist = list(fieldlist)
         if (len(fieldlist) < 1 and len(removeField)>0):
             print("Too many negated fields -- there are no fields left to plot.")
             return(vm)
     else:
         # The field name (or list of names, or wildcard) was specified
         tokens = field.split(',')
         if (msFound):
             fieldlist = []
             removeField = []
             for token in tokens:
                 myloc = token.find('*')
                 if (myloc > 0):   # saw wildcard in name
                     for u in uniqueFields:
                         # Sept 2014: will crash if ms is not found
                         myFieldName = GetFieldNamesForFieldId(u,vm,mymsmd,msFields)[:myloc]
                         if (token[:myloc]==myFieldName[:myloc]):
                             if (DEBUG):
                                 print("Found wildcard match = %s" % GetFieldNamesForFieldId(u,vm,mymsmd,msFields))
                             fieldlist.append(u)
                         else:
                             if (DEBUG):
                                 print("No wildcard match of %s with = %s" % (token[:myloc], GetFieldNamesForFieldId(u,vm,mymsmd,msFields)))
                 elif (myloc==0):  # saw wildcard at start of name
                     for u in uniqueFields:
                         fieldlist.append(u)
                 elif (token in msFields):
                     fieldlist = list(fieldlist)  # needed in case preceding field had ! modifier
                     fieldlist.append(GetFieldIdsForFieldName(token,vm,mymsmd,msFields))
                 elif (token[0] == '!'):
                     if (fieldlist == []):
                         for u in uniqueFields: # uniqueFields is a list of IDs
                             fieldlist.append(u)
                     if (token[1:] in msFields):
                         removeField.append(GetFieldIdsForFieldName(token[1:],vm,mymsmd,msFields))
                     else:
                         print("1) Field %s is not in the ms. It contains: " % (token), uniqueFields, msFieldsUnique)
                         return(vm)
                 else:
                     print("2) Field %s is not in the ms. It contains: " % (token), uniqueFields, msFieldsUnique)
                     return(vm)
             fieldlist = np.array(fieldlist)
             for rm in removeField:
                 fieldlist = fieldlist[np.where(fieldlist != rm)[0]]
             fieldlist = list(fieldlist)
             if (len(fieldlist) < 1 and len(removeField)>0):
                 print("Too many negated fields -- there are no fields left to plot.")
                 return(vm)
         else:
             print("Fields cannot be specified my name if the ms is not found.")
             return(vm)
  elif (type(field) == list):
      # it's a list of integers
      fieldlist = field
  else:
      # It's a single, integer entry
      fieldlist = [field]

  if (len(fieldlist) > 0):
      if (DEBUG):
          print("Finding intersection of ", uniqueFields, fieldlist)
      fieldsToPlot = np.intersect1d(uniqueFields,np.array(fieldlist))
      if (bOverlay):
          fieldsToPlot = np.intersect1d(np.union1d(uniqueFields,uniqueFields2),np.array(fieldlist))
      if (len(fieldsToPlot) < 1):
          print("Requested field not found in solution")
          return(vm)
  else:
      fieldsToPlot = uniqueFields  # use all fields if none are specified
      if (bOverlay):
          fieldsToPlot = np.union1d(uniqueFields,uniqueFields2)
      if (DEBUG):
          print("bOverlay = ", bOverlay)
          print("set fieldsToPlot to uniqueFields = ", fieldsToPlot)
  fieldIndicesToPlot = []

  if (showatmfield == ''):
      showatmfield = fieldsToPlot[0]
  else:
      if (str.isdigit(str(showatmfield))):
          showatmfield = int(str(showatmfield))
          if (showatmfield not in fieldsToPlot):
              print("The showatmfield (%d) is not in the list of fields to plot: " %(showatmfield), fieldsToPlot)
              return(vm)
      else:
          showatmfieldName = showatmfield
          showatmfield = GetFieldIdsForFieldName(showatmfield,vm,mymsmd,msFields)
          if (list(showatmfield) == []):
              print("The showatmfield (%s) is not in the ms." %(showatmfieldName))
              return(vm)
          if (type(showatmfield) == np.ndarray):
              # more than one field IDs exist for this source name, so pick the first
              showatmfield = showatmfield[0]
          if (showatmfield not in fieldsToPlot):
              print("The showatmfield (%d=%s) is not in the list of fields to plot: " %(showatmfield, showatmfieldName), fieldsToPlot)
              return(vm)

  for i in fieldsToPlot:
      match = np.where(i==uniqueFields)[0]
      if (len(match) < 1 and bOverlay):
          match = np.where(i==uniqueFields2)[0]
      fieldIndicesToPlot.append(match[0])
      
  print("Antennas to plot = ", antennasToPlot)
  print("spws to plot = ", spwsToPlot)

# I moved the following line upward to line 1483 on 06-Aug-2013.
#  originalSpwsToPlot = computeOriginalSpwsToPlot(spwsToPlot, originalSpw, tableFormat, debug) 

  if (debug): print("original spws to plot = ", originalSpwsToPlot)
  print("Field IDs to plot: ", fieldsToPlot)
#  print "Field indices to plot: ", fieldIndicesToPlot

  redisplay = False
  myap = 0  # this variable is necessary to make the 'b' option work for 
            # subplot=11, yaxis=both.  It keeps track of whether 'amp' or 
            # 'phase' was the first plot on the page.

  # I added the call to pb.ion() because Remy suggested it.
  if (pb.isinteractive()):
#      print "pylab is set to interactive"
      wasInteractive = True
  else:
#      print "pylab is set to non-interactive"
      wasInteractive = False
  if (interactive):
      if (wasInteractive==False):
          print("Turning on interactive pylab for the duration of this task.")
          pb.ion()
  else:
      if (wasInteractive):
          print("Turning off interactive pylab for the duration of this task.")
          pb.ioff()
# But the call to pb.figure() causes another window to open everytime.
#  pb.figure()

  newylimits = [LARGE_POSITIVE, LARGE_NEGATIVE]
  if casaVersion > '5.8':
      pb.style.use('classic')
  if (bpoly):
    # The number of polarizations cannot be reliably inferred from the shape of
    # the GAIN column in the caltable.  Must use the shape of the DATA column 
    # in the ms.
    if (msFound):
#        (corr_type, corr_type_string, nPolarizations) = getCorrType(msName,spwsToPlot,debug)
        (corr_type, corr_type_string, nPolarizations) = getCorrType(msName,originalSpwsToPlot,debug)
    else:
        print("With no ms available, I will assume ALMA data: XX, YY, and refFreq=first channel.")
        chanFreqGHz = []
        corr_type_string = ['XX','YY']
        corr_type = [9,12]
        nPolarizations = 2
    nPolarizations2 = nPolarizations
    if (corr_type_string == []):
        return(vm)
    polsToPlot = checkPolsToPlot(polsToPlot, corr_type_string)
    if (polsToPlot == []):
        return(vm)
    pb.clf()

    # Here we are only plotting one BPOLY solution, no overlays implemented.
    overlayAntennas = False
    # rows in the table are: antennas 0..nAnt for first spw, antennas 0..nAnt 
    # for 2nd spw...
    pagectr = 0
    pages = []
    xctr = 0
    newpage = 1
    while (xctr < len(antennasToPlot)):
      xant = antennasToPlot[xctr]
      antstring = buildAntString(xant,msFound,msAnt)
      spwctr = 0
      spwctrFirstToPlot = spwctr
      while (spwctr < len(spwsToPlot)):
       ispw = spwsToPlot[spwctr]
       mytime = 0
       while (mytime < nUniqueTimes):
         if (len(uniqueTimes) > 0 and (mytime not in timerangeList)):
             if (debug):
                 print("@@@@@@@@@@@@@@@  Skipping mytime=%d" % (mytime))
             mytime += 1
             continue
         if (newpage == 1):
            pages.append([xctr,spwctr,mytime,0])
            if (debug):
                print("1) appending [%d,%d,%d,%d]" % (xctr,spwctr,mytime,0))
            newpage = 0
         antennaString = 'Ant%2d: %s,  ' % (xant,antstring)
#         fieldIndex = -1 # use this to tell whether we found any match for CAS-7753
         for index in range(nRows):
            # Find this antenna, spw, and timerange combination in the table
            if (xant==ant[index] and sloppyMatch(uniqueTimes[mytime],times[index],solutionTimeThresholdSeconds,
                                                 mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes,
                                                 myprint=debugSloppyMatch) and
                (ispw == cal_desc_id[index]) and (fields[index] in fieldsToPlot)):
                fieldIndex = np.where(fields[index] == uniqueFields)[0]
                if (type(fieldIndex) == list or type(fieldIndex) == np.ndarray):
                    fieldIndex = fieldIndex[0]
                validDomain = [frequencyLimits[0,index], frequencyLimits[1,index]]
                if (msFound):
                    fieldString = msFields[uniqueFields[fieldIndex]] # [0]
                else:
                    fieldString = str(field)
                timeString = ',  t%d/%d  %s' % (mytime,nUniqueTimes-1,utstring(uniqueTimes[mytime],3))
                if (scansForUniqueTimes != []):
                    if (scansForUniqueTimes[mytime]>=0):
                        timeString = ',  scan%d  %s' % (scansForUniqueTimes[mytime],utstring(uniqueTimes[mytime],3))
                if ((yaxis.find('amp')>=0 or amplitudeWithPhase) and myap==0):
                  xframe += 1
                  myUniqueColor = []
                  if (debug):
                      print("v) incrementing xframe to %d" % xframe)
                  adesc = pb.subplot(xframe)
                  previousSubplot = xframe
                  if (ispw==originalSpw[ispw]):
                      # all this was added mistakenly here.  If it causes a bug, remove it.
                      if (overlayTimes and len(fieldsToPlot) > 1):
                        indices = fstring = ''
                        for f in fieldIndicesToPlot:
                            if (f != fieldIndicesToPlot[0]):
                                indices += ','
                                fstring += ','
                            indices += str(uniqueFields[f])
                            if (msFound):
                                fstring += msFields[uniqueFields[f]]
                        if (len(fstring) > fstringLimit):
                            fstring = fstring[0:fstringLimit] + '...'
                        pb.title("%sspw%2d,  fields %s: %s%s" % (antennaString,ispw,
                                indices, fstring, timeString), size=titlesize)
                      else:
                        pb.title("%sspw%2d,  field %d: %s%s" % (antennaString,ispw,
                                uniqueFields[fieldIndex],fieldString,timeString), size=titlesize)
                  else:
                      if (overlayTimes and len(fieldsToPlot) > 1):
                        indices = fstring = ''
                        for f in fieldIndicesToPlot:
                            if (f != fieldIndicesToPlot[0]):
                                indices += ','
                                fstring += ','
                            indices += str(uniqueFields[f])
                            if (msFound):
                                fstring += msFields[uniqueFields[f]]
                        if (len(fstring) > fstringLimit):
                            fstring = fstring[0:fstringLimit] + '...'
                        pb.title("%sspw%2d (%d),  fields %s: %s%s" % (antennaString,ispw,originalSpw[ispw],
                                indices, fstring, timeString), size=titlesize)
                      else:
                        pb.title("%sspw%2d (%d),  field %d: %s%s" % (antennaString,ispw,originalSpw[ispw],
                                uniqueFields[fieldIndex],fieldString,timeString), size=titlesize)
                  amplitudeSolutionX = np.real(scaleFactor[index])+calcChebyshev(polynomialAmplitude[index][0:nPolyAmp[index]], validDomain, frequenciesGHz[index]*1e+9)
                  amplitudeSolutionY = np.real(scaleFactor[index])+calcChebyshev(polynomialAmplitude[index][nPolyAmp[index]:2*nPolyAmp[index]], validDomain, frequenciesGHz[index]*1e+9)
                  amplitudeSolutionX += 1 - np.mean(amplitudeSolutionX)
                  amplitudeSolutionY += 1 - np.mean(amplitudeSolutionY)
                  if (yaxis.lower().find('db') >= 0):
                      amplitudeSolutionX = 10*np.log10(amplitudeSolutionX)
                      amplitudeSolutionY = 10*np.log10(amplitudeSolutionY)
                  if (nPolarizations == 1):
                      pb.plot(frequenciesGHz[index], amplitudeSolutionX, '%s%s'%(xcolor,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                  else:
                      pb.plot(frequenciesGHz[index], amplitudeSolutionX, '%s%s'%(xcolor,bpolymarkstyle), frequenciesGHz[index], amplitudeSolutionY, '%s%s'%(ycolor,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                  if (plotrange[0] != 0 or plotrange[1] != 0):
                      SetNewXLimits([plotrange[0],plotrange[1]],loc=0)
                  if (plotrange[2] != 0 or plotrange[3] != 0):
                      SetNewYLimits([plotrange[2],plotrange[3]],loc=0)
                  xlim=pb.xlim()
                  ylim=pb.ylim()
                  ResizeFonts(adesc,mysize)
                  adesc.xaxis.grid(True,which='major')
                  adesc.yaxis.grid(True,which='major')
                  if (yaxis.lower().find('db')>=0):
                      pb.ylabel('Amplitude (dB)', size=mysize)
                  else:
                      pb.ylabel('Amplitude', size=mysize)
                  pb.xlabel('Frequency (GHz) (%d channels)'%(len(frequenciesGHz[index])), size=mysize)
                  if (xframe == firstFrame):
                      DrawBottomLegendPageCoords(msName, uniqueTimes[mytime], mysize, figfile)
                      pb.text(xstartTitle, ystartTitle,
                              '%s (degamp=%d, degphase=%d)'%(caltable,nPolyAmp[index]-1,
                              nPolyPhase[index]-1),size=mysize,
                              transform=pb.gcf().transFigure) 
                  # draw polarization labels (BPOLY case)
                  x0 = xstartPolLabel
                  y0 = ystartPolLabel
                  for p in range(nPolarizations):
                      if (corrTypeToString(corr_type[p]) in polsToPlot):
                          pb.text(x0, y0-0.03*subplotRows*p, corrTypeToString(corr_type[p])+'',
                                  color=pcolor[p],size=mysize, transform=pb.gca().transAxes)
                      
                  if (xframe == 111 and amplitudeWithPhase):
                     if (len(figfile) > 0):
                         # We need to make a new figure page
                         plotfiles.append(makeplot(figfile,msFound,msAnt,
                                          overlayAntennas,pages,pagectr,
                                          density,interactive,antennasToPlot,
                                          spwsToPlot,overlayTimes,overlayBasebands,
                                                   0,xant,ispw,subplot,resample,
                                          debug,figfileSequential,figfileNumber))
                         figfileNumber += 1
                     donetime = timeUtilities.time()
                     if (interactive):
                        pb.draw()
#                        myinput = raw_input("(%.1f sec) Press return for next page (b for backwards, q to quit): "%(donetime-mytimestamp))
                        myinput = raw_input("Press return for next page (b for backwards, q to quit): ")
                     else:
                        myinput = ''
                     skippingSpwMessageSent = 0
                     mytimestamp = timeUtilities.time()
                     if (myinput.find('q') >= 0):
#                         showFinalMessage(overlayAntennas, solutionTimeSpread, nUniqueTimes)
                         return(vm)
                     if (myinput.find('b') >= 0):
                         if (pagectr > 0):
                             pagectr -= 1
                         #redisplay the current page by setting ctrs back to the value they had at start of that page
                         xctr = pages[pagectr][PAGE_ANT]
                         spwctr = pages[pagectr][PAGE_SPW]
                         mytime = pages[pagectr][PAGE_TIME]
                         myap = pages[pagectr][PAGE_AP]
                         xant = antennasToPlot[xctr]
                         antstring = buildAntString(xant,msFound,msAnt)
                         ispw = spwsToPlot[spwctr]
                         redisplay = True
                     else:
                         pagectr += 1
                         if (pagectr >= len(pages)):
                               pages.append([xctr,spwctr,mytime,1])
                               if (debug):
                                   print("2) appending [%d,%d,%d,%d]" % (xctr,spwctr,mytime,1))
                               newpage = 0
                     pb.clf()

                if (yaxis.find('phase')>=0 or amplitudeWithPhase):
                  xframe += 1
                  myUniqueColor = []
                  if (debug):
                      print("w) incrementing xframe to %d" % xframe)
                  adesc = pb.subplot(xframe)
                  previousSubplot = xframe
                  if (ispw==originalSpw[ispw]):
                        pb.title("%sspw%2d,  field %d: %s%s" % (antennaString,ispw,
                               uniqueFields[fieldIndex],fieldString,timeString), size=titlesize)
                  else:
                        pb.title("%sspw%2d (%d),  field %d: %s%s" % (antennaString,ispw,originalSpw[ispw],
                               uniqueFields[fieldIndex],fieldString,timeString), size=titlesize)
                  phaseSolutionX = calcChebyshev(polynomialPhase[index][0:nPolyPhase[index]], validDomain, frequenciesGHz[index]*1e+9) * 180/math.pi
                  phaseSolutionY = calcChebyshev(polynomialPhase[index][nPolyPhase[index]:2*nPolyPhase[index]], validDomain, frequenciesGHz[index]*1e+9) * 180/math.pi
                  if (nPolarizations == 1):
                      pb.plot(frequenciesGHz[index], phaseSolutionX, '%s%s'%(xcolor,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                  else:
                      pb.plot(frequenciesGHz[index], phaseSolutionX, '%s%s'%(xcolor,bpolymarkstyle), frequenciesGHz[index], phaseSolutionY, '%s%s'%(ycolor,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                  ResizeFonts(adesc,mysize)
                  adesc.xaxis.grid(True,which='major')
                  adesc.yaxis.grid(True,which='major')
                  pb.ylabel('Phase (deg)', size=mysize)
                  pb.xlabel('Frequency (GHz) (%d channels)'%(len(frequenciesGHz[index])), size=mysize)
                  if (plotrange[0] != 0 or plotrange[1] != 0):
                      SetNewXLimits([plotrange[0],plotrange[1]],loc=1)
                  if (plotrange[2] != 0 or plotrange[3] != 0):
                      SetNewYLimits([plotrange[2],plotrange[3]],loc=1)
                  if (amplitudeWithPhase and phase != ''):
                      if (phase[0] != 0 or phase[1] != 0):
                          SetNewYLimits(phase,loc=2)
                  if (xframe == firstFrame):
                      pb.text(xstartTitle, ystartTitle,
                              '%s (degamp=%d, degphase=%d)'%(caltableTitle,
                              nPolyAmp[index]-1,nPolyPhase[index]-1),
                              size=mysize, transform=pb.gcf().transFigure)
                  # draw polarization labels (BPOLY case)
                  x0 = xstartPolLabel
                  y0 = ystartPolLabel
                  for p in range(nPolarizations):
                      if (corrTypeToString(corr_type[p]) in polsToPlot):
                          pb.text(x0, y0-0.03*p*subplotRows, corrTypeToString(corr_type[p])+'',
                                  color=pcolor[p],size=mysize, transform=pb.gca().transAxes)
                      
         # end of 'for' loop over rows
#         Probable fix for CAS-7753, not installed because bpoly not used much
#         if (fieldIndex == -1 and newpage==0):
#             newpage = 1
#             pages = pages[:len(pages)-1]
         redisplay = False
         pb.subplots_adjust(hspace=myhspace, wspace=mywspace)
         if (xframe == lastFrame):
            if (len(figfile) > 0):
                plotfiles.append(makeplot(figfile,msFound,msAnt,
                                  overlayAntennas,pages,pagectr,
                                  density,interactive,antennasToPlot,
                                  spwsToPlot,overlayTimes,overlayBasebands,
                                          1,xant,ispw,subplot,resample,debug,
                                  figfileSequential,figfileNumber))
                figfileNumber += 1
            donetime = timeUtilities.time()
            if (interactive):
               pb.draw()
#               myinput = raw_input("(%.1f sec) Press return for next page (b for backwards, q to quit): "%(donetime-mytimestamp))
               myinput = raw_input("Press return for next page (b for backwards, q to quit): ")
            else:
               myinput = ''
            skippingSpwMessageSent = 0
            mytimestamp = timeUtilities.time()
            if (myinput.find('q') >= 0):
#                showFinalMessage(overlayAntennas, solutionTimeSpread, nUniqueTimes)
                return(vm)
            if (myinput.find('b') >= 0):
                if (pagectr > 0):
                    pagectr -= 1
                #redisplay the current page by setting ctrs back to the value they had at start of that page
                xctr = pages[pagectr][PAGE_ANT]
                spwctr = pages[pagectr][PAGE_SPW]
                mytime = pages[pagectr][PAGE_TIME]
                myap = pages[pagectr][PAGE_AP]
                xant = antennasToPlot[xctr]
                antstring = buildAntString(xant,msFound,msAnt)
                ispw = spwsToPlot[spwctr]
                redisplay = True
            else:
                pagectr += 1
                if (pagectr >= len(pages)):
                    newpage = 1
                else:
                    newpage = 0
            if (debug):
                print("1)Setting xframe to %d" % xframeStart)
            xframe = xframeStart
            if (xctr+1 < len(antennasToPlot)):
                # don't clear the final plot when finished
                pb.clf()
            if (spwctr+1<len(spwsToPlot) or mytime+1<nUniqueTimes):
                # don't clear the final plot when finished
                pb.clf()
            pb.subplots_adjust(hspace=myhspace, wspace=mywspace)
         if (redisplay == False):
             mytime += 1
             if (debug):
                 print("Incrementing mytime to %d" % (mytime))
       # end while(mytime)
       if (redisplay == False):
           spwctr +=1
           if (debug):
               print("Incrementing spwctr to %d" % (spwctr))
      # end while(spwctr)
      if (redisplay == False):
          xctr += 1
    # end while(xctr) for BPOLY
    if (len(figfile) > 0 and pagectr<len(pages)):
       plotfiles.append(makeplot(figfile,msFound,msAnt,overlayAntennas,pages,
                                 pagectr,density,interactive,antennasToPlot,
                                 spwsToPlot,overlayTimes,overlayBasebands,
                                 2,xant,ispw,subplot,resample,debug,
                                 figfileSequential,figfileNumber))
       figfileNumber += 1
    if (len(plotfiles) > 0 and buildpdf):
        pdfname = figfile+'.pdf'
        filelist = ''
        plotfiles = np.unique(plotfiles)
        for i in range(len(plotfiles)):
            cmd = '%s -density %d %s %s.pdf' % (convert,density,plotfiles[i],plotfiles[i].split('.png')[0])
            print("Running command = %s" % (cmd))
            mystatus = os.system(cmd)
            if (mystatus != 0):
                break
            if (cleanup):
                os.system('rm -f %s' % (plotfiles[i]))
            filelist += plotfiles[i].split('.png')[0] + '.pdf '
        if (mystatus != 0):
            print("ImageMagick is missing, no PDF created")
            buildpdf = False
        if (buildpdf):
            # The following 2 lines reduce the total number of characters on the command line, which
            # was apparently a problem at JAO for Liza.
            filelist = ' '.join(au.pruneFilelist(filelist.split()))
            pdfname = au.pruneFilelist([pdfname])[0]
            cmd = '%s %s cat output %s' % (pdftk, filelist, pdfname)
            print("Running command = %s" % (cmd))
            mystatus = os.system(cmd)
            if (mystatus != 0):
                cmd = '%s -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s' % (gs,pdfname,filelist)
                print("Running command = %s" % (cmd))
                mystatus = os.system(cmd)
            if (mystatus == 0):
                print("PDF left in %s" % (pdfname))
                os.system("rm -f %s" % filelist)
            else:
                print("Both pdftk and ghostscript is missing, no PDF created")
    else:
        print("Not building PDF (buildpdf=%s, len(plotfiles)=%d)" % (str(buildpdf), len(plotfiles)))
    return(vm)

#################
# bpoly == false
#################
  msFound = False
  mytb.open(caltable)
  uniqueScanNumbers = sorted(np.unique(mytb.getcol('SCAN_NUMBER')))
  if (ParType == 'Complex'):  # casa >= 3.4
      gain = {}
      for f in range(len(fields)):
          if not mytb.iscelldefined('CPARAM',f): break
          gain[f] = mytb.getcell('CPARAM',f)
  else: # casa 3.3
      gain = {}
#      gain = mytb.getcol('FPARAM')       # e.g. (2, 128, 576)
      if ('FPARAM' in mytb.colnames()):
          for f in range(len(fields)):
              if not mytb.iscelldefined('FPARAM',f): break
              gain[f] = mytb.getcell('FPARAM',f)
      else:
          for f in range(len(fields)):
              gain[f] = mytb.getcell('GAIN',f)
  nPolarizations =  len(gain[0])
  if (debug):
      print("(1)Set nPolarizations = %d" % nPolarizations)
  ggx = {}
  for g in range(len(gain)):
      ggx[g] =  gain[g][0]
  if (nPolarizations == 2):
      ggy = {}
      for g in range(len(gain)):
          ggy[g] =  gain[g][1]
  mytb.close()

  if (debug):
      print("nPolarizations = ", nPolarizations)
  nRows = len(gain)
  if (bOverlay):
        mytb.open(caltable2)
        gain2 = {}
        if (ParType == 'Complex'):
#            gain2 = mytb.getcol('CPARAM')  # will bomb if shape varies
            if (debug): print("ParType = Complex")
            for f in range(len(fields2)):
#                if (debug): print "getting row %d/%d" % (f,len(fields))
                gain2[f] = mytb.getcell('CPARAM',f)
        else:
#            gain2 = mytb.getcol('FPARAM')  # will bomb if shape varies
            if (debug): print("ParType != Complex, it is %s" % (ParType))
            for f in range(len(fields2)):
                if (tableFormat2 == 34):
                    if mytb.iscelldefined('FPARAM',f):
                        gain2[f] = mytb.getcell('FPARAM',f)
                    else:
                        print("No data in row %d of 'FPARAM'" % (f))
                else:
                    if mytb.iscelldefined('GAIN',f):
                        gain2[f] = mytb.getcell('GAIN',f)
                    else:
                        print("No data in row %d of 'GAIN'" % (f))
        mytb.close()
        ggx2 = {}
        for g in range(len(gain2)):
#            print "Appending to ggx: ", gain2[g][0]
            ggx2[g] = gain2[g][0]
        nPolarizations2 = len(gain2[0])
        if (nPolarizations == 2):
            ggy2 = {}
            for g in range(len(gain2)):
                if (debug): print("len(gain2[%d]) = %d" % (g,len(gain2[g])))
                ggy2[g] = gain2[g][1]
        nRows2 = len(gain2)
        if (debug): print("nRows2 = ", nRows2)

  if (tableFormat == 34):
      # CAS-6801, unfortunately corr_type is not available in the caltable
      mytb.open(caltable)
      spectralWindowTable = mytb.getkeyword('SPECTRAL_WINDOW').split()[1]
      if ('OBSERVATION' in mytb.getkeywords()):
          observationTable = mytb.getkeyword('OBSERVATION').split()[1]
      else:
          observationTable = None
      mytb.open(spectralWindowTable)
      refFreq = mytb.getcol('REF_FREQUENCY')
      net_sideband = mytb.getcol('NET_SIDEBAND')
      measFreqRef = mytb.getcol('MEAS_FREQ_REF')
      mytb.close()
      corr_type = None
      if (os.path.exists(msName)):
          try:
              (corr_type, corr_type_string, nPolarizations) = getCorrType(msName,originalSpwsToPlot,debug)
          except:
              print("4) Could not open the associated measurement set tables (%s)." % (msName))
      if (corr_type is None):
          if (observationTable is None):
              corr_type, corr_type_string, nPolarizations = getCorrTypeByAntennaName(msAnt[0].lower())
          else:
              telescope = getTelescopeNameFromCaltableObservationTable(observationTable)
              if (telescope.find('ALMA') >= 0):
                  print("Using telescope name (%s) to set the polarization type." % (telescope))
                  corr_type_string = ['XX','YY']
                  corr_type = [9,12]
              elif (telescope.find('VLA') >= 0):
                  print("Using telescope name (%s) to set the polarization type." % (telescope))
                  corr_type_string = ['RR','LL']
                  corr_type = [5,8]
              else:
                  corr_type, corr_type_string, noPolarizations = getCorrTypeByAntennaName(msAnt[0].lower())
  else:
      try:
          if (DEBUG):
              print("Trying to open %s" % (msName+'/SPECTRAL_WINDOW'))
          mytb.open(msName+'/SPECTRAL_WINDOW')
          refFreq = mytb.getcol('REF_FREQUENCY')
          net_sideband = mytb.getcol('NET_SIDEBAND')
          measFreqRef = mytb.getcol('MEAS_FREQ_REF')
          mytb.close()
          if (debug): print("Calling getCorrtype('%s',%s,%s)" % (msName,str(originalSpwsToPlot),str(debug)))
          (corr_type, corr_type_string, nPolarizations) = getCorrType(msName,originalSpwsToPlot,debug)
          if (corr_type_string == []):
              return(vm)
      except:
          print("4) Could not open the associated measurement set tables (%s). Will not translate antenna names." % (msName))
          print("I will assume ALMA data: XX, YY, and refFreq=first channel.")
          #      chanFreqGHz = []  # comment out on 2014-04-08
          corr_type_string = ['XX','YY']
          corr_type = [9,12]

  if (len(polsToPlot) > len(corr_type)):
      # Necessary for SMA (single-pol) data
      polsToPlot = corr_type_string
  print("Polarizations to plot = ", polsToPlot)
  polsToPlot = checkPolsToPlot(polsToPlot, corr_type_string)
  if (polsToPlot == []):
      return(vm)

  if (len(msAnt) > 0):
      msFound = True
      if ((showatm == True or showtsky==True) and (vm=='' and mymsmd=='')):
          print("Because I could not open the .ms, I will use default weather for showatm.")
#          showatm = False
#          showtsky = False
  else:
      if (xaxis.find('freq')>=0 and tableFormat==33):
          print("Because I could not open the .ms and this is an old caltable, you cannot use xaxis='freq'.")
          return(vm)
      if (showatm == True or showtsky==True):
          print("Because I could not open the .ms, you cannot use showatm or showtsky.")
          return(vm)
  
  if (bpoly == False):
      if (debug):
          print("nPolarizations = ", nPolarizations)
          print("nFields = %d = " % (nFields), uniqueFields)

  if (bOverlay and debug):
        print("nPolarizations2 = ", nPolarizations2)
        print("nFields2 = %d = " % (nFields2), uniqueFields2)
        print("nRows2 = ", nRows2)
  uniqueAntennaIds = np.sort(np.unique(ant))

  yPhaseLabel = 'Phase (deg)'
  tsysPercent = True
  ampPercent = True
  if (VisCal.lower().find('tsys') >= 0):
      if (channeldiff > 0):
          if (tsysPercent):
              yAmplitudeLabel = "Tsys derivative (%_of_median/channel)"
          else:
              yAmplitudeLabel = "Tsys derivative (K/channel)"
      else:
          yAmplitudeLabel = "Tsys (K)"
  elif delay:
      yAmplitudeLabel = 'Delay (nsec)'
  else:
      if (yaxis.lower().find('db')>=0):
          yAmplitudeLabel = "Amplitude (dB)"
      else:
          if (channeldiff > 0):
              if (ampPercent):
                  yAmplitudeLabel = "Amp derivative (%_of_median/channel)"
              else:
                  yAmplitudeLabel = "Amplitude derivative"
              yPhaseLabel = 'Phase derivative (deg/channel)'
          else:
              yAmplitudeLabel = "Amplitude"

  madsigma = channeldiff # for option channeldiff>0, sets threshold for finding outliers
  ampMin = LARGE_POSITIVE
  ampMax = LARGE_NEGATIVE
  PHASE_ABS_SUM_THRESHOLD = 2e-3  # in degrees, used to avoid printing MAD statistics for refant

  pb.clf()
  drewAtmosphere = False
  TDMisSecond = False
  pagectr = 0
  newpage = 1
  pages =  []
  xctr = 0
  myap = 0  # determines whether an amp or phase plot starts the page (in the case of 'both')
            # zero means amplitude, 1 means phase
  redisplay = False
  matchctr = 0
  myUniqueColor = []
  # for the overlay=antenna case, start by assuming the first antenna is not flagged
  firstUnflaggedAntennaToPlot = 0
  lastUnflaggedAntennaToPlot = len(antennasToPlot)-1
  computedAtmSpw = -1
  computedAtmTime = -1
  computedAtmField = -1
  skippingSpwMessageSent = 0
  atmString = ''
  if (showimage and lo1 is None):
      # We only need to run this once per execution.
      getLOsReturnValue = au.getLOs(msName, verbose=False)
      if (getLOsReturnValue != []):
          lo1s = au.interpretLOs(msName,parentms,spwsForIntent=spwsToPlot, mymsmd=mymsmd, verbose=debug)
          foundLO1Message = []  # Initialize so that message is only displayed once per spw
      if debug:
          print("getLOsReturnValue = ", getLOsReturnValue)

  if (channeldiff>0):
      # build blank dictionary:  madstats['DV01']['spw']['time']['pol']['amp' or 'phase' or both]
      #                          where spw, time, pol are each integers
      if (len(msAnt) > 0):
          madstats = dict.fromkeys(msAnt)
      else:
          madstats = dict.fromkeys(['Ant '+str(i) for i in range(len(uniqueAntennaIds))])
#      print "msAnt = ", msAnt
#      print "spwsToPlot = ", spwsToPlot
      for i in range(len(madstats)):
          madstats[madstats.keys()[i]] = dict.fromkeys(spwsToPlot)
          for j in range(len(spwsToPlot)):
              madstats[madstats.keys()[i]][spwsToPlot[j]] = dict.fromkeys(timerangeList) # dict.fromkeys(range(len(uniqueTimes)))
              for k in timerangeList: # range(len(uniqueTimes)):
                  madstats[madstats.keys()[i]][spwsToPlot[j]][k] = dict.fromkeys(range(nPolarizations))
                  for l in range(nPolarizations):
                      if (yaxis == 'both'):
                          madstats[madstats.keys()[i]][spwsToPlot[j]][k][l] = {'amp': None, 'phase': None, 'ampstd': None, 'phasestd': None}
                      elif (yaxis == 'phase'):
                          madstats[madstats.keys()[i]][spwsToPlot[j]][k][l] = {'phase': None, 'phasestd': None}
                      else:
                          # this includes tsys and amp
                          madstats[madstats.keys()[i]][spwsToPlot[j]][k][l] = {'amp': None, 'ampstd': None}
      madstats['platforming'] = {}
#      print "madstats = ", madstats
  myinput = ''
  atmEverBeenCalculated = False
  spwsToPlotInBaseband = []
  frequencyRangeToPlotInBaseband = []
  if (basebands == []):
      # MS is too old to have BBC_NO
      spwsToPlotInBaseband = [spwsToPlot]
      frequencyRangeToPlotInBaseband = [callFrequencyRangeForSpws(mymsmd, spwsToPlot, vm, debug, caltable)]
      basebands = [0]
  elif (overlayBasebands):
      if (list(spwsToPlot) != list(uniqueSpwsInCalTable)):
          # then spws were requested, so treat them all as if in the same baseband, and
          # ignore the basebands parameter
          print("Ignoring the basebands parameter because spws were specified = %s" % (str(spwsToPlot)))
      elif (np.array_equal(np.sort(basebands), np.sort(allBasebands)) == False):
          # Allow the basebands parameter to select the spws
          basebandSpwsToPlot = []
          for baseband in basebands:
              myspws = list(getSpwsForBaseband(vis=msName, mymsmd=mymsmd, bb=baseband, caltable=caltable))
              basebandSpwsToPlot += myspws
          spwsToPlot = np.intersect1d(basebandSpwsToPlot, spwsToPlot)
          print("selected basebands %s have spwsToPlot = %s" % (str(basebands),str(spwsToPlot)))
      else:
          if (debug): print("No spws or basebands selected, will overlay all spws from all basebands.")
      spwsToPlotInBaseband = [spwsToPlot]
      frequencyRangeToPlotInBaseband = [callFrequencyRangeForSpws(mymsmd, spwsToPlot, vm, debug, caltable)]
      basebands = [0]
  else:
      for baseband in basebands:
          myspwlist = []
          for spw in spwsToPlot:
              if (casaVersion >= '4.1.0' and msFound):
                  if (mymsmd.baseband(originalSpwsToPlot[list(spwsToPlot).index(spw)]) == baseband):
                      myspwlist.append(spw)
              else:
                  # need to write a function to retrieve baseband (if I ever run casa 4.0 again)
                  # if (spw != 0): 
                  myspwlist.append(spw)
          spwsToPlotInBaseband.append(myspwlist)
          frequencyRangeToPlotInBaseband.append(callFrequencyRangeForSpws(mymsmd, myspwlist,vm, debug, caltable))
      print("basebands to plot = %s" % (str(basebands)))

  firstTimeMatch = -1    # Aug 5, 2013
  if (overlaySpws or overlayBasebands):
      groupByBaseband = True
  if (groupByBaseband and overlaySpws==False and overlayBasebands==False):
      showBasebandNumber = True
  # Basic nested 'while' loop structure is:  
  #   - antennas
  #     - baseband (if necessary)
  #       - spw
  #         - time
  #           - for i in rows
  while (xctr < len(antennasToPlot)):
    if (debug): print("at top of xctr loop: %d" % (xctr))
    xant = antennasToPlot[xctr]
    bbctr = 0
    spwctr = 0
    spwctrFirstToPlot = 0
    antstring = buildAntString(xant,msFound,msAnt)
    alreadyPlottedAmp = False  # needed for (overlay='baseband', yaxis='both')  CAS-6477
    if (antstring.isdigit()):
        Antstring = "Ant%s" % antstring
    else:
        Antstring = antstring
    finalSpwWasFlagged = False   # inserted on 22-Apr-2014 for g25.27
    while ((bbctr < len(spwsToPlotInBaseband) and groupByBaseband) or
           (spwctr < len(spwsToPlot) and groupByBaseband==False)
           ):
     if (debug): print("at top of bbctr/spwctr loop with bbctr=%d, spwctr=%d" % (bbctr,spwctr))
     if (groupByBaseband):
         baseband = basebands[bbctr]
         spwsToPlot = spwsToPlotInBaseband[bbctr]
         if (debug): print("setting spwsToPlot for baseband %d (bbctr=%d) to %s" % (baseband, bbctr, str(spwsToPlot)))
     else:
         baseband = 0
         if (casaVersion >= '4.1.0'):
             if (getBasebandDict(msName,caltable=caltable, mymsmd=mymsmd) != {}):
                 try:
                     baseband = mymsmd.baseband(originalSpwsToPlot[spwctr])
                     if (baseband not in basebands):
                         spwctr += 1  # was bbctr += 1
                         if (debug): print("A)incrementing spwctr")
#                     baseband = mymsmd.baseband(spwsToPlot[spwctr])
#                     if (baseband not in basebands and bbctr+1 < len(spwsToPlotInBaseband)):
#                         bbctr += 1
                         continue
                 except:
                     pass
     if (debug):
         if (overlayBasebands):
             print("Regardless of baseband (%s), plotting all spws: %s" % (basebands,str(spwsToPlot)))
         else:
             print("Showing baseband %d containing spws:" % (baseband))
             print(str(spwsToPlotInBaseband[bbctr]))
     if (bbctr < len(spwsToPlotInBaseband)):
         if (debug):
             print("A) spwctr=%d,  bbctr=%d < len(spwsToPlotInBaseband)=%d" % (spwctr,bbctr,len(spwsToPlotInBaseband)))
         spwctr = 0
         spwctrFirstToPlot = spwctr
     firstSpwMatch = -1
     while (spwctr < len(spwsToPlot)):
      if (debug): print("at top of spwctr loop, spwctr=%d" % (spwctr))
      allTimesFlaggedOnThisSpw = True # used only by overlay='time'
      if (groupByBaseband == False):
          baseband = 0
          if (casaVersion >= '4.1.0'):
              if (getBasebandDict(msName,caltable=caltable,mymsmd=mymsmd) != {}):
                  try:
                      baseband = mymsmd.baseband(originalSpwsToPlot[spwctr])
                      if (baseband not in basebands):
#                          print "spw %d=%d: baseband %d is not in %s" % (spwsToPlot[spwctr],originalSpwsToPlot[spwctr], baseband, basebands)
                          if (debug): print("B)incrementing spwctr to %d" % (spwctr+1))
                          spwctr += 1
                          continue
                  except:
                      pass
      ispw = spwsToPlot[spwctr]
      ispwInCalTable = list(uniqueSpwsInCalTable).index(ispw)
      mytime = 0
      if (debug):
          print("+++++++ set mytime=0 for ispw=%d, len(chanFreqGHz) = %d" % (ispw, len(chanFreqGHz)))
      if (overlayAntennas):
          xctr = -1
      if (overlayTimes):
          # since the times/scans can vary between spws, redefine nUniqueTimes for each spw
          nUniqueTimes = len(uniqueTimesCopy)
          uniqueTimes = uniqueTimesCopy[:]
          uniqueTimesForSpw = []
          testtime = 0
          while (testtime < nUniqueTimes):
              if (ispw in cal_desc_id[np.where(uniqueTimes[testtime] == times)[0]]):
                  uniqueTimesForSpw.append(uniqueTimes[testtime])
              testtime += 1
          uniqueTimes = uniqueTimesForSpw[:]
          if (tableFormat >= 34):
              scansForUniqueTimes, nUniqueTimes = computeScansForUniqueTimes(uniqueTimes, cal_scans, times, unique_cal_scans)
          else:
              nUniqueTimes = len(uniqueTimes)
#      currentSpwctr = spwctr   # commented out on 2014-04-04 to match task for task01 regression
      if (overlaySpws or overlayBasebands):
          if (xctr >= firstUnflaggedAntennaToPlot):
              if (debug):
                  print("xctr=%d >= firstUnflaggedAntennaToPlot=%d" % (xctr, firstUnflaggedAntennaToPlot))
              spwctr -= 1
      
      firstTimeMatch = -1
      finalTimeMatch = -1 #  for CAS-7820
      while (mytime < nUniqueTimes):
        finalTimerangeFlagged = False  # 04-Aug-2014
        if (debug):
            print("at top of mytime loop: mytime = %d < %d" % (mytime,nUniqueTimes))
            print("timerangeList = %s" % (str(timerangeList)))
            print("timerangeListTimes = %s" % (str(timerangeListTimes)))
            print("debugSloppyMatch = %s" % (str(debugSloppyMatch)))
            print("solutionTimeThresholdSeconds = %s" % (str(solutionTimeThresholdSeconds)))
#        if ((scansToPlot == scansToPlotPerSpw[ispw]).all() == False and False):
#            print "          scansToPlot = ", scansToPlot
#            print "scansToPlotPerSpw[%2d] = " % (ispw), scansToPlotPerSpw[ispw]
        if (len(timerangeList) > 0 and
            (sloppyMatch(uniqueTimes[mytime],timerangeListTimes,solutionTimeThresholdSeconds,
                         mytime, scansToPlot, scansForUniqueTimes, myprint=debugSloppyMatch)==False)):  # task version 
#                         mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes, myprint=debugSloppyMatch)==False)): # au version (with this: test 85 has infinite loop)
            if (debug):
                print("Skipping mytime %d because it is not in %s, or scan %d is not in %s" % (mytime, str(timerangeList),scansForUniqueTimes[mytime],scansToPlotPerSpw[ispw]))
            mytime += 1
            if (debug): print("  0006  incrementing mytime to ", mytime)
            if (mytime == nUniqueTimes and overlayTimes and overlayAntennas):
                # added March 14, 2013 to support the case when final timerange is flagged
                doneOverlayTime = False
                if (debug):
                    print("$$$$$$$$$$$$$$$$$$$$$$  Setting doneOverlayTime=False" % (xframe))
            else:
                if (debug):
                    print("Not setting doneOverlayTime=False because either mytime(%d) != nUniqueTimes(%d) or we are not overlaying time and antenna" % (mytime,nUniqueTimes))
            continue
        if (overlayAntennas):
            xctr += 1
            if (xctr >= len(antennasToPlot)):
                xctr = 0
#                if (mytime == 0):
#                    if (debug): print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ setting firstTimeMatch = -1"
#                    firstTimeMatch = -1  # Aug 5, 2013
            xant = antennasToPlot[xctr]
            antstring = buildAntString(xant,msFound,msAnt)
            if (debug):
                print("mytime=%d, Set xant to %d" % (mytime,xant))
            antennaString = ''
        else:
            antennaString = 'Ant%2d: %s,  ' % (xant,antstring)
        if (overlaySpws or overlayBasebands):
            if (debug): print("C)incrementing spwctr to %d" % (spwctr+1))
            spwctr += 1
            if (spwctr >= len(spwsToPlot)):
                spwctr = 0  # added on 2014-04-04 to match the task
                if (xctr < firstUnflaggedAntennaToPlot):
                    xctr += 1
                    if (xctr == len(antennasToPlot)): 
                        break
                    xant = antennasToPlot[xctr]
                    antstring = buildAntString(xant,msFound,msAnt)
                    if (debug):
                        print("mytime=%d, Set xant to %d" % (mytime,xant))
                    antennaString = 'Ant%2d: %s,  ' % (xant,antstring)
#                if (overlaySpws):               # commented out on 2014-04-04 to match task for task01 regression
#                    spwctr = currentSpwctr      # commented out on 2014-04-04 to match task for task01 regression
#                break  # added on Sep 23, 2013 so that mytime gets incremented # commented out on 2014-04-04 to match task for task01 regression
                if (overlayBasebands):
                    # Added on 7/29/2014 to fix infinite loop in uid___A002_X652932_X20fb bandpass
                    if (mytime == nUniqueTimes):   # +0 fixes  regression 1
                        spwctr = len(spwsToPlot)
                        if (debug): 
                            print("Breaking because spwctr=%d == len(spwsToPlot)=%d" % (spwctr,len(spwsToPlot)))
                        break
            ispw = spwsToPlot[spwctr]
            if (ispw not in uniqueSpwsInCalTable):
                print("spw %d is not in caltable=%s" % (ispw,uniqueSpwsInCalTable))
                return
            ispwInCalTable = list(uniqueSpwsInCalTable).index(ispw)
            if (debug):
                print("----------------------------- spwctr=%d, ispw set to %d, xctr=%d" % (spwctr,ispw,xctr))

        if (newpage==1):
            # add the current page (being created here) to the list
            pages.append([xctr,spwctr,mytime,0])
            if (debug):
                print("top: appending [%d,%d,%d,%d]" % (xctr,spwctr,mytime,0))
            newpage = 0
        gplotx = []
        gploty = []
        channels = []
        xchannels = []
        ychannels = []
        frequencies = []
        xfrequencies = []
        yfrequencies = []
        channels2 = []
        xchannels2 = []
        ychannels2 = []
        frequencies2 = []
        xfrequencies2 = []
        yfrequencies2 = []
        gplotx2 = []
        gploty2 = []
        xflag = []
        yflag = []
        xflag2 = []
        yflag2 = []
        matchFound = False
        matchField = -1
        matchRow = -1
        matchTime = -1
        if (debug): print("looping over all nRows = %d" % (nRows))
        for i in range(nRows):
            if (overlayTimes or overlayAntennas or len(fieldsToPlot)>1 or
                (nFields>1 and len(fieldlist)<nFields)):
                # When there are multiple fields, then matching by scan causes the first
                # matching solution to be displayed every time.  So use the original method
                # of matching by time until I think of something better.
                sm = sloppyMatch(uniqueTimes[mytime],times[i],solutionTimeThresholdSeconds,myprint=False)
            else:
                if (overlayBasebands):
                    sTP = scansToPlot
                else:
                    sTP = scansToPlotPerSpw[ispw]
                sm = sloppyMatch(uniqueTimes[mytime],times[i],solutionTimeThresholdSeconds,
#                                 mytime, scansToPlot, scansForUniqueTimes, myprint=False) # task version
                                 mytime, sTP, scansForUniqueTimes, myprint=False) # au version
                
            if ((ant[i]==xant) and (cal_desc_id[i]==ispw) and sm
                and (mytime in timerangeList)   # this test was added to support multiFieldInTimeOverlay
                ):
                if debug: print("len(chanFreqGHz)=%d, ispw=%d" % (len(chanFreqGHz),ispw))
                if (msFound or tableFormat==34):
                    if (len(chanFreqGHz[ispw]) == 1 and delay==False):
                        if ((skippingSpwMessageSent & (1<<ispw)) == 0):
                            print("Skipping spw=%d because it has only 1 channel." % (ispw))
                            skippingSpwMessageSent |= (1<<ispw)
                        break
                if (fields[i] in fieldsToPlot):
                    interval = intervals[i] # used for CalcAtmTransmission
                    myFieldIndex = np.where(fields[i] == uniqueFields)[0]
                    if (type(myFieldIndex) == list or type(myFieldIndex) == np.ndarray):
                        myFieldIndex = myFieldIndex[0]
                    if (debug):
                        print("%d Found match at field,ant,spw,mytime,time = %d(index=%d),%d,%d,%d,%f=%s" % (matchctr,fields[i],myFieldIndex,xant,ispw,mytime,uniqueTimes[mytime],utstring(uniqueTimes[mytime],4)))
                    if (matchFound):
                        if (myFieldIndex == matchField and matchTime==times[i]):
                            print("WARNING: multiple rows for field=%d,ant=%d,spw=%d,scan=%d,time=%d=%.0f=%s,row=%d. Only showing the first one." % (fields[i],xant,ispw,scansForUniqueTimes[mytime],mytime,uniqueTimes[mytime],utstring(uniqueTimes[mytime],3),i))
                    else:
                        matchFound = True
                        fieldIndex = myFieldIndex
                        matchField = myFieldIndex
                        matchTime = times[i]
                        matchRow = i
                        if (msFound or tableFormat==34):
                            nChannels = len(chanFreqGHz[ispw])
                            solnChannels = len(ggx[i])
                            if (nChannels != solnChannels):
                                print("Mismatch in table between number of solution channels (%d) and spw channels (%d) for spw %d." % (solnChannels,nChannels,ispw))
                                return
                        else:
                            nChannels = len(ggx[0])
                        xflag.append(flags[i][0][:])
                        yflag.append(flags[i][1][:])
                        BRowNumber = i
                        for j in range(nChannels):   # len(chanFreqGHz[ispw])):
                            channels.append(j)  # both flagged and unflagged
                            if (msFound or tableFormat==34):
                                frequencies.append(chanFreqGHz[ispw][j])
                                if (j==0 and debug):
                                    print("found match: ispw=%d, j=%d, len(chanFreqGHz)=%d, chanFreqGHz[0]=%f" % (ispw,j, len(chanFreqGHz),chanFreqGHz[ispw][0]))
                            if (showflagged or (showflagged == False and flags[i][0][j]==0)):
                                gplotx.append(ggx[i][j])
                                xchannels.append(j)
                                if (msFound or tableFormat==34):
                                    xfrequencies.append(chanFreqGHz[ispw][j])
                            if (nPolarizations == 2):
                                if (showflagged or (showflagged == False and flags[i][1][j]==0)):
                                    gploty.append(ggy[i][j])
                                    ychannels.append(j)
                                    if (msFound or tableFormat==34):
                                        yfrequencies.append(chanFreqGHz[ispw][j])
            elif (debug and False):
                if (mytime not in timerangeList):
                    print("mytime %f not in %s" % (mytime,str(timerangeList)))
                if (ant[i]!=xant):
                    print("ant: %d != %d" % (ant[i],xant))
                if (cal_desc_id[i]!=ispw):
                    print("cal_desc_id[i]!=ispw  %d != %d" % (cal_desc_id[i],ispw))
                if (sm == False):
                    print("sm failed")
        # end 'for i in range(nRows)'
        # So, fieldIndex needs to be the correct one when we exit this loop
        if (debug):
            print("finished loop over nRows: matchFound = ", matchFound)
        if (not matchFound and newpage==0):
          if (subplot==11 or (subplot!=11 and firstTimeMatch==-1 and firstSpwMatch==-1)):
            # Fix for CAS-7753   (subplot==11)
            #            for subplot=22,  firstTimeMatch and firstSpwMatch gives spw17 and 23 (or gives 17 and 39)
            # the firstTimeMatch part was needed for regression 65: different antennas having different solution times
            # the firstSpwMatch was needed to keep the first spw in the filename rather than the final spw, e.g. regression 14
#            print "Applying fix for CAS-7753"
            newpage = 1
            pages = pages[:len(pages)-1]
        myspw = originalSpw[ispw]
        if (msFound):
            if debug: 
                print("myspw=", myspw)
                print("len(refFreq)=", len(refFreq))
            if (myspw >= len(refFreq)):
                myspw = ispw
        if (msFound and refFreq[myspw]*1e-9 > 60): 
          # Then this cannot be EVLA data.  But I should really check the telescope name!
          # Then again, Band 1 may not have USB/LSB distinction (nor Band 2 for that matter).
          if debug and False: 
              print("frequencies = ", frequencies)
              print("xfrequencies = ", xfrequencies)
              print("chanFreqGHz[%d] = " % (ispw), chanFreqGHz[ispw])
#          if (refFreq[myspw]*1e-9 > np.mean(frequencies)):
          if (refFreq[myspw]*1e-9 > np.mean(chanFreqGHz[ispw])):  # this is safer (since frequencies might be [])
              sideband = -1
              xlabelString = "%s LSB Frequency (GHz) (%d channels)" % (refTypeToString(measFreqRef[myspw]),len(xchannels))
          else:
              sideband = +1
              xlabelString = "%s USB Frequency (GHz) (%d channels)" % (refTypeToString(measFreqRef[myspw]),len(xchannels))
        else:
            sideband = -1
            xlabelString = "Frequency (GHz) (%d channels)" % (len(xchannels))
        if ((len(frequencies)>0) and (chanrange[1] > len(frequencies))):
            print("Invalid chanrange (%d-%d) for spw%d in caltable1. Valid range = 0-%d" % (chanrange[0],chanrange[1],ispw,len(frequencies)-1))
            return(vm)
        pchannels = [xchannels,ychannels]
        pfrequencies = [xfrequencies,yfrequencies]
        gplot = [gplotx,gploty]
        # We only need to compute the atmospheric transmission if:
        #   * we have been asked to show it,
        #   * there is a non-trivial number of channels, 
        #   * the current field is the one for which we should calculate it (if times are being overlaid)
        #       But this will cause no atmcurve to appear if that field is flagged on the first
        #       antenna; so, I added the atmEverBeenCalculated flag to deal with this.
        #   * the previous calculation is not identical to what this one will be
        #
        if ((showatm or showtsky) and (len(xchannels)>1 or len(ychannels)>1) and
                                                        # this insures a plot if first fieldsToPlot is missing
            ((uniqueFields[fieldIndex]==showatmfield or (uniqueFields[fieldIndex] in fieldsToPlot and overlayTimes)) or
              overlayTimes==False or atmEverBeenCalculated==False) and
            ((overlayTimes==False and computedAtmField!=fieldIndex) or (computedAtmSpw!=ispw) or
             (overlayTimes==False and computedAtmTime!=mytime))):
          atmEverBeenCalculated = True
#          print "Calcing Atm: len(xchannels)=%d, len(ychannels=%d), uniqueFields[%d]=%d ?= %d" % (len(xchannels), len(ychannels),fieldIndex,uniqueFields[fieldIndex],showatmfield)
#          print "             computedAtmField=%d, fieldIndex=%d;  computedAtmSpw=%d, ispw=%d;  computedAtmTime=%d,mytime=%d" % (computedAtmField, fieldIndex, computedAtmSpw, ispw, computedAtmTime,mytime)
          if (True):  # support showatm for overlay='antenna,time'
#            print "Running CalcAtm:  CAF=%d, CAS=%d, CAT=%d, xant=%d" % (computedAtmField, computedAtmSpw, computedAtmTime, xant)
            if (type(fieldIndex) == list or type(fieldIndex) == np.ndarray):
                computedAtmField = fieldIndex[0]
            else:
                computedAtmField = fieldIndex
            computedAtmSpw = ispw
            computedAtmTime = mytime
            atmtime = timeUtilities.time()
            asdm = ''
            (atmfreq,atmchan,transmission,pwvmean,atmairmass,TebbSky,missingCalWVRErrorPrinted,Tsys) = \
               CalcAtmTransmission(channels, frequencies, xaxis, pwv,
                                   vm, msName, asdm, xant, uniqueTimes[mytime],
                                   interval, uniqueFields[fieldIndex],
                                   refFreq[originalSpw[ispw]],
                                   net_sideband[originalSpw[ispw]], mytime, 
                                   missingCalWVRErrorPrinted, mymsmd, caltable,
                                   verbose=DEBUG, maxAtmCalcChannels=maxAtmCalcChannels,
                                   maxAltitude=maxAltitude, Feff=Feff, SBGain=SBGain, 
                                   Trx=Trx, showtsys=showtsys)
            if showtsys: 
                TebbSky = Tsys
            if (showimage):
                if (lo1 is not None):
                    LO1 = lo1
                else:
                  if (getLOsReturnValue == []):
                    if (lo1 is None):
                        print("Because you do not have the ASDM_RECEIVER table, if you want the showimage")
                        print("option to work, then you must specify the LO1 frequency with lo1=.")
#                        return(vm)
                    LO1 = lo1
                  else:
                    if (lo1s is None or lo1s == {}):
                        print("Failed to get LO1, disabling showimage.  Alternatively, you can use printLOsFromASDM and supply the lo1 parameter to plotbandpass.")
                        showimage = False
                        LO1 = None
                    else:
                        if (originalSpw[ispw] not in lo1s.keys()):
                            print("There is a problem in reading the LO1 values, cannot showimage for this dataset.")
                            print("originalSpw[%d]=%d > len(lo1s)=%d" % (ispw, originalSpw[ispw], len(lo1s)))
                            showimage = False
                            LO1 = None
                        else:
                            LO1 = lo1s[originalSpw[ispw]]*1e-9  
                            if (ispw not in foundLO1Message):
                                print("For spw %d (%d), found LO1 = %.6f GHz" % (ispw,originalSpw[ispw],LO1))
                                foundLO1Message.append(ispw)
            if (LO1 is not None):
                frequenciesImage = list(2*LO1 - np.array(frequencies))
                xfrequenciesImage = list(2*LO1 - np.array(pfrequencies[0]))
                yfrequenciesImage = list(2*LO1 - np.array(pfrequencies[1]))
                pfrequenciesImage = [xfrequenciesImage, yfrequenciesImage]
                (atmfreqImage,atmchanImage,transmissionImage,pwvmean,atmairmass,TebbSkyImage,missingCalWVRErrorPrinted,TsysImage) = \
                    CalcAtmTransmission(channels, frequenciesImage, xaxis,
                                        pwv, vm, msName, asdm, xant, uniqueTimes[mytime],
                                        interval, uniqueFields[fieldIndex],
                                        refFreq[originalSpw[ispw]],
                                        net_sideband[originalSpw[ispw]], mytime, 
                                        missingCalWVRErrorPrinted, mymsmd, 
                                        caltable, verbose=DEBUG, 
                                        maxAtmCalcChannels=maxAtmCalcChannels,
                                        maxAltitude=maxAltitude, Feff=Feff, SBGain=SBGain,
                                        Trx=Trx, showtsys=showtsys)
                if showtsys: 
                    TebbSkyImage = TsysImage
                atmfreqImage = list(2*LO1 - np.array(atmfreqImage))
#                atmfreqImage.reverse()  
                atmchanImage.reverse()

            if (overlayTimes):
                atmString = 'PWV %.2fmm, airmass %.2f, maxAlt %.0fkm (field %d)' % (pwvmean,atmairmass,maxAltitude,uniqueFields[fieldIndex])
            else:
                atmString = 'PWV %.2fmm, airmass %.3f, maxAlt %.0fkm' % (pwvmean,atmairmass,maxAltitude)
        elif (False):
            print("Skipping CalcAtm: xframe=%d, len(xchannels)=%d, len(ychannels=%d), uniqueFields[%d]=%d ?= %d" % (xframe,len(xchannels), len(ychannels),fieldIndex,uniqueFields[fieldIndex],showatmfield))
            print("                  computedAtmField=%d, fieldIndex=%d, computedAtmSpw=%d, ispw=%d" % (computedAtmField,fieldIndex,computedAtmSpw,ispw))
            print("                  computedAtmTime=%d, mytime=%d, xant=%d" % (computedAtmTime, mytime, xant))
            
        if (bOverlay):
          for i in range(nRows2):
            if (overlayTimes or overlayAntennas or len(fieldsToPlot)>1 or
                (nFields>1 and len(fieldlist)<nFields)):
                # Not having this path causes Tsys table overlays to behave like overlay='antenna,time' 
                # for caltable2.
                sm = sloppyMatch(uniqueTimes2[mytime],times2[i],solutionTimeThresholdSeconds,myprint=False)
            else:
                if debug:
                    if (mytime >= len(uniqueTimes2)):
                        print("mytime=%d is too long for uniqueTimes2=%d" % (mytime, len(uniqueTimes2)))
                    if (i >= len(times2)):
                        print("i=%d is too long for uniqueTimes2=%d" % (i, len(times2)))
                    if (ispw >= len(scansToPlotPerSpw)):
                        print("ispw=%d is too long for uniqueTimes2=%d" % (ispw, len(scansToPlotPerSpw)))
                if (mytime >= len(uniqueTimes2)):
                    # Fix for CAS-9474: avoid calling sloppyMatch because it will crash.  
                    # Setting sm=False will result in an abort: "no amp data found in second solution."
                    sm = False
                else:
                    sm = sloppyMatch(uniqueTimes2[mytime],times2[i],solutionTimeThresholdSeconds,
                                 mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes,  # au version
                                 myprint=False)
            if ((ant2[i]==xant) and (cal_desc_id2[i]==ispw) and sm
                and (mytime in timerangeList)   # added to match first caltable logic on 2014-04-09
                ):
                if (fields2[i] in fieldsToPlot):
                      xflag2.append(flags2[i][0][:])
                      yflag2.append(flags2[i][1][:])
                      # With solint='2ch' or more, the following loop should not be over
                      # chanFreqGHz2 but over the channels in the solution.
#                      print "len(chanFreqGHz2[%d])=%d" % (ispw,len(chanFreqGHz2[ispw]))
                      for j in range(len(chanFreqGHz2[ispw])):
                        channels2.append(j)
                        frequencies2.append(chanFreqGHz2[ispw][j])
#                        print "len(chanFreqGHz2[%d])=%d, i=%d,j=%d, len(ggx2)=%d, len(ggx2[0])=%d, shape(ggx2) = " % (ispw,len(chanFreqGHz2[ispw]),i,j,len(ggx2),len(ggx2[0])), np.shape(np.array(ggx2))
                        if (showflagged or (showflagged == False and flags2[i][0][j]==0)):
                            gplotx2.append(ggx2[i][j])
                            xchannels2.append(j)
                            xfrequencies2.append(chanFreqGHz2[ispw][j])
#                            print "appending %f to xfrequencies2" % (chanFreqGHz2[ispw][j])
                        else:
                            if (debug):
                                print("********* flags2[%d][0][%d] = %d, showflagged=" % (i,j,flags2[i][0][j]), showflagged)
                        if (nPolarizations2 == 2):
                            if (showflagged or (showflagged == False and flags2[i][1][j]==0)):
                                gploty2.append(ggy2[i][j])
                                ychannels2.append(j)
                                yfrequencies2.append(chanFreqGHz2[ispw][j])
          # end 'for i'
          pchannels2 = [xchannels2,ychannels2]
          pfrequencies2 = [xfrequencies2,yfrequencies2]
          gplot2 = [gplotx2,gploty2]
          # Need to rewrite the xlabel to show the total channel numbers from both caltables.  Note that xaxis must be 'freq' for bOverlay
          if (msFound and refFreq[myspw]*1e-9 > 60): 
              # Then this cannot be EVLA data.  But I should really check the telescope name!
              # Then again, Band 1 may not have USB/LSB distinction (nor Band 2 for that matter).
              if (refFreq[myspw]*1e-9 > np.mean(chanFreqGHz2[ispw])):  # this is safer (since frequencies might be [])
                  sideband = -1
                  xlabelString = "%s LSB Frequency (GHz) (%d, %d channels)" % (refTypeToString(measFreqRef[myspw]),len(xchannels),len(xchannels2))
              else:
                  sideband = +1
                  xlabelString = "%s USB Frequency (GHz) (%d, %d channels)" % (refTypeToString(measFreqRef[myspw]),len(xchannels),len(xchannels2))
          else:
              sideband = -1
              xlabelString = "Frequency (GHz) (%d, %d channels)" % (len(xchannels),len(xchannels2))

#
#  Don't check here, because chanrange refers only to caltable, not caltable2.
#          if (len(frequencies2)>0 and (chanrange[1] > len(frequencies2))):
#              print "Invalid chanrange for spw%d in caltable2. Valid range = 0-%d" % (ispw,len(frequencies2))
#              return(vm)
              
# Prevents crash if long value is set for solutionTimeThresholdSeconds, but prints a lot of
# messages for Tsys with overlay='antenna'.
#        if (len(xchannels) < 1):
#            print "No unflagged data found for (ant,spw,mytime,time) = %d,%d,%d,%.1f=%s" % (xant,ispw,mytime,uniqueTimes[mytime],utstring(uniqueTimes[mytime],4))
#            matchFound = False

        if (matchFound==False): 
            if ((overlayAntennas==False and overlaySpws==False and overlayBasebands==False) or
                (overlayAntennas and xctr+1 >= len(antennasToPlot)) or
                ((overlaySpws or overlayBasebands) and spwctr+1 >= len(spwsToPlot))):
                mytime += 1
#                print "  0005  incrementing mytime to ", mytime
                if (debug):
                    print("a) xctr=%d, Incrementing mytime to %d/%d because NO MATCH FOUND ============" % (xctr, mytime,nUniqueTimes))
            continue
        #  The following variable allows color legend of UT times to match line plot
        myUniqueTime = []
        if (True):   # multiFieldsWithOverlayTime
            # support multi-fields with overlay='time'
            uTPFPS = []
            for f in fieldIndicesToPlot:
                for t in uniqueTimesPerFieldPerSpw[ispwInCalTable][f]:
                    if (sloppyMatch(t, timerangeListTimes, solutionTimeThresholdSeconds,
                                    mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes,  # au version
#                                    mytime, scansToPlot, scansForUniqueTimes,  # task version
                                    myprint=debugSloppyMatch
                                    )):
                        uTPFPS.append(t)
            uTPFPS = np.sort(uTPFPS)
            ctr = 0
            for t in uTPFPS:
                if (debug and False):
                    print("1)checking time %d" % (t))
                if (overlayTimes or overlayAntennas):
                    sm = sloppyMatch(uniqueTimes[mytime],times[i],solutionTimeThresholdSeconds,myprint=False)
                else:
                    sm = sloppyMatch(t, uniqueTimes[mytime], solutionTimeThresholdSeconds,
#                                mytime, scansToPlot, scansForUniqueTimes,  # task version
                                mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes,  # au version
                                myprint=debugSloppyMatch
                                )
                if (sm):
                    if (debug):
                        print("1)setting myUniqueTime to %d" % (mytime))
                    myUniqueTime = mytime
                    ctr += 1
            if (debug): print("Found a match")
#            if (ctr > len(fieldIndicesToPlot) and bOverlay==False):
#                print("multi-field time overlay ***************  why are there more matches (%d) than fields (%d)?" % (ctr,len(fieldIndicesToPlot)))
        
#        print "Overlay antenna %d, myUniqueTime=%d" % (xctr, myUniqueTime)
        if (xframe == xframeStart):
            if (debug):
                print("Clearing the page")
            pb.clf()
        xflag = [item for sublist in xflag for item in sublist]
        yflag = [item for sublist in yflag for item in sublist]
#        pflag = [xflag, yflag]
#        flagfrequencies = [frequencies, frequencies2]
        antstring = buildAntString(xant,msFound,msAnt)
        if (msFound):
            fieldString = msFields[uniqueFields[fieldIndex]]  # [0]
        else:
            fieldString = str(field)
        if (overlayTimes):
            timeString =''
        else:
            timeString = ',  t%d/%d  %s' % (mytime,nUniqueTimes-1,utstring(uniqueTimes[mytime],3))
            if (scansForUniqueTimes != []):
                if (scansForUniqueTimes[mytime]>=0):
                    timeString = ',  scan%d  %s' % (scansForUniqueTimes[mytime],utstring(uniqueTimes[mytime],3))
        spwString = buildSpwString(overlaySpws, overlayBasebands, spwsToPlot, 
                                   ispw, originalSpw[ispw], observatoryName, 
                                   baseband, showBasebandNumber)
        titleString = "%sspw%s,  field %d: %s%s" % (antennaString,spwString,uniqueFields[fieldIndex],fieldString,timeString)
        if (sum(xflag)==nChannels and sum(yflag)==nChannels and showflagged==False):
            if (overlayTimes):
                print("Skip %s, xant=%2d (%s) for time%d=%s all solutions flagged" % (antstring, xant, titleString,mytime,utstring(uniqueTimes[mytime],3)))
                # need to set doneOverlayTime = True if this is the final time,
                # otherwise, we get "subplot number exceeds total subplots" at line 2427
                # but we need to draw the labels at the top of the page, else they will not get done
                if (debug):
                    print("########## uniqueTimes[%d]=%d,  timerangeListTimes[-1]=%d" % (mytime,uniqueTimes[mytime],timerangeListTimes[-1]))
                    print("########## scansToPlotPerSpw[%d]=%s, scansForUniqueTimes=%s" % (ispw,str(scansToPlotPerSpw[ispw]),str(scansForUniqueTimes)))
                if (len(scansToPlotPerSpw[ispw]) < 1):
                    sTPPS = []
                else:
#                    sTPPS = [scansToPlot[ispw][-1]]# added [[-1]]  on 2014-04-04   task version
                    sTPPS = [scansToPlotPerSpw[ispw][-1]]# added [[-1]]  on 2014-04-04  au version
                if (sloppyMatch(timerangeListTimes[-1], uniqueTimes[mytime],
                                solutionTimeThresholdSeconds,
                                mytime, sTPPS, scansForUniqueTimes,  
                                myprint=debugSloppyMatch
                                )):
                    if (overlayAntennas == False or xant==antennasToPlot[-1]):  # 11-Mar-2014
                        doneOverlayTime = True  # 08-Nov-2012
                        finalTimerangeFlagged =  True
                    if (debug):
                        print("###### set doneOverlayTime = True, xant=%d, lastUnflaggedAntennaToPlot=%d, drewAtmosphere=%s" % (xant,lastUnflaggedAntennaToPlot,drewAtmosphere))
                    # draw labels
                    # try adding the following 'if' statement on Jun 18, 2013; it works.
#                    if (drewAtmosphere==False or overlayAntennas==False):
                    # Add the 'and not' case to prevent extra atm/fdms shown if one spw's solutions are all flagged
                    if (drewAtmosphere==False or (overlayAntennas==False and not allTimesFlaggedOnThisSpw)):
                        drawOverlayTimeLegends(xframe,firstFrame,xstartTitle,ystartTitle,
                                               caltable,titlesize,fieldIndicesToPlot,
                                               ispwInCalTable,uniqueTimesPerFieldPerSpw,
                                               timerangeListTimes,
                                               solutionTimeThresholdSeconds,
                                               debugSloppyMatch,
                                               ystartOverlayLegend,debug,mysize,
                                               fieldsToPlot,myUniqueColor,
                                               timeHorizontalSpacing,fieldIndex,
                                               overlayColors,
                                               antennaVerticalSpacing, overlayAntennas,
                                               timerangeList, caltableTitle, mytime,
                                               scansToPlotPerSpw[ispw], scansForUniqueTimes) # au version
                        if LO1 is None and type(lo1s) == dict:  # Fix for SCOPS-4877
                            LO1 = lo1s[myspw]                   # Fix for SCOPS-4877
                        # CAS-8655
                        newylimits = drawAtmosphereAndFDM(showatm,showtsky,atmString,subplotRows,mysize,
                                             TebbSky,TebbSkyImage,plotrange, xaxis,atmchan,
                                             atmfreq,transmission,subplotCols,showatmPoints,
                                             xframe, channels,LO1,atmchanImage,atmfreqImage,
                                             transmissionImage,firstFrame,showfdm,nChannels,
                                             tableFormat,originalSpw_casa33,
                                             chanFreqGHz_casa33,originalSpw,chanFreqGHz,
                                             overlayTimes, overlayAntennas, xant,
                                             antennasToPlot, overlaySpws, baseband,
                                             showBasebandNumber, basebandDict,
                                             overlayBasebands, drewAtmosphere, showtsys)
                        drewAtmosphere = True
                    if (xctr == firstUnflaggedAntennaToPlot or overlayAntennas==False): # changed xant->xctr on 11-mar-2014
                        DrawPolarizationLabelsForOverlayTime(xstartPolLabel,ystartPolLabel,corr_type,polsToPlot,
                                                             channeldiff,ystartMadLabel,subplotRows,gamp_mad,mysize,
                                                             ampmarkstyle,markersize,ampmarkstyle2,gamp_std)
                    elif (overlayAntennas):
                        if (debug):
                            print("xant=%d, firstUnflaggedAntennaToPlot=%d" % (xant,firstUnflaggedAntennaToPlot))
            else:  # not overlaying times
                print("Skipping %s, xant=%d, ispw=%d (%s) all solutions flagged" % (antstring, xant, ispw, titleString))
                if ((overlaySpws or overlayBasebands) and spwctr==spwctrFirstToPlot):
                     spwctrFirstToPlot += 1
                if ((overlaySpws or overlayBasebands) and ispw==spwsToPlotInBaseband[bbctr][-1]):
                    if (debug): print("The final spw was flagged!!!!!!!!!!!!!!")
                    finalSpwWasFlagged =  True  # inserted on 22-Apr-2014 for g25.27
                if (myinput == 'b'):
                    redisplay = False # This prevents infinite loop when hitting 'b' on first screen when ant0 flagged. 2013-03-08
            if (overlayAntennas==False and overlayBasebands==False): # 07/30/2014  added  overlayBasebands==False
                if (doneOverlayTime==False or overlayTimes==False):  # added on 08-Nov-2012
                    finalSpwWasFlagged = False # Added on 23-Apr-2014 for regression61
                    mytime += 1
#                    print " 0003  incrementing mytime to ", mytime
                    if (debug):
                        print("F) all solutions flagged --> incrementing mytime to %d" % mytime)
            if (overlayAntennas): 
                if (xctr == firstUnflaggedAntennaToPlot):
                    firstUnflaggedAntennaToPlot += 1
                    if (firstUnflaggedAntennaToPlot >= len(antennasToPlot)):
                        firstUnflaggedAntennaToPlot = 0
                        if not finalSpwWasFlagged: # Added first test on 23-Apr-2014 for regression61
                            mytime += 1
#                            print " 0002  incrementing mytime to ", mytime
                    if (debug):
                        print("                     A) incrementing mytime to ", mytime)
                        print("----- Resetting firstUnflaggedAntennaToPlot from %d to %d" % (firstUnflaggedAntennaToPlot-1, firstUnflaggedAntennaToPlot))
                        print("-----    = antenna %d" % (antennasToPlot[firstUnflaggedAntennaToPlot]))
                    if (mytime < nUniqueTimes):  # Add this 'if' conditional on 9-22-2015 for CAS-7839
                        continue # Try adding this statement on Apr 2, 2012 to fix bug.
                    mytime -= 1  # Add this on 9-22-2015 for CAS-7839
            if (overlaySpws or overlayBasebands):
                if (xctr == firstUnflaggedAntennaToPlot):
                    firstUnflaggedAntennaToPlot += 1
                    if (firstUnflaggedAntennaToPlot >= len(antennasToPlot)):
                        firstUnflaggedAntennaToPlot = 0
                        if not finalSpwWasFlagged: # Added on 22-Apr-2014 for g25.27 dataset antenna='4'
                            if (overlayBasebands == False or spwctr>len(spwsToPlot)):  # Added on 7/30/2014 for regression 96
                                mytime += 1
#                                print " 0001  incrementing mytime to ", mytime
                    if (debug):
                        print("                     B) incrementing mytime to ", mytime)
                        print("----- Resetting firstUnflaggedAntennaToPlot from %d to %d" % (firstUnflaggedAntennaToPlot-1, firstUnflaggedAntennaToPlot))
                        print("-----    = antenna %d" % (antennasToPlot[firstUnflaggedAntennaToPlot]))
                    if (not finalSpwWasFlagged): # add this test on Apr 22, 2014 to prevent crash on g25.27 dataset with antenna='4,5'
                        continue # Try this 'continue' on Apr 2, 2012 to fix bug -- works.
            if (overlayAntennas==False and subplot==11 
                and not finalSpwWasFlagged    # inserted on 22-Apr-2014 for g25.27
                and not finalTimerangeFlagged # inserted on 04-Aug-2014 for CAS-6812
                ): 
                  # added the case (subplot==11) on April 22, 2012 to prevent crash
                  # on multi-antenna subplot=421
                  if (debug):
                      print("#######  removing [%d,%d,%d,%d]  len(pages) was %d" % (pages[len(pages)-1][PAGE_ANT],
                                                        pages[len(pages)-1][PAGE_SPW],
                                                        pages[len(pages)-1][PAGE_TIME],
                                                        pages[len(pages)-1][PAGE_AP],len(pages)))
                  pages = pages[0:len(pages)-1]
                  newpage = 1
            if (overlayAntennas==False):
                if (doneOverlayTime==False  # inserted on 08-Nov-2012
                    and not finalSpwWasFlagged):  # inserted on 22-Apr-2014 for g25.27
                    if (debug): print("========== continuing before plotting mytime=%d" % (mytime))
                    continue
                elif (debug):
                    print("=========== Not continuing because doneOverlayTime=",doneOverlayTime)
        else:
            allTimesFlaggedOnThisSpw = False
            if (debug):
                print("Not all the data are flagged.  doneOverlayTime=%s" % (str(doneOverlayTime)))
            
        if (firstSpwMatch == -1):
            firstSpwMatch = spwctr
        if (firstTimeMatch == -1):
            firstTimeMatch = mytime
            if (debug):
                print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Setting firstTimeMatch from -1 to ", firstTimeMatch)
        # The following was needed to support overlay='antenna,time' for showatm for QA2 report (CAS-7820)
        if (finalTimeMatch == -1 or finalTimeMatch < mytime):
            if (debug):
                print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Setting finalTimeMatch from %d to %d" % (finalTimeMatch, mytime))
            finalTimeMatch = mytime

######### Here is the amplitude plotting ############
        if (yaxis.find('amp')>=0 or yaxis.find('both')>=0 or yaxis.find('ap')>=0) and doneOverlayTime==False:
          if (overlayBasebands and amplitudeWithPhase): # CAS-6477
              if (float(xframe/10) != xframe*0.1 and alreadyPlottedAmp):
                  xframe -= 2
          if (debug):
              if (mytime < len(uniqueTimes)):
                  print("amp: xctr=%d, xant=%d, myap=%d, mytime=%d(%s), firstTimeMatch=%d, bOverlay=" % (xctr, xant, myap, mytime, utstring(uniqueTimes[mytime],3), firstTimeMatch), bOverlay)
          if (myap==1):
            if (overlayTimes == False or mytime==firstTimeMatch):
              if ((overlaySpws == False and overlayBasebands==False) or 
                  spwctr == spwctrFirstToPlot or 
#                  (spwctr==lowestSpwInFirstScan and str(xframe)[-1] == '0') or 
                  spwctr>len(spwsToPlot)): 
                if (overlayAntennas==False or xctr==firstUnflaggedAntennaToPlot
                    or xctr==antennasToPlot[-1]):  # 2012-05-24, to fix the case where all ants flagged on one timerange
                    xframe += 1
                    if (debug):
                        print("y) incrementing xframe to %d" % xframe)
                    if (debug):
                        print("mytime=%d  ==  firstTimeMatch=%d" % (mytime, firstTimeMatch))
                        print("xctr=%d  ==  firstUnflaggedAntennaToPlot=%d,  antennastoPlot[-1]=%d" % (xctr, firstUnflaggedAntennaToPlot,antennasToPlot[-1]))
                        print("spwctr=%d  >? len(spwsToPlot)=%d" % (spwctr, len(spwsToPlot)))
                    myUniqueColor = []
                    newylimits = [LARGE_POSITIVE, LARGE_NEGATIVE]
          else: # (myap == 0)
            if (overlayTimes == False or mytime==firstTimeMatch):
              if ((overlaySpws == False and overlayBasebands==False) or 
                  spwctr==spwctrFirstToPlot or 
#                  (spwctr==lowestSpwInFirstScan and str(xframe)[-1] == '0') or 
                  (overlayBasebands and amplitudeWithPhase) or # CAS-6477
                  spwctr>len(spwsToPlot)):
                if (overlayAntennas==False or xctr==firstUnflaggedAntennaToPlot
                    or xctr>antennasToPlot[-1]):  # 2012-05-24, to fix the case where all ants flagged on one timerange
                    xframe += 1
                    if (debug):
                        print("Y) incrementing xframe to %d" % xframe)
                    if (debug):
                        print("mytime=%d  ==  firstTimeMatch=%d" % (mytime, firstTimeMatch))
                        print("xctr=%d  ==  firstUnflaggedAntennaToPlot=%d,  antennasToPlot[-1]=%d" % (xctr, firstUnflaggedAntennaToPlot,antennasToPlot[-1]))
                        print("spwctr=%d  >? len(spwsToPlot)=%d" % (spwctr, len(spwsToPlot)))
                    myUniqueColor = []
                    newylimits = [LARGE_POSITIVE, LARGE_NEGATIVE]
                    if (debug):
                        print("myap=%d, mytime == firstTimeMatch" % myap, firstTimeMatch)
                else:
                    if (debug): print("4)Not incrementing xframe from %d" % (xframe))
              else:
                 if (debug): print("2)Not incrementing xframe from %d (spwctr=%d >? len(spwsToPlot)=%d) or (spwctr=%d == spwctrFirstToPlot=%d)" % (xframe,spwctr,len(spwsToPlot),spwctr,spwctrFirstToPlot))
            else:
                if (debug): print("1)Not incrementing xframe from %d" % (xframe))
            if (debug):
                print("$$$$$$$$$$$$$$$$$$$$$$$  ready to plot amp on xframe %d" % (xframe))
            adesc = pb.subplot(xframe)
            if (previousSubplot != xframe):
                drewAtmosphere = False
            previousSubplot = xframe
            alreadyPlottedAmp = True  # needed for (overlay='baseband', yaxis='both')  CAS-6477
#            pb.hold(overlayAntennas or overlayTimes or overlaySpws or overlayBasebands) # not available in CASA6, but never needed
            if (delay):
                gampx = gplotx
            else:
                gampx = np.abs(gplotx)
            if (nPolarizations == 2):
                if (delay):
                    gampy = gploty
                else:
                    gampy = np.abs(gploty)
                if (yaxis.lower().find('db') >= 0):
                    gamp = [10*np.log10(gampx), 10*np.log10(gampy)]
                else:
                    if (channeldiff>0):
                        if (xaxis == 'chan'):
                            gamp0, newx0, gamp0res, newx0res = channelDifferences(gampx, pchannels[0], resample)
                            gamp1, newx1, gamp1res, newx1res = channelDifferences(gampy, pchannels[1], resample)
                            pchannels = [newx0, newx1]
                        else:
                            gamp0, newx0, gamp0res, newx0res  = channelDifferences(gampx, pfrequencies[0], resample)
                            gamp1, newx1, gamp1res, newx1res = channelDifferences(gampy, pfrequencies[1], resample)
                            pfrequencies = [newx0, newx1]
                        gamp = [gamp0, gamp1]
                        gampres = [gamp0res, gamp1res]
                        if (VisCal.lower().find('tsys') >= 0 and tsysPercent):
                            gamp = [100*gamp0/np.median(gampx), 100*gamp1/np.median(gampy)]
                            gampres = [100*gamp0res/np.median(gampx), 100*gamp1res/np.median(gampy)]
                        elif (VisCal.lower().find('tsys') < 0 and ampPercent):
                            gamp = [100*gamp0/np.median(gampx), 100*gamp1/np.median(gampy)]
                            gampres = [100*gamp0res/np.median(gampx), 100*gamp1res/np.median(gampy)]
                        gamp_mad = [madInfo(gampres[0],madsigma,edge,ispw,xant,0), madInfo(gampres[1],madsigma,edge,ispw,xant,1)]
                        gamp_std = [stdInfo(gampres[0],madsigma,edge,ispw,xant,0), stdInfo(gampres[1],madsigma,edge,ispw,xant,1)]
                        if (platformingSigma > 0):
                            platformingThresholdX = gamp_mad[0]['mad']*platformingSigma
                            platformingThresholdY = gamp_mad[1]['mad']*platformingSigma
                        else:
                            platformingThresholdX = platformingThreshold
                            platformingThresholdY = platformingThreshold
                        gamp_platforming = [platformingCheck(gamp[0],platformingThresholdX),
                                            platformingCheck(gamp[1],platformingThresholdY)]
                        for p in [0,1]:
#                            print "gamp_mad[%d] = " % (p), gamp_mad[p]
#                            print "madstats[%s][%d] = " % (Antstring,ispw), madstats[Antstring][ispw]
                            madstats[Antstring][ispw][mytime][p]['amp'] = gamp_mad[p]['mad']
                            madstats[Antstring][ispw][mytime][p]['ampstd'] = gamp_std[p]['std']
                            if (gamp_platforming[p]):
                                if (Antstring not in madstats['platforming'].keys()):
                                    madstats['platforming'][Antstring] = {}
                                if (ispw not in madstats['platforming'][Antstring].keys()):
                                    madstats['platforming'][Antstring][ispw] = {}
                                if (p not in madstats['platforming'][Antstring][ispw].keys()):
                                    madstats['platforming'][Antstring][ispw][p] = []
                                madstats['platforming'][Antstring][ispw][p].append(uniqueTimes[mytime])
                            if (gamp_mad[p]['nchan'] > 0):
                                print("%s, Pol %d, spw %2d, %s, amp: %4d points exceed %.1f sigma (worst=%.2f at chan %d)" % (Antstring, p, ispw, utstring(uniqueTimes[mytime],0), gamp_mad[p]['nchan'], madsigma, gamp_mad[p]['outlierValue'], gamp_mad[p]['outlierChannel']+pchannels[p][0]))
                    else:
                        gamp = [gampx,gampy]
            else:
                if (yaxis.lower().find('db') >= 0):
                    gamp = [10*np.log10(gampx)]
                else:
                    if (channeldiff>0):
                        if (xaxis == 'chan'):
                            gamp0, newx0, gamp0res, newx0res  = channelDifferences(gampx, pchannels[0], resample)
                            pchannels = [newx0]
                        else:
                            gamp0, newx0, gamp0res, newx0res  = channelDifferences(gampx, pfrequencies[0], resample)
                            pfrequencies = [newx0]
                        gamp = [gamp0]
                        gampres = [gamp0res]
                        if (VisCal.lower().find('tsys') >= 0 and tsysPercent):
                            gamp = [100*gamp0/np.median(gampx)]
                            gampres = [100*gamp0res/np.median(gampx)]
                        elif (VisCal.lower().find('tsys') < 0 and ampPercent):
                            gamp = [100*gamp0/np.median(gampx)]
                            gampres = [100*gamp0res/np.median(gampx)]
                        p = 0
                        gamp_mad = [madInfo(gampres[p], madsigma,edge,ispw,xant,p)]
                        gamp_std = [stdInfo(gampres[p], madsigma,edge,ispw,xant,p)]
                        if (platformingSigma > 0):
                            platformingThresholdX = gamp_mad[0]['mad']*platformingSigma
                        else:
                            platformingThresholdX = platformingThreshold
                        gamp_platforming = [platformingCheck(gamp[p], platformingThresholdX)]
                        madstats[Antstring][ispw][mytime][p]['amp'] = gamp_mad[p]['mad']
                        madstats[Antstring][ispw][mytime][p]['ampstd'] = gamp_std[p]['std']
                        if (gamp_platforming[p]):
                            if (Antstring not in madstats['platforming'].keys()):
                                madstats['platforming'][Antstring] = {}
                            if (ispw not in madstats['platforming'][Antstring].keys()):
                                madstats['platforming'][Antstring][ispw] = {}
                            if (p not in madstats['platforming'][Antstring][ispw].keys()):
                                madstats['platforming'][Antstring][ispw][p] = []
                            madstats['platforming'][Antstring][ispw][p].append(mytime)
                        if (gamp_mad[p]['nchan'] > 0):
                            print("%s, Pol %d, spw %2d, %s, amp: %4d points exceed %.1f sigma (worst=%.2f at chan %d)" % (Antstring, p, ispw, utstring(uniqueTimes[mytime],0), gamp_mad[p]['nchan'], madsigma, gamp_mad[p]['outlierValue'], gamp_mad[p]['outlierChannel']+pchannels[p][0]))
                    else:
                        gamp = [gampx]
            if (bOverlay):
                if (delay):
                    gampx2 = gplotx2 + caltable2amplitudeOffset
                else:
                    gampx2 = np.abs(gplotx2) + caltable2amplitudeOffset
                if (nPolarizations == 2):
                    if (delay):
                        gampy2 = gploty2 + caltable2amplitudeOffset
                    else:
                        gampy2 = np.abs(gploty2) + caltable2amplitudeOffset
                    if (yaxis.lower().find('db') >= 0):
                        gamp2 = [10*np.log10(gampx2), 10*np.log10(gampy2)]
                    else:
                        if (channeldiff>0):
                            if (xaxis == 'chan'):
                                gamp2_0, newx0, gamp2_0res, newx0res = channelDifferences(gampx2, pchannels2[0], resample)
                                gamp2_1, newx1, gamp2_1res, newx1res = channelDifferences(gampy2, pchannels2[1], resample)
                                pchannels2 = [newx0, newx1]
                            else:
                                gamp2_0, newx0, gamp2_0res, newx0res = channelDifferences(gampx2, pfrequencies2[0], resample)
                                gamp2_1, newx1, gamp2_1res, newx1res = channelDifferences(gampy2, pfrequencies2[1], resample)
                                pfrequencies2 = [newx0, newx1]
                            gamp2 = [gamp2_0, gamp2_1]
                            gamp2res = [gamp2_0res, gamp2_1res]
                            if (VisCal.lower().find('tsys') >= 0 and tsysPercent):
                                gamp2 = [100*gamp2_0/np.median(gampx2), 100*gamp2_1/np.median(gampy2)]
                                gamp2res = [100*gamp2_0res/np.median(gampx2), 100*gamp2_1res/np.median(gampy2)]
                            elif (VisCal.lower().find('tsys') < 0 and ampPercent):
                                gamp2 = [100*gamp2_0/np.median(gampx2), 100*gamp2_1/np.median(gampy2)]
                                gamp2res = [100*gamp2_0res/np.median(gampx2), 100*gamp2_1res/np.median(gampy2)]
                        else:
                            gamp2 = [gampx2, gampy2]
                else:
                    if (yaxis.lower().find('db') >= 0):
                        gamp2 = [10*np.log10(gampx2)]
                    else:
                        if (channeldiff>0):
                            if (xaxis == 'chan'):
                                gamp2_0, newx0, gamp2_0res, newx0res = channelDifferences(gampx2, pchannels[0], resample)
                                pchannels2 = [newx0]
                            else:
                                gamp2_0, newx0, gamp2_0res, newx0res = channelDifferences(gampx2, pfrequencies[0], resample)
                                pfrequencies2 = [newx0]
                            gamp2 = [gamp2_0]
                            gamp2res = [gamp2_0res]
                            if (VisCal.lower().find('tsys') >= 0 and tsysPercent):
                                gamp2 = [100*gamp2_0/np.median(gampx2)]
                                gamp2res = [100*gamp2_0res/np.median(gampx2)]
                            elif (VisCal.lower().find('tsys') < 0 and ampPercent):
                                gamp2 = [100*gamp2_0/np.median(gampx2)]
                                gamp2res = [100*gamp2_0res/np.median(gampx2)]
                        else:
                            gamp2 = [gampx2]
            if (xaxis.find('chan')>=0 or (msFound==False and tableFormat==33)):    #  'amp'
                if (debug):
                    print("amp: plot vs. channel **********************")
#                pb.hold(True) # not available in CASA6, but never needed
                for p in range(nPolarizations):
                    if (overlayAntennas or overlayTimes):
                        if (corr_type_string[p] in polsToPlot):
                              pdesc = pb.plot(pchannels[p],gamp[p],'%s'%ampmarkstyles[p],
                                              markersize=markersize,
                                              markerfacecolor=overlayColors[xctr],markeredgewidth=markeredgewidth)
                              newylimits =  recalcYlimits(plotrange,newylimits,gamp[p])
                              if (overlayAntennas and overlayTimes==False):
                                  pb.setp(pdesc, color=overlayColors[xctr])
                              elif (overlayTimes and overlayAntennas==False):
                                  pb.setp(pdesc, color=overlayColors[mytime])
                              elif (overlayTimes and overlayAntennas): # try to support time,antenna
                                  if (debug):
                                      print("p=%d, len(fieldsToPlot)=%d, len(timerangeList)=%d" % (p,len(fieldsToPlot),len(timerangeList)))
                                  if (len(fieldsToPlot) > 1 or len(timerangeList)>1):
                                      # The third 'or' below is needed if pol='0' is flagged on antenna 0. -- 2012/10/12
                                      if (p==0 or len(polsToPlot)==1 or myUniqueColor==[]):
                                          myUniqueColor.append(overlayColors[len(myUniqueColor)])
                                      pb.setp(pdesc, color=myUniqueColor[-1])
                                      
                    else:
                        if (corr_type_string[p] in polsToPlot):
#                          print "pcolor[%d]=%s" % (p,pcolor)
                          pb.plot(pchannels[p],gamp[p],'%s%s'%(pcolor[p],ampmarkstyle), markersize=markersize,markeredgewidth=markeredgewidth)
                          newylimits =  recalcYlimits(plotrange,newylimits,gamp[p])
                          if asciiFile:
                              asciiDesc = open('pol%d_vs_chan.txt' % p,'w')
                              for i in range(len(pchannels[p])):
                                  asciiDesc.write('%f %f\n'%(pchannels[p][i], gamp[p][i]))
                              asciiDesc.close()
                if (sum(xflag)>0):
                    myxrange = np.max(channels)-np.min(channels)
                    SetNewXLimits([np.min(channels)-myxrange/20, np.max(channels)+myxrange/20],loc=2)
#                    print "amp: Resetting xaxis channel range to counteract flagged data"
                if (xframe in bottomRowFrames or (xctr+1==len(antennasToPlot) and ispw==spwsToPlot[-1])):
                    pb.xlabel("Channels (%d)" % (len(pchannels[p])), size=mysize)
            elif (xaxis.find('freq')>=0):   # amp
                if (bOverlay):
#                      pb.hold(True)# not available in CASA6, but never needed
                      myxrange = np.abs(xfrequencies[0]-xfrequencies[-1])
                      try:
                          xrange2 = np.abs(xfrequencies2[0]-xfrequencies2[-1])
                      except:
                          print("No amp data found in second solution for spw %d.  Try increasing the solutionTimeThresholdSeconds above %.0f." % (ispw,solutionTimeThresholdSeconds))
                          print("If this doesn't work, email the developer (%s)." % (developerEmail))
                          return(vm)

                      if (np.abs(myxrange/xrange2 - 1) > 0.05 + len(xflag)/len(xchannels)):  # 0.0666 is 2000/1875-1
                         # These line widths are optimal for visualizing FDM over TDM
#                         print "***  Solutions differ in frequency width"
                         width1 = 1
                         width2 = 4
                         # solutions differ in frequency width
                         if (myxrange < xrange2):
                            for p in range(nPolarizations):
                                  if (corrTypeToString(corr_type[p]) in polsToPlot):
                                        pdesc = pb.plot(pfrequencies[p], gamp[p], '%s%s'%(pcolor[p],ampmarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                                        newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,1,chanrangePercent)
                            for p in range(nPolarizations):
                                  if (corrTypeToString(corr_type[p]) in polsToPlot):
                                        pdesc = pb.plot(pfrequencies2[p], gamp2[p], '%s%s'%(p2color[p],ampmarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                                        newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp2[p], sideband,plotrange,xchannels2,debug,2,chanrangePercent)
                         else:
                            for p in range(nPolarizations):
                                  if (corrTypeToString(corr_type[p]) in polsToPlot):
                                        pdesc = pb.plot(pfrequencies2[p], gamp2[p], '%s%s'%(p2color[p],ampmarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                                        newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp2[p], sideband,plotrange,xchannels2,debug,3,chanrangePercent)
                            for p in range(nPolarizations):
                                  if (corrTypeToString(corr_type[p]) in polsToPlot):
                                        pdesc = pb.plot(pfrequencies[p], gamp[p], '%s%s'%(pcolor[p],ampmarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                                        newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,4,chanrangePercent)
                      else:
                         width1 = 1
                         width2 = 2  # Just enough to distinguish one plotted line from the other.
                         # solutions may be different level of smoothing, so plot highest rms first
                         if (madOfDiff(gamp[0]) < madOfDiff(gamp2[0]) and firstPlot != 1):
                            # plot second solution first
#                            print "**** MAD of caltable=%f < caltable2=%f" % (madOfDiff(gamp[0]), madOfDiff(gamp2[0]))
                            for p in range(nPolarizations):
                                if (corrTypeToString(corr_type[p]) in polsToPlot):
                                    pdesc = pb.plot(pfrequencies2[p], gamp2[p], '%s%s'%(p2color[p],ampmarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                                    newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp2[p], sideband,plotrange,xchannels2,debug,5,chanrangePercent)
                            for p in range(nPolarizations):
                                if (corrTypeToString(corr_type[p]) in polsToPlot):
                                    pdesc = pb.plot(pfrequencies[p], gamp[p], '%s%s'%(pcolor[p],ampmarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                                    newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,6,chanrangePercent)
                         else:
                            # plot first solution first
#                            print "**** MAD of caltable=%f >= caltable2=%f" % (madOfDiff(gamp[0]), madOfDiff(gamp2[0]))
                            for p in range(nPolarizations):
                                if (corrTypeToString(corr_type[p]) in polsToPlot):
                                    pdesc = pb.plot(pfrequencies[p], gamp[p], '%s%s'%(pcolor[p],ampmarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                                    newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,7,chanrangePercent)
                            for p in range(nPolarizations):
                                if (corrTypeToString(corr_type[p]) in polsToPlot):
                                    pdesc = pb.plot(pfrequencies2[p], gamp2[p], '%s%s'%(p2color[p],ampmarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                                    newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp2[p], sideband,plotrange,xchannels2,debug,8,chanrangePercent)
                      # must set new limits after plotting  'amp'
                      if (zoom=='intersect'):
                          if (myxrange < xrange2):
                              SetNewXLimits([min(xfrequencies[0],xfrequencies[-1])-myxrange*0.1, max(xfrequencies[0],xfrequencies[-1])+myxrange*0.1],loc=3)
                              SetLimits(plotrange, chanrange, newylimits, channels, frequencies,
                                        pfrequencies, ampMin, ampMax, xaxis, pxl, chanrangeSetXrange,
                                        chanrangePercent,loc=101)
                          else:
#                              print "len(xfrequencies2) = ", len(xfrequencies2)
                              SetNewXLimits([min(xfrequencies2[0],xfrequencies2[-1])-xrange2*0.1, max(xfrequencies2[0],xfrequencies2[-1])+xrange2*0.1],loc=4)
                              slstatus = SetLimits(plotrange, chanrange, newylimits, channels, frequencies2,
                                                   pfrequencies2, ampMin, ampMax, xaxis, pxl,
                                                   chanrangeSetXrange, chanrangePercent,loc=102)
                      else:
                          if (myxrange < xrange2):
                              SetLimits(plotrange, chanrange, newylimits, channels, frequencies,
                                        pfrequencies, ampMin, ampMax, xaxis, pxl,
                                        chanrangeSetXrange, chanrangePercent,loc=103)
                          else:
                              SetLimits(plotrange, chanrange, newylimits, channels, frequencies2,
                                        pfrequencies2, ampMin, ampMax, xaxis, pxl,
                                        chanrangeSetXrange, chanrangePercent,loc=104)
                      # draw polarization and spw labels
                      if (xframe == firstFrame):
                          # draw title including caltable name
                          caltableList = 'c1=' + caltable + ', c2=' + caltable2 # + ' (%s)'%(utstring(uniqueTimes2[mytime],3))
                          pb.text(xstartTitle, ystartTitle, caltableList, size=titlesize,
                                  color='k', transform=pb.gcf().transFigure)
                          if (caltable2amplitudeOffset != 0):
                              pb.text(xstartTitle, 0.935, 'c2 amplitude offset = %.3f' % (caltable2amplitudeOffset),
                                      color='k',size=titlesize,transform=pb.gcf().transFigure)
                elif (bpolyOverlay):
                    if (debug):
                        print("in bpolyOverlay **********************************")
                    matches1 = []
                    for tbp in range(len(timesBP)):
                        if (sloppyMatch(uniqueTimes[mytime], timesBP[tbp], solutionTimeThresholdSeconds,
                                        mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes,  # au version
                                        myprint=debugSloppyMatch
                                        )):
                            matches1.append(tbp)
                    matches1 = np.array(matches1)
                    if (len(matches1) < 1):
                        print("No time match found between %.1f and " % (uniqueTimes[mytime]), timesBP)
                        print("If you are sure the solutions correspond to the same data, you can set solutionTimeThresholdSeconds>=%.0f" % (1+np.ceil(np.abs(timesBP[0]-uniqueTimes[mytime]))))
                        return(vm)
                    matches2 = np.where(xant == np.array(antennasBP))[0]
                    if (len(matches2) < 1):
                        print("No antenna match found: ", xant, antennasBP)
                    if (tableFormat == 33):
                        matches3 = np.where(ispw == np.array(cal_desc_idBP))[0]
                        if (len(matches3) < 1):
                            print("No spw match found: %d not in " % (ispw), cal_desc_idBP)
                    else:
                        matches3 = np.where(ispw == np.array(spwBP))[0]
                        if (len(matches3) < 1):
                            print("No spw match found: %d not in " % (ispw), spwBP)
                    matches12 = np.intersect1d(matches1,matches2)
                    if (len(matches12) < 1):
                        print("No time+antenna match between: ", matches1, matches2)
                    matches = np.intersect1d(matches12, matches3)
                    if (len(matches) < 1):
                        print("No time+antenna+spw match between: ", matches12, matches3)
                    try:
                        index = matches[0]
                        if (debug):
                            print("Match = %d ***********************************" % (index))
                    except:
                        print("No match found for time=%.6f, xant=%d, ispw=%d"  % (uniqueTimes[mytime],xant,ispw))
                        print("antennasBP = ", antennasBP)
                        print("cal_desc_idBP = ", cal_desc_idBP)
                        print("timesBP = ")
                        for i in timesBP:
                            print("%.6f, " % i)
                        return(vm)
                    validDomain = [frequencyLimits[0,index], frequencyLimits[1,index]]
                    cc = calcChebyshev(polynomialAmplitude[index][0:nPolyAmp[index]], validDomain, frequenciesGHz[index]*1e+9)
                    fa = np.array(frequenciesGHz[index])
                    if (xfrequencies[0] < xfrequencies[-1]):
                        matches = np.where(fa>xfrequencies[0])[0]
                        matches2 = np.where(fa<xfrequencies[-1])[0]
                    else:
                        matches = np.where(fa>xfrequencies[-1])[0]
                        matches2 = np.where(fa<xfrequencies[0])[0]
                    if (len(matches) < 1):
                        print("looking for %f-%f GHz inside %f-%f" % (xfrequencies[0],xfrequencies[-1],fa[0],fa[-1]))
                    amplitudeSolutionX = np.mean(gampx)*(cc-np.mean(cc)+1)

                    cc = calcChebyshev(polynomialAmplitude[index][nPolyAmp[index]:2*nPolyAmp[index]], validDomain, frequenciesGHz[index]*1e+9)
                    if (debug):
                        print("nPol=%d, len(xfrequencies)=%d, len(yfrequencies)=%d" % (nPolarizations,len(xfrequencies),len(yfrequencies)))
                    if (nPolarizations > 1):
                        if (yfrequencies[0] < yfrequencies[-1]):
                            matches = np.where(fa>yfrequencies[0])[0]
                            matches2 = np.where(fa<yfrequencies[-1])[0]
                        else:
                            matches = np.where(fa>yfrequencies[-1])[0]
                            matches2 = np.where(fa<yfrequencies[0])[0]
                        amplitudeSolutionY = np.mean(gampy)*(cc-np.mean(cc)+1)
                    if (bpolyOverlay2):
                        validDomain = [frequencyLimits2[0,index], frequencyLimits2[1,index]]
                        cc = calcChebyshev(polynomialAmplitude2[index][0:nPolyAmp2[index]],
                                           validDomain, frequenciesGHz2[index]*1e+9)
                        fa = np.array(frequenciesGHz2[index])
                        if (xfrequencies[0] < xfrequencies[-1]):
                            matches = np.where(fa>xfrequencies[0])[0]
                            matches2 = np.where(fa<xfrequencies[-1])[0]
                        else:
                            matches = np.where(fa>xfrequencies[-1])[0]
                            matches2 = np.where(fa<xfrequencies[0])[0]
                        amplitudeSolution2X = np.mean(gampx)*(cc-np.mean(cc)+1)

                        cc = calcChebyshev(polynomialAmplitude2[index][nPolyAmp2[index]:2*nPolyAmp2[index]],
                                           validDomain, frequenciesGHz2[index]*1e+9)
                        fa = np.array(frequenciesGHz2[index])
                        if (debug):
                            print("nPol=%d, len(xfrequencies)=%d, len(yfrequencies)=%d" % (nPolarizations,len(xfrequencies),len(yfrequencies)))
                        if (nPolarizations > 1):
                            if (yfrequencies[0] < yfrequencies[-1]):
                                matches = np.where(fa>yfrequencies[0])[0]
                                matches2 = np.where(fa<yfrequencies[-1])[0]
                            else:
                                matches = np.where(fa>yfrequencies[-1])[0]
                                matches2 = np.where(fa<yfrequencies[0])[0]
                            amplitudeSolution2Y = np.mean(gampy)*(cc-np.mean(cc)+1)

#                        pb.hold(True)# not available in CASA6, but never needed
                        for p in range(nPolarizations):
                            if (corrTypeToString(corr_type[p]) in polsToPlot):
                                pdesc = pb.plot(pfrequencies[p], gamp[p],'%s%s'%(pcolor[p],ampmarkstyle), markersize=markersize,markeredgewidth=markeredgewidth)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,9,chanrangePercent)
                        if (corrTypeToString(corr_type[0]) in polsToPlot):
                            pdesc = pb.plot(frequenciesGHz[index], amplitudeSolutionX,'%s%s'%(p2color[0],bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                            newylimits = recalcYlimitsFreq(chanrange, newylimits, amplitudeSolutionX, sideband,plotrange,xchannels,debug,10,chanrangePercent)
                            pdesc = pb.plot(frequenciesGHz2[index], amplitudeSolution2X, '%s%s'%(p3color[0],bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                            newylimits = recalcYlimitsFreq(chanrange, newylimits, amplitudeSolution2X, sideband,plotrange,xchannels2,debug,11,chanrangePercent)
                        if (nPolarizations == 2):
                           if (corrTypeToString(corr_type[1]) in polsToPlot):
                              pdesc = pb.plot(frequenciesGHz[index], amplitudeSolutionY,'%s%s'%(p2color[1],bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                              newylimits = recalcYlimitsFreq(chanrange, newylimits, amplitudeSolutionY, sideband,plotrange,ychannels,debug,12,chanrangePercent)
                              pdesc = pb.plot(frequenciesGHz2[index], amplitudeSolution2Y, '%s%s'%(p3color[1],bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                              newylimits = recalcYlimitsFreq(chanrange, newylimits, amplitudeSolution2Y, sideband,plotrange,ychannels2,debug,13,chanrangePercent)
                    else:
#                        pb.hold(True)# not available in CASA6, but never needed
                        for p in range(nPolarizations):
                            if (corrTypeToString(corr_type[p]) in polsToPlot):
                                pdesc = pb.plot(pfrequencies[p], gamp[p],'%s%s'%(pcolor[p],ampmarkstyle), markersize=markersize,markeredgewidth=markeredgewidth)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,10,chanrangePercent)
                        if (corrTypeToString(corr_type[0]) in polsToPlot):
                            pdesc = pb.plot(frequenciesGHz[index], amplitudeSolutionX,'%s%s'%(p2color[0],bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                            newylimits = recalcYlimitsFreq(chanrange, newylimits, amplitudeSolutionX, sideband,plotrange,xchannels,debug,14,chanrangePercent)
                        if (nPolarizations == 2):
                           if (corrTypeToString(corr_type[1]) in polsToPlot):
                              pdesc = pb.plot(frequenciesGHz[index], amplitudeSolutionY,'%s%s'%(p2color[1],bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                              newylimits = recalcYlimitsFreq(chanrange, newylimits, amplitudeSolutionY, sideband,plotrange,ychannels,debug,15,chanrangePercent)
                    # endif (bpolyOverlay2)
                else:
                    # we are not overlaying any B or polynomial solutions      'amp vs. freq'
                    if (showflagged):
                        # Also show the flagged data to see where the flags are
#                        pb.hold(True)  # Matches line 2326 for xaxis='chan' # not available in CASA6, but never needed
                        for p in range(nPolarizations):
                          if (corrTypeToString(corr_type[p]) in polsToPlot):
                            if (overlayAntennas or overlayTimes):
                              pdesc1 = pb.plot(pfrequencies[p], gamp[p], '%s'%ampmarkstyles[p], markersize=markersize,markeredgewidth=markeredgewidth)
                              newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,16,chanrangePercent)
                              if (overlayAntennas and overlayTimes==False):
                                  pb.setp(pdesc1, color=overlayColors[xctr])
                              elif (overlayTimes and overlayAntennas==False):
                                  pb.setp(pdesc1, color=overlayColors[mytime])
                              elif (overlayTimes and overlayAntennas): # try to support antenna,time
                                  if (myUniqueTime != []):
                                      pb.setp(pdesc1, color=overlayColors[myUniqueTime])
                                      # The third 'or' below is needed if pol='0' is flagged on antenna 0. -- 2012/10/12 (original spot)
                                      if (p==0 or len(polsToPlot)==1 or myUniqueColor==[]):
                                          myUniqueColor.append(overlayColors[len(myUniqueColor)])
                                      pb.setp(pdesc1, color=myUniqueColor[-1])
                            else:
                                pdesc = pb.plot(pfrequencies[p], gamp[p], '%s%s'%(pcolor[p],ampmarkstyles[p]), markersize=markersize,markeredgewidth=markeredgewidth)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,17,chanrangePercent)
                    else:   # showing only unflagged data    'amp vs. freq'
#                        pb.hold(True) # not available in CASA6, but never needed
                        for p in range(nPolarizations):
                          if (debug):
                              print("p=%d, polsToPlot=%s, len(fieldsToPlot)=%d, len(timerangeList)=%d, myUniqueTime=" % (p,str(polsToPlot),len(fieldsToPlot),len(timerangeList)), myUniqueTime)
                          if (corrTypeToString(corr_type[p]) in polsToPlot):
                            if (len(gamp[p]) == 0):  # Try this on Apr 2, 2012
#                                print "=============== Skipping flagged data on antenna %d = %s" % (xant,antstring)
                                continue
                            if (overlayAntennas or overlayTimes):
                              pdesc = pb.plot(pfrequencies[p], gamp[p], '%s'%ampmarkstyles[p], 
                                              markersize=markersize,markeredgewidth=markeredgewidth)
                              newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,
                                                             plotrange,xchannels,debug,18,chanrangePercent)
#                              print "   antenna=%s, pol=%d, ispw=%d" % (antstring,p,ispw)
                              if (overlayAntennas and overlayTimes==False):
                                  pb.setp(pdesc, color=overlayColors[xctr])
                              elif (overlayTimes and overlayAntennas==False):
                                  pb.setp(pdesc, color=overlayColors[mytime])
                              elif (overlayTimes and overlayAntennas):     #  try to support antenna,time
                                  if (myUniqueTime != []):
                                      pb.setp(pdesc, color=overlayColors[myUniqueTime])
                                      # The third 'or' below is needed if pol='0' is flagged on antenna 0. -- 2012/10/12 (original spot)
                                      if (p==0 or len(polsToPlot)==1 or myUniqueColor==[]):
                                          myUniqueColor.append(overlayColors[len(myUniqueColor)])
                                      if (debug):
                                          print("myUniqueColor = ", myUniqueColor)
                                      pb.setp(pdesc, color=myUniqueColor[-1])
                            elif (overlaySpws):
                                if overlaySpwDistinguish.find('color') >= 0:
                                    mycolor = overlayColors[list(spwsToPlot).index(ispw)+p*len(spwsToPlot)]
                                    if overlaySpwDistinguish.find('width') >= 0:
                                        if overlaySpwDistinguish.find('width2') >= 0:
                                            linewidth = 1+2*list(spwsToPlot).index(ispw)
                                        else:
                                            linewidth = 1+list(spwsToPlot).index(ispw)
                                    else:
                                        linewidth = 1
                                elif overlaySpwDistinguish.find('width') >= 0:
                                    mycolor = [xcolor,ycolor][p]
                                    if overlaySpwDistinguish.find('width2') >= 0:
                                        linewidth = 1+2*list(spwsToPlot).index(ispw)
                                    else:
                                        linewidth = 1+list(spwsToPlot).index(ispw)
                                else:
                                    mycolor = [xcolor,ycolor][p]
                                    linewidth = 1
                                pdesc = pb.plot(pfrequencies[p], gamp[p], '%s'%(ampmarkstyles[0]), lw=linewidth, color=mycolor,
                                                markersize=markersize,markeredgewidth=markeredgewidth)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,
                                                               plotrange,xchannels,debug,-18,chanrangePercent)
#                                print "   antenna=%s, pol=%d, ispw=%d" % (antstring,p,ispw)
                            else: # show unflagged solutions, no overlay
                               if (corrTypeToString(corr_type[p]) in polsToPlot):
                                  # since there is no overlay, don't use dashed line, so zero here ---------v
                                  pdesc = pb.plot(pfrequencies[p], gamp[p], '%s%s'%(pcolor[p],ampmarkstyles[0]),markersize=markersize,markeredgewidth=markeredgewidth)
                                  newylimits = recalcYlimitsFreq(chanrange, newylimits, gamp[p], sideband,plotrange,xchannels,debug,19,chanrangePercent)
                                  if asciiFile:
                                      asciiDesc = open('pol%d_vs_freq.txt' % p,'w')
                                      for i in range(len(pfrequencies[p])):
                                          asciiDesc.write('%f %f\n'%(pfrequencies[p][i], gamp[p][i]))
                                      asciiDesc.close()
#                        print "newylimits for amp = ", newylimits
                        if (sum(xflag)>0):
#                            print "amp: Resetting xaxis frequency range to counteract flagged data"
                            myxrange = np.max(frequencies)-np.min(frequencies)
                            SetNewXLimits([np.min(frequencies)-0.15*myxrange, np.max(frequencies)+0.15*myxrange],loc=5)
                            
                if (1==1 or (xframe in bottomRowFrames) or (xctr+1==len(antennasToPlot) and ispw==spwsToPlot[-1])):
                    # use 1==1 because spw might change between top row and bottom row of frames
                    pb.xlabel(xlabelString, size=mysize)
            # endif (xaxis=='chan' elif xaxis=='freq'  for 'amp')
            if (overlayTimes):
                timeString =''
            else:
                if (len(uniqueTimes) > mytime):
                    timeString = ',  t%d/%d  %s' % (mytime,nUniqueTimes-1,utstring(uniqueTimes[mytime],3))
                    if (scansForUniqueTimes != []):
                        if (scansForUniqueTimes[mytime]>=0):
                            timeString = ',  scan%d  %s' % (scansForUniqueTimes[mytime],utstring(uniqueTimes[mytime],3))
            spwString = buildSpwString(overlaySpws, overlayBasebands,
                                       spwsToPlot, ispw, originalSpw[ispw],
                                       observatoryName, baseband, 
                                       showBasebandNumber)
            if (overlayTimes and len(fieldsToPlot) > 1):
                indices = fstring = ''
                for f in fieldIndicesToPlot:
                    if (f != fieldIndicesToPlot[0]):
                        indices += ','
                        fstring += ','
                    indices += str(uniqueFields[f])
                    if (msFound):
                        fstring += msFields[uniqueFields[f]]
                if (len(fstring) > fstringLimit):
                    fstring = fstring[0:fstringLimit] + '...'
                titleString = "%sspw%s,  fields %s: %s%s" % (antennaString,spwString,
                                                             indices, fstring, timeString)
            else:
                titleString = "%sspw%s,  field %d: %s%s" % (antennaString,spwString,uniqueFields[fieldIndex],
                                                                fieldString,timeString)

            tsize = titlesize-int(len(titleString)/int(maxCharsBeforeReducingTitleFontSize/subplotCols))
            pb.title(titleString, size=tsize)
            if (abs(plotrange[0]) > 0 or abs(plotrange[1]) > 0):
                SetNewXLimits([plotrange[0],plotrange[1]],loc=6)
            else:
                # Here is 1st place where we eliminate white space on right and left edge of the plots: 'amp'
                if (xaxis.find('chan')>=0):
                    SetNewXLimits([channels[0],channels[-1]],loc=7)
                else:
                    if (zoom != 'intersect'):
                        if (overlaySpws or overlayBasebands):
#                            print "frequencyRangeToPlotInBaseband[%d]=" % (bbctr), frequencyRangeToPlotInBaseband[bbctr]
                            SetNewXLimits(frequencyRangeToPlotInBaseband[bbctr],loc=8)
                        else:
                            SetNewXLimits([frequencies[0], frequencies[-1]],loc=9)
                    if (bOverlay):
                        if (xrange2 > myxrange+0.1 and zoom != 'intersect'):
                            TDMisSecond = True
            if (abs(plotrange[2]) > 0 or abs(plotrange[3]) > 0):
                SetNewYLimits([plotrange[2],plotrange[3]],loc=3)

            ResizeFonts(adesc,mysize)
            adesc.xaxis.grid(True,which='major')
            adesc.yaxis.grid(True,which='major')
            pb.ylabel(yAmplitudeLabel, size=mysize)
            pb.subplots_adjust(hspace=myhspace, wspace=mywspace)
            xlim = pb.xlim()
            ylim = pb.ylim()
            myxrange = xlim[1]-xlim[0]
            yrange = ylim[1]-ylim[0]
#            print "amp: ylim, yrange = ",  ylim, yrange
            if (overlayAntennas == False and overlayTimes == False and bOverlay == False and 
                ((overlaySpws == False and overlayBasebands == False) or 
                 (spwctr==spwctrFirstToPlot or overlaySpwDistinguish.find('color')>=0))):
                # draw polarization labels for no overlay, or overlaySpws/overlayBasebands
                x0 = xstartPolLabel
                y0 = ystartPolLabel
                for p in range(nPolarizations):
                   if (corrTypeToString(corr_type[p]) in polsToPlot):
                      if overlaySpwDistinguish.find('color') >= 0:
                          mySpwIndex = list(spwsToPlot).index(ispw)
                          mycolor = overlayColors[mySpwIndex+p*len(spwsToPlot)]
                          pb.text(x0+mySpwIndex*0.04, y0-subplotRows*p*0.03, corrTypeToString(corr_type[p]),
                                  color=mycolor,size=mysize, transform=pb.gca().transAxes)
                      elif spwctr==spwctrFirstToPlot or (not overlaySpws and not overlayBasebands):
                          # no need to plot it more than once in the same position
                          pb.text(x0, y0-subplotRows*p*0.03, corrTypeToString(corr_type[p]),
                                  color=pcolor[p],size=mysize, transform=pb.gca().transAxes)
                      elif debug:
                          print("Not drawing polarization labels because %d != %d" % (spwctr,spwctrFirstToPlot))
                      if (channeldiff > 0):
                          pb.text(x0, ystartMadLabel-0.03*subplotRows*p,
                                  corrTypeToString(corr_type[p])+' MAD = %.4f, St.Dev = %.4f'%(gamp_mad[p]['mad'],gamp_std[p]['std']),
                                  color=pcolor[p],size=mysize, transform=pb.gca().transAxes)
                if (xframe == firstFrame):
                      # draw title including caltable name
                      caltableList = caltableTitle
                      if (bpolyOverlay):
                            caltableList += ', ' + caltable2 + ' (degamp=%d, degphase=%d)'%(nPolyAmp[index]-1,nPolyPhase[index]-1)
                            if (bpolyOverlay2):
                                  caltableList += ', ' + caltable3 + ' (degamp=%d, degphase=%d)'%(nPolyAmp2[index]-1,nPolyPhase2[index]-1)
                      pb.text(xstartTitle, ystartTitle, caltableList, size=titlesize,
                              color='k', transform=pb.gcf().transFigure)

            elif (overlayAntennas==True and xant==antennasToPlot[-1] and bOverlay == False   # ):
                  and overlayTimes==False):  # try to support antenna,time  avoid antenna labels 'amp'
                    # We do this last, because by then, the limits will be stable.
                    x0 = xstartPolLabel
                    y0 = ystartPolLabel
                    # draw polarization labels for overlayAntennas
                    if (corrTypeToString(corr_type[0]) in polsToPlot):
                      if (channeldiff > 0):
                          pb.text(x0, ystartMadLabel-0.03*subplotRows*0,
                                  corrTypeToString(corr_type[0])+' MAD = %.4f, St.Dev = %.4f'%(gamp_mad[0]['mad'],gamp_std[0]['std']),
                                  color=overlayColors[0],size=mysize, transform=pb.gca().transAxes)
                      if (ampmarkstyle.find('-')>=0):
                          pb.text(x0, y0, corrTypeToString(corr_type[0])+' solid', color=overlayColors[0],size=mysize,
                                  transform=pb.gca().transAxes)
                      else:
                          pb.text(x0+0.02, y0, corrTypeToString(corr_type[0]), color=overlayColors[0],size=mysize,
                                  transform=pb.gca().transAxes)
                          pdesc = pb.plot([x0-0.01], [y0], '%sk'%ampmarkstyle, markersize=markersize,
                                          scalex=False,scaley=False, transform=pb.gca().transAxes,markeredgewidth=markeredgewidth)
                    if (len(corr_type) > 1):
                     if (corrTypeToString(corr_type[1]) in polsToPlot):
                      if (channeldiff > 0):
                          pb.text(x0, ystartMadLabel-0.03*subplotRows*1,
                                  corrTypeToString(corr_type[1])+' MAD = %.4f, St.Dev = %.4f'%(gamp_mad[1]['mad'],gamp_std[1]['std']),
                                  color=overlayColors[0],size=mysize, transform=pb.gca().transAxes)
                      if (ampmarkstyle2.find('--')>=0):
                        pb.text(x0, y0-0.03*subplotRows, corrTypeToString(corr_type[1])+' dashed',
                                color=overlayColors[0],size=mysize, transform=pb.gca().transAxes)
                      else:
                        pb.text(x0+0.02, y0-0.03*subplotRows, corrTypeToString(corr_type[1]),
                                color=overlayColors[0],size=mysize, transform=pb.gca().transAxes)
                        pdesc = pb.plot([x0-0.01], [y0-0.03*subplotRows], '%sk'%ampmarkstyle2,
                                        markersize=markersize, scalex=False,scaley=False,markeredgewidth=markeredgewidth)
                    if (xframe == firstFrame):
                        # draw title including caltable name
                        pb.text(xstartTitle, ystartTitle, caltableTitle, size=titlesize, color='k',
                                transform=pb.gcf().transFigure)
                        DrawAntennaNames(msAnt, antennasToPlot, msFound, mysize)
            elif (overlayTimes==True and bOverlay == False
                  and overlayAntennas==False):  # try to support antenna,time
                doneOverlayTime = True  # assumed until proven otherwise in the 'for' loop
                for f in fieldIndicesToPlot:
                  if (len(uniqueTimesPerFieldPerSpw[ispwInCalTable][f]) > 0):
                    if ((uniqueTimes[mytime] < uniqueTimesPerFieldPerSpw[ispwInCalTable][f][-1]-solutionTimeThresholdSeconds) and
                        (uniqueTimes[mytime] < timerangeListTimes[-1])):
                        if (debug):
                            print("-----------Not finished because %.0f < %.0f-%d for fieldIndex=%d and <%.0f" % (uniqueTimes[mytime], uniqueTimesPerFieldPerSpw[ispwInCalTable][f][-1], solutionTimeThresholdSeconds, f, timerangeListTimes[-1]))
                            print("-----------ispwInCalTable=%d, mytime=%d, len(uniqueTimes) = %d" % (ispwInCalTable, mytime, len(uniqueTimes)))
                        doneOverlayTime = False
                if (debug):
                    print("------doneOverlayTime = %s" % (str(doneOverlayTime)))
                if (doneOverlayTime):
                    # either it is the last time of any times in solution, or the last time
                    # in the list of times to plot
                    if (debug):
                        print("****** on last time = %d for last fieldIndex %d  or %d>=%d" % (mytime,fieldIndex,mytime,timerangeList[-1]))
                    mytime = nUniqueTimes-1
                    # We do this last, because by then, the limits will be broad enough and stable.
                    # draw polarization labels
                    DrawPolarizationLabelsForOverlayTime(xstartPolLabel,ystartPolLabel,corr_type,polsToPlot,
                                                         channeldiff,ystartMadLabel,subplotRows,gamp_mad,mysize,
                                                         ampmarkstyle,markersize,ampmarkstyle2, gamp_std)
                    if (xframe == firstFrame):
                        # draw title including caltable name
                        pb.text(xstartTitle, ystartTitle, caltableTitle, size=titlesize,
                                color='k', transform=pb.gcf().transFigure)
                        drawOverlayTimeLegends(xframe,firstFrame,xstartTitle,ystartTitle,caltable,
                                               titlesize,fieldIndicesToPlot,ispwInCalTable,
                                               uniqueTimesPerFieldPerSpw,
                                               timerangeListTimes, solutionTimeThresholdSeconds,
                                               debugSloppyMatch,ystartOverlayLegend,debug,mysize,
                                               fieldsToPlot,myUniqueColor,timeHorizontalSpacing,
                                               fieldIndex,overlayColors, antennaVerticalSpacing,
                                               overlayAntennas, timerangeList, caltableTitle,
                                               mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes)
            elif (overlayAntennas and overlayTimes):  # Oct 23, 2012
                # This will only happen for overlay='antenna,time'
                if (debug):
                    print("In overlay antenna,time case, xframe=%d, firstFrame=%d, mytime=%d, firstTimeMatch=%d, xctr=%d, firstUnflaggedAntennaToPlot=%d" % (xframe,firstFrame,mytime,firstTimeMatch,xctr,firstUnflaggedAntennaToPlot))
#                if (xframe == firstFrame and mytime == 0 and xctr==firstUnflaggedAntennaToPlot and bOverlay==False):  
                if (xframe == firstFrame and mytime == firstTimeMatch and xctr==firstUnflaggedAntennaToPlot and bOverlay==False):  # bug fix on 2015-08-19 for CAS-7820
                    # draw title including caltable name
                    pb.text(xstartTitle, ystartTitle, caltableTitle, size=titlesize, color='k',
                            transform=pb.gcf().transFigure)
                    DrawBottomLegendPageCoords(msName, uniqueTimes[mytime], mysize, figfile)
                # Adding the following 'for' loop on Mar 13, 2013 to support the case of
                # single time range with overlay='antenna,time'
                if (xant==antennasToPlot[-1]):
                    doneOverlayTime = True  # assumed until proven otherwise in the 'for' loop
                    for f in fieldIndicesToPlot:
                        if (len(uniqueTimesPerFieldPerSpw[ispwInCalTable][f]) > 0):
                            if ((uniqueTimes[mytime] < uniqueTimesPerFieldPerSpw[ispwInCalTable][f][-1]-solutionTimeThresholdSeconds) and
                                (uniqueTimes[mytime] < timerangeListTimes[-1])):
                                if (debug):
                                    print("-----------Not finished because mytime=%d, uniqueTimes[%d]=%.0f < %.0f-%d for fieldIndex=%d and <%.0f" % (mytime,mytime,uniqueTimes[mytime], uniqueTimesPerFieldPerSpw[ispwInCalTable][f][-1], solutionTimeThresholdSeconds, f, timerangeListTimes[-1]))
                                    print("-----------ispwInCalTable=%d, mytime=%d, len(uniqueTimes) = %d" % (ispwInCalTable, mytime, len(uniqueTimes)))
                                doneOverlayTime = False
                    if (debug):
                        print("------doneOverlayTime = %s" % (str(doneOverlayTime)))
                    if (doneOverlayTime):
                        if (debug):
                            print("3412: doneOverlayTime=True, drawOverlayTimeLegends()")
                        # This is necessary for the case that no antennas were flagged for the single timerange selected
                        drawOverlayTimeLegends(xframe,firstFrame,xstartTitle,ystartTitle,caltable,titlesize,
                                               fieldIndicesToPlot,ispwInCalTable,uniqueTimesPerFieldPerSpw,
                                               timerangeListTimes, solutionTimeThresholdSeconds,
                                               debugSloppyMatch,ystartOverlayLegend,debug,mysize,
                                               fieldsToPlot,myUniqueColor,timeHorizontalSpacing,
                                               fieldIndex,overlayColors, antennaVerticalSpacing,
                                               overlayAntennas, timerangeList, caltableTitle,
                                               mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes)
                else:
                    if (debug):
                        print("xant=%d  !=  antennasToPlot[-1]=%d" % (xant, antennasToPlot[-1]))
                        
            # Here is 2nd place where we eliminate any white space on the right and left edge of the plots: 'amp'
            # 
            if (abs(plotrange[2]) > 0 or abs(plotrange[3]) > 0):
                SetNewYLimits([plotrange[2],plotrange[3]],loc=4)
            if (plotrange[0]==0 and plotrange[1]==0):
                if (xaxis.find('chan')>=0):
                    SetNewXLimits([channels[0],channels[-1]],loc=10)
                else:
                    if (zoom != 'intersect'):
                        if (overlaySpws or overlayBasebands):
                            SetNewXLimits(frequencyRangeToPlotInBaseband[bbctr],loc=11)
                        else:
                            SetNewXLimits([frequencies[0], frequencies[-1]],loc=12)
                    if (bOverlay):
#                        print "Checking if %f >= %f" % (xrange2, myxrange)
                        if (xrange2 >= myxrange and zoom != 'intersect'):
                            # This is necessary if caltable2=TDM and caltable=FDM
                            SetNewXLimits([frequencies2[0], frequencies2[-1]],loc=13)
                        if (xrange2 > myxrange+0.1 and zoom != 'intersect'):
                            TDMisSecond = True
            else:
                SetNewXLimits([plotrange[0], plotrange[1]],loc=14)
            if (debug): print("done SetNewXLimits")

            # I need the following line for chanrange to work
            if (chanrange[0] != 0 or chanrange[1] != 0 or chanrangePercent is not None):
                SetLimits(plotrange, chanrange, newylimits, channels, frequencies,
                          pfrequencies, ampMin, ampMax, xaxis, pxl, chanrangeSetXrange, 
                          chanrangePercent, loc=105)

            # Finally, draw the atmosphere and FDM windows, if requested.   'amp'
            if ((overlayAntennas==False and overlayTimes==False) or 
                (overlayAntennas==True and overlayTimes==False and xant==antennasToPlot[-1]) or
                (overlayTimes==True and overlayAntennas==False and doneOverlayTime) or
#                (xant==antennasToPlot[-1] and doneOverlayTime and mytime==nUniqueTimes-1)  # support showatm with overlay='antenna,time'
                (overlayTimes and overlayAntennas and  # Aug 5, 2013
#                 xant==antennasToPlot[-1] and doneOverlayTime and mytime==nUniqueTimes-1 
                 xant==antennasToPlot[-1] and doneOverlayTime and mytime==finalTimeMatch # 2015-08-19 for CAS-7820
                 and not drewAtmosphere)  # added on 2014-12-04 to support case of a flagged antenna CAS-7187
                ):
                if (overlayTimes and overlayAntennas and debug):
                    print("xant=%d, antennasToPlot[-1]=%d, doneOverlayTime=%s" % (xant, antennasToPlot[-1], str(doneOverlayTime)))
                if ((showatm or showtsky) and len(atmString) > 0): 
                    DrawAtmosphere(showatm, showtsky, subplotRows, atmString,
                                   mysize, TebbSky, plotrange, xaxis, atmchan,
                                   atmfreq, transmission, subplotCols,
                                   showatmPoints=showatmPoints, xframe=xframe,
                                   channels=channels,mylineno=lineNumber(),
                                   overlaySpws=overlaySpws,
                                   overlayBasebands=overlayBasebands,
                                   drewAtmosphere=drewAtmosphere,loc=203,
                                   showtsys=showtsys, Trx=Trx)
                    if (LO1 is not None):
                        # Now draw the image band
                        DrawAtmosphere(showatm,showtsky, subplotRows, atmString,
                                       mysize, TebbSkyImage, plotrange, xaxis,
                                       atmchanImage, atmfreqImage, transmissionImage,
                                       subplotCols, LO1, xframe, firstFrame, showatmPoints,
                                       channels=channels,mylineno=lineNumber(),
                                       overlaySpws=overlaySpws,
                                       overlayBasebands=overlayBasebands,
                                       drewAtmosphere=drewAtmosphere,loc=204,
                                       showtsys=showtsys, Trx=Trx)
                    drewAtmosphere = True
                if (xaxis.find('freq')>=0 and showfdm and nChannels <= 256):
                    if (debug):  # amplitude section
                        print("calling showFDM(amplitude), ispw=%d, overlayAntennas=%s, overlayTimes=%s, xant=%d, antennasToPlot[-1]=%d, doneOverlayTime=%s" % (ispw, str(overlayAntennas), str(overlayTimes), xant, antennasToPlot[-1], str(doneOverlayTime)))
                    if (tableFormat == 33):
                        showFDM(originalSpw_casa33, chanFreqGHz_casa33, baseband, showBasebandNumber, basebandDict)
                    else:
                        showFDM(originalSpw, chanFreqGHz, baseband, showBasebandNumber, basebandDict)
            else:
                if (overlayTimes and overlayAntennas and debug):
                    print("xant=%d, antennasToPlot[-1]=%d, doneOverlayTime=%s, mytime=%d, finalTimeMatch=%d" % (xant, antennasToPlot[-1], str(doneOverlayTime), mytime, finalTimeMatch))
            if (debug): print("done drawAtmosphere/FDM check")
            if (bOverlay):
                # draw polarization labels for bOverlay
                x0 = xstartPolLabel
                y0 = ystartPolLabel
                for p in range(nPolarizations):
                    if (corrTypeToString(corr_type[p]) in polsToPlot):
                        pb.text(x0, y0-p*0.03*subplotRows, corrTypeToString(corr_type[p])+'-c1',
                                color=pcolor[p],size=mysize,transform=pb.gca().transAxes)
                        pb.text(x0, y0-p*0.03*subplotRows-0.06*subplotRows, corrTypeToString(corr_type[p])+'-c2',
                                color=p2color[p],size=mysize,transform=pb.gca().transAxes)
            if (debug): print("done pol labels")
            if (bpolyOverlay and xaxis.find('freq')>=0):
                # draw polarization labels for bpolyOverlay
                x0 = xstartPolLabel
                y0 = ystartPolLabel
                if (x2color != xcolor):
                      for p in range(nPolarizations):
                          if (corrTypeToString(corr_type[0]) in polsToPlot):
                              pb.text(x0+0.1, y0-p*0.03*subplotRows, corrTypeToString(corr_type[p]), color=p2color[p],
                                      size=mysize,transform=pb.gca().transAxes)
                if (bpolyOverlay2):
                      for p in range(nPolarizations):
                            if (corrTypeToString(corr_type[0]) in polsToPlot):
                                  pb.text(x0+0.2, y0-p*0.03*subplotRows, corrTypeToString(corr_type[p]),
                                          color=p3color[p], size=mysize,transform=pb.gca().transAxes)

            myIndexTime = uniqueTimesPerFieldPerSpw[ispwInCalTable][fieldIndex][-1]
            if (debug): print("running sloppyMatch")
            matched,mymatch = sloppyMatch(myIndexTime,uniqueTimes,
                                          solutionTimeThresholdSeconds,
                                          mytime, scansToPlotPerSpw[ispw],   # task version and au version
                                          scansForUniqueTimes,myprint=debug,
                                          whichone=True)
            if (debug):
                print("1)done sloppyMatch, mytime=%d, scansForUniqueTimes=%s" % (mytime,str(scansForUniqueTimes)))
                print("ispw=%d" % (ispw))
                print("len(scansToPlotPerSpw)=%d" % (len(scansToPlotPerSpw)))
            # The latter condition is needed to support the scans/timeranges parameters.
            if (matched==False and scansForUniqueTimes[mytime] in scansToPlotPerSpw[ispw]):
                print("---------- Did not find %f within %.0f seconds of anything in %s" % (myIndexTime,solutionTimeThresholdSeconds,str(uniqueTimes)))
                print("---------- uniqueTimesPerFieldPerSpw = %s" % (str(uniqueTimesPerFieldPerSpw)))
                print("Try re-running with a smaller solutionTimeThresholdSeconds (currently %f)" % (solutionTimeThresholdSeconds))
                return
            else:
                # we are on the final time to be plotted
#                if (debug): print "on the final time = %d (scan=%d)" % (mytime,scansForUniqueTimes[mytime])
                mytimeTest = mytime==nUniqueTimes-1 # mytime==myIndexTime  # mytime==mymatch
            if ((xframe == 111 and amplitudeWithPhase) or
                # Following case is needed to make subplot=11 to work for: try to support overlay='antenna,time':  amp
                (xframe == lastFrame and overlayTimes and overlayAntennas and
                 xctr+1==len(antennasToPlot) and
#                 mytime+1==len(uniqueTimes) and  # this worked for nspw <= 4
                 mytimeTest and
                 spwctr<len(spwsToPlot))):
                     if (debug):
                         print(":::: xframe=%d  ==  lastFrame=%d,  amplitudeWithPhase=" % (xframe, lastFrame), amplitudeWithPhase)
                         print(":::: xctr+1=%d == len(antennasToPlot)=%d"  % (xctr+1,len(antennasToPlot)))
                         print(":::: mytimeTest = %s (%d==%d)"  % (mytimeTest, mytime, mymatch))
                         print(":::: spwctr=%d < len(spwsToPlot)=%d"  % (spwctr,len(spwsToPlot)))
                     if (len(figfile) > 0):
                           plotfiles.append(makeplot(figfile,msFound,msAnt,
                                            overlayAntennas,pages,pagectr,
                                            density,interactive,antennasToPlot,
                                            spwsToPlot,overlayTimes,overlayBasebands,
                                                     3,xant,ispw,subplot,resample,
                                            debug,
                                            figfileSequential,figfileNumber))
                           figfileNumber += 1

                     donetime = timeUtilities.time()
                     drewAtmosphere = False # needed for CAS-7187 (subplot=11)
                     if (interactive):
                        pb.draw()
#                        myinput = raw_input(":(%.1f sec) Press return for next page (b for backwards, q to quit): "%(donetime-mytimestamp))
                        myinput = raw_input("Press return for next page (b for backwards, q to quit): ")
                     else:
                        myinput = ''
                     skippingSpwMessageSent = 0
                     mytimestamp = timeUtilities.time()
                     if (myinput.find('q') >= 0):
                         showFinalMessage(overlayAntennas, solutionTimeSpread, nUniqueTimes)
                         return(vm)
                     if (myinput.find('b') >= 0):
                         if (pagectr > 0):
                             pagectr -= 1
                         #redisplay the current page by setting ctrs back to the value they had at start of that page
                         xctr = pages[pagectr][PAGE_ANT]
                         spwctr = pages[pagectr][PAGE_SPW]
                         mytime = pages[pagectr][PAGE_TIME]
                         myap = pages[pagectr][PAGE_AP]
                         xant = antennasToPlot[xctr]
                         antstring = buildAntString(xant,msFound,msAnt)
                         ispw = spwsToPlot[spwctr]
#                         print "Returning to [%d,%d,%d,%d]" % (xctr,spwctr,mytime,myap)
                         redisplay = True
                         if (xctr==pages[0][PAGE_ANT] and spwctr==pages[0][PAGE_SPW] and mytime==pages[0][PAGE_TIME] and pages[0][PAGE_AP]==myap):
                           pb.clf()
                           if (debug):
                               print("2)Setting xframe to %d" % xframeStart)
                           xframe = xframeStart
                           myUniqueColor = []
                           continue
                     else:
                         pagectr += 1
                         if (pagectr >= len(pages)):
#                           print "spwctr=%d, yaxis=%s, xframe=%d, lastFrame=%d, xctr+1=%d, len(antennasToPlot)=%d" % (spwctr, yaxis, xframe, lastFrame, xctr+1,len(antennasToPlot))
                           if (xframe == lastFrame and overlayTimes and overlayAntennas and xctr+1==len(antennasToPlot) and
                               yaxis=='amp'): 
                               # I'm not sure why this works, but is needed to fix CAS-7154
                               myspwctr = spwctr+1
                           else:
                               myspwctr = spwctr
                           pages.append([xctr,myspwctr,mytime,1])
                           if (debug):
                               print("amp: appending [%d,%d,%d,%d]" % (xctr,myspwctr,mytime,1))
                           newpage = 0
                     pb.clf()
                     if (debug):
                         print("3)Setting xframe to %d" % xframeStart)
                     xframe = xframeStart
                     myUniqueColor = []
            else:
                if (debug):
                    print("::: Not done page: Not checking whether we need to set xframe=xframeStart")
                    print("::: xframe=%d  ?=  lastFrame=%d,  amplitudeWithPhase=" % (xframe, lastFrame), amplitudeWithPhase)
                    print("::: xctr+1=%d ?= len(antennasToPlot)=%d"  % (xctr+1,len(antennasToPlot)))
                    print(":::: mytimeTest = %s"  % (mytimeTest))
                    print("::: spwctr=%d ?< len(spwsToPlot)=%d"  % (spwctr,len(spwsToPlot)))
#################################################
######### Here is the phase plotting ############
#################################################
        if (yaxis.find('phase')>=0 or amplitudeWithPhase) and doneOverlayTime==False:
            if (channeldiff > 0):
                pchannels = [xchannels,ychannels]  # this is necessary because np.diff reduces nchan by 1
                pfrequencies = [xfrequencies,yfrequencies]  # this is necessary because np.diff reduces nchan by 1
                if (bOverlay):
                    pchannels2 = [xchannels2,ychannels2]  # this is necessary because np.diff reduces nchan by 1
                    pfrequencies2 = [xfrequencies2,yfrequencies2]  # this is necessary because np.diff reduces nchan by 1
            if (overlayTimes == False or mytime==firstTimeMatch):  
              if ((overlaySpws == False and overlayBasebands==False) or spwctr==spwctrFirstToPlot or 
                  spwctr>spwsToPlot[-1] or 
                  (overlayBasebands and amplitudeWithPhase)) :# CAS-6477
                if (overlayAntennas==False or xctr==firstUnflaggedAntennaToPlot
                    or xctr>antennasToPlot[-1]):  # 2012-05-24, to fix the case where all ants flagged on one timerange
                    xframe += 1
                    if (debug): 
                        print("u) incrementing xframe to %d" % xframe)
                    myUniqueColor = []
                    newylimits = [LARGE_POSITIVE, LARGE_NEGATIVE]
                    if (phase != ''):
                        if ((phase[0] != 0 or phase[1] != 0) and amplitudeWithPhase):
                            newylimits = phase
            if (debug):
                print("$$$$$$$$$$$$$$$$$$$$$$$  ready to plot phase on xframe %d" % (xframe))
            adesc = pb.subplot(xframe)
            if (previousSubplot != xframe):
                drewAtmosphere = False
            previousSubplot = xframe
#            pb.hold(overlayAntennas or overlayTimes) # not available in CASA6, but never needed
            gphsx = np.arctan2(np.imag(gplotx),np.real(gplotx))*180.0/math.pi
            if (nPolarizations == 2):
                gphsy = np.arctan2(np.imag(gploty),np.real(gploty))*180.0/math.pi
                if (debug):
                    print("np.shape(gplotx), (gphsx) = ", np.shape(gplotx), np.shape(gphsx))
                    print("np.shape(gplotxy, (gphsy) = ", np.shape(gploty), np.shape(gphsy))
                if (channeldiff>0):
                    if (xaxis == 'chan'):
                        gphs0, newx0, gphs0res, newx0res = channelDifferences(gphsx, pchannels[0], resample)
                        gphs1, newx1, gphs1res, newx1res = channelDifferences(gphsy, pchannels[1], resample)
                        pchannels = [newx0,newx1]
                    else:
                        gphs0, newx0, gphs0res, newx0res = channelDifferences(gphsx, pfrequencies[0], resample)
                        gphs1, newx1, gphs1res, newx1res  = channelDifferences(gphsy, pfrequencies[1], resample)
                        pfrequencies = [newx0,newx1]
                    gphs = [gphs0, gphs1]
                    gphsres = [gphs0res, gphs1res]
                    gphs_mad = [madInfo(gphsres[0],madsigma,edge,ispw,xant,0), madInfo(gphsres[1],madsigma,edge,ispw,xant,1)]
                    gphs_std = [stdInfo(gphsres[0],madsigma,edge,ispw,xant,0), stdInfo(gphsres[1],madsigma,edge,ispw,xant,1)]
                    for p in [0,1]:
                        madstats[Antstring][ispw][mytime][p]['phase'] = gphs_mad[p]['mad']
                        madstats[Antstring][ispw][mytime][p]['phasestd'] = gphs_std[p]['std']
                        if (gphs_mad[p]['nchan'] > 0):
                            checkAbsSum = np.sum(np.abs(gphs[p]))
                            if (checkAbsSum < PHASE_ABS_SUM_THRESHOLD):
                                if (debug): print("%s, Pol %d, spw %d, %s, phs: not printing because abs sum of all values near zero (%f)" % (Antstring, p, ispw, utstring(uniqueTimes[mytime],0), checkAbsSum))
                            else:
                                print("%s, Pol %d, spw %2d, %s, phs: %4d points exceed %.1f sigma (worst=%.2f at chan %d)" % (Antstring, p, ispw, utstring(uniqueTimes[mytime],0), gphs_mad[p]['nchan'], madsigma, gphs_mad[p]['outlierValue'], gphs_mad[p]['outlierChannel']+pchannels[p][0]))
                else:
                    gphs = [gphsx,gphsy]
            else:  # 1-pol
                if (channeldiff>0):
                    if (xaxis == 'chan'):
                        gphs0, newx0, gphs0res, newx0res = channelDifferences(gphsx, pchannels[0], resample)
                        pchannels = [newx0]
                    else:
                        gphs0, newx0, gphs0res, newx0res = channelDifferences(gphsx, pfrequencies[0], resample)
                        pfrequencies = [newx0]
                    gphs = [gphs0]
                    gphsres = [gphs0res]
                    p = 0
                    gphs_mad = [madInfo(gphsres[p], madsigma, edge, ispw,xant,p)]
                    gphs_std = [stdInfo(gphsres[p], madsigma, edge, ispw,xant,p)]
                    madstats[Antstring][ispw][mytime][p]['phase'] = gphs_mad[p]['mad']
                    madstats[Antstring][ispw][mytime][p]['phasestd'] = gphs_mad[p]['std']
                    if (gphs_mad[p]['nchan'] > 0):
                        checkAbsSum = np.sum(np.abs(gphs[p]))
                        if (checkAbsSum < PHASE_ABS_SUM_THRESHOLD):
                            if (debug): print("%s, Pol %d, spw %d, %s, phs: not printing because all values near zero (%f)" % (Antstring, p, ispw, utstring(uniqueTimes[mytime],0), checkAbsSum))
                        else:
                            print("%s, Pol %d, spw %2d, %s, phs: %4d points exceed %.1f sigma (worst=%.2f at chan %d)" % (Antstring, p, ispw, utstring(uniqueTimes[mytime],0), gphs_mad[p]['nchan'], madsigma, gphs_mad[p]['outlierValue'], gphs_mad[p]['outlierChannel']+pchannels[p][0]))
                else:
                    gphs = [gphsx]
            if (bOverlay):
                  if (debug):
                      print("computing phase for second table")
                  gphsx2 = np.arctan2(np.imag(gplotx2),np.real(gplotx2))*180.0/math.pi
                  if (nPolarizations == 2):
                      gphsy2 = np.arctan2(np.imag(gploty2),np.real(gploty2))*180.0/math.pi
                      if (channeldiff>0):
                          if (xaxis == 'chan'):
                              gphs2_0, newx0, gphs2_0res, newx0res = channelDifferences(gphsx2, pchannels2[0], resample)
                              gphs2_1, newx1, gphs2_1res, newx1res  = channelDifferences(gphsy2, pchannels2[1], resample)
                              pchannels2 = [newx0, newx1]
                          else:
                              gphs2_0, newx0, gphs2_0res, newx0res = channelDifferences(gphsx2, pfrequencies2[0], resample)
                              gphs2_1, newx1, gphs2_1res, newx1res = channelDifferences(gphsy2, pfrequencies2[1], resample)
                              pfrequencies2 = [newx0, newx1]
                          gphs2 = [gphs2_0, gphs2_1]
                          gphs2res = [gphs2_0res, gphs2_1res]
                      else:
                          gphs2 = [gphsx2, gphsy2]
                  else:
                      if (channeldiff>0):
                          if (xaxis == 'chan'):
                              gphs2_0, newx0, gphs2_0res, newx0res = channelDifferences(gphsx2, pchannels2[0], resample)
                              pchannels2 = [newx0]
                          else:
                              gphs2_0, newx0, gphs2_0res, newx0res = channelDifferences(gphsx2, pfrequencies2[0], resample)
                              pfrequencies2 = [newx0]
                          gphs2 = [gphs2_0]
                          gphs2res = [gphs2_0res]
                      else:
                          gphs2 = [gphsx2]

                  if (debug):
                      print("bOverlay is FALSE ===========================")
                
            if (xaxis.find('chan')>=0 or len(xfrequencies) < 1):    # 'phase'
#                pb.hold(True) # not available in CASA6, but never needed
                for p in range(nPolarizations):
                  if (corrTypeToString(corr_type[p]) in polsToPlot):
                    if (overlayAntennas or overlayTimes):
                        pdesc = pb.plot(pchannels[p],gphs[p],'%s'%(phasemarkstyles[p]),markersize=markersize,markeredgewidth=markeredgewidth)
                        newylimits =  recalcYlimits(plotrange,newylimits,gphs[p])  # 10/27/2011
                        if (newylimits[1]-newylimits[0] < minPhaseRange):
                            newylimits = [-minPhaseRange,minPhaseRange]
                        if (phase != ''):
                            if ((phase[0] != 0 or phase[1] != 0) and amplitudeWithPhase):
                                newylimits = phase

                        if (overlayAntennas and overlayTimes==False):
                            pb.setp(pdesc, color=overlayColors[xctr])
                        elif (overlayTimes and overlayAntennas==False):
                            pb.setp(pdesc, color=overlayColors[mytime])
                        elif (overlayTimes):   # try to support antenna,time
                            if (myUniqueTime != []):
                                pb.setp(pdesc, color=overlayColors[myUniqueTime])
                                # The third 'or' below is needed if pol='0' is flagged on antenna 0. -- 2012/10/12 (original spot)
                                if (p==0 or len(polsToPlot)==1 or myUniqueColor==[]):
                                    myUniqueColor.append(overlayColors[len(myUniqueColor)])
                                pb.setp(pdesc, color=myUniqueColor[-1])
                    else:
                        pb.plot(pchannels[p],gphs[p],'%s%s'%(pcolor[p],phasemarkstyles[0]), markersize=markersize,markeredgewidth=markeredgewidth)
                        newylimits =  recalcYlimits(plotrange,newylimits,gphs[p]) # 10/27/2011
                        if (newylimits[1]-newylimits[0] < minPhaseRange):
                            newylimits = [-minPhaseRange,minPhaseRange]
                        if (phase != ''):
                            if ((phase[0] != 0 or phase[1] != 0) and amplitudeWithPhase):
                                newylimits = phase
                if (sum(xflag)>0):
#                    print "phase: Resetting xaxis channel range to counteract flagged data"
                    myxrange = np.max(channels)-np.min(channels)
                    SetNewXLimits([np.min(channels)-myxrange/20, np.max(channels)+myxrange/20],loc=15)
                if (xframe in bottomRowFrames or (xctr+1==len(antennasToPlot) and ispw==spwsToPlot[-1])):
                    pb.xlabel("Channels (%d)"%(len(pchannels[p])), size=mysize)
            elif (xaxis.find('freq')>=0):     # 'phase'
                if (bOverlay):
#                      pb.hold(True) # not available in CASA6, but never needed
                      if (debug):
                          print("Preparing to plot phase from %f-%f for pols:" % (xfrequencies[0],xfrequencies[-1]),polsToPlot)
                          print("Preparing to plot phase from %f-%f for pols:" % (pfrequencies[p][0],pfrequencies[p][-1]),polsToPlot)
                          print("Preparing to plot phase from %f-%f for pols:" % (pfrequencies2[p][0],pfrequencies2[p][-1]),polsToPlot)
                      myxrange = np.abs(xfrequencies[0]-xfrequencies[-1])
                      try:
                          xrange2 = np.abs(xfrequencies2[0]-xfrequencies2[-1])
                      except:
                          print("No phase data found in second solution.  Try increasing the solutionTimeThresholdSeconds above %.0f." % (solutionTimeThresholdSeconds))
                          print("If this doesn't work, email the developer (%s)." % (developerEmail))
                          return(vm)
                      if (np.abs(myxrange/xrange2 - 1) > 0.05 + len(xflag)/len(xchannels)):  # 0.0666 is 2000/1875-1
                         # These line widths are optimal for visualizing FDM over TDM
                         width1 = 1
                         width2 = 4
                         # solutions differ in frequency width, so show the narrower one first
                         if (myxrange < xrange2):
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                                if (debug): print("pb.plot 1")
                                pb.plot(pfrequencies[p], gphs[p], '%s%s'%(pcolor[p],phasemarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband,plotrange,xchannels,debug,20,chanrangePercent)
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                                if (debug): print("pb.plot 2")
                                pb.plot(pfrequencies2[p], gphs2[p], '%s%s'%(p2color[p],phasemarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs2[p], sideband,plotrange,xchannels2,debug,21,chanrangePercent)
                         else:
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                                 if (debug): print("pb.plot 3")
                                 pb.plot(pfrequencies2[p], gphs2[p], '%s%s'%(p2color[p],phasemarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                                 newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs2[p], sideband,plotrange,xchannels2,debug,22,chanrangePercent)
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                                 if (debug): print("pb.plot 4")
                                 pb.plot(pfrequencies[p], gphs[p], '%s%s'%(pcolor[p],phasemarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                                 newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband,plotrange,xchannels,debug,23,chanrangePercent)
                      else:
                         width1 = 1
                         width2 = 1
                         # solutions may be different level of smoothing, so plot highest rms first
#                         pb.hold(True) # not available in CASA6, but never needed
                         if (madOfDiff(gphsx) < madOfDiff(gphsx2)):
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                               if (debug): print("pb.plot 5")
                               pb.plot(pfrequencies2[p], gphs2[p], '%s%s'%(p2color[p],phasemarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                               newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs2[p], sideband,plotrange,xchannels2,debug,24,chanrangePercent)
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                               if (debug): print("pb.plot 6")
                               pb.plot(pfrequencies[p], gphs[p], '%s%s'%(pcolor[p],phasemarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                               newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband,plotrange,xchannels,debug,25,chanrangePercent)
                         else:
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                               if (debug): print("pb.plot 7")
                               pb.plot(pfrequencies[p], gphs[p], '%s%s'%(pcolor[p],phasemarkstyle), linewidth=width2, markersize=markersize,markeredgewidth=markeredgewidth)
                               newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband,plotrange,xchannels,debug,26,chanrangePercent)
                           for p in range(nPolarizations):
                             if (corrTypeToString(corr_type[p]) in polsToPlot):
                               if (debug): print("pb.plot 9")
                               pb.plot(pfrequencies2[p], gphs2[p], '%s%s'%(p2color[p],phasemarkstyle), linewidth=width1, markersize=markersize,markeredgewidth=markeredgewidth)
                               newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs2[p], sideband,plotrange,xchannels2,debug,27,chanrangePercent)
                      # must set new limits after plotting  'phase'
                      (y0,y1) = pb.ylim()
                      if (y1-y0 < minPhaseRange):
                            # this must come before defining ticks 
                            SetNewYLimits([-minPhaseRange,minPhaseRange],loc=5)
                      if (zoom=='intersect'):
                          if (myxrange < xrange2):
                              SetNewXLimits([min(xfrequencies[0],xfrequencies[-1])-myxrange*0.1, max(xfrequencies[0],xfrequencies[-1])+myxrange*0.1],loc=16)
                              SetLimits(plotrange, chanrange, newylimits, channels, frequencies,
                                        pfrequencies, ampMin, ampMax, xaxis,pxl,
                                        chanrangeSetXrange, chanrangePercent,loc=106)
                          else:
                              SetNewXLimits([min(xfrequencies2[0],xfrequencies2[-1])-xrange2*0.1, max(xfrequencies2[0],xfrequencies2[-1])+xrange2*0.1],loc=17)
                              SetLimits(plotrange, chanrange, newylimits, channels, frequencies2,
                                        pfrequencies2, ampMin, ampMax, xaxis,pxl,
                                        chanrangeSetXrange, chanrangePercent,loc=107)
                      else:
                          if (myxrange < xrange2):
                              SetLimits(plotrange, chanrange, newylimits, channels, frequencies,
                                        pfrequencies, ampMin, ampMax, xaxis,pxl,
                                        chanrangeSetXrange, chanrangePercent,loc=108)
                          else:
                              SetLimits(plotrange, chanrange, newylimits, channels, frequencies2,
                                        pfrequencies2, ampMin, ampMax, xaxis,pxl,
                                        chanrangeSetXrange, chanrangePercent,loc=109)
                      # draw polarization and spw labels
                      if (xframe == firstFrame):
                          # draw title including caltable name
                          caltableList = 'c1=' + caltable + ', c2=' + caltable2 # + ' (%s)'%(utstring(uniqueTimes2[mytime],3))
                          pb.text(xstartTitle, ystartTitle, caltableList, size=titlesize,
                                  color='k', transform=pb.gcf().transFigure)
                          if (caltable2amplitudeOffset != 0):
                              pb.text(xstartTitle, 0.935, 'c2 amplitude offset = %.3f' % (caltable2amplitudeOffset),
                                      color='k',size=titlesize,transform=pb.gcf().transFigure)
                elif (bpolyOverlay):
                        matches1 = []
                        for tbp in range(len(timesBP)):
                            if (sloppyMatch(uniqueTimes[mytime], timesBP[tbp], solutionTimeThresholdSeconds,
                                            mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes,
                                            myprint=debugSloppyMatch)):
                                matches1.append(tbp)
                        matches1 = np.array(matches1)
#                        print "time matches: matches1 = ", matches1
                        if (len(matches1) < 1):
                            print("No time match found")
                            print("If you are sure the solutions correspond to the same data, you can set solutionTimeThresholdSeconds=%.0f" % (1+np.ceil(np.abs(timesBP[0]-uniqueTimes[mytime]))))
                            return(vm)
#                        matches1 = np.where(np.floor(uniqueTimes[mytime]) == np.floor(np.array(timesBP)))[0]
                        matches2 = np.where(xant == np.array(antennasBP))[0]
                        if (len(matches2) < 1):
                            print("No antenna match found: ", xant, antennasBP)
#                        print "antenna matches: matches2 = ", matches2

                        if (tableFormat == 33):
                            matches3 = np.where(ispw == np.array(cal_desc_idBP))[0]
                            if (len(matches3) < 1):
                                print("No spw match found: %d not in "% (ispw), cal_desc_idBP)
                        else:
                            matches3 = np.where(ispw == np.array(spwBP))[0]
                            if (len(matches3) < 1):
                                print("No spw match found: %d not in " % (ispw), spwBP)
#                        print "spw matches: matches3 = ", matches3

                        matches12 = np.intersect1d(matches1,matches2)
                        if (len(matches12) < 1):
                            print("No match between: ", matches1, matches2)
#                        print "antenna&time matches: matches12 = ", matches12

                        matches = np.intersect1d(matches12, matches3)
                        if (len(matches) < 1):
                            print("No match between: ", matches12, matches3)
#                        print "antenna&time&spw matches: matches = ", matches

                        try:
                            index = matches[0]  # holds the row number of the matching solution in the BPOLY table
                        except:
                            print("No match found for time=%.6f, xant=%d, ispw=%d"  % (uniqueTimes[mytime],xant,ispw))
                            print("antennasBP = ", antennasBP)
                            print("cal_desc_idBP = ", cal_desc_idBP)
                            print("timesBP = ")
                            for i in timesBP:
                                print("%.6f, " % i)
                            return(vm)
#                        print "phase: Using index = %d/%d (mytime=%d), domain=%.3f,%.3f" % (index,len(polynomialPhase),mytime,frequencyLimits[0,index]*1e-9,frequencyLimits[1,index]*1e-9)
                        if (debug): print("BRowNumber = %d, BPolyRowNumber = %d"  % (BRowNumber, index))
                        validDomain = [frequencyLimits[0,index], frequencyLimits[1,index]]
                        cc = calcChebyshev(polynomialPhase[index][0:nPolyPhase[index]], validDomain, frequenciesGHz[index]*1e+9) * 180/math.pi
                        fa = np.array(frequenciesGHz[index])
                        if (xfrequencies[0] < xfrequencies[-1]):
                            matches = np.where(fa>xfrequencies[0])[0]
                            matches2 = np.where(fa<xfrequencies[-1])[0]
                        else:
                            matches = np.where(fa>xfrequencies[-1])[0]
                            matches2 = np.where(fa<xfrequencies[0])[0]
#                        print "xfrequencies[0] = %f, xfrequencies[-1] = %f" % (xfrequencies[0], xfrequencies[-1])
#                        print "len(matches)=%d, len(matches2)=%d" % (len(matches), len(matches2))
#                        print "fa = ", fa
                        mymean = complexMeanDeg(np.array(cc)[matches[0]:matches2[-1]+1])
                        phaseSolutionX = np.mean(gphsx) - mymean + cc

                        cc = calcChebyshev(polynomialPhase[index][nPolyPhase[index]:2*nPolyPhase[index]], validDomain, frequenciesGHz[index]*1e+9) * 180/math.pi
                        if (nPolarizations > 1):
                            if (yfrequencies[0] < yfrequencies[-1]):
                                matches = np.where(fa>yfrequencies[0])[0]
                                matches2 = np.where(fa<yfrequencies[-1])[0]
                            else:
                                matches = np.where(fa>yfrequencies[-1])[0]
                                matches2 = np.where(fa<yfrequencies[0])[0]
                            mymean = complexMeanDeg(np.array(cc)[matches[0]:matches2[-1]+1])
                            phaseSolutionY = np.mean(gphsy) - mymean + cc
                        if (bpolyOverlay2):
                            validDomain = [frequencyLimits2[0,index], frequencyLimits2[1,index]]
                            cc = calcChebyshev(polynomialPhase2[index][0:nPolyPhase2[index]], validDomain,
                                               frequenciesGHz2[index]*1e+9) * 180/math.pi
                            fa = np.array(frequenciesGHz2[index])
                            if (xfrequencies[0] < xfrequencies[-1]):
                                matches = np.where(fa>xfrequencies[0])[0]
                                matches2 = np.where(fa<xfrequencies[-1])[0]
                            else:
                                matches = np.where(fa>xfrequencies[-1])[0]
                                matches2 = np.where(fa<xfrequencies[0])[0]
                            mymean = complexMeanDeg(np.array(cc)[matches[0]:matches2[-1]+1])
                            phaseSolution2X = np.mean(gphsx) + cc - mymean

                            cc = calcChebyshev(polynomialPhase2[index][nPolyPhase2[index]:2*nPolyPhase2[index]],
                                               validDomain, frequenciesGHz2[index]*1e+9) * 180/math.pi
                            if (yfrequencies[0] < yfrequencies[-1]):
                                matches = np.where(fa>yfrequencies[0])[0]
                                matches2 = np.where(fa<yfrequencies[-1])[0]
                            else:
                                matches = np.where(fa>yfrequencies[-1])[0]
                                matches2 = np.where(fa<yfrequencies[0])[0]
                            mymean = complexMeanDeg(np.array(cc)[matches[0]:matches2[-1]+1])
                            phaseSolution2Y = np.mean(gphsy) + cc - mymean
#                            pb.hold(True) # not available in CASA6, but never needed
                            for p in range(nPolarizations):
                                  if (corrTypeToString(corr_type[p]) in polsToPlot):
                                        pb.plot(pfrequencies[p], gphs[p],'%s%s'%(pcolor[p],phasemarkstyle), markersize=markersize,markeredgewidth=markeredgewidth)
                                        newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband,plotrange,xchannels,debug,28,chanrangePercent)
                            if (corrTypeToString(corr_type[0]) in polsToPlot):
                                pb.plot(frequenciesGHz[index],phaseSolutionX,'%s%s'%(x2color,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, phaseSolutionX, sideband,plotrange,xchannels,debug,29,chanrangePercent)
                                pb.plot(frequenciesGHz2[index],phaseSolution2X,'%s%s'%(x3color,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                                newylimits = recalcYlimitsFreq(chanrange, newylimits, phaseSolution2X, sideband,plotrange,xchannels2,debug,30,chanrangePercent)
                            if (nPolarizations == 2):
                               if (corrTypeToString(corr_type[1]) in polsToPlot):
                                  pb.plot(frequenciesGHz[index],phaseSolutionY,'%s%s'%(y2color,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                                  newylimits = recalcYlimitsFreq(chanrange, newylimits, phaseSolutionY, sideband,plotrange,xchannels,debug,31,chanrangePercent)
                                  pb.plot(frequenciesGHz2[index],phaseSolution2Y,'%s%s'%(y3color,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                                  newylimits = recalcYlimitsFreq(chanrange, newylimits, phaseSolution2Y, sideband,plotrange,xchannels2,debug,32,chanrangePercent)
                        else:
#                            pb.hold(True) # not available in CASA6, but never needed
                            for p in range(nPolarizations):
                                if (corrTypeToString(corr_type[p]) in polsToPlot):
                                    pb.plot(pfrequencies[p], gphs[p],'%s%s'%(pcolor[p],phasemarkstyle), markersize=markersize,markeredgewidth=markeredgewidth)
                                    newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband,plotrange,xchannels,debug,33,chanrangePercent)
                            if (corrTypeToString(corr_type[0]) in polsToPlot):
                               pb.plot(frequenciesGHz[index],phaseSolutionX,'%s%s'%(x2color,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                               newylimits = recalcYlimitsFreq(chanrange, newylimits, phaseSolutionX, sideband,plotrange,xchannels,debug,34,chanrangePercent)
                            if (nPolarizations == 2):
                               if (corrTypeToString(corr_type[1]) in polsToPlot):
                                  pb.plot(frequenciesGHz[index],phaseSolutionY,'%s%s'%(y2color,bpolymarkstyle),markeredgewidth=markeredgewidth,markersize=markersize)
                                  newylimits = recalcYlimitsFreq(chanrange, newylimits, phaseSolutionY, sideband,plotrange,xchannels,debug,35,chanrangePercent)
                        # endif (bpolyOverlay2)
                        # Adding the following 4 lines on March 14, 2013
                        (y0,y1) = pb.ylim()
                        if (y1-y0 < minPhaseRange):
                            # this must come before defining ticks 
                            SetNewYLimits([-minPhaseRange,minPhaseRange],loc=6)
                else:
                    # we are not overlaying any B or polynomial solutions   'phase vs. freq'
#                    pb.hold(True) # not available in CASA6, but never needed
                    for p in range(nPolarizations):
                        if (corrTypeToString(corr_type[p]) in polsToPlot):
                            if (overlayAntennas or overlayTimes):
                              pdesc = pb.plot(pfrequencies[p], gphs[p],'%s'%(phasemarkstyles[p]), markersize=markersize,markeredgewidth=markeredgewidth)
                              newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband,plotrange,xchannels,debug,36,chanrangePercent) # Apr 2, 2012
#                              print "Got newylimits = ", newylimits
                              if (overlayAntennas and overlayTimes==False):
                                  pb.setp(pdesc, color=overlayColors[xctr])
                              elif (overlayTimes and overlayAntennas==False):
                                  pb.setp(pdesc, color=overlayColors[mytime])
                              elif (overlayTimes): # try to support antenna,time
                                  if (myUniqueTime != []):
                                      pb.setp(pdesc, color=overlayColors[myUniqueTime])
                                      # The third 'or' below is needed if pol='0' is flagged on antenna 0. -- 2012/10/12 (original spot)
                                      if (p==0 or len(polsToPlot)==1 or myUniqueColor==[]):
                                          myUniqueColor.append(overlayColors[len(myUniqueColor)])
                                      pb.setp(pdesc, color=myUniqueColor[-1])
                            else:
                              pb.plot(pfrequencies[p], gphs[p],'%s%s'%(pcolor[p],phasemarkstyles[0]), markersize=markersize,markeredgewidth=markeredgewidth)
                              newylimits = recalcYlimitsFreq(chanrange, newylimits, gphs[p], sideband, plotrange,xchannels,debug,37,chanrangePercent)
                        if (sum(xflag)>0):
#                            print "phase frame %d: Resetting xaxis frequency range to counteract flagged data" % (xframe)
                            myxrange = np.max(frequencies)-np.min(frequencies)
                            SetNewXLimits([np.min(frequencies)-0.15*myxrange, np.max(frequencies)+0.15*myxrange],loc=18)
                        if (len(gphs[p]) > 0):
                            if (np.max(gphs[p]) < minPhaseRange and np.min(gphs[p]) > -minPhaseRange):
                                SetNewYLimits([-minPhaseRange,minPhaseRange],loc=7)
                #endif bOverlay

                pb.xlabel(xlabelString, size=mysize)
            #endif xaxis='chan'/freq  for 'phase'
            if (overlayTimes):
                timeString =''
            else:
                timeString = ',  t%d/%d  %s' % (mytime,nUniqueTimes-1,utstring(uniqueTimes[mytime],3))
                if (scansForUniqueTimes != []):
                    if (scansForUniqueTimes[mytime]>=0):
                        timeString = ',  scan%d  %s' % (scansForUniqueTimes[mytime],utstring(uniqueTimes[mytime],3))
            spwString = buildSpwString(overlaySpws, overlayBasebands,
                                       spwsToPlot, ispw, originalSpw[ispw],
                                       observatoryName, baseband, 
                                       showBasebandNumber)
            titleString = "%sspw%s,  field %d: %s%s" % (antennaString,
                      spwString,uniqueFields[fieldIndex],fieldString,timeString)
            tsize = titlesize-int(len(titleString)/int(maxCharsBeforeReducingTitleFontSize/subplotCols))
            pb.title(titleString, size=tsize)
            if (abs(plotrange[0]) > 0 or abs(plotrange[1]) > 0):
                SetNewXLimits([plotrange[0],plotrange[1]],loc=19)

            # Here is 1st place where we eliminate any white space on the right and left edge of the plots: 'phase'
            else:
                if (xaxis.find('chan')>=0):
                    SetNewXLimits([channels[0],channels[-1]],loc=20)
                else:
                    if (zoom != 'intersect'):
                        if (overlaySpws or overlayBasebands):
                            SetNewXLimits(frequencyRangeToPlotInBaseband[bbctr],loc=21)
                        else:
                            SetNewXLimits([frequencies[0], frequencies[-1]],loc=22)
                    if (bOverlay):
                        if (xrange2 > myxrange+0.1 and zoom != 'intersect'):
                            TDMisSecond = True

            if (abs(plotrange[2]) > 0 or abs(plotrange[3]) > 0):
                if (amplitudeWithPhase == False or phase == ''):
                    SetNewYLimits([plotrange[2],plotrange[3]],loc=8)
            if (amplitudeWithPhase and phase != ''):
                if (phase[0] != 0 or phase[1] != 0):
                    SetNewYLimits(phase,loc=9)
                

            (y0,y1) = pb.ylim()
            if (y1-y0 < minPhaseRange):
                  # this must come before defining ticks
                  SetNewYLimits([-minPhaseRange,minPhaseRange],loc=10)
                  SetNewYLimits(newylimits,loc=11)  # added 10/2/2012 for the case of only 1 data point
            if (amplitudeWithPhase and phase != ''):
                if (phase[0] != 0 or phase[1] != 0):
                    SetNewYLimits(phase,loc=12)
            (y0,y1) = pb.ylim()
            ResizeFonts(adesc,mysize)
            adesc.xaxis.grid(True,which='major')
            adesc.yaxis.grid(True,which='major')
            pb.ylabel(yPhaseLabel, size=mysize)
            pb.subplots_adjust(hspace=myhspace, wspace=mywspace)
            ylim = pb.ylim()
            xlim = pb.xlim()
            myxrange = xlim[1]-xlim[0]
            yrange = ylim[1]-ylim[0]
#            print "phase: ylim, yrange = ",  ylim, yrange
            myap = 0
            if (overlayAntennas == False and overlayTimes == False and bOverlay == False and
                ((overlaySpws == False and overlayBasebands == False) or spwctr==spwctrFirstToPlot)):
                # draw polarization labels for no overlay
                x0 = xstartPolLabel
                y0 = ystartPolLabel
                for p in range(nPolarizations):
                      if (corrTypeToString(corr_type[p]) in polsToPlot):
                          pb.text(x0, y0-0.03*subplotRows*p, corrTypeToString(corr_type[p]), color=pcolor[p],
                                  size=mysize, transform=pb.gca().transAxes)
                          if (channeldiff > 0):
                              pb.text(x0, ystartMadLabel-0.03*subplotRows*p,
                                      corrTypeToString(corr_type[p])+' MAD = %.4f, St.Dev = %.4f'%(gphs_mad[p]['mad'],gphs_std[p]['std']),
                                      color=pcolor[p],size=mysize, transform=pb.gca().transAxes)
                if (xframe == firstFrame):
                      # draw title including caltable name
                      caltableList = caltableTitle
                      if (bpolyOverlay):
                            caltableList += ', ' + caltable2 + ' (degamp=%d, degphase=%d)'%(nPolyAmp[index]-1,nPolyPhase[index]-1)
                            if (bpolyOverlay2):
                                  caltableList += ', ' + caltable3 + ' (degamp=%d, degphase=%d)'%(nPolyAmp2[index]-1,nPolyPhase2[index]-1)
                      pb.text(xstartTitle, ystartTitle, caltableList, size=titlesize,
                              color='k', transform=pb.gcf().transFigure)
            elif (overlayAntennas==True and xant==antennasToPlot[-1] and bOverlay==False  # ):
                  and overlayTimes==False):  # try to support antenna,time   avoid antenna labels 'phase'
                # We do this last, because by then, the limits will be stable.
                # draw polarization labels for overlayAntennas
                x0 = xstartPolLabel
                y0 = ystartPolLabel
                if (corrTypeToString(corr_type[0]) in polsToPlot):
                    if (channeldiff > 0):
                        pb.text(x0, ystartMadLabel-0.03*subplotRows*p,
                                corrTypeToString(corr_type[p])+' MAD = %.4f, St.Dev = %.4f'%(gphs_mad[p]['mad'],gphs_std[p]['std']),
                                color=overlayColors[0], size=mysize, transform=pb.gca().transAxes)
                    if (phasemarkstyle.find('-')>=0):
                        pb.text(x0, y0-0.03*subplotRows*0, corrTypeToString(corr_type[0])+' solid', color=overlayColors[0],
                                fontsize=mysize, transform=pb.gca().transAxes)
                    else:
                        pb.text(x0+0.02, y0-0.03*subplotRows*0, corrTypeToString(corr_type[0]), color=overlayColors[0],
                                fontsize=mysize, transform=pb.gca().transAxes)
                        pdesc = pb.plot([x0], [y0+0.015-0*0.03*subplotRows], '%sk'%phasemarkstyle, markersize=markersize,
                                        scalex=False,scaley=False, transform=pb.gca().transAxes,markeredgewidth=markeredgewidth)
                if (len(corr_type) > 1):
                  if (corrTypeToString(corr_type[1]) in polsToPlot):
                    if (channeldiff > 0):
                        pb.text(x0, ystartMadLabel-0.03*subplotRows*p,
                                corrTypeToString(corr_type[p])+' MAD = %.4f, St.Dev = %.4f'%(gphs_mad[p]['mad'],gphs_std[p]['std']),
                                color=overlayColors[0], size=mysize, transform=pb.gca().transAxes)
                    if (phasemarkstyle2.find('--')>=0):
                        pb.text(x0, y0-0.03*subplotRows*1, corrTypeToString(corr_type[1])+' dashed', color=overlayColors[0],
                                fontsize=mysize, transform=pb.gca().transAxes)
                    else:
                        pb.text(x0+0.02, y0-0.03*subplotRows*1, corrTypeToString(corr_type[1]), color=overlayColors[0],
                                fontsize=mysize, transform=pb.gca().transAxes)
                        pdesc = pb.plot([x0], [y0+0.015*subplotRows-0.03*subplotRows*1],'%sk'%phasemarkstyle2, markersize=markersize,
                                        scalex=False,scaley=False, transform=pb.gca().transAxes,markeredgewidth=markeredgewidth)
                if (xframe == firstFrame):
                    # draw title including caltable name
                    pb.text(xstartTitle, ystartTitle, caltableTitle, size=titlesize, color='k',
                            transform=pb.gcf().transFigure)
                    DrawAntennaNames(msAnt, antennasToPlot, msFound, mysize)
            elif (overlayTimes==True and bOverlay == False 
                  and overlayAntennas==False):  # try to support antenna,time
                doneOverlayTime = True # assumed until proven otherwise in the 'for' loop
                for f in fieldIndicesToPlot:
                    if (uniqueTimes[mytime] < uniqueTimesPerFieldPerSpw[ispwInCalTable][f][-1]-solutionTimeThresholdSeconds and
                        uniqueTimes[mytime] < timerangeListTimes[-1]):
                        doneOverlayTime = False
                if (debug):
                    print("------doneOverlayTime = %s" % (str(doneOverlayTime)))
                if (doneOverlayTime):
                    # either it is the last time of any times in solution, or the last time
                    # in the list of times to plot
                    mytime = nUniqueTimes-1
                    # draw polarization labels for overlayTimes
                    # We do this last, because by then, the limits will be broad enough and stable.
                    x0 = xstartPolLabel
                    y0 = ystartPolLabel
                    if (corrTypeToString(corr_type[0]) in polsToPlot):
                      if (channeldiff > 0):
                          p = 0
                          pb.text(x0, ystartMadLabel-0.03*subplotRows*p,
                                  corrTypeToString(corr_type[p])+' MAD = %.4f, St.Dev = %.4f'%(gphs_mad[p]['mad'],gphs_std[p]['std']),
                                  color='k', size=mysize, transform=pb.gca().transAxes)
                      if (phasemarkstyle.find('-')>=0):
                          pb.text(x0, y0, corrTypeToString(corr_type[0])+' solid', color='k',
                                  fontsize=mysize, transform=pb.gca().transAxes)
                      else:
                          pb.text(x0+0.02, y0, corrTypeToString(corr_type[0]), color='k',
                                  fontsize=mysize, transform=pb.gca().transAxes)
                          pdesc = pb.plot([x0], [y0+0.015*subplotRows], '%sk'%phasemarkstyle, markersize=markersize,
                                          scalex=False,scaley=False, transform=pb.gca().transAxes,markeredgewidth=markeredgewidth)
                    if (len(corr_type) > 1):
                      if (corrTypeToString(corr_type[1]) in polsToPlot):
                        if (channeldiff > 0):
                            p = 1
                            pb.text(x0, ystartMadLabel-0.03*subplotRows*p,
                                    corrTypeToString(corr_type[p])+' MAD = %.4f, St.Dev = %.4f'%(gphs_mad[p]['mad'],gphs_std[p]['std']),
                                    color='k', size=mysize, transform=pb.gca().transAxes)
                        if (phasemarkstyle2.find('--')>=0):
                            pb.text(x0, y0-0.03*subplotRows, corrTypeToString(corr_type[1])+' dashed',
                                    color='k',fontsize=mysize, transform=pb.gca().transAxes)
                        else:
                            pb.text(x0+0.02, y0-0.03*subplotRows, corrTypeToString(corr_type[1]),
                                    color='k', fontsize=mysize, transform=pb.gca().transAxes)
                            pdesc = pb.plot([x0], [y0+0.015*subplotRows-0.03*subplotRows], '%sk'%phasemarkstyle2,
                                            markersize=markersize, scalex=False,scaley=False, transform=pb.gca().transAxes,markeredgewidth=markeredgewidth)
                    if (xframe == firstFrame):
                        # draw title including caltable name
                        pb.text(xstartTitle, ystartTitle, caltableTitle, size=titlesize, color='k',
                                transform=pb.gcf().transFigure)
                        drawOverlayTimeLegends(xframe,firstFrame,xstartTitle,ystartTitle,caltable,
                                               titlesize,fieldIndicesToPlot,ispwInCalTable,
                                               uniqueTimesPerFieldPerSpw,
                                               timerangeListTimes, solutionTimeThresholdSeconds,
                                               debugSloppyMatch,ystartOverlayLegend,debug,mysize,
                                               fieldsToPlot,myUniqueColor,timeHorizontalSpacing,
                                               fieldIndex,overlayColors, antennaVerticalSpacing,
                                               overlayAntennas, timerangeList, caltableTitle,
                                               mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes)

            elif (overlayAntennas and overlayTimes):  # Oct 23, 2012
                # This will only happen for: try to support overlay='antenna,time' for 'phase'
                if (xframe == firstFrame and mytime==0 and xctr==firstUnflaggedAntennaToPlot and bOverlay==False):
                    # draw title including caltable name
                    pb.text(xstartTitle, ystartTitle, caltableTitle, size=titlesize, color='k',
                            transform=pb.gcf().transFigure)
                    DrawBottomLegendPageCoords(msName, uniqueTimes[mytime], mysize, figfile)
                
            #endif (overlayAntennas == False and overlayTimes == False and bOverlay == False)
            
            # Here is 2nd place where we eliminate any white space on the right and left edge of the plots: 'phase'
            if (abs(plotrange[2]) > 0 or abs(plotrange[3]) > 0):
                if (amplitudeWithPhase == False or phase == ''):
                    SetNewYLimits([plotrange[2],plotrange[3]],loc=13)
            if (phase != '' and amplitudeWithPhase):
                if (phase[0] != 0 or phase[1] != 0):
                    SetNewYLimits(phase,loc=14)
            if (plotrange[0]==0 and plotrange[1]==0):
                if (xaxis.find('chan')>=0):
                    SetNewXLimits([channels[0],channels[-1]],loc=23)
                else:
                    if (zoom != 'intersect'):
                        if (overlaySpws or overlayBasebands):
                            SetNewXLimits(frequencyRangeToPlotInBaseband[bbctr],loc=24)
                        else:
                            SetNewXLimits([frequencies[0], frequencies[-1]],loc=25)
                    if (bOverlay):
                        if (xrange2 >= myxrange and zoom != 'intersect'):
                            # This is necessary if caltable2=TDM and caltable=FDM
                            SetNewXLimits([frequencies2[0], frequencies2[-1]],loc=26)
                        if (xrange2 > myxrange+0.1 and zoom != 'intersect'):
                            TDMisSecond = True
            else:
                SetNewXLimits([plotrange[0], plotrange[1]],loc=27)

            # I need the following line for chanrange to work
            if (chanrange[0] != 0 or chanrange[1] != 0):
                SetLimits(plotrange, chanrange, newylimits, channels, frequencies,
                          pfrequencies, ampMin, ampMax, xaxis,pxl, chanrangeSetXrange, chanrangePercent,loc=110)

            # Finally, draw the atmosphere and FDM windows, if requested. ' phase'
            if ((overlayAntennas==False and overlayTimes==False) or 
                (overlayAntennas==True and overlayTimes==False and xant==antennasToPlot[-1]) or
                (overlayTimes==True and overlayAntennas==False and doneOverlayTime) or
                (xant==antennasToPlot[-1] and doneOverlayTime)
                ):
                if ((showatm or showtsky) and len(atmString)>0):
                    DrawAtmosphere(showatm, showtsky, subplotRows, atmString,
                                   mysize, TebbSky, plotrange, xaxis, atmchan,
                                   atmfreq, transmission, subplotCols,
                                   showatmPoints=showatmPoints, xframe=xframe,
                                   channels=channels, mylineno=lineNumber(),
                                   overlaySpws=overlaySpws,
                                   overlayBasebands=overlayBasebands,
                                   drewAtmosphere=drewAtmosphere,loc=205,
                                   showtsys=showtsys, Trx=Trx)
                    if (LO1 is not None):
                        DrawAtmosphere(showatm,showtsky, subplotRows, atmString,
                                mysize, TebbSky, plotrange, xaxis, atmchanImage,
                                atmfreqImage, transmissionImage, subplotCols,
                                LO1, xframe, firstFrame, showatmPoints,
                                channels=channels, mylineno=lineNumber(),
                                       overlaySpws=overlaySpws,
                                       overlayBasebands=overlayBasebands,
                                       drewAtmosphere=True, loc=206,showtsys=showtsys, Trx=Trx)
                    drewAtmosphere = True
            
                if (xaxis.find('freq')>=0 and showfdm and nChannels <= 256):
                    if debug or True: # phase section
                        print("calling showFDM(phase), ispw=%d, overlayAntennas=%s, overlayTimes=%s, xant=%d, antennasToPlot[-1]=%d, doneOverlayTime=%s" % (ispw, str(overlayAntennas), str(overlayTimes), xant, antennasToPlot[-1], str(doneOverlayTime)))

                    if (tableFormat == 33):
                        showFDM(originalSpw_casa33, chanFreqGHz_casa33, baseband, showBasebandNumber, basebandDict)              # 'phase'
                    else:
                        showFDM(originalSpw, chanFreqGHz, baseband, showBasebandNumber, basebandDict)

            if (bOverlay):
                # draw polarization labels for bOverlay
                x0 = xstartPolLabel
                y0 = ystartPolLabel
                for p in range(nPolarizations):
                    if (corrTypeToString(corr_type[p]) in polsToPlot):
                        pb.text(x0, y0-p*0.03*subplotRows, corrTypeToString(corr_type[p])+'-c1',
                                color=pcolor[p],size=mysize,transform=pb.gca().transAxes)
                        pb.text(x0, y0-(p*0.03+0.06)*subplotRows, corrTypeToString(corr_type[p])+'-c2',
                                color=p2color[p],size=mysize, transform=pb.gca().transAxes)
            if (bpolyOverlay and xaxis.find('freq')>=0):
                # draw polarization labels for bpolyOverlay
                x0 = xstartPolLabel
                y0 = ystartPolLabel
                if (xcolor != x2color):
                    for p in range(nPolarizations):
                        if (corrTypeToString(corr_type[p]) in polsToPlot):
                            pb.text(x0+0.1, y0-p*0.03*subplotRows, corrTypeToString(corr_type[p]), color=p2color[p],
                                    size=mysize, transform=pb.gca().transAxes)
                if (bpolyOverlay2):
                    for p in range(nPolarizations):
                          if (corrTypeToString(corr_type[p]) in polsToPlot):
                                pb.text(x0+0.2, y0-p*0.03*subplotRows, corrTypeToString(corr_type[p]), color=p3color[p],
                                        size=mysize, transform=pb.gca().transAxes)

        # endif (yaxis='phase')

        redisplay = False

        if (xframe == lastFrame):
          if (debug):
              print("*** mytime+1=%d,  nUniqueTimes=%d, timerangeList[-1]=%d, doneOverlayTime=%s" % (mytime+1, nUniqueTimes,timerangeList[-1],doneOverlayTime))
              print("*** xant=%d, antennasToPlot[-1]=%d, overlayAntennas=%s, overlayTimes=%s" % (xant,antennasToPlot[-1],overlayAntennas,overlayTimes))
              print("*** xframe=%d, lastFrame=%d, xctr=%d, pagectr=%d, spwctr=%d, len(antennasToPlot)=%d, len(spwsToPlot)=%d" % (xframe,lastFrame,xctr,pagectr,spwctr,len(antennasToPlot), len(spwsToPlot)))
          myIndexTime = uniqueTimesPerFieldPerSpw[ispwInCalTable][fieldIndex][-1]
          matched,mymatch = sloppyMatch(myIndexTime,uniqueTimes,solutionTimeThresholdSeconds,
                                        mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes, # au version
#                                        mytime, scansToPlot, scansForUniqueTimes, # task version
                                        whichone=True)
          # The latter condition is needed to support the scans/timeranges parameters.
          if (matched==False and scansForUniqueTimes[mytime] in scansToPlotPerSpw[ispw]):
              print("-------- Did not find %f in %s" % (myIndexTime,str(uniqueTimes)))
              print("Try re-running with a smaller solutionTimeThresholdSeconds (currently %f)" % (solutionTimeThresholdSeconds))
              return
          else:
              # we are on the final time to be plotted
              if (debug):
#                  print "xframe==lastFrame: on the final time = %d (scan=%d)" % (mytime,scansForUniqueTimes[mytime])
                  print("spwctr=%d  len(spwsToPlot)-1=%d, spwsToPlot=" % (spwctr,len(spwsToPlot)-1), spwsToPlot)
                  
              mytimeTest = mytime==nUniqueTimes-1 # mytime==myIndexTime  # mytime==mymatch
          if (debug):
              print("mytimeTest = %s" % (mytimeTest))
#              print "len(scansToPlotPerSpw)=%d, ispw=%d, len(scansForUniqueTimes)=%d, mytime=%d" % (len(scansToPlotPerSpw),ispw,len(scansForUniqueTimes),mytime)
#              if (len(scansToPlotPerSpw)>ispw):
#                  print "len(scansToPlotPerSpw[ispw]) = %d" % (len(scansToPlotPerSpw[ispw]))
          if (scansForUniqueTimes == []):
              # old 3.3 cal tables will land here
              scanTest = False
              scanTest2 = False
          else:
              if (debug):
                  print("ispw=%d len(scansToPlotPerSpw[ispw])=%d   mytime=%d, len(scansForUniqueTimes)=%d, scansForUniqueTimes[%d]=%d" % (ispw,len(scansToPlotPerSpw[ispw]),mytime,len(scansForUniqueTimes),mytime,scansForUniqueTimes[mytime]))
                  print("scansToPlotPerSpw = ", scansToPlotPerSpw)
              if (len(scansToPlotPerSpw[ispw]) == 0):
                  scanTest = False
              else:
                  scanTest = (scansToPlotPerSpw[ispw][-1]==scansForUniqueTimes[mytime])
              highestSpwIndexInSpwsToPlotThatHasCurrentScan = \
                  computeHighestSpwIndexInSpwsToPlotThatHasCurrentScan(spwsToPlot, scansToPlotPerSpw, scansForUniqueTimes[mytime])
              if (highestSpwIndexInSpwsToPlotThatHasCurrentScan == -1):
                  scanTest2 = False
              else:
                  scanTest2 = (spwctr == highestSpwIndexInSpwsToPlotThatHasCurrentScan)
          if ((overlayAntennas==False and overlayTimes==False and overlaySpws==False and overlayBasebands==False)
              # either it is the last time of any, or the last time in the list of times to plot
              or (overlayAntennas==False and overlaySpws==False and overlayBasebands==False and (mytime+1==nUniqueTimes or mytime == timerangeList[-1])) # or mytimeTest)) # mytimeTest removed on July 25,2013
              or (xant==antennasToPlot[-1] and overlayAntennas==True and overlayTimes==False and overlaySpws==False and overlayBasebands==False)
              # inserted mytimeTest below (on Sep 16, 2013 for spectral scan dataset) but breaks VLA
#              or ((mytimeTest or spwctr==len(spwsToPlot)-1) and (overlaySpws or overlayBasebands) and overlayAntennas==False and overlayTimes==False)
              # The following case is needed to prevent frame=225 in test86 (spectral scan dataset with overlay='spw') 
              #   and the lack of showing of 7 of 8 of the spws in final frame of test61.  scanTest2 matches both cases.
              or (scanTest and scanTest2 and overlaySpws and overlayAntennas==False and overlayTimes==False)
              or ((spwctr==len(spwsToPlot)-1) and (overlayBasebands or overlaySpws) and overlayAntennas==False and overlayTimes==False)
              # following case is needed for scans parameter with overlay='time'
              or (overlayTimes and scanTest and overlayAntennas==False)
              # Following case is needed to make subplot=11 to work for: try to support overlay='antenna,time' :  'phase'
              or (xframe == lastFrame and overlayTimes and overlayAntennas and
                  xctr+1==len(antennasToPlot) and
                  mytimeTest and
                  spwctr<len(spwsToPlot))
              or (doneOverlayTime and overlayTimes==True
                  and overlayAntennas==False 
                  )
              ):
            if (debug): 
                print("^^^^^^^^^^^^^ ispw=%d, mytime=%d, len(uniqueTimes)=%d, nUniqueTimes=%d, mytimeTest=%s" % (ispw, mytime,len(uniqueTimes),nUniqueTimes,mytimeTest))
                print("              xctr+1=%d, len(antennasToPlot)=%d, doneOverlayTime=%s, scanTest=%s, scanTest2=%s" % (xctr+1, len(antennasToPlot), doneOverlayTime, scanTest, scanTest2))
            DrawBottomLegendPageCoords(msName, uniqueTimes[mytime], mysize, figfile)

            # added len(pages)>0 on July 30, 2013 to prevent crash when called with single
            # antenna and subplot=11 and all solutions flagged.
            if (len(figfile) > 0 and len(pages) > 0):
                plotfiles.append(makeplot(figfile,msFound,msAnt,
                                          overlayAntennas,pages,pagectr,
                                          density,interactive,antennasToPlot,
                                          spwsToPlot,overlayTimes,overlayBasebands,
                                          4,xant,ispw,subplot,resample,
                                          debug,
                                          figfileSequential,figfileNumber))
                figfileNumber += 1

            myinput = ''
            donetime = timeUtilities.time()
            drewAtmosphere = False # needed for CAS-7187 (subplot=11)
            if (interactive):
                pb.draw()
#                myinput = raw_input("(%.1f sec) Press return for next screen (b for backwards, q to quit): "%(donetime-mytimestamp))
                myinput = raw_input("*Press return for next page (b for backwards, q to quit): ")
            else:
                myinput = ''
            skippingSpwMessageSent = 0
            mytimestamp = timeUtilities.time()
            if (myinput.find('q') >= 0):
                mytime = len(uniqueTimes)
                spwctr = len(spwsToPlot)
                xctr = len(antennasToPlot)
                bbctr = len(spwsToPlotInBaseband)
                break
            if (debug):
                print("4)Setting xframe to %d" % xframeStart)
            xframe = xframeStart
            myUniqueColor = []
            pb.subplots_adjust(hspace=myhspace, wspace=mywspace)
            if (myinput.find('b') >= 0):
                if (pagectr > 0):
                    pagectr -= 1
                #redisplay the current page by setting ctrs back to the value they had at start of that page
                xctr = pages[pagectr][PAGE_ANT]
                spwctr = pages[pagectr][PAGE_SPW]
                mytime = pages[pagectr][PAGE_TIME]
                myap = pages[pagectr][PAGE_AP]
                xant = antennasToPlot[xctr]
                antstring = buildAntString(xant,msFound,msAnt)
                ispw = spwsToPlot[spwctr]
#                print "Returning to [%d,%d,%d,%d]" % (xctr,spwctr,mytime,myap)
                redisplay = True
            else:
                pagectr += 1
                if (pagectr >= len(pages)):
                    newpage = 1
                else:
                    newpage = 0
            if (overlayTimes==True and
                sloppyMatch(uniqueTimesPerFieldPerSpw[ispwInCalTable][fieldIndex][-1],
                            uniqueTimes[mytime],solutionTimeThresholdSeconds,
                            mytime, scansToPlotPerSpw[ispw], scansForUniqueTimes,
                            myprint=debugSloppyMatch)):
                # be sure to avoid any more loops through mytime which will cause 'b' button to fail
                mytime = nUniqueTimes
          else:
              if (debug):
                  print("Not going to new page, uniqueTimes[mytime]=%.8f, uniqueTimesPerFieldPerSpw[ispwInCalTable][fieldIndex=%d][-1]=%.8f" % (uniqueTimes[mytime], fieldIndex, uniqueTimesPerFieldPerSpw[ispwInCalTable][fieldIndex][-1]))
                  print("spwctr=%d ?== (len(spwsToPlot)-1)=%d" % (spwctr,len(spwsToPlot)-1))
        if (redisplay == False):
            if ((overlayAntennas and xctr+1 >= len(antennasToPlot)) or
                ((overlaySpws or overlayBasebands) and spwctr+1 >= len(spwsToPlot)) or                
                (overlayAntennas==False and overlaySpws==False and overlayBasebands==False)):
                mytime += 1
#                print " 0004  incrementing mytime to ", mytime
                if (debug):
                    print("AT BOTTOM OF LOOP: Incrementing mytime to %d (nUniqueTimes=%d), setting firstUnflaggedAntennaToPlot to 0" % (mytime,nUniqueTimes))
                firstUnflaggedAntennaToPlot = 0  # try this
                if (debug):
                    print("AT BOTTOM OF LOOP: Setting firstUnflaggedAntennatoPlot = 0")
                doneOverlayTime = False  # added on 08-nov-2012
                if (overlayBasebands and (uniqueScanNumbers == sorted(scansToPlot))): 
                    if (debug): print("Breaking because scans not specified")
                    break
      # end of while(mytime) loop  endwhile mytime
      if (redisplay == False):
          spwctr += 1
          if (debug):
              print("---------------------------------------- Incrementing spwctr to %d, spwsToPlot=" % (spwctr), spwsToPlot)
              if (spwctr < len(spwsToPlot)):
                  print("---------------------------------------- ispw = %d" % (spwsToPlot[spwctr]))
              else:
                  print("---------------------------------------- done the spws in this baseband (%d)" % (baseband))
#      else:
#          print "redisplay = True"
     # end of while(spwctr) loop
     if (debug): print("B)incrementing bbctr to %d" % (bbctr+1))
     bbctr += 1
    # end of while(bbctr or spwctr) loop  endwhile(bbctr or spwctr)
    if (debug): print("B)finalSpwWasFlagged = %s, len(figfile)=%d" % (finalSpwWasFlagged,len(figfile)))
    if (xant >= antennasToPlot[-1] and xframe != xframeStart):
        # this is the last antenna, so make a final plot
        if (len(figfile) > 0):
              plotfiles.append(makeplot(figfile,msFound,msAnt,overlayAntennas,
                                        pages,pagectr,density,interactive,antennasToPlot,
                                        spwsToPlot,overlayTimes,overlayBasebands,
                                        5,xant,ispw,subplot,resample,debug,
                                        figfileSequential,figfileNumber))
              figfileNumber += 1
    if (redisplay == False):
        xctr += 1
        if (debug):
            print("---------------------------------------- Incrementing xctr to %d (xframe=%d)" % (xctr,xframe))
    if (overlayAntennas):
        if (debug):
            print("Breaking out of antenna loop because we are done -------------------")
        break
  # end of while(xant) loop
  pb.draw()
  if (len(plotfiles) == 1 and figfileSequential):
      # rename the single file to remove ".000"
      newplotfiles = [plotfiles[0].split('.000.png')[0]+'.png']
      print("renaming %s to %s" % (plotfiles[0],newplotfiles[0]))
      os.system('mv %s %s' % (plotfiles[0],newplotfiles[0]))
      plotfiles = newplotfiles
  if (len(plotfiles) > 0 and buildpdf):
    pdfname = figfile+'.pdf'
    filelist = ''
    plotfiles = np.unique(plotfiles)
    for i in range(len(plotfiles)):
      cmd = '%s -density %d %s %s.pdf' % (convert,density,plotfiles[i],plotfiles[i].split('.png')[0])
      print("Running command = %s" % (cmd))
      mystatus = os.system(cmd)
      if (mystatus != 0):
          # Try MacOS typical location
          convert = '/opt/local/bin/convert'
          cmd = '%s -density %d %s %s.pdf' % (convert,density,plotfiles[i],plotfiles[i].split('.png')[0])
          print("Running command = %s" % (cmd))
          mystatus = os.system(cmd)
          if (mystatus != 0):
              print("ImageMagick's convert command not found, no PDF built")
              buildpdf = False
              break
      if (cleanup):
          os.system('rm -f %s' % (plotfiles[i]))
      filelist += plotfiles[i].split('.png')[0] + '.pdf '
    if (buildpdf and (len(plotfiles)>1 or not figfileSequential)):
        # The following 2 lines reduce the total number of characters on the command line, which
        # was apparently a problem at JAO for Liza.
        filelist = ' '.join(au.pruneFilelist(filelist.split()))
        pdfname = au.pruneFilelist([pdfname])[0]
        cmd = '%s %s cat output %s' % (pdftk, filelist, pdfname)
        print("Running command = %s" % (cmd))
        mystatus = os.system(cmd)
        if (mystatus != 0):
            cmd = '%s -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s' % (gs,pdfname,filelist)
            print("Running command = %s" % (cmd))
            mystatus = os.system(cmd)
            gs = '/opt/local/bin/gs'
            if (mystatus != 0):
                cmd = '%s -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s' % (gs,pdfname,filelist)
                print("Running command = %s" % (cmd))
                mystatus = os.system(cmd)
        if (mystatus == 0):
            print("PDF left in %s" % (pdfname))
            os.system("rm -f %s" % filelist)
        else:
            print("Both pdftk and ghostscript are missing, so no PDF built.")
    else:
        print("Not building PDF (buildpdf=%s, len(plotfiles)=%d, figfileSequential=%s)" % (str(buildpdf), len(plotfiles), figfileSequential))
  else:
      print("Not building PDF (buildpdf=%s, len(plotfiles)=%d)" % (str(buildpdf), len(plotfiles)))
  if (wasInteractive and interactive==False):
      print("Restoring interactive pylab.")
      pb.ion()
  if (wasInteractive==False and interactive):
      print("Restoring non-interactive pylab.")
      pb.ioff()
  if (debug):
      print("lastUnflaggedAntennaToPlot = ", lastUnflaggedAntennaToPlot)

  showFinalMessage(overlayAntennas, solutionTimeSpread, nUniqueTimes)

  if (channeldiff>0):
      # Compute median over all antennas, or at least those completed before 'q' was hit
      madstats['median'] = dict.fromkeys(spwsToPlot)
      spwvalue = {}
      spwvalue['amp'] = []
      spwvalue['phase'] = []
      for j in spwsToPlot:
          madstats['median'][j] = dict.fromkeys(timerangeList) # dict.fromkeys(range(len(uniqueTimes)))
          for k in timerangeList: # range(len(uniqueTimes)):
              madstats['median'][j][k] = dict.fromkeys(range(nPolarizations))
              for l in range(nPolarizations):
                  if (yaxis == 'both'):
                      madstats['median'][j][k][l] = {'amp': None, 'phase': None}
                  elif (yaxis == 'phase'):
                      madstats['median'][j][k][l] = {'phase': None}
                  else:
                      # this includes tsys and amp
                      madstats['median'][j][k][l] = {'amp': None}
                  for m in madstats['median'][j][k][l].keys():
                      value = []
                      for i in madstats.keys():  # loop over antennas
                          if (i != 'median' and i != 'platforming'):
                              if (madstats[i][j][k][l][m] is not None):
#                                  print "madstats[%s][%d][%d][%d][%s] = " % (i,j,k,l,m), madstats[i][j][k][l][m]
                                  value.append(madstats[i][j][k][l][m])
                                  spwvalue[m].append(madstats[i][j][k][l][m])
                      madstats['median'][j][k][l][m] = np.median(value)
      # now add another spw which is the median over spw,time,polarization
      if (yaxis == 'both'):
          madstats['median']['median']={'amp': np.median(spwvalue['amp']),
                                        'phase': np.median(spwvalue['phase'])}
      elif (yaxis == 'phase'):
          madstats['median'][j][k][l] = {'phase': np.median(spwvalue['phase'])}
      else:
          madstats['median']['median'] = {'amp': np.median(spwvalue['amp'])}
      
      return(madstats)
  else:
      if (msFound and mymsmd != ''):
          mymsmd.close()
      return(vm)
  # end of plotbandpass3


def getTelescopeNameFromCaltable(caltable):
    mytb = au.createCasaTool(tbtool)
    mytb.open(caltable)
    if ('OBSERVATION' in mytb.getkeywords()):
        observationTable = mytb.getkeyword('OBSERVATION').split()[1]
    else:
        observationTable = None
    mytb.close()
    if (observationTable is None):
        return('')
    else:
        return(getTelescopeNameFromCaltableObservationTable(observationTable))
    

def getTelescopeNameFromCaltableObservationTable(observationTable):
    mytb = au.createCasaTool(tbtool)
    mytb.open(observationTable)
    telescope = mytb.getcell('TELESCOPE_NAME')
    mytb.close()
    return(telescope)

def getCorrTypeByAntennaName(firstAntenna):
    """
    This function is used only if the OBSERVATION table of the caltable is blank and the MS is unavailable.
    """
    print("Using antenna name (%s) to set the polarization type." % (firstAntenna))
    if (firstAntenna.find('ea') >= 0):
        corr_type_string = ['RR','LL']
        corr_type = [5,8]
    elif (firstAntenna.find('dv') >= 0 or firstAntenna.find('da') >= 0 or
          firstAntenna.find('pm') >= 0 or firstAntenna.find('da') >= 0):
        corr_type_string = ['XX','YY']
        corr_type = [9,12]
    else: # SMA
        corr_type_string = ['XX']
        corr_type = [9]
    return(corr_type, corr_type_string, len(corr_type))

def showFinalMessage(overlayAntennas, solutionTimeSpread, nUniqueTimes):
  if (overlayAntennas and solutionTimeSpread > 0 and nUniqueTimes==1):
      print("If not all spws were shown, then try setting solutionTimeThreshold=%.0f seconds" % (solutionTimeSpread+1))

def computeOriginalSpwsToPlot(spwsToPlot, originalSpws, tableFormat, debug):
    if (tableFormat > 33):
        # New caltables use the same numbering as the original ms
        return(spwsToPlot)
    else:
        originalSpwsToPlot = []
        for spw in spwsToPlot:
            originalSpwsToPlot.append(originalSpws[spw])
        return(list(originalSpwsToPlot))

def computeScansForUniqueTimes(uniqueTimes, cal_scans, times, unique_cal_scans,
                               debug=False):
    """
    The original implementation assumed that individual scans map to 
    individual times with no overlap.
    """
    scansForUniqueTimes = []
    nUniqueTimes = len(uniqueTimes)
    for uT in uniqueTimes:
        if (debug): print("-------- Checking uniqueTime = %s" % (str(uT)))
        scansForUniqueTimes.append(cal_scans[list(times).index(uT)])
    if (len(unique_cal_scans) == 1):
        if (unique_cal_scans[0] != -1): 
            if (len(scansForUniqueTimes) != len(np.unique(scansForUniqueTimes))):
                if debug:
                    print("Because there are multiple timestamps per scan, I will not assume there is a one-to-one match.")
            else:
                nUniqueTimes = len(np.unique(scansForUniqueTimes))
        else:
            # This 3.4 table does not have the scan numbers populated
            scansForUniqueTimes = []
            print("Because the scan numbers are either not filled in this table, or the solutions span multiple scans, I will use timestamps instead.")
    else:
        if (len(scansForUniqueTimes) != len(np.unique(scansForUniqueTimes))):
            if debug:
                print("Because there are multiple timestamps per scan, I will not assume there is a one-to-one match.")
        else:
            nUniqueTimes = len(np.unique(scansForUniqueTimes))
    return(scansForUniqueTimes, nUniqueTimes)
    

def calcChebyshev(coeff, validDomain, x):
  """
  Given a set of coefficients,
  this method evaluates a Chebyshev approximation.
  """
  if (type(x) == float or type(x) == int):
       x = [x]
  myxrange = validDomain[1] - validDomain[0]
  x = -1 + 2*(x-validDomain[0])/myxrange
  coeff[0] = 0
  if (True):
      try:
          # python 2.7
          v = np.polynomial.chebyshev.chebval(x,coeff)
      except:
          # python 2.6
          v = np.polynomial.chebval(x,coeff)
  else:
    # manual approach, before I found chebval()
    v = np.zeros(len(x))
    if (len(coeff) > 0):
        v += coeff[0] * 1
    if (len(coeff) > 1):
        v += coeff[1] * (x)
    if (len(coeff) > 2):
        v += coeff[2] * (2*x**2 - 1)
    if (len(coeff) > 3):
        v += coeff[3] * (4*x**3 - 3*x)
    if (len(coeff) > 4):
        v += coeff[4] * (8*x**4 - 8*x**2 + 1)
    if (len(coeff) > 5):
        v += coeff[5] * (16*x**5 - 20*x**3 + 5*x)
    if (len(coeff) > 6):
        v += coeff[6] * (32*x**6 - 48*x**4 + 18*x**2 - 1)
    if (len(coeff) > 7):
        v += coeff[7] * (64*x**7 -112*x**5 + 56*x**3 - 7*x)
    if (len(coeff) > 8):
        v += coeff[8] * (128*x**8 -256*x**6 +160*x**5 - 32*x**2 + 1)
    if (len(coeff) > 9):
        v += coeff[9] * (256*x**9 -576*x**7 +432*x**5 - 120*x**3 + 9*x)
    if (len(coeff) > 10):
        print("Chebyshev polynomials with degree > 10 are not implemented")

  return(v)
    
def ResizeFonts(adesc,fontsize):
#    print "Called ResizeFonts()"
    yFormat = ScalarFormatter(useOffset=False)
    adesc.yaxis.set_major_formatter(yFormat)
    adesc.xaxis.set_major_formatter(yFormat)
    pb.setp(adesc.get_xticklabels(), fontsize=fontsize)
    pb.setp(adesc.get_yticklabels(), fontsize=fontsize)

def complexMeanRad(phases):
    # convert back to real and imaginary, take mean, then convert back to phase
    meanSin = np.mean(np.sin(phases))
    meanCos = np.mean(np.cos(phases))
    return(180*np.arctan2(meanSin, meanCos)/math.pi)

def complexMeanDeg(phases):
    # convert back to real and imaginary, take mean, then convert back to phase
    phases *= math.pi/180
    meanSin = np.mean(np.sin(phases))
    meanCos = np.mean(np.cos(phases))
    return(180*np.arctan2(meanSin, meanCos)/math.pi)

def CalcAtmTransmission(chans,freqs,xaxis,pwv,vm,vis,asdm,antenna,timestamp,
                        interval,field,refFreqInTable, net_sideband=0,
                        mytime=0, missingCalWVRErrorPrinted=False, mymsmd='',
                        caltable='', verbose=False, 
                        maxAtmCalcChannels=MAX_ATM_CALC_CHANNELS,
                        maxAltitude=60.0, atmType=1, h0=1.0, dP=5.0, dPm=1.1,
                        Feff=0.99, SBGain=0.99, Tamb=273.0, Trx=None, showtsys=False):
    """
    Uses the CASA 'at' tool to compute atmospheric model.
    chans: list of all channels, regardless of whether they are flagged
    freqs: list of requested frequencies (in GHz) corresponding to chans
    xaxis: what we are plotting on the xaxis: 'chan' or 'freq'
       if 'chan', then we reverse the direction of the data if necessary (LSB)
    pwv: 'auto', or a value in mm
    vm: '' in the modern era of msmd, or an old valueMapping structure
    vis: measurement set
    asdm: ASDM, passed to getMedianPWV, only if vis directory does not contain any 
          of the following: ASDM_CALWVR, ASDM_CALATMOSPHERE, or CalWVR.xml
    antenna: name string, ID string, or integer ID (passed only to getWeather)
    timestamp: in MJDsec; used with interval to determine timerange to pass to getMedianPWV
    interval: in sec; used with timestamp to determine timerange to pass to getMedianPWV
    field: integer ID; used to find scan observation time to pass to getWeather
    refFreqInTable: in Hz; used along with net_sideband to determine usb/lsb
    mytime: only used to show how the function was called (which timestamp)
    missingCalWVRErrorPrinted: Boolean to prevent multiple printings of same warning
    mymsmd: existing instance of msmd tool
    caltable: only used to get telescope name if mymsmd=='' and vm==''
    maxAtmCalcChannels: set this lower to prevent longer processing times
    atmType:   tropical=1, midLatitudeSummer=2, midLatitudeWinter=3
    Trx: receiver temperature to use when computing Tsys, None -> use au.receiverTrxSpec
    """
    if (casaVersion >= '4.1.0' and mymsmd != ''):
        telescopeName = mymsmd.observatorynames()[0]
    elif (vm != ''):
        telescopeName = au.getObservatoryName(vis)
    elif (os.path.exists(caltable)):
        telescopeName = getTelescopeNameFromCaltable(caltable)
    else:
        telescopeName = 'ALMA'
    if (telescopeName.find('ALMA') >= 0):
        defaultPWV = 1.0   # a reasonable value for ALMA in case it cannot be found
    elif (telescopeName.find('VLA') >= 0):
        defaultPWV = 5.0  
    else:
        defaultPWV = 5.0  
    if (type(pwv) == str):
      if (pwv.find('auto')>=0):
        if (os.path.exists(vis+'/ASDM_CALWVR') or os.path.exists(vis+'/ASDM_CALATMOSPHERE') or
            os.path.exists('CalWVR.xml')):
              if (verbose):
                  print("*** Computing atmospheric transmission using measured PWV, field %d, time %d (%f). ***" % (field,mytime,timestamp))
              timerange = [timestamp-interval/2, timestamp+interval/2]
              if (os.path.exists(vis+'/ASDM_CALWVR') or os.path.exists(vis+'/ASDM_CALATMOSPHERE')):
                  if (verbose):
                      print("Calling au.getMedianPWV('%s',%s,'%s',verbose=False) " % (vis,timerange,asdm))
                  [pwvmean, pwvstd]  = au.getMedianPWV(vis,timerange,asdm,verbose=False)
              else:
                  if (verbose):
                      print("Calling au.getMedianPWV('%s',%s,asdm='',verbose=False) " % (vis,timerange))
                  [pwvmean, pwvstd]  = au.getMedianPWV('.',timerange,asdm='',verbose=False)
              if (verbose):
                  print("retrieved pwvmean = %f" % pwvmean)
              retrievedPWV = pwvmean
              if (pwvmean < 0.00001):
                  pwvmean = defaultPWV
        else:
              pwvmean = defaultPWV
              if (missingCalWVRErrorPrinted == False):
                  missingCalWVRErrorPrinted = True
                  if (telescopeName.find('ALMA')>=0):
                      print("No ASDM_CALWVR, ASDM_CALATMOSPHERE, or CalWVR.xml table found.  Using PWV %.1fmm." % pwvmean)
                  else:
                      print("This telescope has no WVR to provide a PWV measurement. Using PWV %.1fmm." % pwvmean)
      else:
          try:
              pwvmean = float(pwv)
          except:
              pwvmean = defaultPWV
    else:
          try:
              pwvmean = float(pwv)
              print("Using supplied PWV = %f mm" % pwvmean)
          except:
              pwvmean = defaultPWV

    if (verbose):
        print("Using PWV = %.3f mm" % pwvmean)

    # default values in case we can't find them below
    airmass = 1.5
    P = 563.0
    H = 20.0
    T = 273.0


    roundedScanTimes = []
    if (casaVersion >= '4.1.0' and mymsmd != ''):
        if (type(field) == list or type(field) == type(np.ndarray(0))):
            field = field[0]
        if (verbose):
            print("Looking for scans for field integer = %d, type(field)=%s" % (field,str(type(field))))
        if (casaVersion >= '4.4' and casaVersion < '4.6'):
            myscans = []
            myobsids = []
            for obsid in range(mymsmd.nobservations()):
                newScans = list(mymsmd.scansforfield(field, obsid=obsid))
                for myNewScan in newScans:
                    roundedScanTimes.append(np.unique(np.round(mymsmd.timesforscan(myNewScan, obsid=obsid))))
                myscans += newScans
                myobsids += len(newScans)*[obsid]
            myscans = np.array(myscans, dtype=np.int32)
        else:
            myscans = mymsmd.scansforfield(field)
            for myscan in myscans:
                roundedScanTimes.append(np.unique(np.round(mymsmd.timesforscan(myscan))))
#   This method was much slower and not necessary.  Removed for CAS-8065
#        scantimes = mymsmd.timesforscans(myscans) # is often longer than the scans array
#        roundedScanTimes = np.unique(np.round(scantimes,0))
#        if (verbose):
#            print "Running getScansForTimes (%d -> %d)" % (len(scantimes), len(roundedScanTimes))
#        myscans,roundedScanTimes = getScansForTimes(mymsmd,roundedScanTimes) # be sure that each scantime has a scan associated, round to nearest second to save time (esp. for single dish data)
#        if (verbose):
#            print "done"
    else:
        if (verbose):
            print("Looking for scans for field integer = %d, type(field)=%s" % (field,str(type(field))))
        myscans = []
        if (vm != ''):
            myscans = vm.getScansForFieldID(field)
            scantimes = vm.getTimesForScans(myscans)
            roundedScanTimes = scantimes
    if (verbose):
          print("For field %s, Got scans = " % str(field), np.unique(myscans))
    mindiff = 1e20
    bestscan = 1
    for i in range(len(roundedScanTimes)):
        stime = roundedScanTimes[i]
        meantime = np.mean(stime)
        tdiff = np.abs(meantime-timestamp)
        if (tdiff < mindiff):
            bestscan = myscans[i]
            mindiff = tdiff
            if (casaVersion >= '4.4' and casaVersion < '4.6'):
                bestscan_obsid = myobsids[i]
    if (verbose):
          print("For timestamp=%.1f, got closest scan = %d, %.0f sec away" %(timestamp, bestscan,mindiff))
    if (vm != '' or mymsmd != ''):
        if (casaVersion >= '4.4' and casaVersion < '4.6'):
            weatherResult = au.getWeather(vis,bestscan,antenna,verbose,vm,mymsmd,obsid=bestscan_obsid,getSolarDirection=False)
        else:
            weatherResult = au.getWeather(vis,bestscan,antenna,verbose,vm,mymsmd,getSolarDirection=False)
    else:
        weatherResult = None
    if (weatherResult is None):
        conditions = {}
        P = 786  # VLA value, if no ms weather is found
        H = 20
    else:
        [conditions,myTimes,vm] = weatherResult
        P = conditions['pressure']
        H = conditions['humidity']
        T = conditions['temperature']+273.15
    if (P <= 0.0):
        P = 563
    if (H <= 0.0):
        H = 20
    if ((vm!='' or mymsmd!='') and ('elevation' in conditions.keys()) == False):
        # Someone cleared the POINTING table, so calculate elevation from Ra/Dec/MJD
        if (mymsmd != '' and casaVersion >= '4.1.0'):
            if (casaVersion >= '4.4' and casaVersion < '4.6'):
                if (verbose):
                    print("mymsmd.fieldsforscan(bestscan) = ", mymsmd.fieldsforscan(bestscan,obsid=bestscan_obsid))
                myfieldId =  mymsmd.fieldsforscan(bestscan,obsid=bestscan_obsid)[0]
                myscantime = np.mean(mymsmd.timesforscan(bestscan,obsid=bestscan_obsid))
            else:
                if (verbose):
                    print("mymsmd.fieldsforscan(bestscan) = ", mymsmd.fieldsforscan(bestscan))
                myfieldId =  mymsmd.fieldsforscan(bestscan)[0]
                myscantime = np.mean(mymsmd.timesforscan(bestscan))
            telescopeName = mymsmd.observatorynames()[0]
            if (len(telescopeName) < 1):
                telescopeName = 'ALMA'
        elif (vm != ''):
            myfieldId =  vm.getFieldIdsForFieldName(vm.getFieldsForScan(bestscan))
            myscantime = np.mean(vm.getTimesForScans(bestscan))
            telescopeName = au.getObservatoryName(vis)
        else:
            myfieldId = 0
            telescopeName = 'VLA'
        mydirection = au.getRADecForField(vis, myfieldId)
        if (verbose):
            print("Scan =  %d, time = %.1f,  Field = %d, direction = " % (bestscan, myscantime, myfieldId), mydirection)
        if (len(telescopeName) < 1):
            telescopeName = 'ALMA'
        myazel = au.computeAzElFromRADecMJD(mydirection, myscantime/86400., telescopeName, verbose=False)
        conditions['elevation'] = myazel[1] * 180/math.pi
        conditions['azimuth'] = myazel[0] * 180/math.pi
        if (verbose):
            print("Computed elevation = %.1f deg" % (conditions['elevation']))
        
    if (verbose):
        if ('elevation' in conditions.keys()):
            print("CalcAtm: found elevation=%f (airmass=%.3f) for scan:" % (conditions['elevation'],1/np.sin(conditions['elevation']*np.pi/180.)), bestscan)
        print("P,H,T = %f,%f,%f" % (P,H,T))
    if ('elevation' not in conditions.keys()):
        print("Using 45 deg elevation since the actual value is unavailable.")
        airmass = 1.0/math.cos(45*math.pi/180.)
    elif (conditions['elevation'] <= 3):
        print("Using 45 deg elevation instead of %f" % (conditions['elevation']))
        airmass = 1.0/math.cos(45*math.pi/180.)
    else:
        airmass = 1.0/math.cos((90-conditions['elevation'])*math.pi/180.)

    numchan = len(freqs)
    # Set the reference freq to be the middle of the middle two channels
    reffreq=0.5*(freqs[int(numchan/2)-1]+freqs[int(numchan/2)])
    originalnumchan = numchan
    while (numchan > maxAtmCalcChannels):
        numchan //= 2
#        print "Reduce numchan to ", numchan
        chans = range(0,originalnumchan,(originalnumchan//numchan))

    chansep = (freqs[-1]-freqs[0])/(numchan-1)
    nbands = 1
    myqa = au.createCasaTool(qatool)
    fCenter = au.create_casa_quantity(myqa,reffreq,'GHz')
    fResolution = au.create_casa_quantity(myqa,chansep,'GHz')
    fWidth = au.create_casa_quantity(myqa,numchan*chansep,'GHz')
    myat = au.createCasaTool(attool)
    result = myat.initAtmProfile(humidity=H, temperature=au.create_casa_quantity(myqa,T,"K"),
                                 altitude=au.create_casa_quantity(myqa,5059,"m"),
                                 pressure=au.create_casa_quantity(myqa,P,'mbar'),
                                 atmType=atmType,
                                 h0=au.create_casa_quantity(myqa, h0,"km"),
                                 maxAltitude=au.create_casa_quantity(myqa, maxAltitude,"km"),
                                 dP=au.create_casa_quantity(myqa, dP,"mbar"),
                                 dPm=dPm
                                 )
    if type(result) == tuple:
        # CASA 6 (in CASA 5, it is simply one string)
        result = result[0]
    if verbose: 
          au.printNumberOfAtmosphericLayers(result)
    myat.initSpectralWindow(nbands,fCenter,fWidth,fResolution)
    myat.setUserWH2O(au.create_casa_quantity(myqa,pwvmean,'mm'))
#    myat.setAirMass()  # This does not affect the opacity, but it does effect TebbSky, so do it manually.

    n = myat.getNumChan()
    if (casaVersion < '4.0.0'):
        dry = np.array(myat.getDryOpacitySpec(0)['dryOpacity'])
        wet = np.array(myat.getWetOpacitySpec(0)['wetOpacity'].value)
        TebbSky = []
        for chan in range(n):  # do NOT use numchan here, use n
            TebbSky.append(myat.getTebbSky(nc=chan, spwid=0).value)
        TebbSky = np.array(TebbSky)
        # readback the values to be sure they got set
#        rf = myat.getRefFreq().value
#        cs = myat.getChanSep().value # MHz
    else:
        dry = np.array(myat.getDryOpacitySpec(0)[1])
        wet = np.array(myat.getWetOpacitySpec(0)[1]['value'])
        TebbSky = myat.getTebbSkySpec(spwid=0)[1]['value']
        # readback the values to be sure they got set
#        rf = myqa.convert(myat.getRefFreq(),'GHz')['value']
#        cs = myqa.convert(myat.getChanSep(),'MHz')['value'] # MHz
    # for an even-numbered-channel spw (0..127), the center is 63.5, not 128/2=64
    # where 0 = middle of first channel and 127 = middle of final channel

    transmission = np.exp(-airmass*(wet+dry))
    TebbSky *= (1-np.exp(-airmass*(wet+dry)))/(1-np.exp(-wet-dry))

    if (refFreqInTable*1e-9>np.mean(freqs)):
        if ((net_sideband % 2) == 0):
            sense = 1
        else:
            sense = 2
    else:
        if ((net_sideband % 2) == 0):
            sense = 2
        else:
            sense = 1
    if (sense == 1):
        if (xaxis.find('chan')>=0):
            trans = np.zeros(len(transmission))
            Tebb = np.zeros(len(TebbSky))
            for i in range(len(transmission)):
                trans[i] = transmission[len(transmission)-1-i]
                Tebb[i] = TebbSky[len(TebbSky)-1-i]
            transmission = trans
            TebbSky = Tebb
    if showtsys:
        if Trx is None or Trx == 'auto':
            Trx = au.receiverTrxSpec(au.getBand(freqs[0]*1e9))
        Tsys = (Feff*TebbSky + (1.0-Feff)*Tamb + Trx) * ((1.0 + (1.0-SBGain)) / (Feff*np.exp(-airmass*(wet+dry))))
    else:
        Tsys = None

    # Be sure that number of frequencies matched number of transmission values - CAS-10123
    numchan = len(transmission)
    chans = range(len(transmission))
    startFreq = myqa.convert(myat.getChanFreq(0),'GHz')['value']
    endFreq = myqa.convert(myat.getChanFreq(numchan-1),'GHz')['value']
    myat.close()
    myqa.done()
    freq = np.linspace(startFreq, endFreq, numchan)
    # print "startFreq=%f  endFreq=%f " % (startFreq, endFreq)
#    chansepGHz = cs*0.001
#    if sense == 2:
#        freq = np.linspace(rf-((numchan-1)/2.)*chansepGHz, rf+((numchan-1)/2.)*chansepGHz, numchan) 
#    else:
#        freq = np.linspace(rf+((numchan-1)/2.)*chansepGHz, 
#                           rf-((numchan-1)/2.)*chansepGHz, numchan)
    return(freq, chans, transmission, pwvmean, airmass, TebbSky, 
           missingCalWVRErrorPrinted, Tsys)

def checkForNaNs(mylist):
    nans = 0
    for x in mylist:
        if (x != x):
            nans += 1
    return(nans)

def RescaleTrans(trans, lim, subplotRows, lo1=None, xframe=0,
                 overlaySpws=False, overlayBasebands=False):
    # Input: the array of transmission or TebbSky values and current limits
    # Returns: arrays of the rescaled transmission values and the zero point
    #          values in units of the frame, and in amplitude.
    debug = False
    yrange = lim[1]-lim[0]
    if (lo1 is None):
        labelgap = 0.6 # Use this fraction of the margin for the PWV ATM label
    else:
        labelgap = 0.5 # Use this fraction of the margin to separate the top
                       # curve from the upper y-axis
    y2 = lim[1] - labelgap*yrange*TOP_MARGIN/(1.0+TOP_MARGIN)
    y1 = lim[1] - yrange*TOP_MARGIN/(1.0+TOP_MARGIN)
    transmissionRange = np.max(trans)-np.min(trans)
    if (overlaySpws or overlayBasebands) and False:
        if (transmissionRange < 1):
            # Then a transmission range from 0..1 was passed in
            transmissionRange = 1.0
        else:
            # Then a TebbSky was passed in (in Kelvin)
            transmissionRange = 1.0
            transmissiongRange = 290 
    else:
        if (transmissionRange < 0.05):
            # force there to be a minimum range of transmission display 
            # overemphasize tiny ozone lines
            transmissionRange = 0.05

        if (transmissionRange > 1 and transmissionRange < 10):
            # force there to be a minimum range of Tebbsky (10K) to display
            transmissionRange = 10

    # convert transmission to amplitude
    newtrans = y2 - (y2-y1)*(np.max(trans)-trans)/transmissionRange

    # Use edge values
    edgeValueTransmission = trans[-1]
    otherEdgeValueTransmission = trans[0]

    # Now convert the edge channels' transmission values into amplitude
    edgeValueAmplitude = y2 - (y2-y1)*(np.max(trans)-trans[-1])/transmissionRange
    otherEdgeValueAmplitude = y2 - (y2-y1)*(np.max(trans)-trans[0])/transmissionRange

    # Now convert amplitude to frame units, offsetting downward by half
    # the font size
    fontoffset = 0.01*subplotRows
    edgeValueFrame = (edgeValueAmplitude - lim[0])/yrange  - fontoffset
    otherEdgeValueFrame = (otherEdgeValueAmplitude - lim[0])/yrange  - fontoffset

    # scaleFactor is how large the plot is from the bottom x-axis
    # up to the labelgap, in units of the transmissionRange
    scaleFactor = (1+TOP_MARGIN*(1-labelgap)) / (TOP_MARGIN*(1-labelgap))

    # compute the transmission at the bottom of the plot, and label it
    y0transmission = np.max(trans) - transmissionRange*scaleFactor
    y0transmissionFrame = 0
    y0transmissionAmplitude = lim[0]

    if (y0transmission <= 0):
        # If the bottom of the plot is below zero transmission, then label
        # the location of zero transmission instead.
        if (debug):
            print("--------- y0transmission original = %f, (y1,y2)=(%f,%f)" % (y0transmission,y1,y2))
        y0transmissionAmplitude = y1-(y2-y1)*(np.min(trans)/transmissionRange)
        y0transmissionFrame = (y0transmissionAmplitude-lim[0]) / (lim[1]-lim[0])
        y0transmission = 0
    if (debug):
        print("-------- xframe=%d, scaleFactor = " % (xframe), scaleFactor)
        print("edgeValueFrame, other = ", edgeValueFrame, otherEdgeValueFrame)
        print("edgeValueTransmission, other = ", edgeValueTransmission, otherEdgeValueTransmission)
        print("edgeValueAmplitude, otherEdgeValueAmplitude = ", edgeValueAmplitude, otherEdgeValueAmplitude)
        print("y0transmission = %f, y0transmissionFrame = %f" % (y0transmission,y0transmissionFrame))
        print("y0transmissionAmplitude = ", y0transmissionAmplitude)
        print("transmissionRange = ", transmissionRange)
    return(newtrans, edgeValueFrame, y0transmission, y0transmissionFrame,
           otherEdgeValueFrame, edgeValueTransmission,
           otherEdgeValueTransmission, edgeValueAmplitude,
           otherEdgeValueAmplitude, y0transmissionAmplitude)

def RescaleX(chans, lim, plotrange, channels):
    # This function is now only used by DrawAtmosphere when xaxis='chan'.
    # It is only really necessary when len(chans)>MAX_ATM_CALC_CHANNELS.
    #  - September 2012
    # If the user specified a plotrange, then rescale to this range,
    # otherwise rescale to the automatically-determined range.

    # chans = 0..N where N=number of channels in the ATM_CALC
    # channels = 0..X where X=number of channels in the spw, regardless of flagging

    if (len(chans) != len(channels)):
        if (chans[1] > chans[0]):
            atmchanrange = chans[-1]-chans[0]
        else:
            atmchanrange = chans[0]-chans[-1]

        if (channels[1] > channels[0]):
            chanrange = channels[-1]-channels[0]
        else:
            chanrange = channels[0]-channels[-1]
        
        newchans = np.array(chans)*chanrange/atmchanrange
        return(newchans)
    else:
        return(chans)  

def recalcYlimitsFreq(chanrange, ylimits, amp, sideband,plotrange,xchannels,
                      debug=False,location=0, chanrangePercent=None):
  # Used by plots with xaxis='freq'
  # xchannels are the actual channel numbers of unflagged data, i.e. displayed points
  # amp is actual data plotted
  # This function is often overridden by SetLimits which is typically called afterward.
  # It is not clear that it ever takes effect, but it might in some cases. - 2015-03-02
  ylim_debug = False
  amp = np.array(amp)
  if (debug):
      print("recalcYlimitsFreq(loc=%d): len(xchannels) = %d" % (location, len(xchannels)))
  if (len(amp) < 1):
      return(pb.ylim()) # ylimits)
  if (chanrange[0]==0 and chanrange[1] == 0 and plotrange[2] == 0 and plotrange[3]==0
      and chanrangePercent is None):
      if (len(amp) == 1):
          if (ylim_debug):
              print("amp = ", amp)
          ylimits = [amp[0]-0.2, amp[0]+0.2]
      else:
          newmin = np.min(amp)
          newmax = np.max(amp)
          newmin = np.min([ylimits[0],newmin])
          newmax = np.max([ylimits[1],newmax])
          ylimits = [newmin, newmax]
  elif ((abs(chanrange[0]) > 0 or abs(chanrange[1]) > 0)):
      plottedChannels = np.intersect1d(xchannels, range(chanrange[0],chanrange[1]+1))
      if (len(plottedChannels) < 1):
          return(ylimits)
      mylist = np.arange(xchannels.index(plottedChannels[0]), 1+xchannels.index(plottedChannels[-1]))
      if (mylist[-1] >= len(amp)):
          # prevent crash if many channels are flagged
          return(ylimits)
      if (ylim_debug):
          print("Starting with limits = ", ylimits)
          print("Examining channels: ", mylist)
          print("len(amp): %d" % (len(amp)))
          print("type(amp) = %s" % (str(type(amp))))
          print("Examining values: amp[mylist] = %s" % (str(amp[mylist])))
      newmin = np.min(amp[mylist])
      newmax = np.max(amp[mylist])
      newmin = np.min([ylimits[0],newmin])
      newmax = np.max([ylimits[1],newmax])
      #  The following presents a problem with overlays, as it keeps widening forever
#      newmin -= 0.05*(newmax-newmin)
#      newmax += 0.05*(newmax-newmin)
      ylimits = [newmin, newmax]
  elif (chanrangePercent is not None):
      startFraction = (100-chanrangePercent)*0.5*0.01
      stopFraction = 1-(100-chanrangePercent)*0.5*0.01
      if (xchannels == []):
          # prevent crash if many channels are flagged: 2015-04-13
          return(ylimits)
      cr0 = int(np.round(np.max(xchannels)*startFraction))
      cr1 = int(np.round(np.max(xchannels)*stopFraction))
#      print "Plotting %.0f%% of channels %d through %d" % (chanrangePercent, cr0, cr1)
      plottedChannels = np.intersect1d(xchannels, range(cr0, cr1+1))
      if (len(plottedChannels) < 1):
          return(ylimits)
      mylist = np.arange(xchannels.index(plottedChannels[0]), 1+xchannels.index(plottedChannels[-1]))
      if (mylist[-1] >= len(amp)):
          # prevent crash if many channels are flagged
          return(ylimits)
      if (ylim_debug):
          print("Starting with limits = ", ylimits)
          print("Examining channels: ", mylist)
          print("len(amp): %d" % (len(amp)))
          print("type(amp) = %s" % (str(type(amp))))
          print("Examining values: amp[mylist] = %s" % (str(amp[mylist])))
      newmin = np.min(amp[mylist])
      newmax = np.max(amp[mylist])
      newmin = np.min([ylimits[0],newmin])
      newmax = np.max([ylimits[1],newmax])
      ylimits = [newmin, newmax]
  if (ylim_debug):
      print("Returning from loc=%d with limits = %s" % (location, str(ylimits)))
  return ylimits

def recalcYlimits(plotrange, ylimits, amp):
  # Used by plots with xaxis='chan'
  if (len(amp) < 1):
      return(pb.ylim())
  if ((abs(plotrange[0]) > 0 or abs(plotrange[1]) > 0) and (plotrange[2] == 0 and plotrange[3] == 0)):
      x0 = plotrange[0]
      x1 = plotrange[1]
      if (x0 < 0):
          x0 = 0
      if (x1 > len(amp)-1):
          x1 = len(amp)-1
      if (len(amp) > x1 and x0 < x1):
          newmin = np.min(amp[x0:x1])
          newmax = np.max(amp[x0:x1])
          newmin = np.min([ylimits[0],newmin])
          newmax = np.max([ylimits[1],newmax])
          ylimits = [newmin, newmax]
  else:
      ylimits = pb.ylim()  # added on 10/27/2011
#      print "current ylimits = ", ylimits
  return(ylimits)

def SetNewYLimits(newylimits, loc=-1):
    if (False): 
        print("loc=%d, Entered SetNewYLimits with " % (loc), newylimits) 
    newrange = newylimits[1]-newylimits[0]
    if (newrange > 0):
        pb.ylim([newylimits[0], newylimits[1]])

def SetNewXLimits(newxlimits, loc=-1):
#    print "Entered SetNewXLimits from location=%d with range = %.3f" % (loc,np.max(newxlimits)-np.min(newxlimits))
    myxrange = np.abs(newxlimits[1]-newxlimits[0])
    if (myxrange == 0):
        myxrange = 0.001
    mybuffer = 0.01
    if (newxlimits[0] < newxlimits[1]):
        pb.xlim([newxlimits[0]-myxrange*mybuffer,newxlimits[1]+myxrange*mybuffer] )
    else:
#        print "Swapping xlimits order"
        pb.xlim(newxlimits[1]-myxrange*mybuffer, newxlimits[0]+myxrange*mybuffer)

def sloppyMatch(newvalue, mylist, threshold, mytime=None, scansToPlot=[],
                scansForUniqueTimes=[], myprint=False, whichone=False):
    """
    If scan numbers are present, perform an exact match, otherwise compare the
    time stamps of the solutions.
    """
    debug = myprint
    if (debug):
        print("sloppyMatch: scansToPlot = %s" % (str(scansToPlot)))
    mymatch = None
    if (len(scansToPlot) > 0):
        if (mytime >= len(scansForUniqueTimes)):
            print("sloppyMatch() mytime is too large:  mytime=%d >= len(scansForUniqueTimes)=%d: " % (mytime, len(scansForUniqueTimes)), scansForUniqueTimes)
        matched = scansForUniqueTimes[mytime] in scansToPlot
        if (whichone or myprint):
            myscan = scansForUniqueTimes[mytime]
            if (myscan in scansToPlot):
                mymatch = list(scansToPlot).index(myscan)
        if (matched == False and myprint==True):
            print("sloppyMatch: %d is not in %s" % (myscan, list(scansToPlot)))
        elif (myprint==True):
            print("sloppyMatch: %d is in %s" % (myscan, list(scansToPlot)))
    else:
        matched = False
        if (type(mylist) != list and type(mylist)!=np.ndarray):
            mylist = [mylist]
        mymatch = -1
        for i in range(len(mylist)):
            v = mylist[i]
            if (abs(newvalue-v) < threshold):
                matched = True
                mymatch = i
        if (matched == False and myprint==True):
            print("sloppyMatch: %.0f is not within %.0f of anything in %s" % (newvalue,threshold, str([int(round(b)) for b in mylist])))
        elif (myprint==True):
            print("sloppyMatch: %.0f is within %.0f of something in %s" % (newvalue,threshold, str([int(round(b)) for b in mylist])))
    if (whichone ==  False):
        return(matched)
    else:
        return(matched,mymatch)

def sloppyUnique(t, thresholdSeconds):
    """
    Takes a list of numbers and returns a list of unique values, subject to a threshold difference.
    """
    # start with the first entry, and only add a new entry if it is more than the threshold from prior
    sloppyList = [t[0]]
    for i in range(1,len(t)):
        keepit = True
        for j in range(0,i):
            if (abs(t[i]-t[j]) < thresholdSeconds):
                keepit = False
        if (keepit):
            sloppyList.append(t[i])
#    print "sloppyUnique returns %d values from the original %d" % (len(sloppyList), len(t))
    return(sloppyList)

def SetLimits(plotrange, chanrange, newylimits, channels, frequencies, pfrequencies,
              ampMin, ampMax, xaxis, pxl, chanrangeSetXrange, chanrangePercent=None, 
              debug=False, loc=-1):
    """
    This is the place where chanrange actually takes effect.  
    """
    if debug:
        print("SetLimits(%d): freq[0, -1]=" % (loc), frequencies[0], frequencies[-1])
        print("SetLimits(%d): newylimits=" % (loc), newylimits)
    if (abs(plotrange[0]) > 0 or abs(plotrange[1]) > 0):
        SetNewXLimits([plotrange[0],plotrange[1]],loc=30)
        if (plotrange[2] == 0 and plotrange[3] == 0):
            # reset the ylimits based on the channel range shown (selected via plotrange)
            SetNewYLimits(newylimits,loc=15)
    else: # set xlimits to full range
        if (xaxis.find('chan')>=0):
            SetNewXLimits([channels[0],channels[-1]],loc=31)
        else:
            if (chanrangeSetXrange or (chanrange[0]==0 and chanrange[1]==0 and chanrangePercent is None)): # CAS-7965
                SetNewXLimits([frequencies[0], frequencies[-1]],loc=32)
    if (chanrange[0] != 0 or chanrange[1] != 0 or chanrangePercent is not None):
        # reset the ylimits based on the channel range specified (selected via chanrange)
        if (newylimits != [LARGE_POSITIVE, LARGE_NEGATIVE]):
            SetNewYLimits(newylimits,loc=16)
#        print "pxl=%d, chanrange[0]=%d, chanrange[1]=%d, shape(pfreq), shape(freq)=" % (pxl, chanrange[0], chanrange[1]), np.shape(pfrequencies),np.shape(frequencies)
        # Use frequencies instead of pfrequencies, because frequencies are not flagged and
        # will continue to work if chanranze is specified and data are flagged.
        if chanrangeSetXrange:
            if (chanrangePercent is None):
                try:
                    SetNewXLimits([frequencies[chanrange[0]], frequencies[chanrange[1]]],loc=33)  # Apr 3, 2012
#                    print "chanrangePercent=None, Set new xlimits from channels %d thru %d" % (chanrange[0],chanrange[1]-1)
                except:
                    print("a)Invalid chanrange (%d-%d). Valid range = 0-%d" % (chanrange[0],chanrange[1],len(frequencies)-1))
                    return(-1)
            else:
                startFraction = (100-chanrangePercent)*0.5*0.01
                stopFraction = 1-(100-chanrangePercent)*0.5*0.01
                cr0 = int(np.round(np.max(channels)*startFraction))
                cr1 = int(np.round(np.max(channels)*stopFraction))
                try:
                    SetNewXLimits([frequencies[cr0], frequencies[cr1]],loc=34)
#                    print "Set new xlimits from channels %d thru %d" % (cr0,cr1)
                except:
                    print("b)Invalid chanrange (%d-%d). Valid range = 0-%d" % (cr0,cr1,len(frequencies)-1))
                    return(-1)
    if (abs(plotrange[2]) > 0 or abs(plotrange[3]) > 0):
        SetNewYLimits([plotrange[2],plotrange[3]],loc=17)
    return(0)
    
def showFDM(originalSpw, chanFreqGHz, baseband, showBasebandNumber, basebandDict):
    """
    Draws a horizontal bar indicating the location of FDM spws in the dataset.

    Still need to limit based on the baseband -- need dictionary passed in.
    originalSpw: should contain all spws in the dataset, not just the ones
                in the caltable
    baseband: the baseband of the current spw
    showBasebandNumber: force the display of all FDM spws, and their baseband number
    basebandDict: {1:[17,19], 2:[21,23], etc.}  or {} for really old datasets  
    """
    
    # add some space at the bottom -- Apr 25, 2012
    ylim = pb.ylim()
    yrange = ylim[1]-ylim[0]
    pb.ylim([ylim[0]-BOTTOM_MARGIN*yrange, ylim[1]])

    sdebug = False
    if (sdebug):
        print("Showing FDM (%d)" % (len(originalSpw)), originalSpw)
        print("baseband = %d, basebandDict = %s" % (baseband, str(basebandDict)))
    fdmctr = -1
    x0,x1 = pb.xlim()
    y0,y1 = pb.ylim()
    yrange = y1 - y0
    myxrange = x1 - x0
#    pb.hold(True) # not available in CASA6, but never needed
    labelAbove = False  # False means label to the right
    for i in range(len(originalSpw)):
        nchan = len(chanFreqGHz[i])
        # latter 3 values are for ACA with FPS enabled
        if (nchan >= 15 and nchan not in [256,128,64,32,16,248,124,62]):
          # To fix in Sept 2014
          # crashes if showfdm=True but no ms was found
          if (originalSpw[i] in basebandDict[baseband] or showBasebandNumber):
            fdmctr += 1
            verticalOffset = fdmctr*0.04*yrange
            y1a = y0 + 0.03*yrange + verticalOffset
            if (labelAbove):
                y2 = y1a + 0.01*yrange
            else:
                y2 = y1a - 0.016*yrange
#            print "chan=%d: Drawing line at y=%f (y0=%f) from x=%f to %f" % (len(chanFreqGHz[i]),
#                                              y1a,y0,chanFreqGHz[i][0], chanFreqGHz[i][-1])
            f0 = chanFreqGHz[i][0]
            f1 = chanFreqGHz[i][-1]
            if (f1 < f0):
                swap = f1
                f1 = f0
                f0 = swap
            v0 = np.max([f0,x0])
            v1 = np.min([f1,x1])
            if (v1 > v0):
                if (labelAbove):
                    xlabel = 0.5*(v0+v1)
                    if (xlabel < x0):
                        xlabel = x0
                    if (xlabel > x1):
                        xlabel = x1
                else:
                    xlabel = v1+0.02*myxrange
                pb.plot([v0,v1], [y1a,y1a], '-',
                        linewidth=4, color=overlayColors[fdmctr],markeredgewidth=markeredgewidth)
                if (showBasebandNumber):
                    mybaseband = [key for key in basebandDict if i in basebandDict[key]]
                    if (len(mybaseband) > 0):
                        pb.text(xlabel, y2, "spw%d(bb%d)"%(i,mybaseband[0]), size=7)
                    else:
                        pb.text(xlabel, y2, "spw%d(bb?)"%(i), size=7)
                else:
                    pb.text(xlabel, y2, "spw%d"%(i), size=7)
                if (sdebug): print("Plotting spw %d (%d), xlabel=%f, y2=%f, y1a=%f, y0=%f, y1=%f" % (i, originalSpw[i], xlabel, y2, y1a, y0, y1))
            else:
                if (sdebug): print("Not plotting spw %d (%d) because %f < %f" % (i,originalSpw[i],v0,v1))
          else:
              if (sdebug): print("Not plotting spw %d (%d) because it is not in baseband %d (%s)" % (i,originalSpw[i],baseband,basebandDict[baseband]))
        else:
            if (sdebug and False): print("Not plotting spw %d (%d) because fewer than 256 channels (%d)" % (i,originalSpw[i],nchan))
    if (fdmctr > -1):
        pb.ylim([y0,y1])
        pb.xlim([x0,x1])

def DrawAtmosphere(showatm, showtsky, subplotRows, atmString, mysize,
                   TebbSky, plotrange, xaxis, atmchan, atmfreq, transmission,
                   subplotCols, lo1=None, xframe=0, firstFrame=0,
                   showatmPoints=False, channels=[0], mylineno=-1,xant=-1,
                   overlaySpws=False, overlayBasebands=False, 
                   drewAtmosphere=False, loc=-1, showtsys=False, Trx=None):
    """
    Draws atmospheric transmission or Tsky on an amplitude vs. chan or freq plot.
    """
    xlim = pb.xlim()
    ylim = pb.ylim()
    myxrange = xlim[1]-xlim[0]
    yrange = ylim[1]-ylim[0]
    if (not drewAtmosphere and not overlayBasebands):  # CAS-8489 final
        if (lo1 is None):
            # add some space at the top -- Apr 16, 2012
            pb.ylim([ylim[0], ylim[1]+TOP_MARGIN*yrange])
        else:
            pb.ylim([ylim[0], ylim[1]+TOP_MARGIN*yrange*0.5])
        ylim = pb.ylim()
    yrange = ylim[1]-ylim[0]
    #
    ystartPolLabel = 1.0-0.04*subplotRows
    if (lo1 is None):
        transmissionColor = 'm'
        tskyColor = 'm'
    else:
        transmissionColor = 'k'
        tskyColor = 'k'
    if (showatmPoints):
        atmline = '.'
    else:
        atmline = '-'
    if (showatm or showtsky):
        if (showatm):
            atmcolor = transmissionColor
        else:
            atmcolor = tskyColor
        if (lo1 is None and not drewAtmosphere):
            pb.text(0.20, ystartPolLabel, atmString, color=atmcolor, size=mysize, transform=pb.gca().transAxes)
        if (showtsky):
            if showtsys:
                rescaledY = TebbSky
                if Trx == 'auto':
                    scaleFactor = np.mean([ylim[1]-np.max(TebbSky), ylim[0]-np.min(TebbSky)])
                    scaleFactor = 0.5*(ylim[1]+ylim[0]) - np.mean(TebbSky)
                    rescaledY += scaleFactor
            else:
                rescaledY, edgeYvalue, zeroValue, zeroYValue, otherEdgeYvalue, edgeT, otherEdgeT, edgeValueAmplitude, otherEdgeValueAmplitude, zeroValueAmplitude = RescaleTrans(TebbSky, ylim, subplotRows, lo1, xframe, overlaySpws, overlayBasebands)
        else:
            rescaledY, edgeYvalue, zeroValue, zeroYValue, otherEdgeYvalue, edgeT, otherEdgeT, edgeValueAmplitude, otherEdgeValueAmplitude, zeroValueAmplitude = RescaleTrans(transmission, ylim, subplotRows, lo1, xframe, overlaySpws, overlayBasebands)
        if (overlayBasebands and xaxis.find('freq')>=0):
            # use axis coordinates for y-axis only so that transmission can be on common scale
            trans = matplotlib.transforms.blended_transform_factory(pb.gca().transData, pb.gca().transAxes)
            if showtsky:
                pb.plot(atmfreq, TebbSky/300., '%s%s'%(atmcolor,atmline),
                        markeredgewidth=markeredgewidth, transform=trans)
            else:
                pb.plot(atmfreq, transmission, '%s%s'%(atmcolor,atmline),
                        markeredgewidth=markeredgewidth, transform=trans)
        else:
            # use user coordinates
            if (xaxis.find('chan')>=0):
                rescaledX = RescaleX(atmchan, xlim, plotrange, channels)
    #            rescaledX = atmchan  
                pb.plot(rescaledX, rescaledY,'%s%s'%(atmcolor,atmline),markeredgewidth=markeredgewidth)
            elif (xaxis.find('freq')>=0):
    #            print "Calling plot(xmean=%f, ymean=%f)" % (np.mean(atmfreq), np.mean(rescaledY))
                pb.plot(atmfreq, rescaledY, '%s%s'%(atmcolor,atmline),markeredgewidth=markeredgewidth)
        if (lo1 is None):
              xEdgeLabel = 1.01
        else:
            if (xframe == firstFrame):
                xEdgeLabel = -0.10*subplotCols # avoids overwriting y-axis label
            else:
                xEdgeLabel = -0.10*subplotCols
        SetNewXLimits(xlim,loc=35)  # necessary for zoom='intersect'
        if (not overlayBasebands):   # CAS-8489 final
            SetNewYLimits(ylim,loc=18)
        # Now draw the percentage on right edge of plot
        if (not drewAtmosphere):
          if (overlayBasebands and xaxis.find('freq')>=0):   # CAS-8489 final
            trans = matplotlib.transforms.blended_transform_factory(pb.gca().transData, pb.gca().transAxes)
            zeroValue = 0
            zeroValueAmplitude = 0
            edgeValueAmplitude = 1
            if (showtsky):
              edgeT = 300
              if (lo1 is None):
                  pb.text(xlim[1]+0.06*myxrange/subplotCols, edgeValueAmplitude,
                          '%.0fK'%(edgeT), color=atmcolor, size=mysize, transform=trans)
                  pb.text(xlim[1]+0.06*myxrange/subplotCols, zeroValueAmplitude,
                          '%.0fK'%(zeroValue), color=atmcolor, transform=trans,
                          size=mysize)
              else:
                  pb.text(xEdgeLabel, edgeValueAmplitude,'%.0fK'%(edgeT),
                          color=atmcolor,
                          size=mysize, transform=pb.gca().transAxes)
                  pb.text(xEdgeLabel, zeroValueAmplitude,'%.0fK'%(zeroValue),
                          color=atmcolor,
                          size=mysize, transform=pb.gca().transAxes)
            else:
              # showatm=True
              edgeT = 1
              if (lo1 is None):
                  pb.text(xlim[1]+0.05*myxrange/subplotCols, edgeValueAmplitude, 
                          '%.0f%%'%(edgeT*100), color=atmcolor, size=mysize, 
                          transform=trans, va='center')
                  pb.text(xlim[1]+0.05*myxrange/subplotCols, zeroValueAmplitude,
                          '%.0f%%'%(zeroValue*100), color=atmcolor, transform=trans,
                          size=mysize, va='center')
              else:
                  pb.text(xEdgeLabel, edgeValueAmplitude,'%.0f%%'%(edgeT*100),
                          color=atmcolor, va='center',
                          size=mysize, transform=pb.gca().transAxes)
                  pb.text(xEdgeLabel, zeroValueAmplitude,'%.0f%%'%(zeroValue*100),
                          color=atmcolor, va='center',
                          size=mysize, transform=pb.gca().transAxes)
          elif not showtsys:
            if (showtsky):
              if (lo1 is None):
                  # This must be done in user coordinates since another curve
                  # is plotted following this one.
                  pb.text(xlim[1]+0.06*myxrange/subplotCols, edgeValueAmplitude,
                          '%.0fK'%(edgeT), color=atmcolor, size=mysize)
                  pb.text(xlim[1]+0.06*myxrange/subplotCols, zeroValueAmplitude,
                          '%.0fK'%(zeroValue), color=atmcolor,
                          size=mysize)
              else:
                  # This can remain in axes units since it is the final plot.
                  pb.text(xEdgeLabel, otherEdgeYvalue,'%.0fK'%(otherEdgeT),
                          color=atmcolor,
                          size=mysize, transform=pb.gca().transAxes)
                  pb.text(xEdgeLabel, zeroYValue,'%.0fK'%(zeroValue),
                          color=atmcolor,
                          size=mysize, transform=pb.gca().transAxes)
            else:
              # showatm=True
              if (lo1 is None):
                  # This must be done in user coordinates since another curve
                  # is plotted following this one.
                  pb.text(xlim[1]+0.05*myxrange/subplotCols, edgeValueAmplitude,
                          '%.0f%%'%(edgeT*100), color=atmcolor, size=mysize)
                  pb.text(xlim[1]+0.05*myxrange/subplotCols, zeroValueAmplitude,
                          '%.0f%%'%(zeroValue*100), color=atmcolor,
                          size=mysize)
              else:
                  # This can remain in axes units since it is the final plot.
                  pb.text(xEdgeLabel, otherEdgeYvalue,'%.0f%%'%(otherEdgeT*100),
                          color=atmcolor,
                          size=mysize, transform=pb.gca().transAxes)
                  pb.text(xEdgeLabel, zeroYValue,'%.0f%%'%(zeroValue*100),
                          color=atmcolor,
                          size=mysize, transform=pb.gca().transAxes)
        if (lo1 is not None):
            if (xframe == firstFrame):
                pb.text(+1.04-0.04*subplotCols, -0.07*subplotRows,
                        'Signal SB', color='m', size=mysize,
                        transform=pb.gca().transAxes)
                pb.text(-0.03-0.08*subplotCols, -0.07*subplotRows,
                        'Image SB', color='k', size=mysize,
                        transform=pb.gca().transAxes)
    return ylim # CAS-8655

def DrawBottomLegendPageCoords(msName, uniqueTimesMytime, mysize, figfile):
    msName = msName.split('/')[-1]
    bottomLegend = msName + '  ObsDate=' + utdatestring(uniqueTimesMytime)
    if (os.path.basename(figfile).find('regression') == 0):
        regression = True
    else:
        regression = False
    if (regression == False):
        # strip off the seconds from the time to make space for casaVersion
        bottomLegend += '  plotbandpass3 v' \
                  + PLOTBANDPASS_REVISION_STRING.split()[2] + ': ' \
                  + PLOTBANDPASS_REVISION_STRING.split()[3] + ' ' \
                  + PLOTBANDPASS_REVISION_STRING.split()[4][:-3] + ', C' + casaVersion
#    The following should be used going forward, as it is better for long VLA names        
    pb.text(0.04, 0.02, bottomLegend, size=mysize, transform=pb.gcf().transFigure)
#    pb.text(0.1, 0.02, bottomLegend, size=mysize, transform=pb.gcf().transFigure)

def DrawAntennaNames(msAnt, antennasToPlot, msFound, mysize):
    for a in range(len(antennasToPlot)):
        if (msFound):
            legendString = msAnt[antennasToPlot[a]]
        else:
            legendString = str(antennasToPlot[a])
        if (a<maxAntennaNamesAcrossTheTop):
            x0 = xstartTitle+(a*antennaHorizontalSpacing)
            y0 = ystartOverlayLegend
        else:
            # start going down the righthand side
            x0 = xstartTitle+(maxAntennaNamesAcrossTheTop*antennaHorizontalSpacing)
            y0 = ystartOverlayLegend-(a-maxAntennaNamesAcrossTheTop)*antennaVerticalSpacing
        pb.text(x0, y0, legendString,color=overlayColors[a],fontsize=mysize,
                transform=pb.gcf().transFigure)
    
def stdInfo(a, sigma=3, edge=0, spw=-1, xant=-1, pol=-1):
    """
    Computes the standard deviation of a list, then returns the value, plus the
    number and list of channels that exceed sigma*std, and the worst outlier.
    """
    info = {}
    if (edge >= len(a)/2):  # protect against too large of an edge value
        originalEdge = edge
        if (len(a) == 2*(len(a)/2)):
            edge = len(a)/2 - 1 # use middle 2 points
        else:
            edge = len(a)/2  # use central point
        if (edge < 0):
            edge = 0
        print("stdInfo: WARNING edge value is too large for spw%d xant%d pol%d, reducing it from %d to %d." % (spw, xant, pol, originalEdge, edge))
    info['std'] = np.std(a[edge:len(a)-edge])
    chan = []
    outlierValue = 0
    outlierChannel = None
    for i in range(edge,len(a)-edge):
        if (np.abs(a[i]) > sigma*info['std']):
            chan.append(i)
        if (np.abs(a[i]) > np.abs(outlierValue)):
            outlierValue = a[i]
            outlierChannel = i
    info['nchan'] = len(chan)
    info['chan'] = chan
    info['outlierValue'] = outlierValue/info['std']
    info['outlierChannel'] = outlierChannel
    return(info)
    
def madInfo(a, madsigma=3, edge=0, spw=-1, xant=-1, pol=-1):
    """
    Computes the MAD of a list, then returns the value, plus the number and
    list of channels that exceed madsigma*MAD, and the worst outlier.
    """
    info = {}
    if (edge >= len(a)/2):  # protect against too large of an edge value
        originalEdge = edge
        if (len(a) == 2*(len(a)/2)):
            edge = len(a)/2 - 1 # use middle 2 points
        else:
            edge = len(a)/2  # use central point
        if (edge < 0):
            edge = 0
        print("madInfo: WARNING edge value is too large for spw%d xant%d pol%d, reducing it from %d to %d." % (spw, xant, pol, originalEdge, edge))
    info['mad'] = mad(a[edge:len(a)-edge])
    chan = []
    outlierValue = 0
    outlierChannel = None
    for i in range(edge,len(a)-edge):
        if (np.abs(a[i]) > madsigma*info['mad']):
            chan.append(i)
        if (np.abs(a[i]) > np.abs(outlierValue)):
            outlierValue = a[i]
            outlierChannel = i
    info['nchan'] = len(chan)
    info['chan'] = chan
    info['outlierValue'] = outlierValue/info['mad']
    info['outlierChannel'] = outlierChannel
    return(info)
    
def platformingCheck(a, threshold=DEFAULT_PLATFORMING_THRESHOLD):
    """
    Checks for values outside the range of +-threshold.
    Meant to be passed an amplitude spectrum.
    """
    info = {}
    startChan = len(a)/32. - 1
    endChan = len(a)*31/32. + 1
#    print "Checking channels %d-%d for platforming" % (startChan,endChan)
    if (startChan <= 0 or endChan >= len(a)):
        return
    middleChan = (startChan+endChan)/2
    channelRange1 = np.arange(startChan,middleChan+1)
    channelRange2 = np.arange(endChan,middleChan,-1)
    platforming = False
    awayFromEdge = False
    for i in channelRange1:
        if (np.abs(a[i]) > threshold):
            if (awayFromEdge):
#                print "a[%d]=%f" % (i,a[i])
                platforming = True
                return(platforming)
        else:
            awayFromEdge = True
    awayFromEdge = False
    for i in channelRange2:
        if (np.abs(a[i]) > threshold):
            if (awayFromEdge):
                platforming = True
                return(platforming)
        else:
            awayFromEdge = True
    return(platforming)
    
def mad(a, c=0.6745, axis=0):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """
    a = np.array(a)
    good = (a==a)
    a = np.asarray(a, np.float64)
    if a.ndim == 1:
        d = np.median(a[good])
        m = np.median(np.fabs(a[good] - d)) / c
#        print  "mad = %f" % (m)
    else:
        d = np.median(a[good], axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = swapaxes(a[good],0,axis)
        else:
            aswp = a[good]
        m = np.median(np.fabs(aswp - d), axis=0) / c

    return m

def plotbandpassStats(caltable, chanavg=[], channeldiff=5, title='', usetask=False, resample=True, edge=0):
    """
    Calls plotbandpass on a list of caltables with the following naming scheme:
        caltable+'_smoothXch'
    with the channeldiff option (to compute derivative statistics) and plots
    the resulting MADs vs. the level of channel averaging.
    chanavg: an integer list -- if not specified, then it will search for
             the caltables and build it automatically
    usetask: if True, use the casa task plotbandpass rather than analysisUtils version
    resample: if True (default), then resample solution appropriately before computing statistics
    """
    if (usetask):
        from plotbandpass_cli import plotbandpass_cli as plotbandpass
    caltableBase = os.path.basename(caltable)
    caltableDirectory = os.path.dirname(caltable)
    if (chanavg == []):
        if (caltableDirectory == ''):
            caltableDirectory = '.'
        myfiles = os.listdir(caltableDirectory)
        for f in myfiles:
            if (f.find(caltable) == 0 and f.find('.plots') < 0):
                tokens = f.split('_smooth')
                if (len(tokens) < 2):
                    continue
                chavg = int(tokens[1].split('ch')[0])
                chanavg.append(chavg)
        print("found %d files: %s" % (len(chanavg), str(chanavg)))
    chanavg = np.sort(chanavg)  # necessary for plotting connect-the-dots lines
    print("using chanavg = %s" % (str(chanavg)))
    stats = {}
    mytb = au.createCasaTool(tbtool)
    mytb.open(caltable+'_smooth%dch'%(chanavg[0]))
    try:
        vis = mytb.getkeyword('MSName')
        vm = au.ValueMapping(vis)
    except:
        vis = caltable.split('.')[0]
        vm = ''
    mytb.close()
    res = 1
    for ch in chanavg:
        if (resample):
            res = ch
        if (usetask):
            print("Calling plotbandpass('%s_smooth%dch',yaxis='both',interactive=False,channeldiff='%d',resample='%d',edge='%d')" % (caltable, ch, channeldiff, res, edge))
            stats[ch] = plotbandpass(caltable+'_smooth%dch'%(ch),yaxis='both',interactive=False,channeldiff=channeldiff,resample=res,edge=edge)
        else:
            print("Calling plotbandpass3('%s_smooth%dch',yaxis='both',interactive=False,channeldiff='%d',vm=vm,resample='%d',edge='%d')" % (caltable, ch, channeldiff, res, edge))
            stats[ch] = plotbandpass3(caltable+'_smooth%dch'%(ch),yaxis='both',interactive=False,channeldiff=channeldiff,vm=vm,resample=res,edge=edge)
    pb.clf()
#    pb.hold(True) # not available in CASA6, but never needed
    c=['b','c','g','r']
    showdata = True
    showmedian = True
    plots = ['phase','amp']
    for p in range(len(plots)):
        mystats = []
        pb.subplot(2,1,p+1)
        variable = plots[p]
        for ch in chanavg:
            if (showdata):
                for spw in stats[ch]['median'].keys():
                    if (spw != 'median'):
                        for pol in stats[ch]['median'][spw][0].keys():
                            pb.plot([ch],[stats[ch]['median'][spw][0][pol][variable]],'.',color=c[spw],markeredgewidth=markeredgewidth)
                    else:
                        pb.plot([ch],[stats[ch]['median'][spw][variable]],'+',color='k',markeredgewidth=markeredgewidth)
                        mystats.append(stats[ch]['median'][spw][variable])
        if (showmedian):
            pb.plot(chanavg,mystats,'o-',color='k',markeredgewidth=markeredgewidth)
        pb.xlabel('channel average')
        pb.ylabel('MAD of %s derivative' % (variable))
        spw = stats[ch]['median'].keys()[0]  # get first spw
        if (spw == 'median'):
            spw = stats[ch]['median'].keys()[1]
        if (p==0):
            if (vm == ''):
                pb.title('%s, %s' % (vis.split('.')[0], title))
            else:
                pb.title('%s (spw %d = %.0f MHz, %d channels)' % (vis.split('.')[0], spw,
                                                                  vm.spwInfo[spw]['bandwidth']*1e-6,
                                                                  vm.spwInfo[spw]['numChannels']))
    pb.draw()
    if (resample):
        pngname = caltableBase + '.resample.stats.png'
    else:
        pngname = caltableBase + '.stats.png'
    pb.savefig(pngname)
    print("Left plot in: %s" % (pngname))

def callFrequencyRangeForSpws(mymsmd, spwlist, vm, debug=False, caltable=None):
    """
    Returns the min and max frequency of a list of spws.
    """
    if (mymsmd != '' and casaVersion >= '4.1.0'):
        return(frequencyRangeForSpws(mymsmd,spwlist))
    else:
        freqs = []
        if (type(vm) != str):
            if (debug):
                print("vm.spwInfo.keys() = ", vm.spwInfo.keys())
            for spw in spwlist:
                freqs += list(vm.spwInfo[spw]["chanFreqs"])
        else:
            mytb = au.createCasaTool(tbtool)
            try:
                mytb.open(caltable+'/SPECTRAL_WINDOW')
                chanfreq = []
                if (debug):
                    print("getting frequencies of spws = ", originalSpws)
                if (len(spwlist) == 0): # CAS-8489b
                    originalSpws = range(len(mytb.getcol('MEAS_FREQ_REF')))
                    spwlist = originalSpws
                for i in spwlist:   # CAS-8489b
                    # The array shapes can vary.
                    chanfreq.append(mytb.getcell('CHAN_FREQ',i))
                for cf in chanfreq:
                    freqs += list(cf)
                mytb.close()
            except:
                pass
        if (freqs == []):
            return(0,0)
        else:
            return(np.min(freqs)*1e-9, np.max(freqs)*1e-9)

def frequencyRangeForSpws(mymsmd, spwlist):
    """
    Returns the min and max frequency of a list of spws.
    """
    allfreqs = []
    for spw in spwlist:
        allfreqs += list(mymsmd.chanfreqs(spw))
    if (len(allfreqs) == 0):
        print("len allfreqs = zero, spwlist = %s" % (str(spwlist)))
        return(0,0)
    return(np.min(allfreqs)*1e-9, np.max(allfreqs)*1e-9)

def buildSpwString(overlaySpws, overlayBasebands, spwsToPlot, ispw, originalSpw,
                   observatoryName, baseband, showBasebandNumber):
    if (overlayBasebands):
        spwString = ' all'
    elif (overlaySpws and len(spwsToPlot)>1):
        if (observatoryName.find('ALMA') >= 0 or observatoryName.find('ACA') >= 0):
            # show a list of all spws
            spwString = str(spwsToPlot).replace(' ','').strip('[').strip(']')
        else:
            # show the range of spw numbers
            spwString = '%2d-%2d' % (np.min(spwsToPlot),np.max(spwsToPlot))
    elif (ispw==originalSpw):
        spwString = '%2d' % (ispw)
    else:
        spwString = '%2d (%d)' % (ispw,originalSpw)
    if (overlayBasebands==False):
        spwString=appendBasebandNumber(spwString, baseband, showBasebandNumber)
    return(spwString)

def appendBasebandNumber(spwString, baseband, showBasebandNumber):
    if (showBasebandNumber):
        spwString += ', bb%d' % (baseband)
    return(spwString)

def getSpwsForBaseband(bb, vis=None, mymsmd=None, nonchanavg=True, caltable=None):
    needToClose = False
    if (casadef.subversion_revision >= '25753' and vis is not None):
        if (os.path.exists(vis)):
            if (mymsmd is None or mymsmd == ''):
                needToClose = True
                mymsmd = au.createCasaTool(msmdtool)
                mymsmd.open(vis)
                s = mymsmd.spwsforbaseband(bb)
            else:
                s = mymsmd.spwsforbaseband(bb)
        else:
            s = getBasebandDict(vis=vis, caltable=caltable, mymsmd=mymsmd)[bb]
    else:
        s = getBasebandDict(vis=vis, caltable=caltable, mymsmd=mymsmd)[bb]
    spws = []
    for spw in s:
        if (mymsmd.nchan(spw) > 1 or nonchanavg==False):
            spws.append(spw)   
    if needToClose:
        mymsmd.close()
    return(spws)
        
def getBasebandDict(vis=None, spwlist=[], caltable=None, mymsmd=None):
    """
    Builds a dictionary with baseband numbers as the keys and the
    associated spws as the values.  The optional parameter spwlist can
    be used to restrict the contents of the dictionary.
    Note: This is obsoleted by msmd.spwsforbaseband(-1)
    """
    bbdict = {}
    if (vis is not None):
        if (os.path.exists(vis)):
            bbs = au.getBasebandNumbers(vis)
        elif (caltable is not None):
            bbs = au.getBasebandNumbersFromCaltable(caltable)
        else:
            print("Must specify either vis or caltable")
            return
    elif (caltable is not None):
        bbs = au.getBasebandNumbersFromCaltable(caltable)
    else:
        print("Must specify either vis or caltable")
        return
    if (type(bbs) == int):  # old datasets will bomb on msmd.baseband()
        return(bbdict)
    needToClose = False
    if (casaVersion >= '4.1.0' and vis is not None):
        if (os.path.exists(vis)):
            if mymsmd is None or mymsmd == '':
                needToClose = True
                mymsmd = au.createCasaTool(msmdtool)
                mymsmd.open(vis)
            if (spwlist == []):
                nspws = mymsmd.nspw()
                spwlist = range(nspws)
            for spw in spwlist:
                bbc_no = mymsmd.baseband(spw)
                if (bbc_no not in bbdict.keys()):
                    bbdict[bbc_no] = [spw]
                else:
                    bbdict[bbc_no].append(spw)
            if needToClose:
                mymsmd.close()
    if (bbdict == {}):
        # read from spw table
        ubbs = np.unique(bbs)
        for bb in ubbs:
            bbdict[bb] = []
        for i in range(len(bbs)):
            bbdict[bbs[i]].append(i)
    return(bbdict)

#def getScansForTimes(mymsmd, scantimes):
#    myscans = []
#    myscantimes = []
#    scantimes = au.splitListIntoContiguousLists(scantimes)
#    for t in scantimes:
#        mean_t = np.mean(t)
#        range_t = (1+np.max(t)-np.min(t))*0.5
#        scans_t = mymsmd.scansfortimes(mean_t, range_t)
#        if (len(scans_t) > 0):
#            scan = scans_t[0]
#            myscans.append(scan)
#            myscantimes.append(t)
#    return(myscans, myscantimes)
