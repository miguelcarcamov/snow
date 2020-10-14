# Script to image and assess the properties of long baseline calibrators.
# Runs in CASA 4.5.0
# Expects to find the *.split.cal measurement set and the .fluxscale file.
# Unless you want to limit spws (i.e. exclude very narrow ones for
# speed) nothing should need to be changed.

# If the analysis fails (usually only on check source) it's an
# indication that the source is non-point like. The image and png
# should be created regardless.
#
# C. Brogan (Nov 2015)
# T. Hunter (Jan 2016)

###################################
# DATA PROPERTIES
###################################
from __future__ import print_function  # prevents adding old-style print statements
import numpy as np
import pylab as pl
import analysisUtils as au
import glob 
from taskinit import *
from imfit_cli import imfit_cli as imfit
from matplotlib.ticker import MultipleLocator # used by plotPointingResults
try:
    from tclean_cli import tclean_cli as tclean
except:
    print("checksource.py: Cannot import tclean")

def writeOut(f,line):
    print(line)
    f.write(line+'\n')

def version(short=False):
    """
    Returns the CVS revision number as a string.
    """
    myversion = "$Id: checksource.py,v 1.22 2019/03/05 15:26:27 thunter Exp $"
    if (short):
        myversion = myversion.split()[2]
    return myversion

def checksource(overwrite=True, verbose=False, subdir='', splitcal_vis=''):
    """
    Images the phasecal and check source in a manually-calibrated dataset and 
    reports statistics.  Expects to find the *.split.cal measurement set and 
    the corresponding .fluxscale file for it.
    Inputs:
    overwrite: if True, overwrite any existing image from a previous execution
    splitcal_vis: defaults to *.cal, but can be specified as list of strings, 
                  or a comma-delimited string
    Outputs:
    png image plots of each calibrator, and an ASCII file for each dataset
    The name of the ASCII file, and a list of pngs are returned.
    """
    # Read the dataset(s) and get properties
    if (splitcal_vis == ''):
        vislist = glob.glob('*.cal')
    else:
        if (type(splitcal_vis) == str):
            vislist = splitcal_vis.split(',')
        else:
            vislist = splitcal_vis
    print("Checking datasets: ", vislist)
    mymsmd = au.createCasaTool(msmdtool)
    if (len(subdir) > 0):
        if (os.path.exists(subdir)):
            if (subdir[-1] != '/'): 
                subdir += '/'
        else:
            os.mkdir(subdir)
            if (subdir[-1] != '/'): 
                subdir += '/'
    pnglist = []
    textfiles = []
    for vis in vislist:
        mymsmd.open(vis)
        freq=mymsmd.meanfreq(0,unit='GHz')
        # Check Source
        check=mymsmd.fieldsforintent('OBSERVE_CHECK_SOURCE*',True)[0]
        checkid=mymsmd.fieldsforintent('OBSERVE_CHECK_SOURCE*',False)[0]
        checkpos=mymsmd.phasecenter(checkid)
        # Phase calibrator
        phase=mymsmd.fieldsforintent('CALIBRATE_PHASE*',True)[0]
        phaseid=mymsmd.fieldsforintent('CALIBRATE_PHASE*',False)[0]
        phasepos=mymsmd.phasecenter(phaseid)
        if ('OBSERVE_TARGET#ON_SOURCE' in mymsmd.intents()):
            nScienceFields= len(mymsmd.fieldsforintent('OBSERVE_TARGET*',False))
            science = mymsmd.fieldsforintent('OBSERVE_TARGET*',True)[0]
            scienceid = mymsmd.fieldsforintent('OBSERVE_TARGET*',False)[0]
        else:
            nScienceFields = 0
        mymsmd.done()

        floatcell = au.pickCellSize(vis, maxBaselinePercentile=99, 
                                    verbose=verbose)
        cell = au.pickCellSize(vis, maxBaselinePercentile=99, cellstring=True, 
                               verbose=verbose)
#        imsize = int(au.nextValidImsize(int(5.0/floatcell))) # valid when we only had checksources for synthBeam < 0.25
        imsize = int(au.nextValidImsize(int(np.max([5.0,5.0*au.estimateSynthesizedBeam(vis)])/floatcell))) 
        print("imsize = ", imsize)
        region='circle[[%dpix , %dpix], 15pix ]' % (int(imsize/2),int(imsize/2))

        if False:
            # original method (for bands 3-6 only)
            cell = str(np.round(0.015*(100/freq),3))+'arcsec'
            if freq < 116.0:
                imsize = [320,320]
                region='circle[[160pix , 160pix] ,15pix ]'
            else:
                imsize = [680,680]
                region='circle[[340pix , 340pix] ,15pix ]'

        ###################################
        # IMAGE 
        ###################################
        weighting = 'briggs'
        robust = 0.5
        niter = 50
        threshold = '0.0mJy'
        spw=''
        separation = au.angularSeparationOfTwoFields(vis,checkid,phaseid)
        if (nScienceFields > 0):
            separation_pcal_science = au.angularSeparationOfTwoFields(vis,scienceid,phaseid)
            separation_check_science = au.angularSeparationOfTwoFields(vis,scienceid,checkid)

        fieldtype = ['checksource','phasecal']
        field = [check,phase]
        for i,cal in enumerate(field):
            if (not os.path.exists(cal+'_'+vis+'.image') or overwrite):
                os.system('rm -rf '+cal+'_'+vis+'.*')
                if verbose:
                    print("Running tclean('%s', field='%s', cell=%s, imsize=%s, ...)" % (vis, cal, str(cell), str(imsize)))
                tclean(vis=vis,
                       imagename=cal+'_'+vis,
                       field=cal,spw=spw,
                       specmode='mfs',
                       deconvolver='hogbom',
                       imsize = imsize, 
                       cell= cell, 
                       weighting = weighting, 
                       robust = robust,
                       niter = niter, 
                       threshold = threshold, 
                       interactive = False,
                       mask = region,
                       gridder = 'standard')
            png = subdir+fieldtype[i]+'_'+cal+'_'+vis+'.image.png'
            pnglist.append(png)
            au.imviewField(cal+'_'+vis+'.image',radius=30*floatcell,
                           contourImage=cal+'_'+vis+'.mask',levels=[1],
                           plotfile=png)


        ###################################
        # ANALYZE
        ###################################
        ###########
        # PHASE
        ###########
        imagename=phase+'_'+vis
        if verbose:
            print("Running imfit('%s', region='%s')" % (imagename+'.image', region))
        # Fit the phase source to get position and flux
        imagefit=imfit(imagename=imagename+'.image',
                       region=region)      
        fitresults=au.imfitparse(imagefit)

        # Compare the Positions
        phasepos_obs=au.direction2radec(phasepos)
        if fitresults is not None:
            phasepos_fit=','.join(fitresults.split()[:2])
            phasepos_diff=au.angularSeparationOfStrings(phasepos_obs,phasepos_fit,verbose=False)*3600.

        # Compare the Flux densities
        peakIntensity = au.imagePeak(imagename+'.image')
        selffluxfile=glob.glob('*.fluxscale')[0]
        fluxscaleResult = au.fluxscaleParseLog(selffluxfile,field=phase)
        if fluxscaleResult is not None:
            selfflux = fluxscaleResult[0][0]
            phaseflux_fit=float(fitresults.split()[2])
            phaseCoherence = 100*peakIntensity/phaseflux_fit
            phaseflux_diff=100*(selfflux-phaseflux_fit)/selfflux

        # Print the final results and save to file
        textfile = subdir+'calimage_results_'+vis+'.txt'
        textfiles.append(textfile)
        f = open(textfile,'w')
        f.write('\n*************************************************************\n\n')
        line = 'CHECK_SOURCE IMAGE ANALYSIS REPORT (version %s)\n' % version(short=True)
        writeOut(f,line)
        info = au.getFitsBeam(imagename+'.image')
        synthBeam = (info[0]*info[1])**0.5
        if fitresults is None:
            line = "Phasecal %s: imfit failed" % (phase)
        elif fluxscaleResult is not None:
            line= "Phasecal %s: Position difference = %s arcsec = %s synth.beam, Flux %% difference = %s"%(phase,au.roundFiguresToString(phasepos_diff,3), au.roundFiguresToString(phasepos_diff/synthBeam,3), au.roundFiguresToString(phaseflux_diff,3))
            writeOut(f,line)
            line = "    coherence = peakIntensity/fittedFluxDensity = %s%%" % (au.roundFiguresToString(phaseCoherence,3))
        else:
            line = "Phasecal %s: Position difference = %s arcsec = %s synth.beam" % (phase,au.roundFiguresToString(phasepos_diff,3), au.roundFiguresToString(phasepos_diff/synthBeam,3))
        writeOut(f,line)
        f.close()
        if fluxscaleResult is None:
            print("Full checksource analysis is not supported if there is no flux calibrator")
            return textfiles, pnglist

        ###########
        # CHECK
        ###########
        imagename=check+'_'+vis
        # Fit the check source to get position and flux
        if verbose:
            print("Running imfit('%s', region='%s')" % (imagename+'.image', region))
        imagefit=imfit(imagename=imagename+'.image',
                       region=region)      
        fitresults=au.imfitparse(imagefit, deconvolved=True)
        info = au.getFitsBeam(imagename+'.image')
        synthMajor, synthMinor = info[0:2]
        synthBeam = (info[0]*info[1])**0.5

        # Compare the Positions
        checkpos_obs=au.direction2radec(checkpos)
        if fitresults is not None:
            checkpos_fit=','.join(fitresults.split()[:2])
            checkpos_diff=au.angularSeparationOfStrings(checkpos_obs,checkpos_fit,
                                                        verbose=False)*3600.

        # Compare the Flux densities
        selffluxfile=glob.glob('*.fluxscale')[0]
        results = au.fluxscaleParseLog(selffluxfile,field=check)
        peakIntensity = au.imagePeak(imagename+'.image')
        if (results is not None and fitresults is not None):
            selfflux=results[0][0] 
            checkflux_fit=float(fitresults.split()[2])

            checkflux_diff=100*(selfflux-checkflux_fit)/selfflux
            checkCoherence = 100*peakIntensity/checkflux_fit
        if fitresults is not None:
            if verbose: 
                print("Checksource fitresults: ", fitresults)
            deconvolvedMajor = float(fitresults.split()[5])
            deconvolvedMinor = float(fitresults.split()[7])

        # Print the final results and save to file
        f=open(textfile,'a')
        if fitresults is None:
            line = "Checksource %s: imfit failed" % (phase)
        else:
            if (results is not None):
                line= "\nChecksource %s: Position difference = %s arcsec = %s synth.beam, Flux %% difference = %s"%(check ,au.roundFiguresToString(checkpos_diff,3),au.roundFiguresToString(checkpos_diff/synthBeam,3),au.roundFiguresToString(checkflux_diff,3))
                writeOut(f,line)
                line = "    coherence = peakIntensity/fittedFluxDensity = %s%%" % (au.roundFiguresToString(checkCoherence,3))
            else:
                line= "\nChecksource %s: Position difference = %s arcsec = %s synth.beam" % (check ,au.roundFiguresToString(checkpos_diff,3),au.roundFiguresToString(checkpos_diff/synthBeam,3))
            writeOut(f,line)
            line = "    beam size = %s x %s arcsec" % (au.roundFiguresToString(synthMajor,3), au.roundFiguresToString(synthMinor,3))
            writeOut(f,line)
            line = "    apparent deconvolved size = %s x %s arcsec = %s synth.beam area" % (au.roundFiguresToString(deconvolvedMajor,2), au.roundFiguresToString(deconvolvedMinor,2), au.roundFiguresToString(deconvolvedMajor*deconvolvedMinor/(synthBeam**2),2))
        writeOut(f,line)
        line = "    angular separation of phasecal to checksource = %s degree" % (au.roundFiguresToString(separation,3))
        writeOut(f,line)
        if (nScienceFields > 0):
            if (nScienceFields > 1):
                modifier = 'first'
            else:
                modifier = 'only'
            line = "    angular separation of phasecal to %s science field (%d) = %s degree" % (modifier,scienceid,au.roundFiguresToString(separation_pcal_science,3))
            writeOut(f,line)
            line = "    angular separation of checksource to %s science field (%d) = %s degree" % (modifier,scienceid,au.roundFiguresToString(separation_check_science,3))
            writeOut(f,line)
        f.close()
    # end 'for' loop over vislist
    return textfiles, pnglist

def offset(workingdir, vis='', plotfile='', imfitlog=False, spw='', verbose=False):
    """
    Takes a pipeline working directory and find all images of the checksource 
    and produces a plot showing the relative directions of the first two science 
    targets, the phase calibrator, and the checksource, and a vector
    showing the offset of the checksource from its catalog position (computed
    using the results of the CASA task imfit), and a
    text label showing the RAO and DECO offsets.
    workingdir: path to pipeline working directory
    vis: alternate location for a measurement set to consult (ignores *_target.ms)
    Looks first for *chk*iter2.image; if not found, then *chk*iter1.image
    plotfile: default = img+'_offset.png'
    imfitlog: if True, then request imfit to generate log files (*.imfit)
    spw: int or comma-delimited string, if specified, limit to this or these spws
    verbose: print more messages explaining what images it is operating on
    """
    mymsmd = au.createCasaTool(msmdtool)
    if verbose:
        print("workingdir: ", workingdir)
    imglist = sorted(glob.glob(os.path.join(workingdir,'*_chk.spw*image')))
    if len(imglist) == 0:
        print("No check source images found in this directory.")
        return
    # If iter2.image is found, then drop the iter1 version from the list
    for i in imglist: 
        if i.find('iter2') > 0:
            imglist.remove(i.replace('iter2','iter1'))
    if verbose:
        print("Processing %d images:" % (len(imglist)))
        for i in imglist: 
            print(i)
    if vis == '':
        searchpath = os.path.join(workingdir,'*.ms')
        if verbose:
            print("searchpath: ", searchpath)
        allvislist = sorted(glob.glob(searchpath))
        if verbose:
            print("all vis found: " , allvislist)
        vislist = []
        for vis in allvislist:
            if vis.find('_target') < 0:
                vislist.append(vis)
    else:
        vislist = [vis]

    raos = []
    decos = []
    totals = []
    sourcenames = []
    spws = au.parseSpw(vis, spw)
    scienceSpws = au.getScienceSpws(vis, returnString=False)
    spws = np.intersect1d(scienceSpws,spws)
    if verbose:
        print("using spws: ", spws)
    newimglist = []
    for img in imglist: # there will be an image for each spw
        if img.find('spw') > 0 and spw != '':
            myspw = int(img.split('spw')[1].split('.')[0])
            if myspw in spws:
                sourcenames.append(au.imageSource(img))
                newimglist.append(img)
                if verbose:
                    print("Using %s" % (img))
            elif verbose:
                print("Skipping %s" % (img))
        else:
            sourcenames.append(au.imageSource(img))
            newimglist.append(img)
    sourcenames = np.unique(sourcenames)
    pngs = []
    print("vislist = ", vislist)
    imglist = newimglist
    for sourcename in sourcenames:
        for ispw, img in enumerate(imglist): # there will be an image for each spw
            if 'spw' not in img:
                print("No spw in the image name: ", img)
                continue
            spw = int(img.split('spw')[1].split('.')[0])
            # find the first vis that observed this target as check source
            checkid = -1
            for vis in vislist:
#                print "Checking ", vis
                mymsmd.open(vis)
                if spw >= mymsmd.nspw():
                    print("Guessing that spw %d is spw %d in the split ms." % (spw,ispw))
                    spw = ispw
                if 'OBSERVE_CHECK_SOURCE#ON_SOURCE' in mymsmd.intents():
                    checksources = mymsmd.fieldsforintent('OBSERVE_CHECK_SOURCE*',True)
                else:
                    checksources = mymsmd.fieldsforintent('CALIBRATE_DELAY*',True)
                if sourcename in checksources:
                    check = checksources[0]
                    checkid = mymsmd.fieldsforname(sourcename)[0]
                    checkpos = mymsmd.phasecenter(checkid)
                    # Phase calibrator
                    phase = mymsmd.fieldsforintent('CALIBRATE_PHASE*',True)[0]
                    phaseid = mymsmd.fieldsforintent('CALIBRATE_PHASE*',False)[0]
                    phasepos = mymsmd.phasecenter(phaseid)
                    if ('OBSERVE_TARGET#ON_SOURCE' in mymsmd.intents()):
                        nScienceFields = len(mymsmd.fieldsforintent('OBSERVE_TARGET*',False))
                        science = mymsmd.fieldsforintent('OBSERVE_TARGET*',True)[0]
                        scienceid = mymsmd.fieldsforintent('OBSERVE_TARGET*',False)[0]
                        sciencepos = mymsmd.phasecenter(scienceid)
                        if nScienceFields > 1:
                            science2 = mymsmd.fieldsforintent('OBSERVE_TARGET*',True)[1]
                            science2id = mymsmd.fieldsforintent('OBSERVE_TARGET*',False)[1]
                            science2pos = mymsmd.phasecenter(science2id)
                    else:
                        nScienceFields = 0
                    rxBand = mymsmd.namesforspws(spw)[0].split('#')[1].split('_')[-1].lstrip('0') # string
                    break
                else:
                    mymsmd.close()
            if checkid < 0:
                print("Could not find an ms that observed this check source: %s. Try including the vis parameter." % (sourcename))
                continue
            info = au.getFitsBeam(img)
            imsize = info[5]  # size in RA direction
            region = 'circle[[%dpix , %dpix], 15pix ]' % (int(imsize/2),int(imsize/2))
            freq = mymsmd.meanfreq(spw,unit='GHz')
            if imfitlog:
                logfile = img + '.imfit'
            else:
                logfile = ''
            imagefit = imfit(imagename=img, region=region, logfile=logfile)      
            fitresults = au.imfitparse(imagefit, deconvolved=True)
            synthMajor, synthMinor = info[0:2]
            synthBeam = (info[0]*info[1])**0.5
            # Compare the Positions
            checkpos_obs = au.direction2radec(checkpos)
            if fitresults is not None:
                checkpos_fit = ','.join(fitresults.split()[:2])
                print("spw %d: checksource fitted position: " % (spw), checkpos_fit)
                result = au.angularSeparationOfStrings(checkpos_fit, checkpos_obs, True, verbose=False)
                checkpos_diff, deltaLong, deltaLat, deltaLongCosDec, pa = result
            total = checkpos_diff*3600.
            rao = deltaLongCosDec*3600.
            deco = deltaLat*3600.
            print("spw %d: %s offset=%.4f arcsec, RAO=%+.4f, DECO=%+.4f, PA=%.1fdeg" % (spw, sourcename, total, rao, deco, pa))
            totals.append(total)
            raos.append(rao)
            decos.append(deco)
            mymsmd.close()
            if nScienceFields > 1:
                scienceDeg = np.degrees(au.angularSeparationOfDirections(science2pos,sciencepos,True))
            phaseDeg = np.degrees(au.angularSeparationOfDirections(phasepos,sciencepos,True))
            checkDeg = np.degrees(au.angularSeparationOfDirections(checkpos,sciencepos,True))
            if len(raos) == 1:
                pl.clf()
                desc = pl.subplot(111)
                if nScienceFields > 1:
                    pl.plot([0, scienceDeg[3], phaseDeg[3], checkDeg[3]], 
                            [0, scienceDeg[2], phaseDeg[2], checkDeg[2]], 'b+', ms=10, mew=2)
                else:
                    pl.plot([0, phaseDeg[3], checkDeg[3]], [0,phaseDeg[2],checkDeg[2]], 'b+', ms=10, mew=2)
                pl.hold(True)
                pl.axis('equal')
                yrange = np.diff(pl.ylim())[0]
                # reverse RA axis
                x0,x1 = pl.xlim()
                xoffset = 0.15*(x1-x0)
                # Keep a fixed scale among the spws/images
                xscale = 0.5*xoffset/np.max(np.abs([rao,deco]))
            # draw the arrow for each spw's image
            pl.arrow(checkDeg[3], checkDeg[2], rao*xscale, deco*xscale, lw=1, shape='full',
                     head_width=0.15*xoffset, head_length=0.2*xoffset, fc='b', ec='b')
            if len(raos) == 1:
                pl.xlim([x1+xoffset, x0-xoffset])
                yoffset = yrange*0.025
                pl.text(0,                      0+yoffset, 'science', ha='center',va='bottom')
                if nScienceFields > 1:
                    pl.text(scienceDeg[3], scienceDeg[2]+yoffset, 'science (%.1fdeg)'%scienceDeg[0], ha='center',va='bottom')
                    pl.text(scienceDeg[3], scienceDeg[2]-yoffset, science2, ha='center',va='top')
                pl.text(phaseDeg[3], phaseDeg[2]+yoffset, 'phase (%.1fdeg)'%phaseDeg[0], ha='center',va='bottom')
                pl.text(checkDeg[3], checkDeg[2]+yoffset, 'check (%.1fdeg)'%checkDeg[0], ha='center',va='bottom')
                pl.text(0,                      0-yoffset, science, ha='center',va='top')
                pl.text(phaseDeg[3], phaseDeg[2]-yoffset, phase, ha='center',va='top')
                pl.text(checkDeg[3], checkDeg[2]-yoffset, check, ha='center',va='top')
                pl.xlabel('RA offset (deg)')    
                pl.ylabel('Dec offset (deg)')
                projCode = au.projectCodeFromDataset(vis)
                if type(projCode) == str:
                    if verbose:
                        print("Did not find project code")
                    projCode = ''
                else:
                    projCode = projCode[0] + ', Band %s, ' % (rxBand)
                pl.title(projCode + os.path.basename(img).split('.spw')[0] + ', spws=%s'%spws, size=12)
                pl.ylim([pl.ylim()[0]-yoffset*8, pl.ylim()[1]+yoffset*8]) 
                minorLocator = MultipleLocator(0.5) # degrees
                desc.xaxis.set_minor_locator(minorLocator)
                desc.yaxis.set_minor_locator(minorLocator)
        # end 'for' loop over spws/images
        if len(raos) < 1:
            return
        pl.ylim([pl.ylim()[0]-yoffset*7, pl.ylim()[1]+yoffset*15]) 
        rao = np.median(raos)
        raostd = np.std(raos)
        deco = np.median(decos)
        decostd = np.std(decos)
        total = np.median(totals)
        totalstd = np.std(totals)
        raoBeams = rao / synthBeam
        raostdBeams = raostd / synthBeam
        decoBeams = deco / synthBeam
        decostdBeams = decostd / synthBeam
        # draw the median arrow in thick black line
        pl.arrow(checkDeg[3], checkDeg[2], rao*xscale, deco*xscale, lw=2, 
                 shape='full', head_width=0.12*xoffset, 
                 head_length=0.18*xoffset, ec='k', fc='k')
        print("median +- std: offset=%.4f+-%.4f, RAO=%.4f+-%.4f, DECO=%.4f+-%.4f" % (total,totalstd,rao,raostd,deco,decostd))
#        pl.text(checkDeg[3], checkDeg[2]-0.6*xoffset, '$\Delta\\alpha$: %+.4f"+-%.4f"' % (rao,raostd), ha='center')
#        pl.text(checkDeg[3], checkDeg[2]-0.85*xoffset, '$\Delta\\delta$: %+.4f"+-%.4f"' % (deco,decostd), ha='center')
        pl.text(0.05,0.95, '$\Delta\\alpha$: %+.4f"+-%.4f" = %+.2f+-%.2f beams' % (rao,raostd,raoBeams,raostdBeams), ha='left', transform=desc.transAxes)
        pl.text(0.05,0.91, '$\Delta\\delta$: %+.4f"+-%.4f" = %+.2f+-%.2f beams' % (deco,decostd,decoBeams,decostdBeams), ha='left', transform=desc.transAxes)
        if plotfile == '':
            png = img + '_offset.png'
        else:
            png = plotfile
        pl.savefig(png, bbox_inches='tight')
        pl.draw()
        pngs.append(png)
        print("Wrote ", png)
