#!/usr/bin/env python
#
# tsysNormalize: function for normalizing a Tsys cal table to a common elevation
# (the elevation of the first Tsys scan). 
# Based on Neil Phillips tsystransfer.py (originally posted in CAS-4636)
# Todd Hunter (July 12, 2016)
from __future__ import print_function  # prevents adding old-style print statements
import os
import taskinit
import time
import copy
import numpy as np
import shutil
import pylab as pb
#from rmtables import rmtables  # fails in CASA6, but not needed
from analysisUtils import getChannelAveragedScienceSpws, getBasebands

def TsysAfterPowerChange(Vorig, Vnew, TsysOrig=200.0, Tatm=270.0):
    # The maths to derive a new Tsys
    x = (Vnew / Vorig) * (TsysOrig / (TsysOrig + Tatm))
    return Tatm * x / (1.0 - x)

def getPower(vis, scan, spw, duration, fromEnd=False, skipStartSecs=1.0, 
             skipEndSecs=1.0, verbose=True):
    """
    Return a per-antenna list of total power values for the two polarizations of the specified scan and spw.
    duration: number of samples to use starting from the start of the scan
    I think the idea here is that the sky subscan is always the first subscan of a Tsys scan.  If this ever
    changes to one of the loads, then the result will be less than optimal. It would likely result in very
    small changes from the original tsys table (i.e. it will not get normalized).   - Todd
    """
    myms = taskinit.mstool()
    myms.open(vis)
    myms.selectinit(datadescid=spw)
    myms.selecttaql('SCAN_NUMBER==%d AND DATA_DESC_ID==%d AND ANTENNA1==ANTENNA2'%(scan, spw))
#    nrows = myms.nrow()
    print("    Working on spw: %d" % (spw))
    d = myms.getdata(['real','axis_info'], ifraxis=True)
    myms.close()
    if verbose:
        print("keys = ", d.keys())
    if 'real' not in d.keys():
        return
    powers = d['real']
    ants = d['axis_info']['ifr_axis']['ifr_name']
    pols = list(d['axis_info']['corr_axis'])
    idxPol0 = pols.index("XX")
    idxPol1 = pols.index("YY")
    if verbose:
        print("Pol 0,1 indexes: %d, %d" % (idxPol0, idxPol1))
    ts = d['axis_info']['time_axis']['MJDseconds']
    t0 = ts[0]
    tf = ts[-1]
    ts -= t0
    if verbose:
        print("times:", ts)
        print("pols:", pols)
    # choose the time samples we want
    sampStart = -99
    sampEnd = -99
    for i in range(len(ts)):
        if ts[i] > skipStartSecs:
            sampStart = i
            break
    if sampStart >= len(ts):
        sampStart = len(ts) - 1
    for i in range(len(ts)-1, sampStart, -1):
        if tf - ts[i] > skipEndSecs:
            sampEnd = i
            break
    if sampEnd <= sampStart:
        sampEnd = sampStart + 1
    if not fromEnd:
        # take duration samples from start
        for i in range(sampStart+1, sampEnd, 1):
            if ts[i] - ts[sampStart] > duration:
                sampEnd = i
                break
    else:
        # instead from end
        for i in range(sampEnd-1, sampStart, -1):
            if ts[sampEnd] - ts[i] > duration:
                sampStart = i
                break
    if verbose:
        print("chosen sample range: %d to %d" % (sampStart, sampEnd))
    # indexing is pol, baseline(=0), ant, sample
    if verbose:
        print("number of antennas to produce powers for:", len(ants))
    result = []
    for ant in range(len(ants)):
        powersPol0 = powers[idxPol0][0][ant][sampStart:sampEnd]
        powersPol1 = powers[idxPol1][0][ant][sampStart:sampEnd]
        #print "Ant %d powers pol 0: %s, pol 1: %s" % (ant, str(powersPol0), str(powersPol1))
        medianP0 = np.median(powersPol0)
        medianP1 = np.median(powersPol1)
        result.append([medianP0,medianP1])
        #print "Ant %2d (%s) median powers for pols 0,1: %12.6f, %12.6f (nSamples = %d, %d)" % (ant, ants[ant], medianP0, medianP1, len(powersPol0), len(powersPol1))
    return result

def scienceSpwForTsysSpw(mymsmd, tsysSpw):
    """
    Automatically pick one science spw for the specified Tsys spw.  
    Ideally it would be the one with the widest bandwidth, but it will likely make very little
    difference, so just take the first one if there are more than one.
    """
    baseband = mymsmd.baseband(tsysSpw)
    spws = np.intersect1d(mymsmd.almaspws(tdm=True,fdm=True),mymsmd.spwsforbaseband(baseband))
    if 'OBSERVE_TARGET#ON_SOURCE' in mymsmd.intents():
        spws = np.intersect1d(mymsmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE'), spws)
    elif 'CALIBRATE_FLUX#ON_SOURCE' in mymsmd.intents():
        spws = np.intersect1d(mymsmd.spwsforintent('CALIBRATE_FLUX#ON_SOURCE'), spws)
    elif 'CALIBRATE_AMPLI#ON_SOURCE' in mymsmd.intents():
        spws = np.intersect1d(mymsmd.spwsforintent('CALIBRATE_AMPLI#ON_SOURCE'), spws)
        
    if (len(spws) > 1):
        print("Multiple science spws for this tsys spw (%d).  Picking the first one (%d)" % (tsysSpw, spws[0]))
        return(spws[0])
    elif (len(spws) < 1):
        print("No science spws for this tsys spw (%d)" % (tsysSpw))
    else:
        return(spws[0])
    return -1

def tsysNormalize(vis, tsysTable='', newTsysTable='', scaleSpws=[], 
                  verbose=False):
    """
    Generate Tsys entries for one field from other fields, using autocorr
    (linear!) or SQLD data to determine the change in Tsys.
    Inputs:
     vis          the MS
     tsysTable:  the tsys caltable (default = <vis>.tsys)
     newTsysTable:  the new tsys caltable to create (default = <tsysTable>_normalized)
    """

    # intents likely to imply different attenuations or tuning to science-like
    # scans that we are applying Tsys to.
    print("Entered")
    badIntents = ['CALIBRATE_POINTING', 'CALIBRATE_FOCUS', 'CALIBRATE_SIDEBAND_RATIO', 'CALIBRATE_ATMOSPHERE']
    if (tsysTable == ''):
        tsysTable = vis + '.tsys'
    if (not os.path.exists(tsysTable)):
        print("Cannot find Tsys table: ", tsysTable)
        return
    if (not os.path.exists(vis)):
        print("Cannot find measurement set: ", vis)
        return

    t = time.time()
    mytb = taskinit.tbtool()
    mymsmd = taskinit.msmdtool()
    mytb.open(tsysTable,nomodify=False)
    mymsmd.open(vis)
    print("tsysNormalize: initial setup took %.3f seconds" % (time.time() - t))

    # For convenience squish the useful columns into unique lists
    t = time.time()
    tsysSpws     = pb.unique(mytb.getcol("SPECTRAL_WINDOW_ID"))
    tsysScans    = pb.unique(mytb.getcol("SCAN_NUMBER"))
    tsysTimes    = pb.unique(mytb.getcol("TIME"))
    tsysFields   = pb.unique(mytb.getcol("FIELD_ID"))
    tsysAntennas = pb.unique(mytb.getcol("ANTENNA1"))
    if type(scaleSpws) == str:
        scaleSpws = [int(i) for i in scaleSpws.split(',')]
    if len(scaleSpws) < len(tsysSpws):
        scaleSpws = []
        for tsysSpw in tsysSpws:
            scaleSpws.append(scienceSpwForTsysSpw(mymsmd, tsysSpw))
        print("Identified autocorrelation spws to use: ", scaleSpws)
    print("Tsys Spws (%d):"%len(tsysSpws), tsysSpws)
    print("Tsys Scans (%d):"%len(tsysScans), tsysScans)
    print("Tsys Times (%d):"%len(tsysTimes), tsysTimes)
    print("Tsys Fields (%d):"%len(tsysFields), tsysFields)
    print("Tsys Antennas (%d):"%len(tsysAntennas), tsysAntennas)

    # Gather the power levels to use in the normalization process
    refPowers = {}
    refScans = {}
    for f in tsysFields:
        scanFieldsTab = mytb.query('FIELD_ID==%d'%f)
        fieldTsysScans = pb.unique(scanFieldsTab.getcol("SCAN_NUMBER"))
        scanFieldsTab.close()
        fieldAllScans = mymsmd.scansforfield(f)
        fieldNonTsysScans = [x for x in fieldAllScans if x not in fieldTsysScans]
        fieldName = mymsmd.namesforfields(f)[0]
        if (len(fieldNonTsysScans) < 1):
            # Then there is no non-tsys scan for this field, e.g. which can happen in a mosaic where the Tsys scan has a different field ID,
            # but in this case the field name will have other scans with different field IDs, so revert to using field names.  Using field
            # names might work from the outset, but I have not tried it.
            fieldAllScans = mymsmd.scansforfield(fieldName)
            fieldNonTsysScans = [x for x in fieldAllScans if x not in fieldTsysScans]
            if (len(fieldNonTsysScans) < 1):
                print("****** This field (id=%d, name=%s) appears to have no non-Tsys-like-scans, and thus cannot be normalized." % (f,fieldName))
                return -1
        scienceLikeScans = []
        for s in fieldNonTsysScans:
            intents = mymsmd.intentsforscan(s)
            good = True
            for i in intents:
                for b in badIntents:
                    if i.startswith(b):
                        good = False
                        break
            if good: scienceLikeScans.append(s)
        powerRefScans = []
        for s in fieldTsysScans:
            minDist = 9999999
            refScan = -1
            for r in scienceLikeScans:
                dist = abs(r - s)
                if dist < minDist:
                    minDist = dist
                    refScan = r
            powerRefScans.append(refScan)
        print("Field %d (%s) Tsys scans:"%(f,fieldname), fieldTsysScans, ", All scans:", fieldAllScans, ", Non-Tsys scans:", fieldNonTsysScans, ", Non-Tsys science-like scans:", scienceLikeScans)
        for i in range(len(fieldTsysScans)):
            print("        Tsys scan %3d power reference scan: %3d" % (fieldTsysScans[i], powerRefScans[i]))
            refScans[fieldTsysScans[i]] = powerRefScans[i]
        if verbose:
            print("populating powers corresponding to each Tsys scan on field %d..." % (f))
        for i in range(len(fieldTsysScans)):
            refPowers[fieldTsysScans[i]] = []
            for spw in scaleSpws:
                if verbose:
                    print("calling getPower(vis, %d, %d, 10.0, %s)"%(powerRefScans[i], spw, str(powerRefScans[i] < fieldTsysScans[i])))
                p = getPower(vis, powerRefScans[i], spw, 10.0, powerRefScans[i] < fieldTsysScans[i], verbose=verbose)
                refPowers[fieldTsysScans[i]].append(p)
            if verbose:
                print("powers to use for Tsys scan %d:"%fieldTsysScans[i], refPowers[fieldTsysScans[i]])
    if verbose: print(refPowers)

    print("tsysNormalize: summarising Tsys table took %.3f seconds" % (time.time() - t))
    t = time.time()

    # Now copy the original Tsys caltable and update all the values in the new one.
    if (newTsysTable == ''):
        newTsysTable = tsysTable + '_normalized'
    if (os.path.exists(newTsysTable)):
        shutil.rmtree(newTsysTable)
    mytb.copy(newTsysTable)
    mytb.close()
    mytb.open(newTsysTable, nomodify=False)
    startRefPower = refPowers[tsysScans[0]]
    for i in range(1,len(tsysScans)):
        # need to adjust each successive Tsys
        refPower = refPowers[tsysScans[i]]
        for ispw in range(len(tsysSpws)):
            spw = tsysSpws[ispw]
            for ant in range(len(tsysAntennas)):
                tsysSubTab1 = mytb.query("SCAN_NUMBER==%d AND SPECTRAL_WINDOW_ID==%d AND ANTENNA1==%d"%(tsysScans[i],tsysSpws[ispw],ant))
                tsys1 = tsysSubTab1.getcell('FPARAM', 0)
                newTsys = tsysSubTab1.getcell('FPARAM', 0)
                for pol in range(len(tsys1)):
                    for chan in range(len(tsys1[pol])):
                        a = TsysAfterPowerChange(refPowers[tsysScans[i]][ispw][ant][pol], startRefPower[ispw][ant][pol], tsys1[pol][chan])
                        newTsys[pol][chan] = a
                    print("Scan %2d spw %2d pol %d mean %.1f --> %.1f" % (tsysScans[i], spw, pol, np.mean(tsys1[pol]), np.mean(newTsys[pol])))
                tsysSubTab1.putcell('FPARAM', 0, newTsys)
                tsysSubTab1.close()
    mymsmd.close()
    mytb.close()

def fieldsMatch(f0, f1):
    """
    Checks whether any entry in one list appears in a second list.
    """
    for f in f0:
        if f in f1: return True
    return False

def tsysTransfer(vis, scaleSpws='', tsysTable='', newTsysTable='', 
                 verbose=False, overwrite=True, printAntenna=0, printPol=0):
    """
    Generate a new Tsys table where the entries for one field are propagated to
    other fields which do not have a measured Tsys, using autocorr
    (linear!) or SQLD data to determine the change in Tsys.
    Input:
     vis          the MS
     scaleSpws    the autocorr or SQLD SpWs to use for scaling (integer list or
          comma-delimited string, default is the channel-averaged science spws)
     tsysTable:   if blank, then try vis+'.tsys'
     newTsysTable:   if blank, then try vis+'.newtsys'
     printAntenna: print the before/after values for this antenna ID
     printPol: print the before/after values for this polarization (0 or 1)
    Returns: nothing
    """
    # intents likely to imply different attenuations or tuning to science-like
    # scans that we are applying Tsys to.
    badIntents = ['CALIBRATE_POINTING', 'CALIBRATE_FOCUS', 
                  'CALIBRATE_SIDEBAND_RATIO', 'CALIBRATE_ATMOSPHERE']
    if type(scaleSpws) == str:
        if (len(scaleSpws) > 0):
            scaleSpws = [int(i) for i in scaleSpws.split(',')]
    if (tsysTable == ''):
        tsysTable = vis + '.tsys'
        if not os.path.exists(tsysTable):
            tsysTables = glob.glob(os.path.join(vis,'*tsyscal.tbl'))
            if len(tsysTables) < 1:
                print("Could not find any tsys tables.")
                return
            tsysTable = tsysTables[0]
    if not os.path.exists(tsysTable):
        print("Could not find tsys table: %s" % (tsysTable))
        return
    if (newTsysTable == ''):
        newTsysTable = vis + '.newtsys'
    if overwrite and os.path.exists(newTsysTable):
        print("Removing pre-existing newTsysTable: ", newTsysTable)
        rmtables(newTsysTable) 
        if os.path.exists(newTsysTable):
            shutil.rmtree(newTsysTable)
    if (not os.path.exists(tsysTable)):
        print("Cannot find Tsys table: ", tsysTable)
        return
    if (not os.path.exists(vis)):
        print("Cannot find measurement set: ", vis)
        return

    t = time.time()
    mytb = taskinit.tbtool()
    mymsmd = taskinit.msmdtool()
    mytb.open(tsysTable,nomodify=False)
    mymsmd.open(vis)
    print("tsysTransfer: initial setup took %.3f seconds" % (time.time() - t))

    # For convenience squish the useful columns into unique lists
    t = time.time()
    tsysSpws     = pb.unique(mytb.getcol("SPECTRAL_WINDOW_ID"))
    tsysBasebands = getBasebands(mymsmd, tsysSpws)
    tsysScans    = pb.unique(mytb.getcol("SCAN_NUMBER"))
    tsysTimes    = pb.unique(mytb.getcol("TIME"))
    tsysFields   = pb.unique(mytb.getcol("FIELD_ID"))
    tsysAntennas = pb.unique(mytb.getcol("ANTENNA1"))
    finalScan = np.max(mymsmd.scannumbers())
    print("Tsys SpWs (%d):"%len(tsysSpws), tsysSpws)
    print("Tsys Basebands (%d):"%len(tsysSpws), tsysBasebands)
    print("Tsys Scans (%d):"%len(tsysScans), tsysScans)
    print("Tsys Times (%d):"%len(tsysTimes), tsysTimes)
    print("Tsys Fields (%d):"%len(tsysFields), tsysFields)
    print("Tsys Antennas (%d):"%len(tsysAntennas), tsysAntennas)
    if (len(scaleSpws) == 0):
        # number of scaleSpws should not exceed number of Tsys spws
        scaleSpws = np.unique(getChannelAveragedScienceSpws(vis,mymsmd=mymsmd))
        scaleBasebands = getBasebands(mymsmd, scaleSpws)
        if scaleBasebands != tsysBasebands:
            print("re-ordering scaleSpws to match Tsys basebands")
            newScaleSpws = []
            for baseband in tsysBasebands:
                newScaleSpws.append(scaleSpws[scaleBasebands.index(baseband)])
            scaleSpws = newScaleSpws
            scaleBasebands = tsysBasebands[:]
        print("Getting power from spws: ", scaleSpws)

    tsysScanTimes = {}
    for s in tsysScans:
        st = mytb.query('SCAN_NUMBER==%d'%s)
        ts = st.getcol("TIME")
        st.close()
        tsysScanTimes[s] = sum(ts) / float(len(ts))
        if verbose:
            print("Tsys scan %d assumed time: %.4f" % (s, tsysScanTimes[s]))

    refPowers = {}
    refScans = {}
    tsysScansOnField = {}
    for f in tsysFields:
        scanFieldsTab = mytb.query('FIELD_ID==%d'%f)
        fieldTsysScans = pb.unique(scanFieldsTab.getcol("SCAN_NUMBER"))
        scanFieldsTab.close()
        tsysScansOnField[f] = fieldTsysScans
        fieldAllScans = mymsmd.scansforfield(f)
        fieldName = mymsmd.namesforfields(f)[0]
        fieldNonTsysScans = [x for x in fieldAllScans if x not in fieldTsysScans]
        if (len(fieldNonTsysScans) < 1):
            # Then there is no non-tsys scan for this field, e.g. which can happen in a mosaic where the Tsys scan has a different field ID,
            # but in this case the field name will have other scans with different field IDs, so revert to using field names.  Using field
            # names might work from the outset, but I have not tried it.
            fieldAllScans = mymsmd.scansforfield(fieldName)
            fieldNonTsysScans = [x for x in fieldAllScans if x not in fieldTsysScans]
            if (len(fieldNonTsysScans) < 1):
                print("****** This field (id=%d, name=%s) appears to have no non-Tsys-like-scans, and thus cannot be normalized." % (f,fieldName))
                return -1
        print("Field %d (%s) Tsys scans:"%(f,fieldName), fieldTsysScans, ", All scans:", fieldAllScans, ", Non-Tsys scans:", fieldNonTsysScans)
        scienceLikeScans = []
        for s in fieldNonTsysScans:
            intents = mymsmd.intentsforscan(s)
            good = True
            for i in intents:
                for b in badIntents:
                    if i.startswith(b):
                        good = False
                        break
            if good: scienceLikeScans.append(s)
        powerRefScans = []
        for s in fieldTsysScans:
            minDist = 9999999
            refScan = -1
            for r in scienceLikeScans:
                dist = abs(r - s)
                if dist < minDist:
                    minDist = dist
                    refScan = r
            powerRefScans.append(refScan)
        if verbose:
            print("Field %d (%s) Tsys scans:"%(f,fieldName), fieldTsysScans, ", All scans:", fieldAllScans, ", Non-Tsys scans:", fieldNonTsysScans, ", Non-Tsys science-like scans:", scienceLikeScans)
        for i in range(len(fieldTsysScans)):
            if verbose:
                print("        Tsys scan %3d power reference scan: %3d" % (fieldTsysScans[i], powerRefScans[i]))
            refScans[fieldTsysScans[i]] = powerRefScans[i]
        if verbose:
            print("populating powers corresponding to each Tsys scan on field %d..." % (f))
        for i in range(len(fieldTsysScans)):
            refPowers[fieldTsysScans[i]] = []
            for spw in scaleSpws:
                if verbose:
                    print("powerRefScans: ", powerRefScans)
                    print("calling getPower(vis, %d, %d, 10.0, %s)"%(powerRefScans[i], spw, str(powerRefScans[i] < fieldTsysScans[i])))
                p = getPower(vis, powerRefScans[i], spw, 10.0, powerRefScans[i] < fieldTsysScans[i], verbose=verbose)
                refPowers[fieldTsysScans[i]].append(p)
            #print "powers to use for Tsys scan %d:"%fieldTsysScans[i], refPowers[fieldTsysScans[i]]
#    print refPowers

    print("tsysTransfer: summarising Tsys table took %.3f seconds" % (time.time() - t))

    t = time.time()
    mytb.copy(newTsysTable)
    mytb.close()
    # re-open original table as read-only
    mytb.open(tsysTable)
    mytbNew = taskinit.tbtool()
    mytbNew.open(newTsysTable, nomodify=False)
    print("tsysTransfer: Copying Tsys table from '%s' to '%s' took %.3f seconds" % (tsysTable, newTsysTable, time.time() - t))
    anyProcessingNeeded = False

    # Loop over each Tsys scan
    for i in range(len(tsysScans)-1):
        tsysTime0 = tsysScanTimes[tsysScans[i]]
        tsysTime1 = tsysScanTimes[tsysScans[i+1]]
        tsysTimeGap = tsysTime1 - tsysTime0
        tsysFields0 = mymsmd.fieldsforscan(tsysScans[i])  # current Tsys scan
        tsysFields1 = mymsmd.fieldsforscan(tsysScans[i+1]) # next Tsys scan
        # loop over all scans between the current Tsys scan and the next one
        startScan = tsysScans[i]+1
        stopScan = tsysScans[i+1]
#        if finalScan > stopScan and i==len(tsysScans)-1:
#            print "There are more scans after the final Tsys scan, extending the range of scans accordingly."
#            stopScan = finalScan
        for scan in range(startScan, stopScan):
            if 'CALIBRATE_POINTING#ON_SOURCE' in mymsmd.intentsforscan(scan): continue
            processingNeeded = False
            fields = mymsmd.fieldsforscan(scan)
            times = mymsmd.timesforscan(scan)
            startTime = times[0]
            endTime = times[-1]
            print("Processing scan %d with fields %s, between Tsys scan %d (fields %s) and %d (fields %s)" % (scan, str(fields[0]), tsysScans[i], str(tsysFields0[0]), tsysScans[i+1], str(tsysFields1[0])))
            print("    Scan %d starts %.3f sec after preceding Tsys, and ends %.3f sec before next Tsys" % (scan, startTime - tsysTime0, tsysTime1 - endTime))
            # There are a few possible cases to deal with:
            # 1) this was a power reference scan for a Tsys scan, in which case only produce one extra Tsys, at the opposite end of the scan, or none if there are Tsys scans for the same field at both ends
            fieldMatchesPriorTsysField = fieldsMatch(fields, tsysFields0)
            fieldMatchesNextTsysField = fieldsMatch(fields, tsysFields1)
            priorScanIsTsys = scan == tsysScans[i]+1
            nextScanIsTsys = scan == tsysScans[i+1]-1
            bracketingTsysFieldsMatch = fieldsMatch(tsysFields0, tsysFields1)
            scanIsNotRefScan = scan != refScans[tsysScans[i]] and scan != refScans[tsysScans[i+1]]
            if fieldMatchesPriorTsysField and fieldMatchesNextTsysField and priorScanIsTsys and nextScanIsTsys:
                print("    Nothing needed for scan %d as bracketed immediately by two Tsys scans of same field" % scan)
            # The most common case for wanting to do the transfer: science field bracketed by phase cal, or phase cal without Tsys immediately before/after
            elif bracketingTsysFieldsMatch and (not fieldMatchesPriorTsysField or scanIsNotRefScan):
                # The two Tsys scans that bracket this scan are taken on the same field; 
                # and either this scan is not on the field of the prior Tsys scan, or
                # this scan is not a reference scan
                processingNeeded = True
                priorScanToUse = tsysScans[i]
                nextScanToUse = tsysScans[i+1]
            elif (not bracketingTsysFieldsMatch and fields[0] in tsysScansOnField.keys()):
                candidateScans = np.array(tsysScansOnField[fields[0]])
                if (scan < candidateScans[0] or scan > candidateScans[-1]):
                    print("    The bracketing Tsys fields do not match, and there are not two scans to interpolate between.")
                else:
                    processingNeeded = True
                    priorScanToUse = np.max(candidateScans[np.where(candidateScans < scan)])
                    nextScanToUse = np.min(candidateScans[np.where(candidateScans > scan)])
                    print("    The bracketing Tsys fields do not match, but there are two scans to interpolate between: %d and %d." % (priorScanToUse,nextScanToUse))
            elif (not bracketingTsysFieldsMatch):
                # This section added by Todd for initial phase calibrator scans when Tsys taken on science target only.
                # Not sure what to do yet, though.
                print("    The bracketing Tsys fields do not match, and Tsys was never taken on this field. No processing will be done.")
                if False:
                    processingNeeded = True
                    if i+1 < len(tsysScans):
                        print("    Extrapolating from subsequent Tsys scan: %d" % (tsysScans[i+1]))
                        priorScanToUse = tsysScans[i+1]
                        nextScanToUse = tsysScans[i+1]
                    else:
                        print("    Extrapolating from prior Tsys scan: %d" % (tsysScans[i+1]))
                        priorScanToUse = tsysScans[i]
                        nextScanToUse = tsysScans[i]
            else:
                print("    This scan arrangement is unexpected.  No processing will be done.")
                print("      bracketingTsysFieldsMatch = %s" % bracketingTsysFieldsMatch)
                print("      fieldMatchesPriorTsysField = %s" % fieldMatchesPriorTsysField)
                print("      fieldMatchesNextTsysField = %s" % fieldMatchesNextTsysField)
                print("      priorScanIsTsys = %s" % priorScanIsTsys)
                print("      nextScanIsTsys = %s" % nextScanIsTsys)
                print("      scanIsNotRefScan = %s" % scanIsNotRefScan)
                print("      %s in tsysScansOnField(%s) = %s" % (fields[0], tsysScansOnField.keys(), fields[0] in tsysScansOnField.keys()))            
            if processingNeeded:
                anyProcessingNeeded = True
                print("    For scan %d will generate two Tsys entries for beginning and end of scan, interpolating reference from scans %d and %d" % (scan, priorScanToUse, nextScanToUse))
                for ispw in range(len(scaleSpws)):
                    spw = scaleSpws[ispw]
                    startPower = getPower(vis, scan, spw, 10.0, False, verbose=verbose)
                    endPower   = getPower(vis, scan, spw, 10.0, True, verbose=verbose)
                    for ant in range(len(tsysAntennas)):
                        tsysSubTab0 = mytb.query("SCAN_NUMBER==%d AND SPECTRAL_WINDOW_ID==%d AND ANTENNA1==%d"%(priorScanToUse, tsysSpws[ispw], ant))
                        tsysSubTab1 = mytb.query("SCAN_NUMBER==%d AND SPECTRAL_WINDOW_ID==%d AND ANTENNA1==%d"%(nextScanToUse, tsysSpws[ispw], ant))
                        # sanity check for duplicate entries
                        if tsysSubTab0.nrows() != 1 or tsysSubTab1.nrows() != 1:
                            print("WARNING!!! not one result row for (scan,ant,spw) query in Tsys table. Scan %d: %d rows, Scan %d: %d rows." % (priorScanToUse, tsysSubTab0.nrows(), nextScanToUse, tsysSubTab1.nrows()))
                        tsys0 = tsysSubTab0.getcell('FPARAM', 0)
                        tsys1 = tsysSubTab1.getcell('FPARAM', 0)
                        tsysSubTab1.close()
                        startTsys = copy.copy(tsys0)
                        endTsys = copy.copy(tsys0)  # just a placeholder, new values will be filled in below
                        startRefPower = refPowers[priorScanToUse]
                        endRefPower = refPowers[nextScanToUse]
                        tsysTime0 = tsysScanTimes[priorScanToUse]
                        tsysTime1 = tsysScanTimes[nextScanToUse]
                        tsysTimeGap = tsysTime1 - tsysTime0
                        for pol in range(len(tsys0)):
                            for chan in range(len(tsys0[pol])):
                                startTsys0 = TsysAfterPowerChange(startRefPower[ispw][ant][pol], startPower[ant][0], tsys0[pol][chan])
                                startTsys1 = TsysAfterPowerChange(  endRefPower[ispw][ant][pol], startPower[ant][0], tsys1[pol][chan])
                                endTsys0   = TsysAfterPowerChange(startRefPower[ispw][ant][pol],   endPower[ant][0], tsys0[pol][chan])
                                endTsys1   = TsysAfterPowerChange(  endRefPower[ispw][ant][pol],   endPower[ant][0], tsys1[pol][chan])
                                if tsysTimeGap == 0:
                                    startTsys[pol][chan] = startTsys0
                                    endTsys[pol][chan] = endTsys0
                                else:
                                    startTsys[pol][chan] = ((startTime - tsysTime0) * startTsys1 + (tsysTime1 - startTime) * startTsys0) / tsysTimeGap
                                    endTsys[pol][chan]   = ((endTime - tsysTime0) * endTsys1 + (tsysTime1 - endTime) * endTsys0) / tsysTimeGap
                                if chan == len(tsys0[pol]) / 2 and ant==printAntenna and pol==printPol:
                                    print("    ispw=%d spw=%d ant=%d pol=%d chan=%d: TsysBefore: %.1f K, TsysScanStart: %.1f K (interp %.1f,%.1f), TsysScanEnd: %.1f K (interp %.1f,%.1f), TsysAfter: %.1f K" % (ispw, spw, ant, pol, chan, tsys0[pol][chan], startTsys[pol][chan], startTsys0, startTsys1, endTsys[pol][chan], endTsys0, endTsys1, tsys1[pol][chan]))
                        for f in fields:
                            nr = mytbNew.nrows()
                            tsysSubTab0.copyrows(newTsysTable, nrow=1)
                            if verbose: 
                                print("setting tsys at row %d" % nr)
                            mytbNew.putcell('FPARAM', nr, startTsys)
                            mytbNew.putcell('TIME', nr, startTime)
                            mytbNew.putcell('FIELD_ID', nr, f)
                            mytbNew.putcell('SCAN_NUMBER', nr, scan)
                            nr = mytbNew.nrows()
                            tsysSubTab0.copyrows(newTsysTable, nrow=1)
                            if verbose:
                                print("setting tsys at row %d" % nr)
                            mytbNew.putcell('FPARAM', nr, endTsys)
                            mytbNew.putcell('TIME', nr, endTime)
                            mytbNew.putcell('FIELD_ID', nr, f)
                            mytbNew.putcell('SCAN_NUMBER', nr, scan)
                        tsysSubTab0.close()
                        # end loop over fields (f)
                    # end loop over antennas (ant)
                # end loop over spws (ispw)
            # end if processingNeeded
        # end loop over scans between tsysScans (scan)
        mytbNew.flush()
    # end loop over Tsys scans
    if not anyProcessingNeeded:
        print("Because no processing was needed the new Tsys table is identical to the original.")
    # TODO: These cleanups should be done also on an exception too
    print("Closing tables...")
    mytbNew.unlock()
    mytbNew.close()
    mytbNew.done()
    mymsmd.close()
    mytb.close()
    mytb.done()

