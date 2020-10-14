# compare focus results of two asdms
# First version imported by R. Lucas.  All subsequent edits by T. Hunter
from __future__ import print_function  # prevents adding old-style print statements
import pylab as pl
from TelCal.AsdmReader import AsdmPlotter
from asdm import *
import numpy.ma as ma
import numpy as np
from scipy.optimize import leastsq
from math import *

# note I have changed asdmReader to exclude the last point of WVR in all subscans.
# TODO do this only if it overflows the time for the correlator data.


class WvrCorrection(AsdmPlotter):

    def getWvrData(self, doFlag = True):
        '''
        Make a dictionary with the WVR data
        '''
        
        #self.setX('time')
        self.setAntennas('')
        self.setX('time') ##  elevation')
        ## self.setX('time elevation')
        self.setY('wvr1 wvr2 wvr3 wvr4')
        self.configure()
        self.fillData()
        wvrs = ['WVR1','WVR2','WVR3','WVR4']
        wvrData = {}
        timeData = {}
        elevationData = {}
        for i in range(self.getNumReaderTables()):
            t = self.getReaderTable(i)
            f = np.array(t.getFlag())
            if not ma.all(f):
                key = t.getAntenna()
                # if key=='DV17': print f
                x = np.array(t.getX())
                y = np.array(t.getY())
                if t.getXQuantity() == t.ELEVATION:
                    elevationData[key] = ma.masked_array(x, f, copy=True)
                elif t.getXQuantity() == t.TIME:
                    channel = t.getYQuantity()-t.WVR_1
                    if key not in list(wvrData.keys()):
                        wvrData[key] = ma.zeros((4,len(y)))
                    wvrData[key][channel,:] = ma.masked_array(y, f, copy=True)
                    timeData[key] = ma.masked_array(x,f, copy=True)
        self.wvrData = wvrData
        self.timeData = timeData
        self.elevationData = elevationData
        self.airMassData = {}
        for key in list(self.elevationData.keys()):
            self.airMassData[key] = 1./np.sin(self.elevationData[key])

    def getAsdm(self):
        self.asdm = ASDM()
        po = ASDMParseOptions()
        # The calibration result is an ASDM without ExecBlock table. It must
        # be defined "asALMA" to be readable by the method setFromFile().
        po.asALMA()
        po.loadTablesOnDemand(True)
        self.asdm.setFromFile(self.getInputFileName(), po)
        

    def getFrequency(self):
        self.getAsdm()
        # find a full resolution row...
        frequency = {}
        for row in self.asdm.configDescriptionTable().get():
            if str(row.spectralType())  ==  'CHANNEL_AVERAGE':
                ddIds = row.dataDescriptionId()
                break
        for ddId in ddIds:
            spwId = self.asdm.dataDescriptionTable().getRowByKey(ddId).spectralWindowId()
            row = self.asdm.spectralWindowTable().getRowByKey(spwId)
            frequency[row.basebandName()] = row.chanFreqArray()[0].get()
        self.frequency = frequency

    def getPhaseData(self):
        phaseData = {}
        self.setX('time')
        self.setY('phaseant')
        self.configure()
        self.fillData()
        timeData = {}
        for i in range(self.getNumReaderTables()):
            t = self.getReaderTable(i)
            f = np.array(t.getFlag())
            if not ma.all(f):
                key = (t.getAntenna(), t.getBasebandNames()[0], t.getPolarizations()[0])
                x = np.array(t.getX())
                y = np.array(t.getY())
                phaseData[key] = ma.masked_array(y, f, copy=True)
                timeData[key] = ma.masked_array(x,f, copy=True)
        self.phaseData = phaseData
        self.timeData = timeData
        

        

    def getWvrCoefficients(self, calDataSet=None):
        if calDataSet is None:
            self.cdm = self.asdm
            # calDataSet=self.getInputFileName()
        else:
            self.cdm = ASDM()
            po = ASDMParseOptions()
            # The calibration result is an ASDM without ExecBlock table. It must
            # be defined "asALMA" to be readable by the method setFromFile().
            po.asALMA()
            po.loadTablesOnDemand(True)
            self.cdm.setFromFile(calDataSet, po)

        cdRows = self.cdm.calDataTable().get()
        self.wvrPathCoeff = {}
        self.wvrRefTemp = {}
        self.water = {}
        self.crId = {}
        for row in cdRows:
            scan = row.scanSet()[0]
            cid = row.calDataId()
            for row in self.cdm.calWVRTable().get():
                if row.calDataId() == cid:
                    antenna = row.antennaName()
                    key = (scan, antenna)
                    # simplified version ignoring antenna-based interpolation...
                    # if antenna == 'DA50':
                    #    print row
                    self.wvrPathCoeff[key] = np.array(row.pathCoeff())[0,:,0]
                    # print row.refTemp()
                    self.wvrRefTemp[key] = np.array([row.refTemp()[0][i].get() for i in range(4)])
                    self.water[key] = row.water().get()
                    self.crId[scan] = row.calReductionId()
       
        
        

    def getWvrCorrection(self, calDataSet=None, scan=None, removeAnts=[]):
        if calDataSet is None:
            self.cdm = self.asdm
            # calDataSet=self.getInputFileName()
        else:
            self.cdm = ASDM()
            po = ASDMParseOptions()
            # The calibration result is an ASDM without ExecBlock table. It must
            # be defined "asALMA" to be readable by the method setFromFile().
            po.asALMA()
            po.loadTablesOnDemand(True)
            print('calDataSet: ', calDataSet)
            self.cdm.setFromFile(calDataSet, po)
            
        cdRows = self.cdm.calDataTable().get()
        cid = None
        for row in cdRows:
            if scan is None:
                cid = row.calDataId()
                break
            elif scan in row.scanSet():
                cid = row.calDataId()
                break
            
        if cid is None:
            return
        wvrPathCoeff = {}
        wvrRefTemp = {}
        water = {}
        for row in self.cdm.calWVRTable().get():
            if row.calDataId() == cid:
                antenna = row.antennaName()
                # simplified version ignoring antenna-based interpolation...
                wvrPathCoeff[antenna] = np.array(row.pathCoeff())[0,:,0]
                wvrRefTemp[antenna] = np.array([row.refTemp()[0][i].get() for i in range(4)])
                water[antenna] = row.water().get()
                cdId = row.calReductionId()

        self.water = water
        self.wvrPathCoeff = wvrPathCoeff
        self.wvrRefTemp = wvrRefTemp
        # get the wvrCouplings in the CalReduction table
        print('cdId ',cdId)
        row = self.cdm.calReductionTable().getRowByKey(cdId)
        paramSet = row.paramSet()
        self.skyCoupling = {}
        for par in paramSet:
            if par.split('=')[0].split('[')[0] == 'skyCoupling':
                ant = par.split('=')[0].split('[')[1].strip(']')
                value = par.split('=')[1]
                self.skyCoupling[ant] = float(value)
        
        wvrCorrection = {}
        wvrChanCorrection = {}
        for antenna in list(self.wvrData.keys()):
            wvrCorrection[antenna] = 0.
            wvrChanCorrection[antenna] = []
            for i in range(4):
                wvrChanCorrection[antenna].append((self.wvrData[antenna][i]-self.wvrRefTemp[antenna][i])*self.wvrPathCoeff[antenna][i])
                wvrCorrection[antenna] += (self.wvrData[antenna][i]-self.wvrRefTemp[antenna][i])*self.wvrPathCoeff[antenna][i]
        self.wvrCorrection = wvrCorrection
        self.wvrChanCorrection = wvrChanCorrection
        self.wvrAntennas = []
        for ant in list(self.wvrCorrection.keys()):
            if ant[0:2] != 'CM' and ant not in removeAnts:
                self.wvrAntennas.append(ant)
        self.wvrAntennas.sort()
        self.averageWvrCorrection = self.wvrCorrection[list(self.wvrCorrection.keys())[0]] 
        self.averageWvrCorrection = 0
        self.averageWvrChanCorrection = [0,0,0,0]
        numAveraged = 0
        for ant in self.wvrAntennas:
            self.averageWvrCorrection +=  self.wvrCorrection[ant]
            for i in range(4):
                self.averageWvrChanCorrection[i] += self.wvrChanCorrection[ant][i]
            numAveraged += 1
        self.averageWvrCorrection /= numAveraged
        for i in range(4):
            self.averageWvrChanCorrection[i] /= numAveraged


    def getCoeffFromWater(self):
        self.tweak = {}
        self.tweakError = {}
        averageWater = 0
        numAveraged = 0
        for ant in self.wvrAntennas:
            averageWater += self.water[ant]
            numAveraged += 1
        averageWater /= numAveraged
        
        for ant in self.wvrAntennas:
            self.tweak[ant] = averageWater/self.water[ant]
            self.tweakError[ant] = 0

        self.averageTweakedWvrCorrection = 0
        numAveraged = 0
        for ant in self.wvrAntennas:
            self.averageTweakedWvrCorrection +=  self.wvrCorrection[ant]*self.tweak[ant]
            numAveraged += 1
        self.averageTweakedWvrCorrection /= numAveraged
        self.raw_rms = {}
        self.tweaked_rms = {}
        for ant in self.wvrAntennas:
            yy = (self.wvrCorrection[ant]-self.averageWvrCorrection)
            self.raw_rms[ant] = sqrt(yy.var())
            yy = (self.wvrCorrection[ant]*self.tweak[ant]-self.averageTweakedWvrCorrection)
            self.tweaked_rms[ant] = sqrt(yy.var())
            
    def tweakCoeffs(self, airmassMax=1.5, np=20):

        self.tweak = {}
        self.tweakError = {}
        self.averageTweakedWvrCorrection =  self.averageWvrCorrection
        self.averageTweakedWvrCorrection = 0

        def errorfunc(p, y, x):
            a,b = p
            return y - x/a -b
        

    
        for ant in self.wvrAntennas:
            yy = self.wvrCorrection[ant]
            xx = self.averageWvrCorrection
            aa = self.airMassData[ant]
            tt = self.timeData[ant]
            np = 20 #
            dx = xx[np:]-xx[0:-np]
            dy = yy[np:]-yy[0:-np]
            da = aa[np:]-aa[0:-np]
            dt = tt[np:]-tt[0:-np]
            # mask
            mask = abs(da) < 0.01 # typ 0.1 degree elev
            mask = ma.logical_or(mask, (abs(aa[np:]) > airmassMax)) 
            mask = ma.logical_or(mask, (abs(aa[0:-np]) > airmassMax)) 
            da = ma.masked_array(da, mask, copy=True)
            dx = ma.masked_array(dx, mask, copy=True)
            dy = ma.masked_array(dy, mask, copy=True)
            output = leastsq(errorfunc,  (1., 0.), (dy,dx), full_output=True)
            # print ant, output
            p = output[0]
            res = errorfunc(p, dy, dx)
            try:
                print(ant, p[0], sqrt(res.var()*output[1][0,0]))
                self.tweak[ant] = p[0]
                self.tweakError[ant] =  sqrt(res.var()*output[1][0,0])
            except:
                print('error fitting %s data' % ant)

        numAveraged = 0
        for ant in self.wvrAntennas:
            ## self.tweak[ant] = sqrt(self.averageWvrCorrection.var()/self.wvrCorrection[ant].var())
            if ant in list(self.tweak.keys()):
                self.averageTweakedWvrCorrection +=  self.wvrCorrection[ant]*self.tweak[ant]
                numAveraged += 1
        self.averageTweakedWvrCorrection /= numAveraged
        self.raw_rms = {}
        self.tweaked_rms = {}
        for ant in self.wvrAntennas:
            if ant in list(self.tweak.keys()):
                yy = (self.wvrCorrection[ant]-self.averageWvrCorrection)
                self.raw_rms[ant] = sqrt(yy.var())
                yy = (self.wvrCorrection[ant]*self.tweak[ant]-self.averageTweakedWvrCorrection)
                self.tweaked_rms[ant] = sqrt(yy.var())
        
            
    def process(self, calDataSet=None, removeAnts=[]):
        self.setSubscans('1')  # ignore the loads...
        self.setAntennas('')
        self.setBasebands('BB_1')
        self.setPolarizations('XX')
        self.setPhaseCorrection('AP_UNCORRECTED')
        self.setScans('')
        self.setAntennas('')
        self.getWvrData()
        self.getAsdm()
        print("... getWvrCorrection")
        self.getWvrCorrection(calDataSet=calDataSet, removeAnts=removeAnts)
        print("... tweakCoeffs")
        self.tweakCoeffs()
        #self.getCoeffFromWater()
        print("... get Efficiencies")
        self.getEfficiencies()

        
    def summaryPlot(self, calDataSet=None, refant=None):
        pl.clf()
        ax = pl.subplot(411)
        ax.set_xticklabels([])
        numAntennas = len(self.wvrAntennas)
        pl.xlim(-0.5, numAntennas-0.5)
        for ant in self.wvrAntennas:
            pl.plot(self.wvrAntennas.index(ant), self.water[ant], 'ob')
            # pl.plot(self.wvrAntennas.index(ant), self.wvrCoupling[ant]*self.tweak[ant], 'og')
        pl.ylabel('Water (mm)')

        ax = pl.subplot(412)
        ax.set_xticklabels([])
        numAntennas = len(self.wvrAntennas)
        pl.xlim(-0.5, numAntennas-0.5)
        for ant in self.wvrAntennas:
            pl.plot(self.wvrAntennas.index(ant), self.skyCoupling[ant], 'ob')
            if ant in list(self.tweak.keys()):
                pl.plot(self.wvrAntennas.index(ant), self.skyCoupling[ant]*self.tweak[ant], 'og')
        pl.ylabel('Sky Coupling')

        
        ax = pl.subplot(413)
        ax.set_xticklabels([])
        numAntennas = len(self.wvrAntennas)
        pl.xlim(-0.5, numAntennas-0.5)
        for ant in self.wvrAntennas:
            if ant in list(self.tweak.keys()):
                pl.plot(self.wvrAntennas.index(ant), self.tweak[ant], 'ob')
                pl.errorbar(self.wvrAntennas.index(ant), self.tweak[ant], yerr=self.tweakError[ant], color='b')
        # pl.ylim(0.95, 1.05)
        pl.ylabel('Scale Factor')
        
        ax = pl.subplot(414)
        ax.set_xticklabels([])
        numAntennas = len(self.wvrAntennas)
        pl.xlim(-0.5, numAntennas-0.5)
        for ant in self.wvrAntennas:
            if ant in list(self.tweak.keys()):
                pl.plot(self.wvrAntennas.index(ant), self.raw_rms[ant], 'or')
                pl.plot(self.wvrAntennas.index(ant), self.tweaked_rms[ant], 'og')
        # pl.ylim(0.0, 0.002)
        pl.xticks(np.arange(numAntennas), self.wvrAntennas,
                     horizontalalignment='left',
                     size='small',
                     rotation='-45')
        
        pl.ylabel('rms')
        pl.savefig('%s-corrPlot.png'%self.getInputFileName())

    def antennaPlot(self, ant):
        pl.clf()
        ax = pl.subplot(311)
        y1 = self.wvrCorrection[ant]
        y2 = self.averageWvrCorrection
        # y3 = self.wvrCorrection[ant]-self.averageWvrCorrection
        raw_rms = sqrt(y1.var())
        pl.plot(self.timeData[ant], y1,'.r')
        pl.plot(self.timeData[ant], y2,'.b')
        # pl.plot(self.timeData[ant], y3,'.g')
        pl.xlabel('Time')
        pl.ylabel('Wvr Correction')
        pl.xlim(-100.)
        pl.title('%s - %s'%(self.getInputFileName(), ant))
                    
        ax = pl.subplot(312)
        y1 = self.wvrCorrection[ant]-self.averageWvrCorrection
        y3 = y1[1:]-y1[0:-1]
        tt = self.timeData[ant]
        dt = tt[1:]-tt[0:-1]
        ym = ma.masked_array(y3, mask=abs(dt<2.))
        #y2 = self.averageTweakedWvrCorrection
        #y3 = self.wvrCorrection[ant]*self.tweak[ant]-self.averageTweakedWvrCorrection
        tweaked_rms = sqrt(y1.var())
        pl.plot(self.timeData[ant], y1,',r')
        pl.plot(self.timeData[ant][1:], ym,'.m')
        #pl.plot(self.timeData[ant], y2,'.b')
        #pl.plot(self.timeData[ant], y3,'.g')
        pl.xlabel('Time')
        pl.ylabel('Wvr C-Aver')
        pl.xlim(-100.)
        pl.title('%s - %s raw'%(self.getInputFileName(), ant))
        
        ax = pl.subplot(313)
        #y1 = self.wvrCorrection[ant]*self.tweak[ant]
        #y2 = self.averageTweakedWvrCorrection
        y1 = self.wvrCorrection[ant]*self.tweak[ant]-self.averageTweakedWvrCorrection
        y3 = y1[1:]-y1[0:-1]
        tt = self.timeData[ant]
        dt = tt[1:]-tt[0:-1]
        ym = ma.masked_array(y3, mask=abs(dt<2.))
        tweaked_rms = sqrt(y1.var())
        #pl.plot(self.timeData[ant], y1,'.r')
        #pl.plot(self.timeData[ant], y2,'.b')
        pl.plot(self.timeData[ant], y1,',g')
        pl.plot(self.timeData[ant][1:], ym,'.m')
        pl.xlabel('Time')
        pl.ylabel('Wvr C-Av (tweaked)')
        pl.xlim(-100.)
        pl.title('%s - %s tweaked'%(self.getInputFileName(), ant))
        pl.savefig('%s-%s-antennaPlot.png'%(self.getInputFileName(), ant))

    def getEfficiencies(self):
        #wvr.skyCoupling = {}
        #for ant in self.wvrAntennas:
        #    wvr.skyCoupling[ant] = 0.975*self.tweak[ant]
        self.skyCouplingString = ''
        self.skyCouplingTweaked = {}
        for ant in self.wvrAntennas:
            try:
                self.skyCouplingTweaked[ant] = self.skyCoupling[ant]*self.tweak[ant]
                self.skyCouplingString += '%s:%6.4f '%(ant, self.skyCouplingTweaked[ant])
            except:
                print('problem with %s' %ant)


    def getInputData(self):
        '''
        Open the skyCoupling.dat file that contains the used altitudes, temperatures,
        pressures for each antenna
        '''
        sc = open('skyCoupling.dat')
        self.dist = {}
        self.alti = {}
        self.temp = {}
        self.pres = {}
        self.water = {}
        self.coef = {}
        for l in sc.readlines():
            ls = l.split()
            ant = ls[1]
            self.alti[ant] = float(ls[5])
            self.dist[ant] = float(ls[3])
            self.temp[ant] = float(ls[7])
            self.pres[ant] = float(ls[9])
            self.water[ant] =float(ls[13])
            self.coef[ant] = pl.array([float(ls[14+i]) for i in range(4)])

        antennas =  list(self.alti.keys())
        

        # altitude
        pl.clf()
        ax = pl.subplot(3,1,1)
        pl.ylabel('Altitude [m]')
        for i in range(len(antennas)):
            pl.plot(self.dist[antennas[i]], self.alti[antennas[i]],'or')
        pl.title(asdm)
        # temp
        #ax = pl.subplot(4,1,2)
        #for i in range(len(antennas)):
        #    pl.plot(self.dist[antennas[i]], self.temp[antennas[i]],'og')
        # pres
        #ax = pl.subplot(4,1,3)
        #for i in range(len(antennas)):
        #    pl.plot(self.dist[antennas[i]], self.pres[antennas[i]],'ob')
        # water
        ax = pl.subplot(3,1,2)
        pl.ylabel('Water [mm]')
        for i in range(len(antennas)):
            pl.plot(self.dist[antennas[i]], self.water[antennas[i]]*1000.,'oc')
        # coefs
        ax = pl.subplot(3,1,3)
        pl.ylabel('Coefficients [$\mu$m/K]')
        colors = ['r','g','b','c']
        for j in range(4):
            for i in range(len(antennas)):
                pl.plot(self.dist[antennas[i]], self.coef[antennas[i]][j],'.'+colors[j])

        pl.xlabel('Distance [m]')
            


    def getLoadTemps(self):

        self.ambientLoadTemp = {}
        self.hotLoadTemp = {}

        for row in self.getDataset().calDeviceTable().get():

            ant =  self.getDataset().antennaTable().getRowByKey(row.antennaId()).name()
            #  print ant
            if ant not in list(self.ambientLoadTemp.keys()):
                for i in range(len(row.calLoadNames())):
                    if str(row.calLoadNames()[i]) == 'AMBIENT_LOAD':
                        self.ambientLoadTemp[ant] = row.temperatureLoad()[i].get()
                    if str(row.calLoadNames()[i]) == 'HOT_LOAD':
                        self.hotLoadTemp[ant] = row.temperatureLoad()[i].get()


    def getLoadMeasurements(self):

        self.ambientLoadMeasurements = {}
        self.hotLoadMeasurements = {}
        for ant in list(self.wvrData.keys()):
            self.ambientLoadMeasurements[ant] = [None,None,None,None]
            self.hotLoadMeasurements[ant] = [None,None,None,None]
            
            for i in range(4):
                amb = []
                hot = []
                for v in self.wvrData[ant][i]:
                    if abs(v - self.ambientLoadTemp[ant]) < 4.:
                        amb.append(v)
                    if abs(v - self.hotLoadTemp[ant]) < 4.:
                        hot.append(v)
                self.ambientLoadMeasurements[ant][i] = pl.array(amb).mean()
                self.hotLoadMeasurements[ant][i] = pl.array(hot).mean()


## asdm='uid___A002_X901da2_X10'
## w = WvrCorrection(asdm)
## w.getWvrData()
## w.getLoadTemps()
## w.getLoadMeasurements()


## antennas = w.wvrData.keys()
## antennas.sort()
## for a in antennas:
##     if a[0:2] != 'CM':
##         ma = pl.array(w.ambientLoadMeasurements[a])
##         ta = float(w.ambientLoadTemp[a])
##         mh = pl.array(w.hotLoadMeasurements[a])
##         th = float(w.hotLoadTemp[a])
##         print "%5s Ambient measured %6.2f %6.2f %6.2f %6.2f K Load %6.2f K" %\
##               (a, ma[0], ma[1], ma[2], ma[3], ta)
##         ma -= ta
##         print "            diff     %6.2f %6.2f %6.2f %6.2f K" %\
##               (ma[0], ma[1], ma[2], ma[3])
##         ma /= ta /100.
##         print "            ratio     %6.2f %6.2f %6.2f %6.2f %%" %\
##               (ma[0], ma[1], ma[2], ma[3])
##         print "%5s Hot     measured %6.2f %6.2f %6.2f %6.2f K Load %6.2f K" %\
##               (a, mh[0], mh[1], mh[2], mh[3], th)
##         mh -= th
##         print "            diff     %6.2f %6.2f %6.2f %6.2f K" %\
##               (mh[0], mh[1], mh[2], mh[3])
##         mh /= th /100.
##         print "            diff     %6.2f %6.2f %6.2f %6.2f %%" %\
##               (mh[0], mh[1], mh[2], mh[3])
##         print ""
