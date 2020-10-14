# plot Antennas position results relative to one antenna
# First version imported by R. Lucas.  All subsequent edits by T. Hunter
#
from __future__ import print_function  # prevents adding old-style print statements
from asdm import *
import pylab as pl
from scipy.optimize import leastsq
import tmUtils as tm
import datetime
from math import *

def stationToAntenna(stationPosition, geo):
    '''
    return vector in antenna coordinates of geocentric coordinates (geo) relative to geocentric station position
    '''
    xx,yy,zz = stationPosition[0],stationPosition[1],stationPosition[2]
    lat = asin(zz / sqrt(xx ** 2 + yy ** 2 + zz ** 2))
    lon = atan2(yy, xx)
    alma = [0, 0, 0]
    alma[0] = -geo[0] * sin(lon) + geo[1] * cos(lon)
    alma[1] = -geo[0] * cos(lon) * sin(lat) - geo[1] * sin(lon) * sin(lat) + geo[2] * cos(lat)
    alma[2] = geo[0] * cos(lon) * cos(lat) + geo[1] * sin(lon) * cos(lat) + geo[2] * sin(lat)
    return alma

def antennaToStation(padPosition, position):
    """
    Computes the corrections to make to absolute geocentric coordinates
    XYZ (pad position) from an offset in local coordinates (antenna position).
    padPosition: vector [X,Y,Z] geocentric coords
    position: antenna position in local coords (e.g. from ASDM_ANTENNA table)
    Returns: the corrections to apply in geocentric frame as dX, dY, dZ
    """
    xx,yy,zz = padPosition[0],padPosition[1],padPosition[2]
    lat = asin(zz / sqrt(xx ** 2 + yy ** 2 + zz ** 2))
    lon = atan2(yy, xx)
    itrf_correction = []
    itrf_correction.append(-sin(lon)*position[0] \
                           -sin(lat)*cos(lon)*position[1] + \
                            cos(lat)*cos(lon)*position[2])
    itrf_correction.append(+cos(lon)*position[0] \
                           -sin(lat)*sin(lon)*position[1] + \
                            cos(lat)*sin(lon)*position[2])
    itrf_correction.append(+cos(lat)*position[1] + \
                            sin(lat)*position[2])
    return pl.array(itrf_correction)
    
def ellipsoidHeight(X,Y,Z):
    """
    Return the height above WGS84 earth ellipsoid of geocentric location X,Y,Z
    """
    a = 6378137.
    f = 1./298.2572236
    e2 = 2*f-f*f
    epsilon = e2/(1.-e2)
    b = a*(1.-f)
    p = sqrt(X*X+Y*Y)
    q = atan((Z*a)/(p*b))
    phi = atan( (Z+epsilon*b*sin(q)**3)/(p-e2*a*cos(q)**3) )
    nu = a / sqrt( 1-e2*sin(phi)**2)
    h = p/cos(phi) - nu
    return h

class AntPosResult(ASDM):
    '''
    class to analyse the results on an antenna position measurement.
    '''

    def getData(self, caldm, scanList=None):
        """
        """
        
        # self._cdm = ASDM()
        self._caldm = caldm
        # try:
        po = ASDMParseOptions()
        # The calibration result is an ASDM without ExecBlock table. It must
        # be defined "asALMA" to be readable by the method setFromFile().
        po.asALMA()
        po.loadTablesOnDemand(True)
        print('setFromFile '+caldm)
        self.setFromFile(caldm, po)
        self._calData = self.calDataTable().get()
        if scanList==None:
            self.scanList = self.getScanList()
        else:
            self.scanList = scanList
        # open as well the original asdm, for weather table...
        asdm = caldm[0:caldm.find('_delays')]    
        self._asdm = asdm
        po.asALMA()
        po.loadTablesOnDemand(True)
        self.asdm = ASDM()
        print('setFromFile '+asdm)
        self.asdm.setFromFile(asdm, po)
                
        self.scanTimes = {}
        self.scanArrayTimes = {}
        self.antennaTemps = {}
        self.antennaPressures = {}
        # long integer -> to MJD
        self.observedTime = ((self._calData[0].startTimeObserved().get()
                              +self._calData[0].endTimeObserved().get())/2)/1e9/86400.
       #except:
        #    self._calData = None
        #    print "Could not read the result"

    def getPressureErrors(self):
        """
        """
        self._station = self.asdm.stationTable().get()
        self._antenna = self.asdm.antennaTable().get()

        self.heightError = {}
        self.centralPressure = {}
        self.deltaPressures = {}
        self.minPressure = 1e10
        self.maxPressure = -1e10

        antennas = []
        for  r in self._station:
            if str(r.name()) == "MeteoTB2":
                centralStationId = r.stationId()
            
        if centralStationId == Tag(0):
            print("== no central station")
            return
        refPos = self.asdm.stationTable().getRowByKey(centralStationId).position()
        refVector = pl.array([refPos[0].get(),refPos[1].get(),refPos[2].get()])
        for row in self.asdm.stationTable().get():
            ant = row.name()
            antennas.append(ant)
            stationId = row.stationId()
            r0 = self.asdm.stationTable().getRowByKey(stationId)
            pos =  r0.position()
            vector = pl.array([pos[0].get(), pos[1].get(), pos[2].get()])
            h1 = sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
            h1 -=  sqrt(refVector[0]**2+refVector[1]**2+refVector[2]**2)
            h0 = ellipsoidHeight(vector[0], vector[1], vector[2])
            h0 -= ellipsoidHeight(refVector[0], refVector[1], refVector[2])
            # print '%s / %s h1 %f h0 %f h1-h0 %f' %( ant, r0.name(),  h1, h0, h1-h0 )
            self.heightError[ant] = h0 - h1
        
        for row in self.asdm.calDataTable().get():
            scan = row.scanSet()[0]
            if scan not in list(self.scanArrayTimes.keys()):
                start = row.startTimeObserved().get()
                end = row.endTimeObserved().get()
            rows = self.asdm.weatherTable().getByContext(centralStationId)
            for r in rows:
                ttt = r.timeInterval().start().get()
                if (ttt > start) and (ttt < end):
                    found = True
                    self.centralPressure[scan] = r.pressure().get()
            self.deltaPressures[scan] = {}
            self.minPressure = min(self.minPressure, self.centralPressure[scan])
            self.maxPressure = max(self.minPressure, self.centralPressure[scan])

            for ant in antennas:
                self.deltaPressures[scan][ant] = self.centralPressure[scan]*6.5e-3/293.5*self.heightError[ant]*5.26

                    


        

    def getDelaysFromAntposFile(self, fileName):
        """
        """

        f = open(fileName)
        correctedName = fileName[0:fileName.find('baseline.dat')]+'pressure.dat' 
        ff = open(correctedName, 'w')
        refant = ''
        for l in f.readlines():
            # print l
            w = l.split()
            ll = ''
            if w[0] == 'SCAN':
                # first try to get the reference antenna
                if refant == '':
                    k = 8
                    while k < len(w):
                        if float(w[k+2]) == 0:
                            refant = w[k]
                        k += 4
                    
                scan  = int(w[1])
                deltaDir = pl.array([float(w[2]),float(w[3]),float(w[4])]) 
                dir = pl.array([float(w[5]),float(w[6]),float(w[7])])
                oldDir = dir - deltaDir
                for k in range(8):
                    ll += w[k] + ' '
                k = 8
                while k < len(w):
                    ant = w[k]
                    delay = float(w[k+1])
                    dp = self.deltaPressures[scan-1][ant]-self.deltaPressures[scan-1][refant]
                    oldDeltaDelay = 0.00227*dp/100/299792458.*oldDir[2]
                    dp = self.deltaPressures[scan][ant]-self.deltaPressures[scan][refant]
                    deltaDelay = 0.00227*dp/100/299792458. *dir[2] 
                    if ant == 'DA65': print(scan, ant, (deltaDelay-oldDeltaDelay)*1e12)
                    newDelay = delay - (deltaDelay - oldDeltaDelay)
                    ll += '%s  %10.4e %s %s '% (w[k], newDelay, w[k+2], w[k+3])
                    k += 4
                ll += '\n'
            else:
                ll = l
            ff.write(ll)
        ff.close()
        f.close()
        return correctedName

    def getPressures(self, flaggedmeteo, useWeatherStations=True, scaleHeight=500.):
        
        """
        Get for each antenna the pressure of the closest weather station
        """
        self._weather = self.asdm.weatherTable().get()
        self._station = self.asdm.stationTable().get()
        self._antenna = self.asdm.antennaTable().get()
        antennas = []
        wStationId = {}
        wStationName = {}
        wStationDistance = {}
        flagged_meteo = flaggedmeteo.split()
        count = {}
        self.meanDeltaPressure = {}
        
        centralStationId = Tag(0)
        #for r in self._station:
        #    if str(r.name()) == "MeteoCentral":
        #        centralStationId = r.stationId()
        for  r in self._station:
            if str(r.name()) == "MeteoTB2":
                centralStationId = r.stationId()
            
        if centralStationId == Tag(0):
            print("== no central station")
            return
        refPos = self.asdm.stationTable().getRowByKey(centralStationId).position()
        refVector = pl.array([refPos[0].get(),refPos[1].get(),refPos[2].get()])
        for row in self._antenna:
            ant = row.name()
            antennas.append(ant)
            count[ant] = 0
            self.meanDeltaPressure[ant] = 0
            if useWeatherStations:
                stationId = row.stationId()
                r0 = self.asdm.stationTable().getRowByKey(stationId)

                d2min = 1e12
                for r in self._station:
                    if (str(r.type()) == 'WEATHER_STATION') and (str(r.name()) not in flagged_meteo):
                        d2 = 0
                        for i in range(3):
                            d2 += (r0.position()[i].get()-r.position()[i].get())**2
                        if d2 < d2min: 
                            rows = self.asdm.weatherTable().getByContext(r.stationId())
                            # test th epressure
                            if rows[0].pressure().get() > 1000:
                                # 
                                wStationName[ant] = r.name()
                                wStationId[ant] = r.stationId()
                                wStationDistance[ant] = sqrt(d2)
                                d2min = d2
                print('%s/%s : Weather station %15s   distance %10.2f m' \
                    %(ant, r0.name(), wStationName[ant], wStationDistance[ant]))            
        
        self.deltaPressures = {}
        self.centralPressure = {}
        self.centralWaterPressure = {}
        self.centralTemperature = {}
        self.minPressure = 1e10
        self.maxPressure = -1e10
        
        for row in self.asdm.calDataTable().get():
            if str(row.calType()) == "CAL_WVR":
                scan = row.scanSet()[0]
                if scan not in list(self.scanArrayTimes.keys()):
                    start = row.startTimeObserved().get()
                    end = row.endTimeObserved().get()

                self.deltaPressures[scan] = {}
                rows = self.asdm.weatherTable().getByContext(centralStationId)
                for r in rows:
                    ttt = r.timeInterval().start().get()
                    if (ttt > start) and (ttt < end):
                        found = True
                        self.centralPressure[scan] = r.pressure().get()
                        self.centralTemperature[scan] = r.temperature().get()
                for wvrrow in self.asdm.calWVRTable().get():
                    #print wvrrow.calDataId(), row.calDataId()
                    if wvrrow.antennaName() == self.refAntenna:
                        if wvrrow.calDataId() == row.calDataId():
                            water = wvrrow.water().get() # meters
                            break
                # assuming scale height of 1000m
                scaleHeight = 1000.
                self.centralWaterPressure[scan] = self.centralTemperature[scan]*water*1000./217.*100*(1000./scaleHeight)  ## in pascals.
                print("=== scan %2s pres %7.3f mb temp %7.3f K w %6.3f mm ppH2O %7.3f mb" %\
                    (scan, self.centralPressure[scan]/100., self.centralTemperature[scan], water*1000, self.centralWaterPressure[scan]/100.))
                self.minPressure = min(self.minPressure, self.centralPressure[scan])
                self.maxPressure = max(self.minPressure, self.centralPressure[scan])

                for ant in antennas:
                    # print "antenna ", ant 
                    water = 0
                    for wvrrow in self.asdm.calWVRTable().get():
                        if wvrrow.antennaName() == ant:
                            if wvrrow.calDataId() == row.calDataId():
                                water = wvrrow.water().get() # meters
                                break
                    temp = self.centralTemperature[scan]
                    water_pressure = temp*water*1000./217.*100.*(1000./scaleHeight)  # pascals
                    self.deltaPressures[scan][ant] = \
                        - (water_pressure-self.centralWaterPressure[scan] )            
                    if useWeatherStations:
                        rows = self.asdm.weatherTable().getByContext(wStationId[ant])
                        sRow = self.asdm.stationTable().getRowByKey(wStationId[ant])
                        pos = sRow.position()
                        padVector = pl.array([pos[0].get(),pos[1].get(),pos[2].get()]) 
                        diffVector = padVector - refVector
                        diffHeight = sqrt(padVector[0]**2+padVector[1]**2+padVector[2]**2)
                        diffHeight -= sqrt(refVector[0]**2+refVector[1]**2+refVector[2]**2)
                        found = False
                        pres = 0
                        for r in rows:
                            ttt = r.timeInterval().start().get()
                            if (ttt > start) and (ttt < end):
                                found = True
                                pres = r.pressure().get()
                                temp = r.temperature().get()
                        if found:
                            self.deltaPressures[scan][ant] += \
                                pres - self.centralPressure[scan]*(1.-6.5e-3/293.5*diffHeight)**5.26 
                    # if scan>1:
                    self.meanDeltaPressure[ant] += self.deltaPressures[scan][ant]
                    count[ant] += 1

        for ant in list(count.keys()):
            self.meanDeltaPressure[ant] /= count[ant]
                
                 


    def getScanList(self):
        """
        Get the scan list from calData
        """
        
        scanList = []
        for row in self._calData:
            if str(row.calType()) == 'CAL_DELAY':
                scanList.append(row.scanSet()[0])
        return scanList

    def getCalPositions(self):
        """
        """
        self._calPositions = self.calPositionTable().get()
        self.positionOffset = {}
        self.positionError = {}
        self.refAntenna = None
        self.station = {}
        self.antennaPosition = {}
        self.stationPosition = {}
        self.absolutePosition = {}
        self.delayRms = {}
        self.residualOffset = {}
        for row in self._calPositions:
            # print row
            ant = row.antennaName()
            self.positionOffset[ant] = pl.array([row.positionOffset()[i].get() for i in range(3)])
            self.positionError[ant] = pl.array([row.positionErr()[i].get() for i in range(3)])
            
            self.station[ant] = row.stationName()
            self.stationPosition[ant] = pl.array([row.stationPosition()[i].get() for i in range(3)])
            self.antennaPosition[ant] = pl.array([row.antennaPosition()[i].get() for i in range(3)])
            self.absolutePosition[ant] = self.stationPosition[ant] + antennaToStation(self.stationPosition[ant], self.antennaPosition[ant]+self.positionOffset[ant])
            self.refAntenna = row.refAntennaNames()[0]
            # print self.refAntenna
            if self.refAntenna != ant:
                self.delayRms[ant] = row.delayRms()
            else:
                self.delayRms[ant] = 0
            for i in range(3):
                if self.positionError[ant][i] == 0:
                     self.positionError[ant][i] = 999.
        # print self.delayRms
        self.antennas = list(self.positionOffset.keys())
        # not in the data anyway.
        self.refPosition = pl.array([2225061.869, -5440061.953, -2481682.085])

    def getAntennaTemperatures(self):
        """
        """

        # Get scan times
        self.scanTimes = {}
        self.antennaTemps = {}
        for row in self.calDataTable().get():
            scan = row.scanSet()[0]
            if scan not in list(self.scanTimes.keys()):
                start = row.startTimeObserved().toFITS()
                end = row.endTimeObserved().toFITS()
                self.scanTimes[scan] = [start, end]
            # print self.scanTimes[scan]
            self.antennaTemps[scan] = {}
            for ant in self.antennas:
                if ant[0:2]=='CM':
                    # CM antennas: quadrapod, 4 total:
                    # GET_METR_TEMPS_04[0]...[3]
                    mpts = {'METR_TEMPS_04': [0,1,2,3]}
            
                elif  ant[0:2]=='DA':
                    # DA antennas: apex, 3 total:
                    # GET_METR_TEMPS_14[3] + GET_METR_TEMPS_15[0] + GET_METR_TEMPS_15[1]
                    mpts = {'METR_TEMPS_14': [3], 'METR_TEMPS_15':[0,1]}
            
                elif  ant[0:2]=='DV':
                    # DV antennas: hexapod-subreflector interface cylinder, 8 total:
                    # GET_METR_TEMPS_00[0]...[3] + GET_METR_TEMPS_01[0]...[3]
                    mpts = {'METR_TEMPS_00': [0,1,2,3], 'METR_TEMPS_01': [0,1,2,3]}
                elif  ant[0:2]=='PM':
                    # PM antennas: 'ambient' temperature, 1 total:
                    # GET_METR_TEMPS_1B[0]
                    mpts = {'METR_TEMPS_1B': [0]}
                temps = pl.array([])
                for m in list(mpts.keys()):
                    t = tm.get_tmc_data(ant,'Mount',m,
                                        self.scanTimes[scan][0],self.scanTimes[scan][1],
                                        removefile=True, verbose=False)
                    start = self.scanTimes[scan][0][0:self.scanTimes[scan][0].find('.')]
                    s = datetime.datetime.strptime(start, '%Y-%m-%dT%H:%M:%S')
                    end = self.scanTimes[scan][1][0:self.scanTimes[scan][1].find('.')]
                    e = datetime.datetime.strptime(end, '%Y-%m-%dT%H:%M:%S')
                    q = pl.find((pl.array(t['datetime']) > s) & (pl.array(t['datetime']) < e))
                    for k in mpts[m]:
                        temps = pl.append(temps, pl.array(t['value'])[q,k]/100.)
                print("== ", ant, temps.shape, temps.mean(), sqrt(temps.var()))
                self.antennaTemps[scan][ant] = pl.array(temps).mean()
                    

        



    def setRefAntenna(self, newRefAnt):
        """
        """
        print("setRefAnt: ",  self.refAntenna, newRefAnt)
        ra = self.refAntenna
        print(ra, self.positionError[ra], newRefAnt, self.positionError[newRefAnt])
        for ant in self.antennas:
            self.positionOffset[ant] -= self.positionOffset[newRefAnt]
            self.positionError[ant] = pl.sqrt(self.positionError[ant]**2+self.positionError[newRefAnt]**2)
        self.refAntenna = newRefAnt
        self.positionError[newRefAnt] = pl.zeros(3) 
        print(ra, self.positionError[ra], newRefAnt, self.positionError[newRefAnt])
        
    def getDelays(self, antennaList=None, delayResult=None):
        '''
        read the delays from on-line reduction or from off-line delayResult
        '''

        self.delayOffset = {}
        self.delayError = {}
        if antennaList == None:
            self.antennaList = []

        apc={'C': 'AP_CORRECTED', 'U': 'AP_UNCORRECTED'}
        corr={'AP_CORRECTED':'C', 'AP_UNCORRECTED':'U'}
        if delayResult != None:
            delay_cdm = ASDM()
            try:
                po = ASDMParseOptions()
                po.asALMA()
                po.loadTablesOnDemand(True)
                delay_cdm.setFromFile(delayResult, po)
            except:
                print("problem reading ", delayResult)
                return
        else:
            delay_cdm = self
            
        for row in delay_cdm.calDelayTable().get():
            calDataRow = delay_cdm.calDataTable().getRowByKey(row.calDataId())
            scan = calDataRow.scanSet()[0]
            antenna = row.antennaName()
            self._rant = row.refAntennaName()
            if (antennaList == None) and (antenna not in self.antennaList):
                self.antennaList.append(antenna)
            if (scan in self.scanList) and (antenna in self.antennaList):
                key = (antenna, corr[str(row.atmPhaseCorrection())], str(row.basebandName()), scan)
                self.delayOffset[key] = row.delayOffset()
                self.delayError[key] = row.delayError()
        
    def getDirections(self):
        '''
        get the source directions prom the phase cal results in the file
        (shoudl be able to get them from a off-line result as well)
        '''
        
        self.direction = {}
        for row in self.calPhaseTable().get():
            calDataRow = self.calDataTable().getRowByKey(row.calDataId())
            if calDataRow != None:
                scan = calDataRow.scanSet()[0]
                if (scan in self.scanList) and (scan not in list(self.direction.keys())):
                    az = row.direction()[0].get()
                    el = row.direction()[1].get()
                    self.direction[scan] = (sin(az)*cos(el), cos(az)*cos(el), sin(el))


    def getWetPath(self):
        '''
        Get the wet path from the on-line CalWVR table;
        should be able to get them from an off-line reduction as well.
        '''
        
        self.wetPath = {}
        for row in self.calWVRTable().get():
            calDataRow = self.calDataTable().getRowByKey(row.calDataId())
            antenna = row.antennaName()
            if calDataRow != None:
                scan = calDataRow.scanSet()[0]
                if (scan in self.scanList) and ((antenna,scan) not in list(self.wetPath.keys())):
                    self.wetPath[(antenna,scan)] = row.wetPath()[0]
    
                
                
    def solveAntenna(self, antenna, nc=8, correction='C', flag='', wetpath=0., diff = False):
        '''
        solve for antenna positions using baseband/polarization delays in on-line result
        optionally the fit uses scan-to-scan differences (diff=True)
        subtract the wet path from WVR (wetpath=1.) 
        '''

        delays = [[] for i in range(nc)] # 8 delays for each scan
        delayerrors = [[] for i in range(nc)] # 8 delays for each scan
        coords = [[] for i in range(3)] # 3 directions
        scans=[]
        paths=[]
        # print self.scanList
        flagList = [int(fl) for fl in flag.split()]
        print("flagging scans ", flagList)
        for scan in self.scanList:
            if scan not in flagList:
                # print scan
                for i in range(3):
                    coords[i].append(self.direction[scan][i])
                idelay = 0
                for ipolar in range(2):
                    ibaseband = 0
                    for baseband in ['BB_1','BB_2', 'BB_3', 'BB_4']:
                        if idelay<nc:
                            delays[idelay].append(self.delayOffset[(antenna, correction, baseband, scan)][ipolar]*1e12)
                            delayerrors[idelay].append(self.delayError[(antenna, correction, baseband, scan)][ipolar]*1e12)
                        ibaseband += 1
                        idelay += 1
                scans.append(int(scan))
                paths.append(self.wetPath[(antenna,scan)]-self.wetPath[(self._rant,scan)])

        # for a differentiaal fit:
        for i in range(8):
            delays[i] = pylab.array(delays[i])
            delayerrors[i] = pylab.array(delayerrors[i])
            if diff:
                delays[i] = delays[i][1:]-delays[i][0:-1] 
                delayerrors[i] = pl.sqrt(delayerrors[i][1:]**2+delayerrors[i][0:-1]**2)
        for i in range(3):
            coords[i] = pl.array(coords[i])
            if diff:
                coords[i] = coords[i][1:]-coords[i][0:-1]
        paths = -pylab.array(paths)/299792458.*1e12
        if diff:
            paths = paths[1:] - paths[0:-1]
            
        delays = pylab.array(delays)
        coords = pylab.array(coords)
        scans = pylab.array(scans)
        if diff:
            scans = scans[1:]
        delayerrors = pylab.array(delayerrors)
        delays = delays - paths*wetpath

        def errorfunc(pars, delays, coords):
            '''
            coords[0:2] = pointing direction cosines
            data[0:8] delays for 8 baseband / polars
            
            '''
            errs = []
            for i in range(nc):
                err = delays[i].copy()
                for j in range(3):
                    err -= pars[j]/1000. * coords[j] *1e12/299792458. # pars in mm, result in ps
                err -= pars[3+i]  # the delay offset, in ps
                errs.extend(err)
            errs= pl.array(errs)
            return errs

        p0 = pylab.zeros(nc+3)
        rawrms = sqrt(errorfunc(p0, delays, coords).var())
        pars, cov, info, mesg, ier = leastsq(errorfunc, p0, (delays, coords), full_output=1)
        rms = sqrt(errorfunc(pars, delays, coords).var())
        pars = pl.array(pars)
        epars = rms * pl.sqrt(cov)
            
        print("")
        print("Antenna %s: raw rms %10.3f ps -- fit rms %10.3f ps " % (antenna, rawrms, rms))
        print(" X= %10.3f (%5.3f) mm  Y= %10.3f (%5.3f) mm  Z=%10.3f (%5.3f) mm "% \
              (pars[0], epars[0,0], pars[1], epars[1,1], pars[2], epars[2,2]))
        offs = pars[3:]/1000.
        print("bb delays offsets X  %10.3f ns %10.3f ns %10.3f ns %10.3f ns" % (offs[0], offs[1], offs[2], offs[3]))
        print("bb delays offsets Y  %10.3f ns %10.3f ns %10.3f ns %10.3f ns" % (offs[4], offs[5], offs[6], offs[7]))

        #
        # compute baseband average and residuals
        np = len(delays[0])
        averDelays = pylab.zeros(np)
        averErrors = pylab.zeros(np)
        averResiduals = pylab.zeros(np)
        for ic in range(nc):
            averDelays += delays[ic]-pars[3+ic]
            averErrors += delayerrors[ic]**2
            averResiduals += delays[ic]-pars[3+ic]
        averDelays /= nc
        averResiduals /= nc
        averErrors = pl.sqrt(averErrors/nc)
        
        for i in range(3):
            averResiduals -= coords[i]*pars[i]/1000. *1e12/299792458.
            


        pl.clf()

        axes = ['X', 'Y', 'Z']

        for iaxis in range(3):
            ax = pl.subplot(2,2,iaxis+1)
            pl.errorbar(coords[iaxis], averDelays, fmt = '.', color='r')
            pl.errorbar(coords[iaxis], averResiduals, fmt = '.', color='g')
            pl.errorbar(coords[iaxis], paths, fmt = '.', color='c')
            if diff:
                pl.xlabel('$\Delta '+axes[iaxis])
            else:
                pl.xlabel(axes[iaxis])
            pl.ylabel('Delay [ps]')

        ax = pl.subplot(2,2,4)
        pl.errorbar(scans, averDelays, fmt = '.', color='r')
        pl.errorbar(scans, averResiduals, fmt = '.', color='g')
        pl.errorbar(scans, paths, fmt = '.', color='c')

        ax = pl.subplot(2,2,1)
        pl.text(0., 1.08, '%s %s ref %s'%(self._asdm, antenna, self._rant), transform = ax.transAxes) 
        pl.text(0., 1.01, 'X= %8.3f (%5.3f) mm Y= %8.3f(%5.3f) mm Z=%8.3f(%5.3f) mm rms=%8.3f ps'%\
                (pars[0], epars[0,0], pars[1], epars[1,1], pars[2], epars[2,2], rms),
                transform = ax.transAxes) 
        
        pl.savefig('%s-%s.png'%(self._asdm,a))
        # return pars

    def plotArray(self, coord='Z', xmax = 8., zmax=20., fit=False):
        """
        """

        pl.clf()
        xlabel = ['W (X) [km]', 'SE (-Y+X ) [km]','NE(X+Y) [km]']
        color = {'A': 'oc', 'J': '.k', 'N': 'or', 'S': 'og', 'P': 'ob', 'W': 'ro', 'T': '.k'}
        hscale=1000.
        # self.refPosition = pl.array([2225061.869, -5440061.953, -2481682.085])
        # print 'ref position ', refPosition
        icoord = ['X','Y','Z'].index(coord)

        # get horizontal positions
        self.stationVector = {}
        for a in self.antennas:
            self.stationVector[a] = stationToAntenna(self.refPosition,
                                                     self.stationPosition[a]-self.refPosition)
            self.stationVector[a] += self.antennaPosition[a]-self.antennaPosition[self.refAntenna]

        iplot = 0
        
        for k in range(3):
            iplot += 1
            ax = pl.subplot(3,1,iplot)
            for a in self.antennas:
                s = self.station[a]
                if k==0:
                    hscale = 1000.
                    pl.plot(self.stationVector[a][0]/hscale, self.positionOffset[a][icoord]*1000., color[s[0]])
                elif k==1:
                    hscale = 1000.*sqrt(2.)
                    pl.plot((-self.stationVector[a][0]+self.stationVector[a][1])/hscale,
                            self.positionOffset[a][icoord]*1000, color[s[0]])
                elif k==2:
                    hscale = 1000.*sqrt(2.)
                    pl.plot((self.stationVector[a][0]+self.stationVector[a][1])/hscale,
                            self.positionOffset[a][icoord]*1000, color[s[0]])
            pl.ylim(-zmax, zmax)
            pl.xlim(-xmax, xmax)
            pl.xlabel(xlabel[k])

            pl.ylabel('$\delta Z$  [mm]')

        pl.figtext(0.1, 0.95, 'Z offset vs horizontal antenna position  %s' % type)
        pl.savefig('Z_vs_Arms_%s.png' % type)

        if not fit:
            return

        X = [[],[],[]]
        Z = []
        W = []
        for a in self.antennas:
            if a != self.refAntenna:
                Z.append(self.positionOffset[a][icoord])
                W.append(1./self.positionError[a][icoord])
                for k in range(3):
                    X[k].append(self.stationVector[a][k])
        Z = pl.array(Z)
        X = pl.array(X)
        W = pl.array(W)

        def func(p, X):
            return p[0] +p[1]*X[0] +p[2]*X[1] +p[3]*X[0]**2+p[4]*X[0]*X[1]+p[5]*X[1]**2 \
                  +X[0]**3*p[6] +X[0]**2*X[1]*p[7] +X[0]*X[1]**2*p[8] + X[1]**3*p[9]

        def errorfunc(p, Z, X, W):
            '''
            subtract a 3rd order polynomial
            '''
            return W * (Z - func(p,X))
        
        npars=10
        p0 = pl.zeros(npars)
        val =  errorfunc(p0, Z, X, W)
        self.rawRms = sqrt(val.var())/sqrt((W*W).sum())

        pars, cov, info, mesg, ier = leastsq(errorfunc, p0, (Z, X, W), full_output=1)

        # print pars
        self.fitPolynomial = pars
        res = errorfunc(pars, Z, X, W)
        self.residualRms = sqrt(res.var())/sqrt((W*W).sum())
        print("rms raw %6.3f mm residual %6.3f mm "%( self.rawRms*1000., self.residualRms*1000.))
        ax = pl.subplot(3,1,1)
        hscale = 1000.
        pl.plot(X[0]/hscale, res*1000,'ow')
        pl.ylim(-zmax, zmax)
        pl.xlim(-xmax, xmax)
        ax = pl.subplot(3,1,2)
        hscale = 1000.*sqrt(2.)
        pl.plot((-X[0]+X[1])/hscale, res*1000,'ow')
        pl.ylim(-zmax, zmax)
        pl.xlim(-xmax, xmax)
        ax = pl.subplot(3,1,3)
        hscale = 1000.*sqrt(2.)
        pl.plot((X[0]+X[1])/hscale, res*1000,'ow')
        pl.ylim(-zmax, zmax)
        pl.xlim(-xmax, xmax)
        self.residualOffset = {}
        ra = self.refAntenna
        for a in self.antennas:
            self.residualOffset[a] = self.positionOffset[a][icoord] \
                                     - (func(pars, self.stationVector[a])-func(pars, self.stationVector[ra]))
            
        

        
    
        
