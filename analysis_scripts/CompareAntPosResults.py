# plot Antennas position results relative to one antenna
# First version imported by R. Lucas.  All subsequent edits by T. Hunter
#
from __future__ import print_function  # prevents adding old-style print statements
import math
from asdm import *
from AntPosResult import *
import pylab as pl
from WvrCorrection import *

def stationToAntenna(stationPosition, geo):
    xx,yy,zz = stationPosition[0],stationPosition[1],stationPosition[2]
    lat = asin(zz / sqrt(xx ** 2 + yy ** 2 + zz ** 2))
    lon = atan2(yy, xx)
    alma = [0, 0, 0]
    alma[0] = -geo[0] * sin(lon) + geo[1] * cos(lon)
    alma[1] = -geo[0] * cos(lon) * sin(lat) - geo[1] * sin(lon) * sin(lat) + geo[2] * cos(lat)
    alma[2] = geo[0] * cos(lon) * cos(lat) + geo[1] * sin(lon) * cos(lat) + geo[2] * sin(lat)
    return alma
        
    
class CompareAntPosResults:

    def __init__(self):

        self.refAntenna = None
        self.results = []
        self.resultColors = []
        self.wvrResults = {}
        self.antennas=[]
        self.stationPosition = {}
        self.antennaPosition = {}
        self.station = {}

    def addResult(self, result, wvrResult = None):

        if len(self.results) == 0:
            self.refAntenna = result.refAntenna
        elif self.refAntenna != result.refAntenna:
            if self.refAntenna in result.antennas:
                ra = result.refAntenna
                print(ra, result.positionError[ra])
                result.setRefAntenna(self.refAntenna)
                print(self.refAntenna, result.positionError[ra])
            else:
                print('Current reference antenna not available in result')
                return
        for a in result.antennas:
            if a not in self.antennas:
                # Store the first position encountered
                self.antennas.append(a)
                self.stationPosition[a] = result.stationPosition[a]      

                self.antennaPosition[a] = result.antennaPosition[a]
                self.station[a] = result.station[a]
            else:
                # check for changes
                stationChange = result.stationPosition[a]-self.stationPosition[a] 
                antennaChange = result.antennaPosition[a]-self.antennaPosition[a]
                padChange = self.station[a] != result.station[a]
                dp = pl.sqrt((stationChange**2).sum())
                da = pl.sqrt((antennaChange**2).sum())
                if padChange:
                    print("%s: %s had been moved from pad %s to pad %s" %\
                          (result._asdm, a, result.station[a], self.station[a]))
                    print("=> ignoring the data set for %s" % a)
                    result.positionError[a] = [999.,999.,999.]
                else:
                    if da>1.e-6:
                        print("%s: %s pad %s antenna vector moved by %.6f m " %\
                            (result._asdm, a, result.station[a], da))
                        result.positionOffset[a] += antennaChange
                        result.antennaPosition[a] = self.antennaPosition[a]
                        print("=> change added to position offset")
                    if dp>1.e-6:
                        print("%s: %s pad %s pad vector moved by %.6f m" %\
                            (result._asdm, a, result.station[a], dp))
                        antennaChange = stationToAntenna(result.stationPosition[a], stationChange)
                        result.positionOffset[a] += antennaChange
                        result.stationPosition[a] = self.stationPosition[a]
                        print("=> change rotated and added to position offset")
               #  print '== ', a, result.positionError[a]
            
        self.antennas.sort()
        self.results.append(result)
        if wvrResult != None:
            self.wvrResults[result._asdm] = wvrResult
        self.refPosition = result.refPosition
        

    def plotAntenna(self, antenna, errmax=2.0, doResiduals = False, greenThreshold=1.0,
                    minSnr=5.0, doWvrCorrection=False):

        pl.clf()
        axis=['X [mm]','Y [mm]','Z [mm]']
        k = 0
        wmean = {}
        station = self.station[antenna]
        wmax = 0.
        for r in self.results:
            k += 1
            wv = []
            water = self.wvrResults[r._asdm].water
            for key in list(water.keys()):
                if antenna==key[1]:
                    wv.append(water[key])
            wv = pl.array(wv)
            wmean[r] = wv.mean()
            if wmax < wv.mean():
                wmax = wv.mean()
        rms = pl.zeros(3)
        mean = pl.zeros(3)
        meanAbsolute = pl.zeros(3)
        median = pl.zeros(3)
        mad = pl.zeros(3)
        ax = [None,None,None,None]
        distance = sqrt(((self.stationPosition[antenna]-self.refPosition)**2).sum())
        height = stationToAntenna(self.refPosition, self.stationPosition[antenna]-self.refPosition)[2]
        
        items = []
        texts = []
        if doWvrCorrection:
            rows = 4
        else:
            rows = 3
        for i in range(rows): # dX, dY, dZ, (dP)
            error = []
            weight = []
            resError = []
            absolutePosition = []
            for j in range(2):
                ax[i] = pl.subplot(rows, 3, 1+3*i+j)

                if j==0: 
                    pl.xlabel('Water [mm]')
                if j==1:
#                    pl.xlabel('Date [days since %s]'%(getObservationStartDateFromASDM(r._asdm)))
                    pl.xlabel('Date [days]')
                mjds = []
                for r in self.results:
                    nAntennas = len(list(r.positionOffset.keys()))
                    text = r._asdm[r._asdm.find('_X')+1:] + ' (%d)'%(nAntennas)
                    item = []
                    # plot in time order
                    k = len(self.results) - self.results.index(r)
                    # rather, use time (mjd)
                    mjd = r.observedTime - self.results[0].observedTime
                    # only one scan in the tc_antpos result...
                    deltaT = 0.
                    if antenna in list(r.positionOffset.keys()):
                        if r.antennaTemps != {}:
                            scan = list(r.antennaTemps.keys())[0]
                            if antenna in list(r.antennaTemps[scan].keys()):
                                deltaT = r.antennaTemps[scan][antenna] - r.antennaTemps[scan][self.refAntenna]
                                print(antenna, deltaT)
                        if j==0:
                            if i<3:
                                absolutePosition.append(r.absolutePosition[antenna][i])
                                error.append(r.positionOffset[antenna][i])
                                weight.append(1./r.positionError[antenna][i]**2)
                                if i==2 and doResiduals:
                                    resError.append(r.residualOffset[antenna])
                        if (rows > 3):
                            dp = r.meanDeltaPressure[antenna] - r.meanDeltaPressure[r.refAntenna]
                        # subtract the delta pressure effect!! 
                        # if (j==0) and (i==2):
                        #    r.positionOffset[antenna][i] -= dp/100*0.0027
                        # print "=== ", antenna, i, j, r.positionError[antenna]
                        #if antenna == 'DA65': 
                        #    print '== %s dp %10.3f mbar' % (antenna, dp/100.)
                        if i<3:
                            if r.positionError[antenna][i] <= 0.001*errmax:
                                if j==0:
                                    item = pl.errorbar(wmean[r]*1000., r.positionOffset[antenna][i]*1000,
                                                       yerr =  r.positionError[antenna][i]*1000, fmt=r.color,
                                                       markersize=5)
                                    
                                elif j==1:
                                    # item =  pl.errorbar(deltaT, r.positionOffset[antenna][i]*1000,
                                    #                    yerr =  r.positionError[antenna][i]*1000, fmt=r.color)
                                    mjds.append(mjd)                                    
                                    item =  pl.errorbar(mjd, r.positionOffset[antenna][i]*1000,
                                                        yerr =  r.positionError[antenna][i]*1000, fmt=r.color,
                                                        markersize=5)
                            else:
                                print("Position error of %s (%.1fmm) is > threshold to include (%.1fmm) on the vs. Date plot." % (antenna,r.positionError[antenna][i]*1000, errmax))
                            if i==2 and doResiduals:
                                if j==0:
                                    item = pl.plot(wmean[r]*1000., r.residualOffset[antenna]*1000, 'ow',
                                                   markersize=5)
                                elif j==1:
                                    mjds.append(mjd)                                    
                                    item = pl.plot(mjd, r.residualOffset[antenna]*1000, 'ow',
                                                   markersize=5)
                        else:
                            # The bottom right corner plot
                            if j==0:
                                item = pl.errorbar(wmean[r]*1000., dp/100,
                                                   fmt=r.color,
                                                   markersize=5)
                            elif j==1:
                                mjds.append(mjd)                                    
                                item =  pl.errorbar(mjd, dp/100,
                                                    fmt=r.color,
                                                    markersize=5)

                        if (text not in texts) and (len(item)>0):
                            items.append(item[0])
                            texts.append(text)
                            
                            
                if j==0:
                    if i<3:
                        pl.ylabel(axis[i])
                    else:
                        pl.ylabel('Delta p [mb]')
                    pl.xlim(0, wmax*1000.+0.5)
                elif j==1:
                    n = len(self.results)-1
                    # The previous logic assumed the datasets were processed in order of observation.
                    # mjdmin = self.results[len(self.results)-1].observedTime - self.results[0].observedTime
                    # mjdmax = 0.5
                    # The new logic does not assume this ordering.
                    mjdmin = min(mjds) 
                    mjdmax = max(mjds)
                    pl.xlim(mjdmin-0.5, mjdmax+0.5)

            error = pl.array(error)
            absolutePosition = pl.array(absolutePosition)
            resError = pl.array(resError)
            weight = pl.array(weight)

            if i<3:
                meanAbsolute[i] = pl.average(absolutePosition, weights=weight)
                mean[i] = pl.average(error, weights=weight)
                median[i] = pl.median(error)
                mad[i] = self.computeMAD(error)

                
                # contributions from intrisic error and dispersion of measurements...
                rms[i] = sqrt(1./pl.sum(weight) + pl.average((error-mean[i])**2, weights=weight))
                # rms[i] = sqrt(1./pl.sum(weight))
                # rms[i] = sqrt(pl.average((error-mean[i])**2, weights=weight))
            if i==2 and doResiduals:
                resMean = pl.average(resError, weights=weight)
                rms[i] = sqrt(1./pl.sum(weight) + pl.average((error-mean[i])**2, weights=weight))
                # resRms = sqrt(1./pl.sum(weight))
                # resRms = sqrt(pl.average((resError-resMean)**2, weights=weight))

        if len(items) > 0:
            pl.figlegend(items, texts, loc='upper right', labelspacing=0.1,
                         prop={'size':'small'})
        coord=['X','Y','Z']
        recommendUpdate = 0
        for i in range(3):
            snr = abs(mean[i])/rms[i]
            if (snr > minSnr):
                if (mean[i]*1000 < greenThreshold):
                    mycolor = 'g' # it is small but worth updating
                else:
                    mycolor='r' # it is large and worth updating
                recommendUpdate += 1
            else:
                mycolor='k'
            pl.figtext(0.96, 0.04+(3-i)*0.03,
                       '%s wmean %6.3f rms %6.3f mm'%(coord[i],mean[i]*1000, rms[i]*1000),
                       ha='right', va='bottom', color=mycolor)
        # median and MADs
        for i in range(3):
            pl.figtext(0.96, 0.17+(3-i)*0.03,
                       '%s med. %6.3f mad %6.3f mm'%(coord[i],median[i]*1000, mad[i]*1000),
                       ha='right', va='bottom')
        # total correction
        pl.figtext(0.96, 0.16, 'total correction %6.3f mm' % (1000*(mean[0]**2+mean[1]**2+mean[2]**2)**0.5),
                   ha='right', va='bottom')
        if doResiduals:
            pl.figtext(0.96, 0.04,
                       '%s resMean %6.3f rms %6.3f mm'%(coord[2],resMean*1000, resRms*1000),
                       ha='right', va='bottom')
        return pl.array(mean)*1000, pl.array(rms)*1000, distance, height, meanAbsolute, recommendUpdate

    def computeMAD(self, a, c=0.6745, axis=0):
        a = np.array(a)
        if a.ndim == 1:
            d = pl.median(a)
            m = pl.median(np.fabs(a - d) / c)
        elif (a.ndim > 1):
            d = pl.median(a, axis=axis)
            if axis > 0:
                aswp = swapaxes(a,0,axis)
            else:
                aswp = a
            m = pl.median(np.fabs(aswp - d) / c, axis=0)
        else:
            m = 0
        return m
    
    def plotAll(self, type=None, errmax=2., doWvrCorrection=False,
                returnRms=False, drawFits=[]):
        aa = []
        xx = []
        yy_x = []
        yy_y = []
        yy_z = []
        pngs = []
        for a in self.antennas:
            if a != self.refAntenna:
                mean, rms, distance, height, meanAbsolute, recommendUpdate = self.plotAntenna(a, errmax=errmax, doWvrCorrection=doWvrCorrection)
                station = self.station[a]
                pl.figtext(0.1, 0.95, '%s/%s d=%6.0fm h=%3.0fm %s'%\
                           (a, station, distance, height, type), ha='left', va='top')
                png = 'XYZ_%s_%s_%s.png'%(type, a, station)
                pl.savefig(png)
                pngs.append(png)
                print('%s/%s %s %6.0f %3.0f X=%6.3f(%5.3f)  Y=%6.3f(%5.3f)  Z=%6.3f(%5.3f) [mm]' %\
                      (a, station, type, distance, height, mean[0], rms[0], mean[1], rms[1], mean[2], rms[2]))
                xx.append(distance)
                yy_z.append(rms[2]) # Z component
                yy_y.append(rms[1])
                yy_x.append(rms[0])
                aa.append(a)
        pl.clf()
        # Just get the plot limits, since Z will have the largest scatter
        pl.plot(xx,yy_z,'or')
        ylims = pl.ylim()

        pl.clf()
        pl.plot(xx,yy_x,'or')
        slope = ''
        if (len(drawFits) > 0):
            pl.plot([np.min(xx),np.max(xx)], [drawFits[0][0]+drawFits[1][0]*np.min(xx),
                                     drawFits[0][0]+drawFits[1][0]*np.max(xx)], 'k-')
            slope = '(slope=%.3fmm/km)' % (drawFits[1][0]*1000)
        pl.ylim(ylims)
        pl.xlabel('Distance (m)')
        pl.ylabel('rms of X component (mm)')
        pl.title('%d datasets from %s to %s %s' % (len(self.results), self.mjdToDate(self.results[-1].observedTime),
                                                self.mjdToDate(self.results[0].observedTime), slope))
        for i,a in enumerate(aa):
            if (xx[i] > 2000 or yy_x[i] > 1.0):
                pl.text(xx[i],yy_x[i],a,size=8,weight='bold',va='top',ha='left')
        png = 'x_rms_vs_distance.png'
        pl.savefig(png)
        pngs.append(png)

        pl.clf()
        pl.plot(xx,yy_y,'or')
        if (len(drawFits) > 0):
            pl.plot([np.min(xx),np.max(xx)], [drawFits[0][1]+drawFits[1][1]*np.min(xx),
                                     drawFits[0][1]+drawFits[1][1]*np.max(xx)], 'k-')
            slope = '(slope=%.3fmm/km)' % (drawFits[1][1]*1000)
        pl.ylim(ylims)
        pl.xlabel('Distance (m)')
        pl.ylabel('rms of Y component (mm)')
        pl.title('%d datasets from %s to %s %s' % (len(self.results), self.mjdToDate(self.results[-1].observedTime),
                                                self.mjdToDate(self.results[0].observedTime),slope))
        for i,a in enumerate(aa):
            if (xx[i] > 2000 or yy_y[i] > 1.0):
                pl.text(xx[i],yy_y[i],a,size=8,weight='bold',va='top',ha='left')
        png = 'y_rms_vs_distance.png'
        pl.savefig(png)
        pngs.append(png)

        pl.clf()
        pl.plot(xx,yy_z,'or')
        if (len(drawFits) > 0):
            pl.plot([np.min(xx),np.max(xx)], [drawFits[0][2]+drawFits[1][2]*np.min(xx),
                                     drawFits[0][2]+drawFits[1][2]*np.max(xx)], 'k-')
            slope = '(slope=%.3fmm/km)' % (drawFits[1][2]*1000)
        pl.xlabel('Distance (m)')
        pl.ylabel('rms of Z component (mm)')
        pl.title('%d datasets from %s to %s %s' % (len(self.results), self.mjdToDate(self.results[-1].observedTime),
                                                self.mjdToDate(self.results[0].observedTime),slope))
        for i,a in enumerate(aa):
            if (xx[i] > 2000 or yy_z[i] > 1.0):
                pl.text(xx[i],yy_z[i],a,size=8,weight='bold',va='top',ha='left')
        png = 'z_rms_vs_distance.png'
        pl.savefig(png)
        pngs.append(png)
        if returnRms:
            return(pngs, [aa, xx, yy_x, yy_y, yy_z])
        else:
            return(pngs)

    def addResults(self, antposlist, type=None, errmax=2., fit=False,
                   useWeatherStations=False, doWvrCorrection=False ):

        colors=['ro','go','bo','co','yo','mo','ko',
                'rs','gs','bs','cs','ys','ms','ks',
                'rv','gv','bv','cv','yv','mv','kv',
                'r^','g^','b^','c^','y^','m^','k^',
                'r<','g<','b<','c<','y<','m<','k<',
                'r>','g>','b>','c>','y>','m>','k>']
        icolor = 0
        for r in antposlist:
            result = r
            if type == 'DRY':
                k  = r.find('_WET')
                result = r[0:k]+'_DRY'+r[k+4:]
        
            asdm = result[:result.find('_delays')]
            # print result
            apr = AntPosResult()
            apr.getData(result)
            apr.getCalPositions()
            apr.color = colors[icolor % 42]
            if doWvrCorrection:
                print("== apr.getPressures('')")
                apr.getPressures('', useWeatherStations=useWeatherStations )
            wvr = WvrCorrection(asdm)
            wvr.getAsdm()
            wvr.getWvrCoefficients()
            
            # apr.plotArray(coord='Z', xmax = 8., zmax=20., fit=fit)
            self.addResult(apr, wvrResult = wvr)
            icolor += 1

    def mjdToDate(self, mjd):
        jd = mjd + 2400000.5
        jd = jd + 0.5
        F, I = math.modf(jd)
        I = int(I)
        A = math.trunc((I - 1867216.25)/36524.25)
        if I > 2299160:
            B = I + 1 + A - math.trunc(A / 4.)
        else:
            B = I
        C = B + 1524
        D = math.trunc((C - 122.1) / 365.25)
        E = math.trunc(365.25 * D)
        G = math.trunc((C - E) / 30.6001)
        day = C - E + F - math.trunc(30.6001 * G)
        if G < 13.5:
            month = G - 1                
        else:
            month = G - 13
        if month > 2.5:
            year = D - 4716            
        else:
            year = D - 4715
        return '%4d-%02d-%02d' % (year, month, day)
