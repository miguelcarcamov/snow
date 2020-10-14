#!/usr/bin/env python
# O2SounderPlayer.py

"""
O2SounderPlayer.py
started by dbarkats

$Id: O2SounderPlayer.py,v 1.1 2014/12/25 20:48:00 dbarkats Exp $
"""

# Execution starts here
import os, sys
import socket
import getpass
import pickle
import glob
username = getpass.getuser()
from pylab import *
from datetime import datetime,timedelta
import analysisUtils as au



class o2SounderPlayer(object):
    """
    
    """
  

    def __init__(self):
        """

        """
     
        

    def readLv2Data(self, filename):
        ##TODO
        #read the profiles at other places than Zenith. RIght now only reading Zenith, not N, S, or A
        ##
        """
        Simple function to read level2 data file
        output is P['t'], P['T'], P['VD'], P['LD'], P['RH'] (Profiles of T: temperature, VD: vapor density, LD: liquid density, RH: relative humidity).
        output is M['t'], M['Tamb'](K),M['Rh'],M['P'](mb),M['Tir'](K),M['R'] (rain)
        output is I['t'],I['IV'] (Int. Vapor(cm),I['IL'](Int. Liquid(mm)), I['CB'] Cloud Base(km) 
        """
        
        #initialize output variables
        P={}
        M={}
        I={}

        # height variable. TODO: should be read from file
        h = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25, 5.50, 5.75, 6.00, 6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25, 8.50, 8.75, 9.00, 9.25, 9.50, 9.75,10.00]

        profileCols = ['t','T','VD','LD','RH']
        metCols = ['t','Tamb','Rh','P','Tir','R']
        intCols = ['t','IV','IL','CB']
        for k in profileCols:
            P[k]=[]
        for k in metCols:
            M[k]=[]
        for k in intCols:
            I[k]=[]
            
        f= open(filename, 'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            
            sline = line.split(',')
            record = sline[0]
            if record == 'Record': continue
            recordType = sline[2]
            if recordType == '101': continue
            time = datetime.strptime(sline[1],"%m/%d/%y %H:%M:%S")
            #time = au.date2num(time)
            if recordType == '201':
                ind = 3
                M['t'].append(time)
                for i in metCols[1:]:
                    M[i].append(float(sline[ind]))
                    ind+=1
            if recordType == '301':
                ind = 3
                I['t'].append(time)
                for i in intCols[1:]:
                    I[i].append(float(sline[ind]))
                    ind+=1
            if recordType[0] == '4':
                angle = sline[3]
                if angle == 'Zenith':
                    P['t'].append(time)
                    if recordType[2]   == '1': ind = 'T'
                    elif recordType[2] == '2': ind = 'VD'
                    elif recordType[2] == '3': ind = 'LD'
                    elif recordType[2] == '4': ind = 'RH'
                    P[ind].append(map(float,sline[4:-1]))

            
        return h, P, M, I

    
            
        
        
