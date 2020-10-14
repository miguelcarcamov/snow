"""
Utilities for getting computing logs data for containers at the AOS

$Id: compUtils.py,v 1.25 2019/08/01 14:47:43 thunter Exp $
"""
from __future__ import print_function  # prevents adding old-style print statements
import os
import sys
import distutils.spawn
import shutil
import glob
import datetime
import numpy as np

def get_host_name():
    """
    Returns  the  hostname
    """
    hostname = 'http://computing-logs.aiv.alma.cl'
    return  hostname

def get_root_url_for_date(date):
    """
    Returns the root URL of the computing logs web I/F for the given date.

    The argument date should be an ISO-8601 date string (YYYY-MM-DD).
    The returned URL already contains the date.
    """
    year = date[:4]
    mm = date[5:7]
    hostname = get_host_name()
    if date <= '2016-10-01':
        aosdir = 'AOS'
    elif date <= '2017-01-27':
        aosdir = 'AOS64'
    else:
        aosdir = 'APE1'
    return "%s/index.php?dir=%s/CONTAINER/%s/" % (hostname, aosdir, date)

def get_root_url_for_abm_container(antenna,date):
    """
    Returns the root URL of the computing logs web I/F for the given date.

    The argument date should be an ISO-8601 date string (YYYY-MM-DD).
    The returned URL already contains the date.
    """
    url = get_root_url_for_date(date)
    url += 'alma/logs/' + antenna.lower() + '-abm/CONTROL/' + antenna.upper() + '/' 
    return url

def get_root_url_for_ccc_container(date):
    """
    Returns the root URL of the computing logs web I/F for the given date.

    The argument date should be an ISO-8601 date string (YYYY-MM-DD).
    The returned URL already contains the date.
    """
    url = get_root_url_for_date(date)
    url += 'alma/logs/cob-cc/CORR/CCC/' 
    return url

def retrieve_abm_container_data_files(antenna, date, time='*', overwrite=False, verbose=True):
    """
    Retrieve abm container data files via HTTP.

    Parameters are something like:
    antenna = 'DV01'
    container = 'acsStartContainer_cppContainer'
    date = '2010-04-24'  # ISO-8601 date or datetime string

    Return the path to the concatenated file if successful, otherwise '_CURL_FAILED_'.
    """

    isodate = get_datetime_from_isodatetime(date).date().strftime('%Y-%m-%d')
    inputdate = datetime.datetime.strptime(date, '%Y-%m-%d')
    
    rooturl = get_root_url_for_abm_container(antenna,date)

    today = datetime.datetime.today()
    twentydaysago = today + datetime.timedelta(days=-20)

    unzip = 1
    extension = 'txt.gz'

    completeurl = '%s' % (rooturl)
    completeurl = completeurl.replace('index.php?dir=','')
    directory = completeurl.replace('http://','')
    print(completeurl)
    if (os.path.exists(directory) == False or overwrite):
        print("Retrieving %s" % (completeurl))
        wget = distutils.spawn.find_executable('wget',path=':'.join(sys.path)+':'+os.environ['PATH'])
        cmd = wget + ' -r -l1 --no-parent -A.gz %s' % (completeurl) # -o %s' % (completeurl, outfile)
        # will write to new subdirectory: computing-logs.aiv.alma.cl/AOS/CONTAINER/2014-04-29/alma/logs/dv25-abm/CONTROL/DV25/
        print("Calling: ", cmd)
        exitcode = os.system(cmd)
        if exitcode == 0:
            if  unzip:
                files = glob.glob(directory+'/*.gz')
                for f in files:
                    os.system('gunzip -f %s' %f)
                files = glob.glob(directory+'/*')
                for f in files:
                    if (f.find('.tar') >= 0):
                        os.system('tar -C %s -xvf %s' % (directory,f))
                        os.remove(f)
            files = sorted(glob.glob(directory+'/*'))
            print(files)
            allfiles = catAllFiles(files, outfilename=files[0][:-13]) # strip off the time string at end
            print("concatenated file = ", allfiles)
            return allfiles
        else:
            print('Retrieval failed. Check permissions on directory and set outpath if necessary')
            return '_CURL_FAILED_'
    else:
        files = sorted(glob.glob(directory+'/*'))
        outfilename = files[0][:-13] # strip off the time string at end
        if not os.path.exists(outfilename):
            allfiles = catAllFiles(files, outfilename=outfilename) 
        if (verbose):
            print("Directory already present, returning name of file: ", outfilename)
        return outfilename

def retrieve_ccc_container_data_files(date, prefix='acsStartContainer_cppContainer', time='*', overwrite=False, verbose=True):
    """
    Retrieve ccc container data files via HTTP.
    Parameters are something like:
    date = '2010-04-24'  # ISO-8601 date or datetime string

    Return the path to the concatenated file if successful, otherwise '_CURL_FAILED_'.
    """

    isodate = get_datetime_from_isodatetime(date).date().strftime('%Y-%m-%d')
    inputdate = datetime.datetime.strptime(date, '%Y-%m-%d')
    
    rooturl = get_root_url_for_ccc_container(date)

    today = datetime.datetime.today()
    twentydaysago = today + datetime.timedelta(days=-20)

    unzip = 1
    extension = 'txt.gz'

    completeurl = '%s' % (rooturl)
    completeurl = completeurl.replace('index.php?dir=','')
    directory = completeurl.replace('http://','')
    print(completeurl)
    if (os.path.exists(directory) == False or overwrite):
        print("Retrieving %s" % (completeurl))
        wget = distutils.spawn.find_executable('wget',path=':'.join(sys.path)+':'+os.environ['PATH'])
        cmd = wget + ' -r -l1 --no-parent -A.gz %s' % (completeurl) # -o %s' % (completeurl, outfile)
        # will write to new subdirectory: computing-logs.aiv.alma.cl/AOS/CONTAINER/2014-04-29/alma/logs/dv25-abm/CONTROL/DV25/
        print("Calling: ", cmd)
        exitcode = os.system(cmd)
        if exitcode == 0:
            if  unzip:
                files = glob.glob(directory+'/%s*.gz' % (prefix))
                for f in files:
                    os.system('gunzip -f %s' %f)
                files = glob.glob(directory+'/*')
                for f in files:
                    if (f.find('.tar') >= 0):
                        os.system('tar -C %s -xvf %s' % (directory,f))
                        os.remove(f)
            files = sorted(glob.glob(directory+'/%s*'%(prefix)))
            print(files)
            allfiles = catAllFiles(files, outfilename=files[0][:-13]) # strip off the time string at end
            print("concatenated file = ", allfiles)
            return allfiles
        else:
            print('Retrieval failed. Check permissions on directory and set outpath if necessary')
            return '_CURL_FAILED_'
    else:
        files = sorted(glob.glob(directory+'/*'))
        outfilename = files[0][:-13] # strip off the time string at end
        if not os.path.exists(outfilename):
            allfiles = catAllFiles(files, outfilename=outfilename) 
        if (verbose):
            print("Directory already present, returning name of file: ", outfilename)
        return outfilename

def catAllFiles(files, outfilename='allfiles', remove=True):
    with open(outfilename, 'wb') as outfile:
        for filename in files:
            with open(filename) as readfile:
                shutil.copyfileobj(readfile, outfile)
            if (remove):
                os.remove(filename)
    return(outfilename)

def get_datetime_from_isodatetime(isodatetime):
    """
    Return a datetime.datetime object for given ISO-8601 date/datetime string.

    The argument isodatetime should be in YYYY-MM-DDThh:mm:ss or YYYY-MM-DD
    (in the latter case, 00:00:00 is assumed).
    Return 0001-01-01T00:00:00 if an invalid string is given.
    """

    datelist = isodatetime.split('T')
    if len(datelist) == 1:  # date only
        timelist = [0, 0, 0]
        datelist = datelist[0].split('-')
    elif len(datelist) == 2:  # date and time
        timelist = datelist[1].split(':')
        datelist = datelist[0].split('-')
    else:
        print("Date %s is invalid." % isodatetime)
        return datetime.date(1, 1, 1)

    if (len(datelist) == 3) and (len(timelist) == 3):
        microsec = int(1e6 * (float(timelist[2]) - int(float(timelist[2]))))
        timelist[2] = int(float(timelist[2]))
        return datetime.datetime( \
            int(datelist[0]), int(datelist[1]), int(datelist[2]), \
            int(timelist[0]), int(timelist[1]), int(timelist[2]), microsec )
    else:
        print("Date '%s' is invalid." % isodatetime)
        return datetime.date(1, 1, 1)

def parseIndexHtml(filename='index.html'):
    if not os.path.exists(filename):
        print("File not found: ", filename)
        return
    print("Trying to open ", filename)
    f = open(filename,'r')
    names = []
    lastTouched = []
    for line in f.readlines():
        if (line.find('a href="./') > 0):
            # this was the original code
            name = line.split('a href="./')[1].split('"')[0]
            print("name = ", name)
            names.append(name)
            datestamp = line.split('<td ')[2].split('>')[1].split('<')[0]
            print("name: %s, datestamp: %s" % (name,datestamp))
            lastTouched.append(convertDDMMMYYYYHHMM(datestamp))
        elif (line.find('a href="log') > 0):
            # this is what works on October 30, 2018
            name = "log" + line.split('a href="log')[1].split('"')[0]
            print("name = ", name)
            names.append(name)
            datestamp = line.split('<td ')[2].split('>')[1].split('<')[0]
            print("name: %s, datestamp: %s" % (name,datestamp))
            lastTouched.append(convertDDMMMYYYYHHMM(datestamp))
    names = np.array(names)
    lastTouched = np.array(lastTouched)
    names = names[np.argsort(lastTouched)]
    lastTouched = np.sort(lastTouched)
    return(names, lastTouched)

def convertDDMMMYYYYHHMM(date_input):
    """
    Converts a date/time string like '12-Oct-2018 11:04' to '2018-10-12 11:04:00'
    """
    date_input = date_input.strip()
    date_object = datetime.datetime.strptime(date_input, '%d-%b-%Y %H:%M')
    return(date_object.strftime('%Y-%m-%d %H:%M:%S'))

def withinOneDay(filedate, startDate):
    if (abs(computeIntervalBetweenTwoDays(filedate,startDate)) > 1):
        return False
    else:
        return True

def retrieve_aos_system_logs(startTime, stopTime, overwrite=False, verbose=False, aosdir=''):
    """
    Retrieve AOS/SYSTEM logs via HTTP.  They are always in .xml.gz format.
    Tries to get the minimum number of files that span the dataset, then unzips them.
    Will work for times that span (up to) two different days.
    Inputs:
    startTime = '2015-05-31 05:42:26 UT' or '2015-05-31T05:42:26'  or  '2015-05-31 05:42:26'
    stopTime =  '2015-05-31 06:42:26 UT'  etc.
    aosdir: typical options are 'AOS', 'AOS64', 'APE1', 'APE2', 'TFINT', 'OFFLINE'
          leaving it blank will try to determine it by the date
    Returns: a list of log files
    -Todd Hunter
    """
    startTime = startTime.rstrip(' UT').replace('T',' ')
    stopTime = stopTime.rstrip(' UT').replace('T',' ')
    startDate = startTime.split()[0]
    stopDate = stopTime.split()[0]
    dates = np.unique([startDate,stopDate])
    overlappingXMLFiles = []
    print("startDate = ", startDate)
    if aosdir == '':
        if startDate < '2016-10-01':
            # AOS ends on 2016-09-30
            aosdir = 'AOS'
        elif startDate < '2017-01-23':
            # AOS64 ends on 2017-01-22
            aosdir = 'AOS64'
        elif startDate < '2017-05-13':
            # APE1 starts 2016-12-02
            aosdir = 'APE1'
        if verbose: print("Using aosdir: ", aosdir)
    for myDate in dates:
        overlappingFiles = []
        rooturl = 'http://computing-logs.aiv.alma.cl/%s/SYSTEM/%s' % (aosdir,myDate)
        completeurl = rooturl
        directory = completeurl.replace('http://','')
        print("Retrieving list of files at %s" % (completeurl))
        wget = distutils.spawn.find_executable('wget',path=':'.join(sys.path)+':'+os.environ['PATH'])
        os.system('rm -f index.html')
        cmd = wget + ' -l1 --no-remove-listing %s' % (completeurl) 
        print("Running: ", cmd)
        exitcode = os.system(cmd)
        if os.path.exists('index.html'):
            indexFile = 'index.html'
        else:
            indexFile = myDate
        myresult = parseIndexHtml(indexFile)
        if myresult is None:
            return
        availableFiles, lastTouched = myresult
        print("Found %d total files" % (len(availableFiles)))
        if verbose: print(str(availableFiles))
        lastfile = -2
        for i,f in enumerate(availableFiles):
            if (lastTouched[i] > startTime):
                validStartTimeInFileName = False
                startTimeInFileName = os.path.basename(f).split('--')
                if (len(startTimeInFileName) == 1):
                    startTimeInFileName = os.path.basename(f).split('_')
                if len(startTimeInFileName) > 1:
                    validStartTimeInFileName = startTimeInFileName[1].replace('T',' ') < stopTime
                if (lastTouched[i] < stopTime or validStartTimeInFileName):
                    overlappingFiles.append(f)
                    lastfile = i
                elif (i-1 == lastfile):
                    # get the one final file that ends after the observing ends
                    overlappingFiles.append(f)
        print("Found %d overlapping files over %s to %s" % (len(overlappingFiles), startTime, stopTime))
        print(str(overlappingFiles))
        if (os.path.exists(directory)):
            allFilesExist = True
            for f in overlappingFiles:
                f2 = f.replace('.tar.gz','.xml')
                f2 = f2.replace('.gz','')
                f2 = f2.replace('__','--').replace('_',':')
                dir1 = os.path.join(directory,f)
                dir2 = os.path.join(directory,f2)
                if not os.path.exists(dir1) and not os.path.exists(dir2):
                    print("Setting allFilesExist=False because %s does not exist and %s does not exist" % (dir1,dir2))
                    allFilesExist = False
                    break
        else:
            print("Setting allFilesExist=False because directory %s does not exist" % (directory))
            allFilesExist = False
        if (not os.path.exists(directory) or overwrite or not allFilesExist):
            for i,f in enumerate(overlappingFiles):
                if (not os.path.exists(directory+'/'+f)):
                    cmd = wget + ' -r -l1 --no-parent %s' % (completeurl+'/'+f)
                    print("Executing: ", cmd)
                    exitcode = os.system(cmd)
        if allFilesExist:
            files = glob.glob(os.path.join(directory,'*.tar'))
            if len(files) > 0:
                target = directory
                for f in files:
                    cmd = 'tar xvf %s -C %s' % (f,target)
                    print("Running: ", cmd)
                    os.system(cmd)
                    cmd = 'rm ' + f
                    print("Running: ", cmd)
                    os.system(cmd)
        else:
            files = glob.glob(os.path.join(directory,'*.gz'))
            if (len(files) > 0):
                print("Gunzipping %d files" % (len(files)))
                for f in files:
                    os.system('gunzip -f %s' % f)
            files = glob.glob(os.path.join(directory,'*.tar'))
            if len(files) > 0:
                print("Untarring %d files" % (len(files)))
                target = directory
                for f in files:
                    cmd = 'tar xvf %s -C %s' % (f,target)
                    print("Running: ", cmd)
                    os.system(cmd)
                    cmd = 'rm ' + f
                    print("Running: ", cmd)
                    os.system(cmd)
        overlappingXMLFiles += glob.glob(os.path.join(directory,'*.xml'))
    return overlappingXMLFiles

def searchForLogFile(ts, loc='AOS'):
    """
    Determines the log file associated to a specific instant. - R. Amestica
    ts: timestamp of format '2015-07-23T17:29:42.378'
    loc: location, e.g. 'AOS'
    """
    import urllib2, re, dateutil.parser
    tsrex = '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}'
    frex = 'log(' + tsrex + ')(?:--|_)(' + tsrex + ')(?:--' + loc + ')?.xml.gz'
    pageurlfmt = 'http://computing-logs.aiv.alma.cl/index.php?page=%d&dir=' + loc + '/SYSTEM/%s/'
    logurlfmt = 'http://computing-logs.aiv.alma.cl/' + loc + '/SYSTEM/%s/%s'
    dts = dateutil.parser.parse(ts)
    url = pageurlfmt % (1, dts.date())
    html = urllib2.urlopen(url).read()
    npages = int(re.search('\tof ([0-9]+).*', html).groups()[0])
    for page in range(1, npages + 1):
        for candidate in [[dateutil.parser.parse(i.group(j)) for j in [1, 2]] + [i.group(0)] for i in re.finditer(frex, html)]:
            if candidate[0] <= dts and candidate[1] >= dts:
                return logurlfmt % ((candidate[0].date(), candidate[2].replace(':', '%3A')))
        html = urllib2.urlopen(pageurlfmt % (page, dts.date())).read()
    return None

def computeIntervalBetweenTwoDays(date1, date2):
    """
    Takes 2 strings of format 'YYYY-MM-DD' and returns the number of
    days between them.
    """
    d1 = date1.replace('-','').replace('/','')
    d2 = date2.replace('-','').replace('/','')
    delta = datetime.date(int(d1[0:4]),int(d1[4:6]),int(d1[6:])) - \
            datetime.date(int(d2[0:4]),int(d2[4:6]),int(d2[6:]))
    return(delta.days)

