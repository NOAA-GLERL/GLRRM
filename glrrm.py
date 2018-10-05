#/bin/python

#----------------------------------------------------------------
#  This is a really basic example of how the overall GLRRM might
#  be structured to use the databank repository.  The actual
#  implementation of GLRRM will be very different, once it is
#  developed.  This is mainly intended as a demonstration of how
#  to use the repository, while also illustrating the basic structure
#  that Tim Hunter had in mind when he developed the design documents.
#-----------------------------------------------------------------

import sys
from copy import copy, deepcopy
import datetime

import databank
import databank_io
import databank_util

#------------------------------------------------------
#  Define lake surface areas for each basin
#  We will use the coordinated values, in sq meters
sup_area = 8.21e10
mhu_area = 1.17e11
stc_area = 1.11e09
eri_area = 2.57e10

#------------------------------------------------------
#  Define a few other constants that I will use.
#  e.g. I am assuming a fixed 30-day month
#
seconds_per_month = 30.0 * 86400.0

#---------------------------------------------------------------------------------
def read_main_config(filename):
    outdir = None
    sdate  = None
    edate  = None
    suplev = -99.9
    mhulev = -99.9
    stclev = -99.9
    erilev = -99.9

    with open(filename, "r") as f:
        for line in f:
            s1 = line_for_parsing(line).strip()
            if len(s1) > 0:
                if s1.find('#') < 0:
                    p=parse1(s1)
                    if (p):
                        if p[0].strip() == 'title1':
                            title1 = p[1]
                        elif p[0].strip() == 'title2':
                            title2 = p[1]
                        elif p[0].strip() == 'title3':
                            title3 = p[1]
                        elif p[0].strip() == 'startdate':
                            print('Found start date: ', p[1])
                            y,m,d = p[1].split(',')
                            sdate = datetime.date(int(y), int(m), int(d))
                        elif p[0].strip() == 'enddate':
                            print('Found end date: ', p[1])
                            y,m,d = p[1].split(',')
                            edate = datetime.date(int(y), int(m), int(d))
                        elif p[0].strip() == 'output directory':
                            outdir = p[1].strip()
                        elif p[0].strip() == 'sup start level':
                            s = p[1].split()
                            suplev = float(s[0])
                        elif p[0].strip() == 'mhu start level':
                            s = p[1].split()
                            mhulev = float(s[0])
                        elif p[0].strip() == 'st. c start level':
                            s = p[1].split()
                            stclev = float(s[0])
                        elif p[0].strip() == 'eri start level':
                            s = p[1].split()
                            erilev = float(s[0])
                    else:
                        print('parse1 failed')
                        print('  line=[' + s1 + ']') 

    return outdir, sdate, edate, suplev, mhulev, stclev, erilev



#---------------------------------------------------------------------------------
#  Given a line of text (i.e. a string)
#  1. strip off any comments (everything at/after the first # character)
#  2. convert to all lowercase
#
def line_for_parsing(line):
    i = line.find('#')
    s = line
    if (i >= 0):
       s = line[0:i-1].rstrip()
    return s.lower()


#---------------------------------------------------------------------------------
#  Given a line of text (i.e. a string)
#  Look for the first colon.  If you find one return a tuple with
#  everything to the left of the colon and everything to the right
#  of the colon.
#  If no colon is found, return None.
#  e.g.
#     'startdate: 2001,01,01'          -> ('startdate', ' 2001,01,01')
#     'starttime: 2001-01-01 18:00:00' -> ('starttime', ' 2001-01-01 18:00:00')
#     'startdate= 2001,01,01'          -> None
#
def parse1(line):
    i=line.find(':')
    if (i > 0):
        a = line[0:i]
        b = line[i+1:]
        return a,b


#---------------------------------------------------------------------------------
#  Given a start and end date, how many months are included?  Note that the 
#  day is ignored.  Wouldn't want to do that in a real model, but this is just
#  a quick and dirty example.
#  Given the following values:
#   sdate          edate          num_months
#   2017-01-15     2017-01-15       1
#   2017-01-15     2017-04-15       4
#   2017-01-31     2017-01-01       1      # note this one well
#   2017-01-15     2018-06-30       18
#
def number_of_months(sdate, edate):
    sdtt = sdate.timetuple()
    edtt = edate.timetuple()
    sm = sdtt.tm_mon
    sy = sdtt.tm_year
    em = edtt.tm_mon
    ey = edtt.tm_year
    nm = (ey-sy)*12 + em - sm + 1
    return nm



#---------------------------------------------------------------------------------
#  A nonsensical MONTHLY Lake Superior regulation
#
def silly_supreg(dvault, sdate, edate, sup_startlev):
    #
    #  Get nbs values for the entire period of interest.
    #  Notice that we are getting the NBS values expressed in
    #  units of meters over the lake surface.  Databank is
    #  doing that conversion for us.
    #
    #  Note that nbs_sup is a full DataSeries object. The actual 
    #  data values are accessible in nbs_sup.dataVals
    #
    nbs_sup = dvault.withdraw(kind='nbs', units='meters', intvl='mon', loc='sup',
                              first=sdate, last=edate)

    num_months = number_of_months(sdate, edate)
    
    #
    #  Create blank lists for the 3 timeseries we will create.
    #  Daily Sup levels
    #  Daily StMarys flows
    #
    suplevd = list(range(num_months+1))
    smrflow = list(range(num_months+1))

    slev_today = sup_startlev
    for i in range(0, num_months):
            
        #
        #  St Marys flow is computed by an extremely unrealistic
        #  rule created by Tim for this demo. 
        #
        #  Assume a basic flow in the St. Marys River of 2100 cms.
        #  If the nbs for Sup is more than that, increase the flow
        #  by 1/2 of the excess.
        #  If the nbs for Sup is less than that, decrease the flow
        #  by 1/2 of the deficiency.
        #  If the resulting level of Sup is > 183.80m or < 183.00m, then
        #  adjust the flow to keep Superior levels within that range.
        #
        nbs = (nbs_sup.dataVals[i] * sup_area) / seconds_per_month     # convert to cms
        smflow = 2100.0 + (nbs-2100)/2.0               # st mary's flow in cms
        cis   = nbs - smflow/sup_area                  # change in storage for Sup, in cms
        slevd = (cis * seconds_per_month) / sup_area   # sup level delta for today, meters
        newlev = slev_today + slevd
        if newlev > 183.80:
            adj = (newlev - 183.8) * sup_area           # cubic meters
            smflow = smflow + adj/seconds_per_month
            newlev = 183.80
        elif newlev < 183.00:
            adj = (183.0 - newlev) * sup_area           # cubic meters
            smflow = smflow - adj/seconds_per_month
            newlev = 183.00
            
        suplevd[i] = newlev
        smrflow[i] = smflow
        slev_today = newlev
        
    #
    #  Now store the suplevd and smrflow sequences in the vault.  Remember that
    #  the deposit_data() method will build the DataSeries object using the specified
    #  metadata and data values. If we had already built a DataSeries object we could 
    #  use dvault.deposit() and pass it the DataSeries object.
    #
    dvault.deposit_data(kind='eomlev', units='meters', intvl='mon', loc='sup', 
                     first=sdate, last=edate, values=suplevd)
    dvault.deposit_data(kind='flow', units='cms', intvl='mon', loc='stmarys', 
                     first=sdate, last=edate, values=smrflow)

    return True

#---------------------------------------------------------------------------------
def silly_midlakes(dvault, sdate, edate, mhu_startlev, stc_startlev, eri_startlev):
    nbs_mhu = dvault.withdraw(kind='nbs', units='meters', intvl='mon', loc='mhu',
                              first=sdate, last=edate)
    nbs_stc = dvault.withdraw(kind='nbs', units='meters', intvl='mon', loc='stc',
                              first=sdate, last=edate)
    nbs_eri = dvault.withdraw(kind='nbs', units='meters', intvl='mon', loc='eri',
                              first=sdate, last=edate)
    smrflow = dvault.withdraw(kind='flow', units='cms', intvl='mon', loc='smr',
                              first=sdate, last=edate)
                              
    num_months = number_of_months(sdate, edate)
    
    #
    #  Create blank lists for the timeseries we will create.
    #  Daily lake levels
    #  Daily river flows
    #
    mhulevd = [None] * (num_months+1)
    stcflow = [None] * (num_months+1)
    
    mlev_today = mhu_startlev
    slev_today = stc_startlev
    for i in range(0, num_months):
        
        #
        #  Now determine what the result would be for Mhu levels, using
        #  that St. Marys flow along with the Mhu nbs.
        #
        smrf = smrflow.dataVals[i]
        mnbs = nbs_mhu.dataVals[i]
        mhu_cms = smrf + (mnbs*mhu_area / seconds_per_month)     # total supply in cms
        
        #  Assume a basic flow in the St. Clair River of 5100 cms.
        #  If mhu_cms is more than that, increase the flow by 1/2 of the excess.
        #  If mhu_cms is less than that, decrease the flow by 1/2 of the deficiency.
        #  If the resulting level of Mhu is > 177.50m or < 175.50m, then
        #  simply adjust the level to stay in range.
        #
        scflow = 5100.0 + (mhu_cms-5100)/2.0
        mlevd = (mhu_cms * seconds_per_month) / mhu_area    # mhu level delta for today, meters
        newlev = mlev_today + mlevd
        if newlev > 177.50:
            adj = (newlev - 177.5) * mhu_area           # cubic meters
            scflow = scflow + adj/seconds_per_month
            newlev = 177.50
        elif newlev < 175.50:
            adj = (175.5 - newlev) * sup_area           # cubic meters
            scflow = scflow - adj/seconds_per_month
            newlev = 183.00
        mhulevd[i] = newlev
        mlev_today = newlev
        stcflow[i] = scflow
        
        
        #  At this point, Midlakes should compute:
        #    Lake StClair level,
        #    Detroit River flow,
        #    Lake Erie level
        #    Niagara River flow
        #  But I'm gonna skip that for now, because this is all just
        #  ridiculous calculation for the sole purpose of illustrating how
        #  databank would be used.  No need to do those for this purpose.
        #  
        
        
    #
    #  Store the MH level and St Clair river flow sequences in the vault.  Remember that
    #  the deposit_data() method will build the DataSeries object using the specified
    #  metadata and data values. If we had already built a DataSeries object we could 
    #  use dvault.deposit() and pass it the DataSeries object.
    #
    dvault.deposit_data(kind='eomlev', units='meters', intvl='mon', loc='mhu', 
                     first=sdate, last=edate, values=mhulevd)
    dvault.deposit_data(kind='flow', units='cms', intvl='mon', loc='stcriver', 
                     first=sdate, last=edate, values=stcflow)

    return True

        
#--------------------------------------------------------------------------------
#
#  Create the data repository object.  This will persist
#  IN MEMORY until the script exits.  It does NOT create any
#  kind of run-to-run persistent file or database.  This is
#  the behavior requested by the committee.
#
the_vault = databank.DataVault()

#
#  Read the main configuration file.  This file will specify
#  which models are being used, dates, etc.  i.e. the overall
#  behavior for this run. The procedure returns a tuple with:
#  (output_dir, start_date, end_date, sup_level, mh_level,
#  stc_level, eri_level).  A real implementation, of course,
#  would need to contain much more.
#
maincfg = read_main_config('./data/example/example_config.txt')
outdir      = maincfg[0]
model_sdate = maincfg[1]
model_edate = maincfg[2]
suplev      = maincfg[3]
mhulev      = maincfg[4]
stclev      = maincfg[5]
erilev      = maincfg[6]

#
#  Read monthly NBS values for each lake, and store each of them
#  into the vault.
#
#  Note that each of these read operations creates/clears the
#  DataSeries object, which I reuse, because I don't need to keep
#  it around persistently.  The data is stored into the vault and I
#  can retrieve a copy at any time.
#
#  Another thing to note is that I don't need to worry about units
#  at this stage, because the databank will normalize things internally
#  and I can specify the units I need when I retrieve the data.
#
#  Here I am reading files that have the data in mm for the first 3 
#  lakes, and cms for Erie. This is largely to illustrate/prove that
#  this works.
#
nbs = databank_io.read_file('data/example/nbs_sup_mm.txt')
the_vault.deposit(nbs)

nbs = databank_io.read_file('data/example/nbs_mhu_mm.txt')
the_vault.deposit(nbs)

nbs = databank_io.read_file('data/example/nbs_stc_mm.txt')
the_vault.deposit(nbs)

nbs = databank_io.read_file('data/example/nbs_eri_cms.txt')
the_vault.deposit(nbs)

#
#  This short code section is not really relevant to the larger goal
#  of this fake model. It is, instead, just a test/proof of how the
#  databank handles bad data and how it can be used for multiple 
#  "sets" of similar data. Here I am reading a file that contains
#  another "set" of monthly Erie NBS values. The file is identical 
#  to nbs_eri_cms.txt, but with some bad data values
#  in the file. This will show how the databank handles that.
#  Notice how I specify set='t1' to differentiate this data from the 
#  previously read monthly Erie NBS data. Both will be stored in the
#  vault.
#
nbs = databank_io.read_file('data/example/nbs_eri_faulty.txt', missing_value=-99999, set='t1')
the_vault.deposit(nbs)
nbs = None     # reset the object
nbs = the_vault.withdraw(kind='nbs', units='mm', intvl='mon', loc='eri', set='t1')
fileout = outdir + 'eri_data_test.txt'
databank_io.write_file(filename=fileout, file_format="table", missing_value=-99999.9,
                       series=nbs, width=10, prec=2, overwrite=True)

#
#  Now I need to know what the overall period of record is for
#  the data that I stored.  I will retrieve each of the monthly
#  NBS data sets, and find the overlapping period of record.
#
#  I don't care about the units because I will not be using the
#  actual data here, but I have to specify something, so I arbitrarily 
#  have chosen 'cfs'.
#
data_start = datetime.date(1900, 1, 1)      # initialize to far past date
data_end   = datetime.date(2099, 1, 1)      # initialize to far future date

nbs = the_vault.withdraw(kind='nbs', units='cfs', intvl='mon', loc='sup')
if nbs.startDate > data_start:
    data_start = nbs.startDate
if nbs.endDate < data_end:
    data_end = nbs.endDate

    
nbs = the_vault.withdraw(kind='nbs', units='cfs', intvl='mon', loc='mhu')
if nbs.startDate > data_start:
    data_start = nbs.startDate
if nbs.endDate < data_end:
    data_end = nbs.endDate


nbs = the_vault.withdraw(kind='nbs', units='cfs', intvl='mon', loc='stc')
if nbs.startDate > data_start:
    data_start = nbs.startDate
if nbs.endDate < data_end:
    data_end = nbs.endDate

    
nbs = the_vault.withdraw(kind='nbs', units='cfs', intvl='mon', loc='eri')
if nbs.startDate > data_start:
    data_start = nbs.startDate
if nbs.endDate < data_end:
    data_end = nbs.endDate
   
#
#  Do I have sufficient data stored to run my model for the period that
#  was specified in the config file?
#
if model_sdate < data_start:
    print('Insufficient data to run the model.')
    print('Data in files starts too late.')
    sys.exit(1)

if model_edate > data_end:
    print('Insufficient data to run the model.')
    print('Data in files ends too soon.')
    sys.exit(1)

#
#  Compute the sequence of lake levels using our
#  super-unrealistic water balance models.
#  Note how these models retrieve their inputs from the
#  vault and store their outputs into the vault.
#  I am simply passing to them the required period of record
#  and starting levels.
#
print('run silly_supreg')
ok = silly_supreg(the_vault, model_sdate, model_edate, suplev)
if not ok:
    print('Error in the Superior regulation model')
    sys.exit(1)

print('run silly_midlakes')
ok = silly_midlakes(the_vault, model_sdate, model_edate, mhulev, stclev, erilev)
if not ok:
    print('Error in the middle lakes model')
    sys.exit(1)

print('output timeseries data into files')
ds = the_vault.withdraw(kind='eomlev', units='meters', intvl='mon', loc='sup', 
                     first=model_sdate, last=model_edate)
fileout = outdir+'suplev_test.txt'
databank_io.write_file(filename=fileout, file_format="table",
                       series=ds, overwrite=True)
                    
ds = the_vault.withdraw(kind='eomlev', units='meters', intvl='mon', loc='mhu', 
                     first=model_sdate, last=model_edate)
fileout = outdir+'mhulev_test.txt'
databank_io.write_file(filename=fileout, file_format="table",
                       series=ds, overwrite=True)
                    
ds = the_vault.withdraw(kind='flow', units='cms', intvl='mon', loc='stmarys', 
                     first=model_sdate, last=model_edate)
fileout = outdir+'smrflow_test.txt'
databank_io.write_file(filename=fileout, file_format="table",
                       series=ds, overwrite=True)
                    
ds = the_vault.withdraw(kind='flow', units='cms', intvl='mon', loc='stcriver', 
                     first=model_sdate, last=model_edate)
fileout = outdir+'stcflow_test.txt'
databank_io.write_file(filename=fileout, file_format="table",
                       series=ds, overwrite=True)
    

