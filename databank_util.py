import datetime as dt
import databank

#==========
#  IMPORTANT:  Please note the above use of "import databank".  This is
#  an inelegant circular import that offends my sensibilities, but I'm kind of
#  stuck with it unless I were to combine databank.py and databank_util.py into
#  a single file, which is also ugly.  I only do it for use in the getLakeArea()
#  function, where I need to be able to look up the normalized data location
#  name. In truth, this seems to imply to me that I should reorganize this
#  whole thing, moving all of those names,etc into a separate .py file, but 
#  I really don't want to do that at this time. So I am gonna let this one bit 
#  of nastiness persist -- for now.
#  Tim Hunter  2018sep17
#===========


#-------------------------
#  Define a "missing value" for dates and other variable types.
#  For dates, I am using datetime.date(9999,9,9) because it's easily
#  recognized as an out-of-range value, and none of our valid data will 
#  ever have year=9999.
#  For numeric values, I am just assigning very large negative numbers.
#-------------------------
MISSING_DATE = dt.date(9999, 9, 9)
MISSING_REAL = -9.9e29
MISSING_INT  = -999999999
MISSING_TEST = -999999998

linear_units = ('mm',  'cm',  'm',  'km',  'in',  'ft',  'yd',  'mi')
areal_units  = ('mm2', 'cm2', 'm2', 'km2', 'in2', 'ft2', 'yd2', 'mi2')
cubic_units  = ('mm3', 'cm3', 'm3', 'km3', 'in3', 'ft3', 'yd3', 'mi3')
rate_units   = ('cms', '10cms', 'cfs', 'tcfs')

_coordLakeArea = (
    ('su', 8.21e10),
    ('mi', 5.78e10),
    ('hu', 5.96e10),
    ('sc', 1.11e09),
    ('er', 2.57e10),
    ('on', 1.90e10),
    ('mh', 1.17e11),
)

#---------------------------
#  define the breaking point for each quarter-month
#---------------------------
qtr_month_start_end_days = (
        ((1,7), (8,14), (15,21), (22,28)),         # qtr-months for 28-day month
        ((1,7), (8,14), (15,21), (22,29)),         # qtr-months for 29-day month
        ((1,8), (9,15), (16,23), (24,30)),         # qtr-months for 30-day month
        ((1,8), (9,15), (16,23), (24,31))          # qtr-months for 31-day month
)


#-------------------------------------------------------------------------------
#  Some useful utility functions
#-------------------------------------------------------------------------------
#
#--------------------------------------------------------------------
def getLakeArea(loc=None):
    if not loc: return None
    try:
        l = databank.DataLocation(loc)
        s = l.primaryName().lower()
        #
        #  First handle the obvious case where the specified location 
        #  is the lake.
        #
        for t in _coordLakeArea:
            if t[0].lower() == s:
                return t[1]
        
        #
        #  What to do when the location is a river? This presents
        #  a bit of a dilemma. I am, for the now anyway, going to
        #  make an assumption that the proper lake area to use is the
        #  "upstream" lake.  So for the St. Mary's river, this will return
        #  the area of Lake Superior. For the StClair River, you get the
        #  area of Michigan-Huron. For Detroit River you get Lake StClair.
        #  For Niagara River, you get Lake Erie.
        #  
        #  This may need to be revisited as actual models are developed.
        #  Tim Hunter -- 2018Sep17
        #
        if s=='smr':  return _coordLakeArea[0][1]      # lake superior
        if s=='scr':  return _coordLakeArea[6][1]      # lake mich-hur
        if s=='det':  return _coordLakeArea[2][1]      # lake stclair
        if s=='nia':  return _coordLakeArea[3][1]      # lake erie
        
        return None
    except:
        print('Error: Bad location (' + loc + ') in getLakeArea')
        return None

#--------------------------------------------------------------------
def days_in_month(year=None, month=None):
    ''' Determine # of days in month.  This is maybe out of place... 
        not intended to be called by user but just a helper function 
    '''
    if year==None:
        raise Exception('No year specified in days_in_month()')
    if month==None:
        raise Exception('No month specified in days_in_month()')
        
    try:
        if month==12:
            return (dt.datetime(year+1,1,1)-dt.datetime(year,12,1)).days
        else:
            return (dt.datetime(year,month+1,1)-dt.datetime(year,month,1)).days
    except:
        raise Exception('Error computing number of days in a month')


#--------------------------------------------------------------------
def getFridayDate(year=None, month=None, day=None):
    '''Get the date of the most recent (i.e. PRECEDING) friday given a yr,mo, day.  
	This is to match the previous convention that the CGLRRM used for weekly data beginning on friday.'''

    in_date=dt.date(year, month, day)
    idx=in_date.weekday() 
	# idx= 0  1  2  3  4  5  6
	#      M  T  W  R  F  Sa Su
    shift=[3, 4, 5, 6, 0, 1, 2]

    friday_date = in_date-dt.timedelta(days=shift[idx])
    return friday_date


#-----------------------------------------------------------------------------
#  Determine the last day of the current month for any datetime.date object.
#  e.g. if any_date = datetime.date(2001,3,15), this will return datetime.date(2001,3,31)
#-----------------------------------------------------------------------------
def last_day_of_month(any_date):
    next_month = any_date.replace(day=28) + dt.timedelta(days=4)
    return next_month - dt.timedelta(days=next_month.day)

#-----------------------------------------------------------------------------
#  Given a passed argument, attempt to translate it into a valid
#  datetime.date object.
#  The item passed in may be datetime object or a string.
#  If it's a datetime object, just pass it right back as the return value.
#  If it's a string, try to parse it by assuming isoformat (yyyy-mm-dd).  If that
#  works, pass back the valid datetime.date object. If not, pass back the defined
#  datetime.date object for MISSING_DATE.
#-----------------------------------------------------------------------------
def date_from_entry(in_date):
    if isinstance(in_date, dt.date):
        return in_date
    elif isinstance(in_date, str):
        try:
            return dt.datetime.strptime(in_date, '%Y-%m-%d').date()
        except:
            return MISSING_DATE
    else:
        return MISSING_DATE

#--------------------------------------------------------------------
def days_in_qtr_mon(year=None, month=None, qtr=None):
    ''' Determine number of days in the qtrmonth'''
    if not year:
        raise Exception('No year specified in days_in_qtr_mon()')
    if not month:
        raise Exception('No month specified in days_in_qtr_mon()')
    if not qtr:
        raise Exception('No quarter specified in days_in_qtr_mon()')

    try:
        sd, ed = getQtrMonthStartEnd(year=year, month=month, qtr=qtr)
        return ed - sd + 1
    except:
        raise Exception('Error computing number of days in a qtr-month')

#--------------------------------------------------------------------
def getQtrMonthStartEnd(year=None, month=None, qtr=None):
    ''' Determine first and last day of the qtrmonth'''
    if not year:
        raise Exception('No year specified in getQtrMonthStartEnd()')
    if not month:
        raise Exception('No month specified in getQtrMonthStartEnd()')
    if not qtr:
        raise Exception('No quarter specified in getQtrMonthStartEnd()')

    try:
        days = days_in_month(year=year, month=month)
        sd = qtr_month_start_end_days[days-28][qtr-1][0]
        ed = qtr_month_start_end_days[days-28][qtr-1][1]
        return sd, ed
    except:
        raise Exception('Error finding start/end of a qtr-month')

#-------------------------------------------------------
#  values = list of data values
#           Any value < MISSING_TEST is considered "missing"
#  oldstart = starting date for the old data list
#  oldend   = ending date for the old data list
#  newstart = starting date for the result data list
#  newend   = ending date for the result data list
#  intvl    = interval of the data.  Must be one of the following:
#             ['dy', 'wk', 'qm', 'mn']
#
#  The dates must have already been verified to be set such that
#  oldstart <= newstart  and  newend <= oldend.
#-------------------------------------------------------
def trimDataValues(values=None, oldstart=None, oldend=None,
                newstart=None, newend=None, intvl=None):
    if not values:   return None
    if not oldstart: return None
    if not oldend:   return None
    if not newstart: return None
    if not newend:   return None
    if not intvl:    return None
    
    #
    #  Transform (if needed) dates into datetime.date objects
    #
    try:
        olds = date_from_entry(oldstart)
        olde = date_from_entry(oldend)
        news = date_from_entry(newstart)
        newe = date_from_entry(newend)
    except:
        raise Exception('Invalid date specification for trimDataValues()')

    ok = True
    if olds == MISSING_DATE:  ok = False
    if olde == MISSING_DATE:  ok = False
    if news == MISSING_DATE:  ok = False
    if newe == MISSING_DATE:  ok = False
    if not ok:
        raise Exception('Invalid date specification for trimDataValues()')

    #
    #  Transform interval spec into the primary name
    #
    ok = False
    if intvl.lower()=='dy':  ok = True
    if intvl.lower()=='wk':  ok = True
    if intvl.lower()=='qm':  ok = True
    if intvl.lower()=='mn':  ok = True
    if not ok:
        raise Exception('Invalid interval specified to trimDataValues()')
    
    #
    if intvl.lower()=='dy':
        # index of start/end days
        d1 = (news-olds).days
        d2 = (newe-olds).days
        return values[d1:d2]
        
    if intvl.lower()=='mn':
        # index of start/end months
        y = news.year - olds.year
        m = news.month - olds.month
        m1 = y*12 + m
        y = newe.year - olds.year
        m = newe.month - olds.month
        m2 = y*12 + m + 1
        return values[m1:m2]

#-------------------------------------------------------
#  value = data value to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm', 'inch' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  area  = area in sq meters
#  secs  = number of seconds over which the linear amount was accumulated
#-------------------------------------------------------
def valueLinearToRate(value=None, oldu=None, newu=None, area=None, secs=None):
    if value==None: return None
    if not oldu:  return None
    if not newu:  return None
    if not area:  return None
    if not secs:  return None

    try:
        v1 = [value]
        m  = linearConvert(v1, oldu, 'm')
        vcms = [MISSING_REAL if v<MISSING_TEST else v*(area/secs) for v in m]
        v2 = rateConvert(vcms, 'cms', newu)
        return v2[0]
    except:
        raise Exception('Unable to convert ' + oldu + '->' + newu)

        
#-------------------------------------------------------
#  value = data value to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm3', 'in3' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  secs  = number of seconds over which the volume was accumulated
#-------------------------------------------------------
def valueCubicToRate(value=None, oldu=None, newu=None, secs=None):
    if value==None: return None
    if not oldu:  return None
    if not newu:  return None
    if not secs:  return None

    try:
        v1 = [value]
        m3 = cubicConvert(v1, oldu, 'm3')
        vcms = [MISSING_REAL if v<MISSING_TEST else v/secs for v in m3]
        v2 = rateConvert(vcms, 'cms', newu)
        return v2[0]
    except:
        raise Exception('Unable to convert ' + oldu + '->' + newu)

        
#-------------------------------------------------------
#  value = data value to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm3', 'in3' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  secs  = number of seconds over which the volume was accumulated
#-------------------------------------------------------
def valueRateToLinear(value=None, oldu=None, newu=None, area=None, secs=None):
    if value==None: return None
    if not oldu:  return None
    if not newu:  return None
    if not area:  return None
    if not secs:  return None

    try:
        v1 = [value]
        vcms = rateConvert(v1, oldu, 'cms')
        vm = [MISSING_REAL if v<MISSING_TEST else v*secs/area for v in vcms]
        v2 = linearConvert(vm, 'm', newu)
        return v2[0]
    except:
        raise Exception('Unable to convert ' + oldu + '->' + newu)

#-------------------------------------------------------
#  value = data value to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm3', 'in3' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  secs  = number of seconds over which the volume was accumulated
#-------------------------------------------------------
def valueRateToCubic(value=None, oldu=None, newu=None, secs=None):
    if value==None: return None
    if not oldu:  return None
    if not newu:  return None
    if not secs:  return None

    try:
        v1 = [value]
        vcms = rateConvert(v1, oldu, 'cms')
        vm3 = [MISSING_REAL if v<MISSING_TEST else v*secs for v in vcms]
        v2 = cubicConvert(vm3, 'm3', newu)
        return v2[0]
    except:
        raise Exception('Unable to convert ' + oldu + '->' + newu)

#-------------------------------------------------------
#  values = list of data values
#           Any value < MISSING_TEST is considered "missing"
#  oldstr = unit string for incoming data (e.g. 'mm', 'm', 'ft')
#  newstr = unit string for outgoing data (e.g. 'mm', 'm', 'ft')
#-------------------------------------------------------
def linearConvert(values=None, oldstr=None, newstr=None):
    if not values: return None
    if not oldstr: return None
    if not newstr: return None
    if not isinstance(oldstr, str):
        raise Exception('Invalid old units specification in linearConvert.')
    if not isinstance(newstr, str):
        raise Exception('Invalid new units specification in linearConvert.')

    try:    
        #
        #  Determine the conversion factor...
        #  m1 = old units -> meters
        #  m2 = meters -> new units
        #
        m1 = 0.0
        if oldstr=='mm': m1 = 0.001 
        if oldstr=='cm': m1 = 0.01 
        if oldstr=='m' : m1 = 1.0 
        if oldstr=='km': m1 = 1000.0 
        if oldstr=='in': m1 = 0.0254
        if oldstr=='ft': m1 = 0.3048
        if oldstr=='yd': m1 = 0.9144
        if oldstr=='mi': m1 = 1609.34
                
        m2 = 0.0
        if newstr=='mm': m2 = 1000.0 
        if newstr=='cm': m2 = 100.0 
        if newstr=='m' : m2 = 1.0 
        if newstr=='km': m2 = 0.001 
        if newstr=='in': m2 = 39.3701
        if newstr=='ft': m2 = 3.28084
        if newstr=='yd': m2 = 1.09361
        if newstr=='mi': m2 = 0.0006214
                
        mult = m1 * m2
        if mult > 0.0:
            return [MISSING_REAL if v<MISSING_TEST else v * mult for v in values]
        else:
            raise Exception('Invalid conversion specified: ' + oldstr 
                          + '->' + newstr)
    except:
        raise Exception('Error converting ', + oldstr + '->' + newstr)
        
#-------------------------------------------------------
#  values = list of data values
#           Any value < MISSING_TEST is considered "missing"
#  oldstr = unit string for incoming data (e.g. 'mm2', 'm2', 'ft2')
#  newstr = unit string for outgoing data (e.g. 'mm2', 'm2', 'ft2')
#-------------------------------------------------------
def arealConvert(values=None, oldstr=None, newstr=None):
    if not values: return None
    if not oldstr: return None
    if not newstr: return None
    if not isinstance(oldstr, str):
        raise Exception('Invalid old units specification in arealConvert.')
    if not isinstance(newstr, str):
        raise Exception('Invalid new units specification in arealConvert.')

    try:    
        #
        #  Determine the conversion factor...
        #  m1 = old units -> square meters
        #  m2 = square meters -> new units
        #
        m1 = 0.0
        if oldstr=='mm2': m1 = 1.0e-6 
        if oldstr=='cm2': m1 = 0.0001
        if oldstr=='m2' : m1 = 1.0 
        if oldstr=='km2': m1 = 1.0e+6 
        if oldstr=='in2': m1 = 0.00064516
        if oldstr=='ft2': m1 = 0.092903
        if oldstr=='yd2': m1 = 0.836127
        if oldstr=='mi2': m1 = 2.59e+6
                
        m2 = 0.0
        if newstr=='mm2': m2 = 1.0e+6 
        if newstr=='cm2': m2 = 10000.0 
        if newstr=='m2' : m2 = 1.0 
        if newstr=='km2': m2 = 1.0e-6 
        if newstr=='in2': m2 = 1550
        if newstr=='ft2': m2 = 10.7639
        if newstr=='yd2': m2 = 1.19599
        if newstr=='mi2': m2 = 3.861e-7
                
        mult = m1 * m2
        if mult > 0.0:
            return [MISSING_REAL if v<MISSING_TEST else v * mult for v in values]
        else:
            raise Exception('Invalid conversion specified: ' + oldstr 
                          + '->' + newstr)
    except:
        raise Exception('Error converting ', + oldstr + '->' + newstr)
        
#-------------------------------------------------------
#  values = list of data values
#           Any value < MISSING_TEST is considered "missing"
#  oldstr = unit string for incoming data (e.g. 'mm3', 'm3', 'ft3')
#  newstr = unit string for outgoing data (e.g. 'mm3', 'm3', 'ft3')
#-------------------------------------------------------
def cubicConvert(values=None, oldstr=None, newstr=None):
    if not values: return None
    if not oldstr: return None
    if not newstr: return None
    if not isinstance(oldstr, str):
        raise Exception('Invalid old units specification in cubicConvert.')
    if not isinstance(newstr, str):
        raise Exception('Invalid new units specification in cubicConvert.')

    try:    
        #
        #  Determine the conversion factor...
        #  m1 = old units -> cubic meters
        #  m2 = cubic meters -> new units
        #
        m1 = 0.0
        if oldstr=='mm3': m1 = 1.0e-9 
        if oldstr=='cm3': m1 = 1.0e-6
        if oldstr=='m3' : m1 = 1.0 
        if oldstr=='km3': m1 = 1.0e+9 
        if oldstr=='in3': m1 = 1.6387e-5
        if oldstr=='ft3': m1 = 0.0283168
        if oldstr=='yd3': m1 = 0.764555
        if oldstr=='mi3': m1 = 4.168e+9
                
        m2 = 0.0
        if newstr=='mm3': m2 = 1.0e+9 
        if newstr=='cm3': m2 = 1.0e+6 
        if newstr=='m3' : m2 = 1.0 
        if newstr=='km3': m2 = 1.0e-9
        if newstr=='in3': m2 = 61023.7
        if newstr=='ft3': m2 = 35.3147
        if newstr=='yd3': m2 = 1.30795
        if newstr=='mi3': m2 = 2.3991e-10
                
        mult = m1 * m2
        if mult > 0.0:
            return [MISSING_REAL if v<MISSING_TEST else v * mult for v in values]
        else:
            raise Exception('Invalid conversion specified: ' + oldstr 
                          + '->' + newstr)
    except:
        raise Exception('Error converting ', + oldstr + '->' + newstr)
        
#-------------------------------------------------------
#  values = list of data values
#           Any value < MISSING_TEST is considered "missing"
#  oldstr = unit string for incoming data (e.g. 'cms', 'tcfs' )
#  newstr = unit string for outgoing data (e.g. 'cms', 'tcfs')
#-------------------------------------------------------
def rateConvert(values=None, oldstr=None, newstr=None):
    if not values: return None
    if not oldstr: return None
    if not newstr: return None
    if not isinstance(oldstr, str):
        raise Exception('Invalid old units specification in rateConvert.')
    if not isinstance(newstr, str):
        raise Exception('Invalid new units specification in rateConvert.')

    try:    
        #
        #  Determine the conversion factor...
        #  m1 = old units ->  meters per second
        #  m2 = meters per second -> new units
        #
        m1 = 0.0
        if oldstr=='cms':   m1 = 1.0
        if oldstr=='10cms': m1 = 10.0
        if oldstr=='cfs' :  m1 = 0.0283168
        if oldstr=='tcfs':  m1 = 28.3168 
                
        m2 = 0.0
        if newstr=='cms':   m2 = 1.0
        if newstr=='10cms': m2 = 0.1
        if newstr=='cfs' :  m2 = 35.3147
        if newstr=='tcfs':  m2 = 0.0353147
                
        mult = m1 * m2
        if mult > 0.0:
            return [MISSING_REAL if v<MISSING_TEST else v * mult for v in values]
        else:
            raise Exception('Invalid conversion: ' + oldstr + '->' + newstr)
    except:
        raise Exception('Error converting ' + oldstr + '->' + newstr)


#-------------------------------------------------------
#  values = list of data values to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm', 'in' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  area  = effective area in square meters
#  intvl = interval of the data ('dy', 'wk', 'qm', 'mn')
#  first = start date (datetime.date)
#  last  = end date (datetime.date)
#-------------------------------------------------------
def linearToRate(values=None, oldu=None, newu=None, area=None, 
                 intvl=None, first=None, last=None):
    if not values: return None
    if not oldu:   return None
    if not newu:   return None
    if not area:   return None
    if not intvl:  return None
    if not first:  return None
    if not last:   return None

    #
    #  daily data is relatively clean/easy because each value
    #  can use the same number of seconds for conversion.
    #  Thus, we can process it using list comprehension.
    #
    if intvl.lower()=='dy':
        secs = 86400
        try:
            m  = linearConvert(values, oldu, 'm')
            vcms = [MISSING_REAL if v<MISSING_TEST else v*(area/secs) for v in m]
            newv = rateConvert(vcms, 'cms', newu)
            return newv
        except:
            raise Exception('Unable to convert ' + oldu + '->' + newu)

    #
    #  monthly we just need to determine how many days 
    #  in each month.
    #
    if intvl.lower()=='mn':
        newv = []
        y = first.year
        m = first.month
        for val in values:
            d = days_in_month(year=y, month=m)
            secs = d * 86400
            v2 = valueLinearToRate(value=val, oldu=oldu, 
                 newu=newu, area=area, secs=secs)
            newv.append(v2)
            m += 1
            if m>12:
                m = 1
                y = y + 1
        return newv
    
#-------------------------------------------------------
#  values = list of data values to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm', 'in' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  area  = effective area in square meters
#  intvl = interval of the data ('dy', 'wk', 'qm', 'mn')
#  first = start date (datetime.date)
#  last  = end date (datetime.date)
#-------------------------------------------------------
def rateToLinear(values=None, oldu=None, newu=None, area=None, 
                 intvl=None, first=None, last=None):
    if not values: return None
    if not oldu:   return None
    if not newu:   return None
    if not area:   return None
    if not intvl:  return None
    if not first:  return None
    if not last:   return None

    #
    #  daily data is relatively clean/easy because each value
    #  can use the same number of seconds for conversion.
    #  Thus, we can process it using list comprehension.
    #
    if intvl.lower()=='dy':
        secs = 86400
        try:
            vcms = rateConvert(values, oldu, 'cms')
            vm = [MISSING_REAL if v<MISSING_TEST else v*secs/area for v in vcms]
            newv = linearConvert(vm, 'm', newu)
            return newv
        except:
            raise Exception('Unable to convert ' + oldu + '->' + newu)

    #
    #  monthly we just need to determine how many days 
    #  in each month.
    #
    if intvl.lower()=='mn':
        newv = []
        y = first.year
        m = first.month
        for val in values:
            d = days_in_month(year=y, month=m)
            secs = d * 86400
            v2 = valueRateToLinear(value=val, oldu=oldu, 
                 newu=newu, area=area, secs=secs)
            newv.append(v2)
            m += 1
            if m>12:
                m = 1
                y = y + 1
        return newv
    
#-------------------------------------------------------
#  values = list of data values to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm3', 'in3' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  intvl = interval of the data ('dy', 'wk', 'qm', 'mn')
#  first = start date (datetime.date)
#  last  = end date (datetime.date)
#-------------------------------------------------------
def cubicToRate(values=None, oldu=None, newu=None, 
                 intvl=None, first=None, last=None):
    if not values: return None
    if not oldu:   return None
    if not newu:   return None
    if not intvl:  return None
    if not first:  return None
    if not last:   return None

    #
    #  daily data is relatively clean/easy because each value
    #  can use the same number of seconds for conversion.
    #  Thus, we can process it using list comprehension.
    #
    if intvl.lower()=='dy':
        secs = 86400
        try:
            m3 = cubicConvert(values, oldu, 'm3')
            vcms = [MISSING_REAL if v<MISSING_TEST else v/secs for v in m3]
            newv = rateConvert(vcms, 'cms', newu)
            return newv
        except:
            raise Exception('Unable to convert ' + oldu + '->' + newu)

    #
    #  monthly we just need to determine how many days 
    #  in each month.
    #
    if intvl.lower()=='mn':
        newv = []
        y = first.year
        m = first.month
        for val in values:
            d = days_in_month(yr=y, mo=m)
            secs = d * 86400
            v2 = valueCubicToRate(value=val, oldu=oldu, 
                 newu=newu, secs=secs)
            newv.append(v2)
            m += 1
            if m>12:
                m = 1
                y = y + 1
        return newv
    
#-------------------------------------------------------
#  values = list of data values to be converted
#            Any value < MISSING_TEST is considered "missing"
#  oldu  = unit string for incoming data (e.g. 'cm3', 'in3' )
#  newu  = unit string for outgoing data (e.g. 'cms', 'tcfs')
#  intvl = interval of the data ('dy', 'wk', 'qm', 'mn')
#  first = start date (datetime.date)
#  last  = end date (datetime.date)
#-------------------------------------------------------
def rateToCubic(values=None, oldu=None, newu=None, 
                intvl=None, first=None, last=None):
    if not values: return None
    if not oldu:   return None
    if not newu:   return None
    if not intvl:  return None
    if not first:  return None
    if not last:   return None

    #
    #  daily data is relatively clean/easy because each value
    #  can use the same number of seconds for conversion.
    #  Thus, we can process it using list comprehension.
    #
    if intvl.lower()=='dy':
        secs = 86400
        try:
            vcms = rateConvert(values, oldu, 'cms')
            vm3 = [MISSING_REAL if v<MISSING_TEST else v*secs for v in vcms]
            newv = cubicConvert(vm3, 'm3', newu)
            return newv
        except:
            raise Exception('Unable to convert ' + oldu + '->' + newu)

    #
    #  monthly we just need to determine how many days 
    #  in each month.
    #
    if intvl.lower()=='mn':
        newv = []
        y = first.year
        m = first.month
        for val in values:
            d = days_in_month(yr=y, mo=m)
            secs = d * 86400
            v2 = valueRateToCubic(value=val, oldu=oldu, 
                 newu=newu, secs=secs)
            newv.append(v2)
            m += 1
            if m>12:
                m = 1
                y = y + 1
        return newv
    

#--------------------------------------------------------------------
#  oldunits, newunits must be specified as strings, and must have a matching
#  entry in the tuples defined at the top.
#--------------------------------------------------------------------
def convertValues(values=None, oldunits=None, newunits=None, 
                  area=None, intvl=None, first=None, last=None):
    
    if not values:   return None
    if not oldunits: return None
    if not newunits: return None
    
    if not isinstance(oldunits, str):
        raise Exception('Invalid oldunits specification in convertValues.')
    if not isinstance(newunits, str):
        raise Exception('Invalid newunits specification in convertValues.')

    #
    #  If the conversion is within the same kind of units, we can
    #  do it directly.
    #
    try:
        if (oldunits in linear_units) and (newunits in linear_units):
            return linearConvert(values, oldunits, newunits)
        elif (oldunits in areal_units) and (newunits in areal_units):
            return arealConvert(values, oldunits, newunits)
        elif (oldunits in cubic_units) and (newunits in cubic_units):
            return cubicConvert(values, oldunits, newunits)
        elif (oldunits in rate_units) and (newunits in rate_units):
            return rateConvert(values, oldunits, newunits)
    except:
        raise Exception('Error converting ' + oldunits
                  + ' to ' + newunits)
    
    #
    #  If the conversion request is cross-group (e.g. cm -> cms)
    #  then we need to do it in two steps.
    #  Note that this is only valid to/from rate units.  Conversion
    #  from, for example, cm -> ft3 makes no sense.
    #     linear  ->  meters
    #     areal   ->  invalid
    #                 area <> rate doesn't work
    #     cubic   ->  cubic meters
    #     rate    ->  cubic meters per second
    #
    #  We require the number of seconds in all cases.
    #  If we are converting linear <-> rate, we also need the area.
    #  If we are doing cubic <-> rate, we do not need area.
    #
    #  If doing daily or weekly we could process every value with
    #  the same conversion factors, but in the general case each 
    #  value must be computed independently because the number
    #  of seconds will vary (e.g. February != June).  That will
    #  be handled in the routines called from this section of code.
    #
    #  We require the start/end dates for data intervals greater than
    #  weekly, because we have to compute the number of days for each
    #  value.
    #
    if not intvl: 
        print('Error: A unit conversion operation needed the time interval, but none was given.')
        raise Exception('Error converting ' + oldunits + ' to ' + newunits)
                  
    if (intvl == 'qm') or (intvl == 'mn') or (intvl == 'yr'):
        if not first: 
            print('Error: Unit conversion needed start date, but none was given.')
            raise Exception('Error converting ' + oldunits + ' to ' + newunits)
        if not last:  
            print('Error: Unit conversion needed end date, but none was given.')
            raise Exception('Error converting ' + oldunits + ' to ' + newunits)

    try:
        if (oldunits in linear_units) and (newunits in rate_units):
            if not area:
                print('Error: A unit conversion operation needed the area, but none was given.')
                raise Exception('Error converting ' + oldunits + ' to ' + newunits)
            return linearToRate(values=values, oldu=oldunits, newu=newunits, 
                   area=area, intvl=intvl, first=first, last=last)
        elif (oldunits in rate_units) and (newunits in linear_units):
            if not area:
                print('Error: A unit conversion operation needed the area, but none was given.')
                raise Exception('Error converting ' + oldunits + ' to ' + newunits)
            return rateToLinear(values=values, oldu=oldunits, newu=newunits,
                   area=area, intvl=intvl, first=first, last=last)
        elif (oldunits in cubic_units) and (newunits in rate_units):
            return cubicToRate(values=values, oldu=oldunits, newu=newunits, 
                   intvl=intvl, first=first, last=last)
        elif (oldunits in rate_units) and (newunits in cubic_units):
            return rateToCubic(values=values, oldu=oldunits, newu=newunits, 
                   intvl=intvl, first=first, last=last)
    except:
        raise Exception('Error converting ' + oldunits + ' to ' + newunits)
         
    raise Exception('Invalid conversion specified; ' + oldunits + ' to ' + newunits)
         
