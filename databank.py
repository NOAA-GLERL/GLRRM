#/bin/python

import sys
from copy import copy, deepcopy
import datetime
import databank_util as util

#--------------------------------------------------------------------
#  Method to get the primary name for a particular metadata
#  item when we do NOT already have an instantiated metadata
#  object. 
#  The caller must specify the type of metadata and a string
#  that they think will be a match.  This routine will try to
#  create a metadata object for that type of metadata and then
#  map the string that was passed in to the set of valid
#  names.  If that is all successful, the primary name for that
#  kind of metadata is returned.
#  For example:
#    (meta='kind',     name='runoff')      ->  'run'
#    (meta='units',    name='centimeters') ->  'cm'
#    (meta='location', name='ogoki')       ->  'og'
#    (meta='kind',     name='length')      ->  'na'   (no match)
#    (meta='kind',     name='meters')      ->  'na'   (no match)
#    (meta='interval', name='sup')         ->  'na'   (no match)
#--------------------------------------------------------------------
def getPrimaryName(meta=None, name=None):
    if not meta:
        raise Exception('Missing "meta=" for getPrimaryName()')
    if not name:
        raise Exception('Missing "name=" for getPrimaryName()')
        
    if not isinstance(meta, str):
        raise Exception('Invalid "meta=" for getPrimaryName()')
    if not isinstance(name, str):
        raise Exception('Invalid "name=" for getPrimaryName()')
        
    obj = None
    try:
        if meta.lower() == 'kind':
            obj = DataKind(name)
        if meta.lower() == 'units':
            obj = DataUnits(name)
        if meta.lower() == 'interval':
            obj = DataInterval(name)
        if meta.lower() == 'location':  
            obj = DataLocation(name)

        if obj:
            return obj.primaryName()
        else:
            return 'na'
    except:
        return 'na'


#-----------------------------------------------------------------------------
#  Define the classes used to store metadata "things" like kind, units, etc.
#--------------------------------------------------------------------------
#  BaseMeta is the base class and all of the actual functionality is
#  provided by methods of BaseMeta, but the subclasses provide the
#  specific character strings for the metadata items.  Two tuples are
#  defined for each kind of metadata. The _inputStrings tuple provides
#  various lookup strings that can be searched for a match in order to
#  identify specifically which item is being referenced.  The _outputStrings
#  tuple is intended for use when writing output to a file.
#
#  When editing these 2 tuples, remember that they are a matched set.
#  The rows must match.  e.g. Whatever quantity is referenced in row 3 of
#  _inputStrings, that same quantity must be referenced by row 3 of the
#  _outputStrings.
#
#  For the _inputStrings tuple, the very first entry should be a short
#  (2 or 3 characters) reference name. This is the name that will be used 
#  when building dictionary lookup keys. It will also be used as the
#  assigned metadata value for DataSeries objects. 
#  The rest of the tuple is really just freeform. There is no prescribed 
#  minimum or maximum number of entries.  I suggest that the set of input 
#  strings be fairly large and complete in order to accomodate all of the
#  reasonable possibilities.  Keep in mind that all comparisons will be
#  done in lowercase. So please keep all input strings 100% lowercase.
#
#  Conversely, the _outputStrings tuple should be kept minimal. For now 
#  I am defining only a short and long entry.  The entries should be 
#  specified with the desired case, because they will be used exactly
#  as specified.  Make sure that an equivalent for each output string is
#  included in the input string tuple so that if/when an output file is
#  later used as an input file, the name read from the file will be a
#  match when comparing input strings.
#---------------------------------------------------------------------------
#
#  BaseMeta should never be used directly. Please always use one of the
#  derived subclasses.
#-----------------------------------------------------------------------------
class BaseMeta():
    def __init__(self, initval):
        if initval:
            if isinstance(initval, str):
                self.myValue = self.intValueFromString(initval)
            else:
                self.myValue = 0
        else:
            self.myValue = 0

    #---------------------------------------------------
    def className(self):
        return self.__class__.__name__

    #---------------------------------------------------
    #  Method to get the primary name for a particular entry
    #  in a subclass when we have an instantiated object.  
    #  The primary name is the first entry in that 
    #  value's _inputStrings tuple.
    #  If unable to find any match, returns 'na'.
    #
    def primaryName(self):
        try:
            return self.inputName(0)
        except:
            return 'na'

    #---------------------------------------------------
    #  For inputName() and outputName(), first collapse the index
    #  into the valid range, then access the tuple entry.
    #  If called with no column specified, the default is to
    #  return the first entry in the tuple (index 0).  If a
    #  value is given for column, then that will be the index
    #  used.  
    #  For example, if self is an object of type DataUnits
    #  and self.myValue is 3 (i.e. the units are meters), then
    #   inputName() will return "m", then first entry in that tuple.
    #   inputName(0) will return "m"
    #   inputName(1) will return "meter"
    #   inputName(2) will return "meters"
    #   inputName(3) [or greater] will return "meters"
    #
    def inputName(self, column=0):
        try:
            the_tuple = type(self)._inputStrings[self.myValue]
            maxcol = len(the_tuple) - 1
            c = max(0, min(maxcol, column))
            return the_tuple[c]
        except:
            return 'n/a'

    #---------------------------------------------------
    def outputName(self, column):
        try:
            the_tuple = type(self)._outputStrings[self.myValue]
            maxcol = len(the_tuple) - 1
            c = max(0, min(maxcol, column))
            return the_tuple[c]
        except:
            return 'n/a'

    #---------------------------------------------------
    def outputNameShort(self):
        try:
            the_tuple = type(self)._outputStrings[self.myValue]
            return the_tuple[0]
        except:
            return 'n/a'

    #---------------------------------------------------
    def outputNameLong(self):
        try:
            the_tuple = type(self)._outputStrings[self.myValue]
            maxcol = len(the_tuple) - 1
            return the_tuple[maxcol]
        except:
            return 'n/a'
            
    #---------------------------------------------------
    #  Method to convert a string into the index for that
    #  metadata value by finding the matching entry in _inputStrings.
    #  If unable to find any match, returns 0 (undefined).
    #
    #  For example, if self is an object of type DataUnits, and this
    #  routine is passed the string "meters", this will return the
    #  index value 3, because the set of names for "meters" is the
    #  third tuple in the _inputStrings list.
    #
    def intValueFromString(self, string):
        j = -1
        s = string.lower().strip()
        for i in range(len(self._inputStrings)):
            if j < 0:
                try:
                    for k in range(len(self._inputStrings[i])):
                        t = self._inputStrings[i][k].lower().strip()
                        if t == s:
                            return i
                except:
                    pass
        return 0

#--------------------------------------------------------------------------------
"""The DataKind class contains the defined data kinds (a.k.a types) and their 
   associated string literals.
"""
class DataKind(BaseMeta):
    _inputStrings = (
        ('na', 'n/a', 'undef', 'undefined'),
        ('prc', 'prec',    'precip',    'precipitation'),
        ('run', 'runf',    'runoff'),
        ('evp', 'evap',    'evaporation'),
        ('nbs', 'net_basin_supply', 'net basin supply'),
        ('mlv', 'meanlev', 'mean_level', 'mean level', 'mean water level'),
        ('blv', 'bomlev',  'beginning_of_month_level', 
                'beginning of month level', 'beginning-of-month water level'),
        ('elv', 'eomlev',  'end_of_month_level',  'end of month level',
                'end-of-month water level'),
        ('flw', 'flow', 'channel flow'),
        ('con', 'cons',    'consumptive_use', 'consuse', 'consumptive use'),
        ('icw', 'iceweed', 'ice_weed_retardation', 'ice & weed retardation')
    )

    _outputStrings = (
        ('undef',      'undefined'),
        ('Precip',     'Precipitation'),
        ('Runoff',     'Runoff'),
        ('Evap',       'Evaporation'),
        ('NBS',        'Net Basin Supply'),
        ('MeanLev',    'Mean Water Level'),
        ('BOMLev',     'Beginning-Of-Month Water Level'),
        ('EOMLev',     'End-Of-Month Water Level'),
        ('Flow',       'Channel Flow'),
        ('ConsUse',    'Consumptive Use'),
        ('IceWeed',    'Ice & Weed Retardation')
    )


#--------------------------------------------------------------------------------
"""The DataUnits class contains the defined data units and their associated string literals.
"""
class DataUnits(BaseMeta):
    _inputStrings = (
        ('na', 'n/a', 'undef', 'undefined'),
        ('mm',  'millimeter',  'millimeters'),
        ('cm',  'centimeter',  'centimeters'),
        ('m',   'meter',       'meters'),
        ('km',  'kilometer',   'kilometers'),
        ('in',  'inch',        'inches'),
        ('ft',  'foot',        'feet'),
        ('yd',  'yard',        'yards'),
        ('mi',  'mile',        'miles'),
        ('mm2', 'sqmm',        'sq_mm',    'square_mm',  
                'square_millimeters', 'square millimeters'),
        ('cm2', 'sqcm',        'sq_cm',    'square_cm',  
                'square_centimeters', 'square centimeters'),
        ('m2',  'sqm',         'sq_m',     'square_m',   
                'square_meters', 'square meters'),
        ('km2', 'sqkm',        'sq_km',    'square_km',  
                'square_kilometers', 'square kilometers'),
        ('in2', 'sqin',        'sq_in',    'square_in',  
                'square_inches', 'square inches'),
        ('ft2', 'sqft',        'sq_ft',    'square_ft',  
                'square_feet', 'square feet'),
        ('yd2', 'sqyd',        'sq_yd',    'square_yd',  
                'square_yards', 'square yards'),
        ('mi2', 'sqmi',        'sq_mi',    'square_mi',  
                'square_miles', 'square miles'),
        ('mm3',                'cu_mm',    'cubic_mm',   
                'cubic_millimeters', 'cubic millimeters'),
        ('cm3',                'cu_cm',    'cubic_cm',
                'cubic_centimeters', 'cubic centimeters'),
        ('m3',                 'cu_m',     'cubic_m',
                'cubic_meters', 'cubic meters'),
        ('km3',                'cu_km',    'cubic_km',
                'cubic_kilometers', 'cubic kilometers'),
        ('in3',                'cu_in',    'cubic_in',
                'cubic_inches', 'cubic inches'),
        ('ft3',                'cu_ft',    'cubic_ft', 
                'cubic_feet', 'cubic feet'),
        ('yd3',                'cu_yd',    'cubic_yd',
                'cubic_yards', 'cubic yards'),
        ('mi3',                'cu_mi',    'cubic_mi',
                'cubic_miles', 'cubic miles'),
        ('cms', 'm3s',         'cu_ms',    'cubic_ms',
                'cubic_meters_per_second', 'cubic meters per second'),
        ('10cms', '10m3s',  '10m3/s', '10*cms',    'tens_cms',
                'tens_of_cubic_meters_per_second', 
                'ten cubic meters per second'),
        ('cfs', 'ft3s',        'cu_fts',   'cubic_fts',  
                'cubic_feet_per_second', 'cubic feet per second'),
        ('tcfs',               'thousand_cfs',
                'thousand_cubic_feet_per_second',
                'thousand cubic feet per second'),
    )

    _outputStrings = (
        ('undef', 'undefined'),
        ('mm',    'Millimeters'),
        ('cm',    'Centimeters'),
        ('m',     'Meters'),
        ('km',    'Kilometers'),
        ('in',    'Inches'),
        ('ft',    'Feet'),
        ('yd',    'Yards'),
        ('mi',    'Miles'),
        ('mm2',   'Square Millimeters'),
        ('cm2',   'Square Centimeters'),
        ('m2',    'Square Meters'),
        ('km2',   'Square Kilometers'),
        ('in2',   'Square Inches'),
        ('ft2',   'Square Feet'),
        ('yd2',   'Square Yards'),
        ('mi2',   'Square Miles'),
        ('mm3',   'Cubic Millimeters'),
        ('cm3',   'Cubic Centimeters'),
        ('m3',    'Cubic Meters'),
        ('km3',   'Cubic Kilometers'),
        ('in3',   'Cubic Inches'),
        ('ft3',   'Cubic Feet'),
        ('yd3',   'Cubic Yards'),
        ('mi3',   'Cubic Miles'),
        ('cms',   'Cubic Meters per Second'),
        ('10cms', 'Ten Cubic Meters per Second'),
        ('cfs',   'Cubic Feet per Second'),
        ('tcfs',  'Thousand Cubic Feet per Second'),
    )


#--------------------------------------------------------------------------------
"""The DataInterval class contains the defined time intervals and their associated string literals.
"""
class DataInterval(BaseMeta):
    _inputStrings = (
        ('na', 'n/a', 'undef', 'undefined'),
        ('dy',  'day',   'dly',   'daily'),
        ('wk',  'week',  'wkly',  'weekly'),
        ('qm',  'qmon',  'qtrm',  'qtrmon', 'qtrmonth', 'qtrmonthly', 
            'qtr-month', 'qtr-monthly', 'quartermon', 'quartermonth', 
            'quartermonthly', 'quarter-month', 'quarter-monthly'),
        ('mn',  'mon',   'month',  'monthly'),
        ('yr',  'year',  'yearly', 'annual')
    )

    _outputStrings = (
        ('undef', 'undefined'),
        ('dy',    'Daily'),
        ('wk',    'Weekly'),
        ('qmon',  'Qtr-monthly'),
        ('mn',    'Monthly'),
        ('yr',    'Annual')
    )


#--------------------------------------------------------------------------------
"""The DataLocation class contains the defined locations and their associated string literals.
"""
class DataLocation(BaseMeta):
    _inputStrings = (
        ('na', 'n/a', 'undef', 'undefined'),
        ('su', 'sup', 'lksup', 'lakesup', 'superior', 'lake superior'),
        ('mi', 'mic', 'lkmic', 'lakemic', 'michigan', 'lake michigan'),
        ('hu', 'hur', 'lkhur', 'lakehur', 'huron',    'lake huron'),
        ('sc', 'stc', 'lkstc', 'lakestc', 'stclair',  'lake stclair'),
        ('er', 'eri', 'lkeri', 'lakeeri', 'erie',     'lake erie'),
        ('on', 'ont', 'lkont', 'lakeont', 'ontario',  'lake ontario'),
        ('mh', 'mih', 'mhu', 'mhn', 'lakemhu', 'lakemhn', 'michhur',
               'michuron', 'mich-huron', 'lake mich-huron', 
               'michigan-huron', 'lake michigan-huron'),
        ('og',  'ogoki',   'odiv',      'ogokidiv',   'ogoki diversion'),
        ('ll',  'longlac', 'long lac',  'longlacdiv', 
                'longlac diversion', 'long lac diversion'),
        ('oll', 'olldiv',  'ogoki-longlac', 'ogoki-longlac diversion'),
        ('chi', 'chidiv',  'chicago',    'chicago diversion'),
        ('wel', 'weldiv',  'welland',    'welland diversion'),
        ('smr', 'stmriv',  'stmarys',    'stmarys river', 'st. marys river'),
        ('scr', 'stcriv',  'stcriver',   'stclair river', 'st. clair river'),
        ('det', 'detriv',  'detriver',   'detroit', 'detroit river'),
        ('nia', 'niariv',  'niariver',   'niagara', 'niagara river'),
        ('stl', 'stlriv',  'stlriver',   'stlawrence', 'stlawrence river',
                'st. lawrence river')
    )

    #
    #  *** NOTE ***
    #  Lake St. Clair and the St. Clair River introduces an issue for us.
    #  The name "St. Clair" could potentially refer to either one.
    #  So what do I do?  I have chosen to define StClair as one of the 
    #  output strings for both. That way it can be used, in context, for
    #  either one.  But I am leaving it as a valid input string for only
    #  the lake, because we cannot have duplications in that set.
    #  Users of this class will need to be cognizant of that situation.
    #
    _outputStrings = (
        ('undef', 'undefined'),
        ('su',    'Superior',       'Lake Superior'),
        ('mi',    'Michigan',       'Lake Michigan'),
        ('hu',    'Huron',          'Lake Huron'),
        ('sc',    'StClair',        'Lake StClair'),
        ('er',    'Erie',           'Lake Erie'),
        ('on',    'Ontario',        'Lake Ontario'),
        ('mh',    'Michigan-Huron', 'Lake Michigan-Huron'),
        ('og',    'Ogoki',          'Ogoki Diversion'),
        ('ll',    'LongLac',        'LongLac Diversion'),
        ('oll',   'Ogoki-LongLac',  'Ogoki-LongLac Diversion'),
        ('chi',   'Chicago',        'Chicago Diversion'),
        ('wel',   'Welland',        'Welland Diversion'),
        ('smr',   'StMarys',        'St. Marys River'),
        ('scr',   'StClair',        'St. Clair River'),
        ('det',   'Detroit',        'Detroit River'),
        ('nia',   'Niagara',        'Niagara River'),
        ('stl',   'StLawrence',     'St. Lawrence River')
    )

#--------------------------------------------------------------------------------
#  Define the dataseries class that stores a single timeseries of data along with
#  its metadata.
#--------------------------------------------------------------------------------
class DataSeries(object):
    """The dataseries class contains a single series of data.

    Metadata may be left undefined ('na') at the time one of these objects is 
    created. Obviously, though, once data is assigned, we will need the 
    metadata values to be valid. 
    
    The start and end dates are always defined as a particular DAY, even for
    the larger timesteps. For the larger timesteps, the start date is the
    first day of the first period and the enddate is the last day of the last
    period.  For example:

    Suppose we have weekly data for the seven "regulation weeks"
    (Friday-Thursday) of January 6, 2017 through February 23, 2017, then the
    start and end dates are exactly those days.

    If we have quarter-monthly data for the first 3 quarters of March, 2017,
    then startdate is 3/1/2017 and enddate is 3/23/2017.

    If we have a dataset with monthly data for April 2010 through November 2016,
    then startdate is 4/1/2010 and enddate is 11/30/2016.

    Annual start/end dates are always January 1 and December 31 of the respective
    years.
    """

    def __init__(self, kind=None, units=None, intvl=None, loc=None, set=None,
                       first=None, last=None, values=None):
        self.dataKind     = 'na'
        self.dataUnits    = 'na'
        self.dataInterval = 'na'
        self.dataLocation = 'na'
        self.dataSet      = 'na'
        self.startDate    = util.MISSING_DATE
        self.endDate      = util.MISSING_DATE
        self.dataVals     = []
                       
        #
        #  Handle metadata initialization
        #  They are all specified with text strings.
        #
        if kind:
            try:
                self.dataKind = getPrimaryName(meta='kind', name=kind)
            except:
                raise Exception('Invalid kind specifier in DataSeries init')
        else:
            self.dataKind = 'na'

        if units:
            try:
                self.dataUnits = getPrimaryName(meta='units', name=units)
            except:
                raise Exception('Invalid units specifier in DataSeries init')
        else:
            self.dataUnits = 'na'

        if intvl:
            try:
                self.dataInterval = getPrimaryName(meta='interval', name=intvl)
                if isinstance(intvl, str):
                    self.dataInterval = DataInterval(intvl).primaryName()
            except:
                raise Exception('Invalid interval specifier in DataSeries init')
        else:
            self.dataInterval = 'na'

        if loc:
            try:
                self.dataLocation = getPrimaryName(meta='location', name=loc)
                if isinstance(loc, str):
                    self.dataLocation = DataLocation(loc).primaryName()
            except:
                raise Exception('Invalid location specifier in DataSeries init')
        else:
            self.dataLocation = 'na'

        if set:
            try:
                self.dataSet = set
            except:
                raise Exception('Invalid set specifier in DataSeries init')
        else:
            self.dataSet = 'na'

        if first:
            self.startDate = util.date_from_entry(first)
        else:
            self.startDate = util.MISSING_DATE
            
        if last:
            self.endDate = util.date_from_entry(last)
        else:
            self.endDate = util.MISSING_DATE

        if values:
            self.dataVals = values
        else:
            self.dataVals = []
            
        #
        #  Future enhancement: compute the required number of values
        #  based on the specified start/end dates and interval.  Then
        #  compare that to the number of entries in the values list.
        #

    #---------------------------------------------------------------------
    def printSummary(self):
        print('Summary of DataSeries...')
        print(' kind  = ', DataKind(self.dataKind).outputNameLong())
        print(' units = ', DataUnits(self.dataUnits).outputNameLong())
        print(' intvl = ', DataInterval(self.dataInterval).outputNameLong())
        print(' loc   = ', DataLocation(self.dataLocation).outputNameLong())
        print(' set   = ', self.dataSet)
        print(' start date = ', str(self.startDate))
        print(' end date   = ', str(self.endDate))
        print(' data values = ', self.dataVals)

    #---------------------------------------------------------------------
    def printOneLineSummary(self):
        print(
            ' kind  = ', DataKind(self.dataKind).outputNameLong(),  ';',
            ' units = ', DataUnits(self.dataUnits).outputNameLong(), ';',
            ' intvl = ', DataInterval(self.dataInterval).outputNameLong(),  ';',
            ' loc = ',   DataLocation(self.dataLocation).outputNameLong(),  ';',
            ' set = ',   self.dataSet, ';',
            ' dates = ', str(self.startDate),
            ' to ',      str(self.endDate)
        )

    #---------------------------------------------------------------------
    def getOneLineSummary(self):
        return(
            'kind='  + DataKind(self.dataKind).outputNameLong() + ';',
            'units=' + DataUnits(self.dataUnits).outputNameLong() + ';',
            'intvl=' + DataInterval(self.dataInterval).outputNameLong() + ';',
            'loc='   + DataLocation(self.dataLocation).outputNameLong() + ';',
            'set='   + self.dataSet + ';',
            'dates=' + str(self.startDate) +
            ' to '   + str(self.endDate)
        )

    #---------------------------------------------------------------------
    def add_data(self, newData):
        """Add a continuous timeseries of data to the stored data.

        On success, it returns True.
        If there is a problem, it returns False.
        """
        #
        #  Compare the metadata values.  Reminder -- these are the
        #  2- or 3-character strings, all lowercase.
        #
        if newData.dataKind != self.dataKind:
            print('Error. Mismatched data kinds in DataSeries.add_data')
            raise TypeError('Invalid attempt to add data to a DataSeries')

        if newData.dataUnits != self.dataUnits:
            print('Error. Mismatched data units in DataSeries.add_data')
            raise TypeError('Invalid attempt to add data to a DataSeries')

        if newData.dataInterval != self.dataInterval:
            print('Error. Mismatched intervals in DataSeries.add_data')
            raise TypeError('Invalid attempt to add data to a DataSeries')

        if newData.dataLocation != self.dataLocation:
            print('Error. Mismatched locations in DataSeries.add_data')
            raise TypeError('Invalid attempt to add data to a DataSeries')

        if newData.dataSet != self.dataSet:
            print('Error. Mismatched set designators in DataSeries.add_data')
            raise TypeError('Invalid attempt to add data to a DataSeries')

        #
        #  Compare dates
        #
        if newData.startDate == util.MISSING_DATE:
            print("Missing start date specification in call to add_data().")
            return False

        if newData.endDate == util.MISSING_DATE:
            print("Missing end date specification in call to add_data().")
            return False

        #
        #  Call the appropriate merge routine
        #
        if newData.dataInterval == getPrimaryName(meta='interval', name='dly'):
            try:
                self.mrg_daily_data(newData)
            except:
                raise Exception('Error attempting to merge daily data.')
        elif newData.dataInterval == getPrimaryName(meta='interval', name='monthly'):
            try:
                self.mrg_monthly_data(newData)
            except:
                raise Exception('Error attempting to merge monthly data.')
        else:
            print('Unable to merge data because interval is invalid.')
            return False

        return True

    #---------------------------------------------------------------------
    def mrg_daily_data(self, newData):
        """Merge an update set of continuous daily data to the stored data.
        NOT intended to be called by user code, but rather for use by other
        class methods.
        It overwrites any old values in the period of the new
        data, but preserves existing values outside that period.
        The dates and data list are modified in place.

        On success, it returns True.
        If there is a problem, it returns False.

        Note that metadata (kind, units, interval) are NOT verified here.
        They are assumed to have been verified prior to calling this function.
        """
        if newData.startDate == util.MISSING_DATE:
            print("Missing start date specification in call to mrg_daily_data().")
            return False

        if newData.endDate == util.MISSING_DATE:
            print("Missing end date specification in call to merge_daily_data().")
            return False

        #
        #  Determine the new date extents.
        #  Just save these values for now
        #
        mrgStart = min(self.startDate, newData.startDate)
        mrgEnd   = max(self.endDate,   newData.endDate)

        #
        #  create a template list of the new size, filled with missing data values
        #
        ndays = mrgEnd.toordinal() - mrgStart.toordinal() + 1
        mrgData = [util.MISSING_REAL] * ndays

        #
        #  replace the appropriate slice of mrgData with the old data
        #  Don't forget that the end index is always 1 past the desired end
        #
        i = self.startDate.toordinal() - mrgStart.toordinal()
        j = self.endDate.toordinal() - mrgStart.toordinal() + 1
        m = 0
        n = self.endDate.toordinal() - self.startDate.toordinal() + 1
        mrgData[i:j] = self.dataVals[m:n]

        #
        #  replace the appropriate slice of mrgData with the new data
        #  Don't forget that the end index is always 1 past the desired end
        #
        i = newData.startDate.toordinal() - mrgStart.toordinal()
        j = newData.endDate.toordinal() - mrgStart.toordinal() + 1
        m = 0
        n = newData.endDate.toordinal() - newData.startDate.toordinal() + 1
        mrgData[i:j] = newData.dataVals[m:n]

        #
        #  replace/overwrite old object values with the new ones
        #
        self.startDate = mrgStart
        self.endDate = mrgEnd
        self.dataVals = mrgData
        return True

    #---------------------------------------------------------------------
    def mrg_monthly_data(self, newData):
        """Merge an update set of continuous monthly data to the stored data.
        NOT intended to be called by user code, but rather for use by other
        class methods.
        It overwrites any old values in the period of the new
        data, but preserves existing values outside that period.
        The dates and data list are modified in place.

        On success, it returns True.
        If there is a problem, it returns False.

        Note that metadata (kind, units, interval) are NOT verified here.
        They are assumed to have been verified prior to calling this function.
        """
        if newData.startDate == util.MISSING_DATE:
            print("Missing start date specification in call to add_data().")
            return False

        if newData.endDate == util.MISSING_DATE:
            print("Missing end date specification in call to add_data().")
            return False

        if newData.startDate.day != 1:
            print("Invalid start date for monthly data.  Must be 1.")
            return False

        d = last_day_of_month(newData.endDate)
        if newData.day != d.day:
            print("Invalid end date for monthly data.  Must be last day of the month.")
            return False

        #
        #  Determine the new date extents.
        #  Just save these values for now
        #
        mrgStart = min(self.startDate, newData.startDate)
        mrgEnd = max(self.endDate,   newData.endDate)

        #
        #  create a template list of the new size, filled with missing data values
        #
        mons = ((mrgEnd.year - mrgStart.year)*12
                + (mrgEnd.month - mrgStart.month) + 1)
        mrgData = [util.MISSING_REAL] * mons

        #
        #  replace the appropriate slice of mrgData with the old data
        #  Don't forget that the end index is always 1 past the desired end
        #
        sd = self.startDate
        ed = self.endDate
        i = (sd.year - mrgStart.year)*12 + (sd.month - mrgStart.month)
        j = (ed.year - mrgEnd.year)*12 + (ed.month - mrgEnd.month) + 1
        m = 0
        n = (ed.year - sd.year)*12 + (ed.month - sd.month) + 1
        mrgData[i:j] = self.dataVals[m:n]

        #
        #  replace the appropriate slice of mrgData with the new data
        #  Don't forget that the end index is always 1 past the desired end
        #
        sd = newData.startDate
        ed = newData.endDate
        i = (sd.year - mrgStart.year)*12 + (sd.month - mrgStart.month)
        j = (ed.year - mrgEnd.year)*12 + (ed.month - mrgEnd.month) + 1
        m = 0
        n = (ed.year - sd.year)*12 + (ed.month - sd.month) + 1
        mrgData[i:j] = newData.dataVals[m:n]

        #
        #  replace/overwrite old object values with the new ones
        #
        self.startdate = mrgStart
        self.enddate = mrgEnd
        self.dataVals = mrgData

        return True

#--------------------------------------------------------------------------------
#  Define the DataVault class that stores a bunch of DataSeries objects and will
#  be used as the repository for GLRRM data.
#
#  When a dataset is stored into the vault, the data values will be "normalized" 
#  to a specific set of units. When data is withdrawn, the values will be 
#  converted into whatever was requested. When this conversion requires a lake
#  area to be used (e.g. when converting from cms -> mm) we will use the 
#  coordinated lake areas by default. The caller will have the option to specify
#  an alternate area.
#  *** MODELERS BEWARE ***
#  If you use this functionality to specify an alternative lake area, please be
#  certain that you are consistent, and ALWAYS use it.  As you can imagine, if
#  you were to store values converted via coordinated areas (no area specified), 
#  then withdraw them using an alternative area, then store them again using 
#  the coordinated area (no area specified) -- You will end up changing the 
#  data values in the vault in a way you don't intend.
#  
#  Normalized units are:
#    precipitation      cubic meters per second
#    runoff             cubic meters per second
#    evaporation        cubic meters per second
#    net basin supply   cubic meters per second
#    mean level         meters
#    bom level          meters
#    eom level          meters
#    flow               cubic meters per second
#    consumptive use    cubic meters per second
#    ice-weed ret       cubic meters per second
#
#  Note that lake levels are stored using a generic elevation without any
#  consideration of IGLD55, IGLD85, etc.  For now, it is the responsibility
#  of the model to account for any IGLD conversion issues.
#--------------------------------------------------------------------------------
class DataVault(object):

    #--------------------------------
    #  Specify the normalized units to use for each kind
    #--------------------------------
    _normalizedUnits = (
        ('prc', 'cms'),
        ('run', 'cms'),
        ('evp', 'cms'),
        ('nbs', 'cms'),
        ('mlv', 'm'),
        ('blv', 'm'),
        ('elv', 'm'),
        ('flw', 'cms'),
        ('con', 'cms'),
        ('icw', 'cms')
    )
    def getNormalizedUnits(self, kind=None):
        if kind:
            if isinstance(kind, str):
                s = kind.lower()
            else:
                return 'na'
        else:
            return 'na'
        try:
            for t in self._normalizedUnits:
                if t[0].lower() == s:
                    return t[1]
            return 'na'
        except:
            return 'na'
    
    
    #--------------------------------------------------------
    def __init__(self):
        self.vault = {}               # the dictionary object

    #-------------------------------------------------------------------
    #  Construct a lookup key for our dictionary from EITHER:
    #    1) The metadata in a DataSeries object, if ds is provided.
    #       OR
    #    2) The specific kind, interval, location specified.
    #  If ds is given, kind, intvl, loc, set will be ignored.
    #
    #  Returns a text string that looks like this:
    #    kind_intvl_loc_set
    #  where:
    #     kind  = the character string defined as the first entry in the
    #             _inputStrings tuple for this data kind. e.g. For runoff data,
    #             the _inputStrings tuple is ('run', 'runf', 'runoff').  We will
    #             use 'run' as the kind string.  Typically, this will be a
    #             3-character string, but it is not required. We will just use
    #             that very first entry in the tuple. Keep in mind that this
    #             implicitly assumes the entry is unique, which is required
    #             anyway.
    #     intvl = the character string defined as the first entry in the
    #             _inputStrings tuple for this data interval. Same
    #             story as the data kind.
    #     loc   = the character string defined as the first entry in the
    #             _inputStrings tuple for this location. Same story as kind/intvl.
    #     set   = character string name of the data "set" associated with this data.
    #             e.g. "fcst1", "fcst2", "fcst3", etc.  Length is unrestricted, but
    #             try not to go too crazy, and no imbedded spaces.
    #             If set is not specified, a default value of "na" will be used 
    #             in the key string.
    #
    #  For example, if the function call looks like:
    #     myVault._construct_vault_key(kind='runoff', intvl='daily', loc='ont')
    #  then the kind/intvl/loc values will turn into 'run'/'dy'/'on', and the
    #  resulting lookup key will be 'run_dy_on_na'.
    #
    #  If the function call looks like:
    #     myVault._construct_vault_key(kind='runoff', intvl='daily', loc='ont', set='abc')
    #  then the resulting lookup key will be 'run_dy_on_abc'.
    #
    #-------------------------------------------------------------------
    @classmethod
    def _construct_vault_key(thisclass, ds=None, kind=None, 
                             intvl=None, loc=None, set=None):
        if ds:
            kstr = ds.dataKind
            istr = ds.dataInterval
            lstr = ds.dataLocation
            sstr = ds.dataSet
            key  = kstr + '_' + istr + '_' + lstr + '_' + sstr
            return key

        #
        #  If all 3 items are specified, e.g.
        #     _construct_vault_key(kind='nbs', intvl='daily', loc='erie')
        #  First create temporary metadata objects for each, then
        #  get the correct lookup name for each.
        #
        if kind and intvl and loc:
            kobj = DataKind(kind)
            iobj = DataInterval(intvl)
            lobj = DataLocation(loc)
            
            kstr = kobj.primaryName()
            istr = iobj.primaryName()
            lstr = lobj.primaryName()
            key = kstr + '_' + istr + '_' + lstr
            if set:
                key = key + '_' + set
            else:
                key = key + '_na'
            return key
        else:
            raise ValueError('Missing DataSeries information in _construct_vault_key')
            


    #-------------------------------------------------------------------
    def printVault(self):
        for key in self.vault:
            print('key=', key, ':', self.vault[key].getOneLineSummary())

    #-------------------------------------------------------------------
    #  The deposit() function is how a user adds data to the vault.
    #  series is a DataSeries object.
    #-------------------------------------------------------------------
    def deposit(self, series, lake_area=None):
        try:
            key = type(self)._construct_vault_key(series)
        except:
            raise Exception('databank.deposit: error getting the key')

        #
        #  If user did not specify a lake area, then assign a value (if needed)
        #
        if not lake_area:
            lake_area = util.getLakeArea(series.dataLocation)
            
        #
        #  Create a "normalized" DataSeries object such that the data
        #  units and values conform to the prescribed data units
        #  for storage in the vault.
        #
        tempvals = copy(series.dataVals)      # default is to use data as-is
        normstr = self.getNormalizedUnits(kind=series.dataKind)
        try:
            #
            #  If needed, convert data units.
            #
            if series.dataUnits != normstr:
                try:
                    tempvals = None
                    if series.dataUnits != normstr:
                        oldstr = series.dataUnits
                        oldvals = copy(series.dataVals)
                        
                        if (oldstr in util.linear_units) and normstr=='m': 
                            tempvals = util.convertValues(values=oldvals, 
                                    oldunits=oldstr, newunits=normstr, intvl=series.dataInterval)
                        elif (oldstr in util.rate_units) and normstr=='cms': 
                            tempvals = util.convertValues(values=oldvals, 
                                    oldunits=oldstr, newunits=normstr, intvl=series.dataInterval)
                        elif (oldstr in util.areal_units):
                            raise Exception('Error: datavault unable to store '
                                          + 'areal datasets.')
                        elif normstr=='m':
                            tempvals = util.convertValues(values=oldvals, 
                                    oldunits=oldstr, newunits=normstr,
                                    area=lake_area, first=series.startDate, 
                                    last=series.endDate, intvl=series.dataInterval)
                        elif normstr=='cms':
                            tempvals = util.convertValues(values=oldvals, 
                                    oldunits=oldstr, newunits=normstr,
                                    area=lake_area, first=series.startDate, 
                                    last=series.endDate, intvl=series.dataInterval)
                        else:
                            print('series.dataUnits=', series.dataUnits)
                            print('normstr=', normstr)
                            raise Exception('Unhandled data units conversion.')
                except:
                    raise Exception('Unable to do required data conversion.')
        except:
            raise Exception('Unable to create dataset for the datavault.')

        #
        #  Create a temporary dataset that contains the data to be added, 
        #  in the correct units.
        #
        tds = DataSeries(kind=series.dataKind, units=normstr, loc=series.dataLocation,
                    intvl=series.dataInterval, set=series.dataSet, 
                    first=series.startDate, last=series.endDate, values=tempvals)
            
        #
        #  Do we already have a data series like this?
        #  If so, we will merge them.
        #
        try:
            #
            #  Get the old data set.  Note that this line will fail
            #  with an exception if there is no matching old dataset.
            #  That will cause us to drop into the except block and skip the
            #  rest of the try block.  Seems a little ugly to me, but I
            #  guess that is the "python way" to do it?  I don't see any
            #  search functionality that works "cleanly".
            #
            old = self.vault[key]

            #
            #  Getting to here means we have an old dataset. Proceed with
            #  checking and merging.
            #
            #  Verify that the data sets have matching metadata.
            #  This should actually never be an issue, but verifying is good.
            #
            if old.dataKind != new.dataKind:
                raise ValueError('Data kind mismatch')
            if old.dataInterval != new.dataInterval:
                raise ValueError('Data interval mismatch')
            if old.dataLocation != new.dataLocation:
                raise ValueError('Data location mismatch')

            #
            #  Merge the two DataSeries objects
            #
            try:
                old.add_data(tds)
            except:
                raise('Error merging the new data into the old.')

        except:
            #
            #  No old dataset, so just add this new one to the vault.
            #
            try:
                self.vault.update({key:tds})
            except:
                print('Error adding new dataset to vault')

    #---------------------------------------------------------------
    #  Equivalent to deposit, but with all fields individually specified.
    #---------------------------------------------------------------
    def deposit_data(self, kind=None, units=None, intvl=None, loc=None, set=None,
                     first=None, last=None, values=None, lake_area=None):
        try:
            dk = DataKind(kind).primaryName()
            du = DataUnits(units).primaryName()
            di = DataInterval(intvl).primaryName()
            dl = DataLocation(loc).primaryName()
            dt = set
            ds = DataSeries(kind=dk, units=du, intvl=di, loc=dl, set=dt,
                 first=first, last=last, values=values)
        except:
            raise Exception(
                'Error creating temporary DataSeries object in deposit_data()')
            return

        try:
            self.deposit(ds, lake_area)
        except:
            raise Exception('Error adding data series object in deposit_data()')
            return

            
    #----------------------------------------------------------------
    #  Caller must specify the kind, interval, location and units.
    #  first/last are optional.
    #  If invalid specifiers are given, returns with a exception.
    #  If all works correctly, it returns a DataSeries object, if the 
    #    data is in the vault.
    #  If specified ok, but data does not exist in the vault, then
    #    this also returns None, but no exception is generated.
    #----------------------------------------------------------------
    def withdraw(self, kind=None, units=None, intvl=None, loc=None, set=None,
                 first=None, last=None, lake_area=None):

        #
        #  Verify that all metadata strings were validly specified.
        #  Do a full translation from string to object and back in
        #  order to verify a valid specification.
        #
        dk = 'na'
        if kind:
            try:
                dk = DataKind(kind).primaryName()
            except:
                dk = 'na'
        if dk=='na':
            raise Exception('Invalid or missing kind specification '
                           + 'to DataVault.withdraw()')

        du = 'na'
        if units:
            try:
                du = DataUnits(units).primaryName()
            except:
                du = 'na'
        if du=='na':
            raise Exception('Invalid or missing units specification '
                           + 'to DataVault.withdraw()')

        di = 'na'
        if intvl:
            try:
                di = DataInterval(intvl).primaryName()
            except:
                di = 'na'
        if di=='na':
            raise Exception('Invalid or missing interval specification '
                           + 'to DataVault.withdraw()')

        dl = 'na'
        if loc:
            try:
                dl = DataLocation(loc).primaryName()
            except:
                dl = 'na'
        if dl=='na':
            raise Exception('Invalid or missing location specification '
                           + 'to DataVault.withdraw()')

        dt = 'na'
        if set:
            dt = set

        #
        #  Construct the vault key from the lookup strings
        #
        key = type(self)._construct_vault_key(kind=dk, intvl=di, loc=dl, set=dt)
        
        #
        #  Get a temporary dataset
        #
        try:
            tds = self.vault[key]
        except:
            raise Exception('Unable to find requested data in the vault')

        #
        #  Determine the period of record for the DataSeries that
        #  will be returned
        #
        if not first:
            newfirst = tds.startDate
        else:
            d = util.date_from_entry(first)
            newfirst = max(tds.startDate, d)

        if not last:
            newlast = tds.endDate
        else:
            d = util.date_from_entry(last)
            newlast = min(tds.endDate, d)

        #
        #  Trim the old dataset to match this new period
        #
        trimvals = util.trimDataValues(values=tds.dataVals,
                oldstart=tds.startDate, oldend=tds.endDate,
                newstart=newfirst, newend=newlast,
                intvl=tds.dataInterval)
            
        #
        #  Use the correct lake area in subsequent data conversion
        #
        if lake_area:
            lkarea = lake_area
        else:
            lkarea = util.getLakeArea(tds.dataLocation)
        
        #
        #  Now build a final dataset that has the correct units and 
        #  period of record.
        #
        try:
            newvals = util.convertValues(values=trimvals,
                    oldunits=tds.dataUnits, newunits=du, 
                    intvl=tds.dataInterval, area=lkarea, 
                    first=newfirst, last=newlast)
            kstr = tds.dataKind;
            istr = tds.dataInterval;
            lstr = tds.dataLocation;
            rds = DataSeries(kind=kstr, units=units, intvl=istr, loc=lstr,
                    first=newfirst, last=newlast, set=dt, values=newvals)
            return rds
        except:
            raise Exception('Error while attempting to convert data units in '
                          + 'DataVault.withdraw()')
