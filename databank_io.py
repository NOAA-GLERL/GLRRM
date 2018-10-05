import os
import re
import datetime as dt
import databank as databank
import databank_util as util


# set fill value for missing data or non-existent values in table
fill = util.MISSING_REAL

month_names = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',  
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')

# =================================================================
# ==============   BEGIN READ SECTION ============================
# =================================================================

def __get_valid_lines(filename):
    '''--------------------------------------------------------------------
    Read the specified file, returning an ordered list of strings.
    Each string is a line in the file that has valid data and all
    comments are removed. Comments are defined as anything on a
    line that comes after a # character.  If a line ends up being
    blank, it is not included in the returned list.  Quick example...
    File content:
      # this is a header comment
      # Year Mon Day Value
        1990  01  01 10.50
        1990  01  02  3.10
        1990  01  03 -9.99    # this is a missing value
        1990  01  04 21.00
        #  a bunch of data is missing here
                            # BLANK LINE
        1990  12  31  4.44

    Returned list:
    ('  1990  01  01 10.50', '  1990  01  02  3.10',
     '  1990  01  03 -9.99', '  1990  01  04 21.00',
     '  1990  12  31  4.44')

    The 9-line file turned into a 5-string list.
    --------------------------------------------------------------------
    '''
    linelist = []
    with open(filename, "r") as f:
        for line in f:
            i = line.find('#')
            if i == -1:
                s1 = line
            else:
                s1 = line[:i]
            s2 = s1.rstrip()
            if len(s2) > 0:
                linelist.append(s2)
    return linelist

#--------------------------------------------------------------------
def __get_meta_data(all_lines):
    '''
    Parse all lines of the file looking for lines that match the format
    of metadata specifiers.  When one of those is encountered, try to
    parse out the metadata in it.

    Note that with old CGLRRM files, the kind (type) and location metadata
    was not part of the header.  It was suggested to include it as part of
    the comments, but they were not part of a defined header entry.

    The old CGLRRM looked at only the first character of interval
    specification, so we will handle that as a special case.
    '''
    mkind  = 'na'
    munits = 'na'
    mintvl = 'na'
    mloc   = 'na'

    for line in all_lines:
        try:
            if line.find(':'):
                items = line.split(':')
                if len(items) == 2:
                    lstr = items[0].strip().lower()
                    rstr = items[1].strip().lower()
                    if (lstr == 'kind') or (lstr == 'type'):
                        mkind = databank.DataKind(rstr).primaryName()
                    if lstr == 'location':
                        mloc = databank.DataLocation(rstr).primaryName()
                    if lstr == 'units':
                        munits = databank.DataUnits(rstr).primaryName()
                    if lstr == 'interval':
                        mintvl = databank.DataInterval(rstr).primaryName()
                        if mintvl == 'na':
                            t = 'na'
                            if rstr[0] == 'd': t = 'dy'
                            if rstr[0] == 'w': t = 'wk'
                            if rstr[0] == 'q': t = 'qm'
                            if rstr[0] == 'm': t = 'mn'
                            if rstr[0] == 'y': t = 'yr'
                            mintvl = databank.DataInterval(t).primaryName()
        except:
            raise Exception('Error in header of the file')

    return mkind, munits, mintvl, mloc


#--------------------------------------------------------------------
def __detect_format(filename, data_lines, intvl):
    ''' Detect the format of filename based on the given data_lines
    and intvl.

    Loop through data_lines and check the number of items on each line,
    check for the type of delimiter, and confirm that all date strings
    can be cast as integers and that all data values can be cast as
    floats.
    Modification on 10/1/2018 -- gracefully accept faulty data values.

    Returns:
    --------
    a string representing format name (cglrrm, column, table or unknown)
    '''
    yy = None
    mm = None
    dd = None
    qq = None
    range_ok = True # are the date values possible (month=[1,12])?
    data_ok = True  # are the dates and data values actually numbers?
    format_ok = True # does this file match and acceptable format?
    if not data_lines: raise Exception('data_lines are empty in '+filename)

    if intvl == 'dy':
        #      ACCEPTED DAILY FORMATS
        #  CGLRRM: YYYY MM QQ VAL1 VAL2 ... VAL7 (VAL8)
        #  TABLE: YYYY-MM, VAL1, VAL2, ... VAL31(,)
        #  COLUMN: YYYY-MM-DD, VAL1(,)
        for num, line in enumerate(data_lines):
            cglrrm = False
            table = False
            column = False
            items = [s.strip() for s in line.split() if s]  # old format (space delim)
            if (len(items) == 10 or len(items) == 11) and ',' not in line:
                cglrrm = True
                first = 3
                yy = items[0]
                mm = items[1]
                qq = items[2]
            else:
                items = [s.strip() for s in line.split(',') if s] # new format (comma delim)
                if len(items) == 32:
                    table = True
                    first = 2
                    yy = items[0].split('-')[0]
                    mm = items[0].split('-')[1]
                if len(items) == 2:
                    column = True
                    first = 2
                    yy = items[0].split('-')[0]
                    mm = items[0].split('-')[1]
                    dd = items[0].split('-')[2]

            # ensure one and only one format was matched
            if sum([cglrrm, table, column]) != 1: format_ok = False

            #  check that dates are numbers
            #  data may be a problem, but that will be handled later.
            try:
                yy = int(yy)
                mm = int(mm)
                if qq: qq = int(qq)
                if dd: dd = int(dd)
            except:
                data_ok = False
                break

            # check date range
            if yy < 1: range_ok = False
            if mm < 1 or mm > 12: range_ok = False
            if qq:
                if qq < 1 or qq > 4: range_ok = False
            if dd:
                if dd < 1 or dd > 31: range_ok = False


            if not all([format_ok, data_ok, range_ok]): break


    if intvl == 'wk':
        #      ACCEPTED WEEKLY FORMATS
        #  CGLRRM: YYYY MM DD VAL1
        #  COLUMN: YYYY-MM-DD, VAL1(,)
        for num, line in enumerate(data_lines):
            cglrrm = False
            table = False
            column = False
            items = [s.strip() for s in line.split() if s]  # old format (space delim)
            if len(items) == 4 and ',' not in line:
                cglrrm = True
                first = 3
                yy = items[0]
                mm = items[1]
                dd = items[2]
            else:
                items = [s.strip() for s in line.split(',') if s] # new format (comma delim)
                if len(items) == 2:
                    column = True
                    first = 1
                    yy = items[0].split('-')[0]
                    mm = items[0].split('-')[1]
                    dd = items[0].split('-')[2]

            # ensure one and only one format was matched
            if sum([cglrrm, column]) != 1: format_ok = False

            # check that dates are numbers
            #  data may be a problem, but that will be handled later.
            try:
                yy = int(yy)
                mm = int(mm)
                dd = int(dd)
            except:
                data_ok = False
                break

            # check that dates are possible
            if yy < 1:             range_ok = False
            if mm < 1 or mm > 12:  range_ok = False
            if dd < 1 or dd > 31: range_ok = False

            if not all([format_ok, data_ok, range_ok]): break

    if intvl == 'qm':
        #    ACCEPTED QTR-MONTHLY FORMATS
        #  CGLRRM: YYYY QQ VAL1 VAL2 ... VAL12
        #  TABLE: YYYY-MM, VAL1, VAL2, VAL3, VAL4(,)
        #  COLUMN: YYYY-MM-QQ, VAL1(,)
        for num, line in enumerate(data_lines):
            cglrrm = False
            table = False
            column = False
            items = [s.strip() for s in line.split() if s]
            if len(items) == 14 and ',' not in line:
                cglrrm = True
                first = 2
                yy = items[0]
                qq = items[1]
            else:
                items = [s.strip() for s in line.split(',') if s]
                if len(items) == 5:
                    table = True
                    first = 1
                    yy = items[0].split('-')[0]
                    mm = items[0].split('-')[1]
                if len(items) == 2:
                    column = True
                    first = 1
                    yy = items[0].split('-')[0]
                    mm = items[0].split('-')[1]
                    qq = int(items[0].split('-')[2])

            # ensure one and only one format was matched
            if sum([cglrrm, table, column]) != 1: format_ok = False

            # check that dates are numbers
            #  data may be a problem, but that will be handled later.
            try:
                yy = int(yy)
                if qq: qq = int(qq)
                if mm: mm = int(mm)
            except:
                data_ok = False
                break

            # check that dates are possible
            if yy < 1: range_ok = False
            if qq:
                if qq < 1 or qq > 4: range_ok = False
            if mm:
                if mm < 1 or mm > 12:  range_ok = False

            if not all([format_ok, data_ok, range_ok]): break


    if intvl == 'mn':
        #      ACCEPTED MONTHLY FORMATS
        #  CGLRRM: YYYY VAL1 VAL2 ... VAL12
        #  TABLE: YYYY, VAL1, VAL2, ... VAL12(,)
        #  COLUMN: YYYY-MM, VAL1(,)
        for num, line in enumerate(data_lines):
            cglrrm = False
            table = False
            column = False
            items = [s.strip() for s in line.split() if s]
            if len(items) == 13 and ',' not in line:
                cglrrm = True
                first = 2
                yy = int(items[0])
            else:
                items = [s.strip() for s in line.split(',') if s]
                if len(items) == 13:
                    table = True
                    first = 2
                    yy = items[0].split('-')[0]
                if len(items) == 2:
                    column = True
                    first = 1
                    yy = items[0].split('-')[0]
                    mm = items[0].split('-')[1]


            # ensure one and only one format was matched
            if sum([cglrrm, table, column]) != 1: format_ok = False

            # check that dates are numbers
            #  data may be a problem, but that will be handled later.
            try:
                yy = int(yy)
                if mm: mm = int(mm)
            except:
                data_ok = False
                break


            # check that dates are possible
            if yy < 1:             range_ok = False
            if mm:
                if mm < 1 or mm > 12:  range_ok = False

            if not all([format_ok, data_ok, range_ok]): break


    # raise exceptions
    if not format_ok:
        raise Exception('Error determining format for ' +filename+' near line '
                        '#:'+str(num)+' (not counting hdrs/comments)')
    if not range_ok:
        raise Exception('Dates exceed acceptable range in '+filename+' near '
                        'line #:'+str(num)+' (not counting hdrs/comments)')
    if not data_ok:
        raise Exception('Data vals or dates could not be cast as floats in '
                        + filename+' near line #:'+str(num)+' (not counting '
                        'hdrs/comments)')

    if cglrrm: return 'cglrrm'
    if table: return 'table'
    if column: return 'column'
    return 'unknown'


#--------------------------------------------------------------------
def read_file(filename, kind=None, units=None, intvl=None, loc=None, 
              set=None, missing_value=None):
    '''
    Read in data and metadata from filename and return a dataseries
    object.

    Parameters
    ----------
    filename: string
        The desired filename to read from.
    kind : string, optional
        The kind of data ('precip') to be read in.  If filename    does
        not contain kind metadata (header info), kind must be set
        on the function call.
    units : string, optional
        The units of data ('mm') to be read in.  If filename does not
        contain units metadata (header info), units must be set.
    intvl : string, optional
        The intierval of data ('weekly') to be read in.  If filename
        does not contain intvl metadata (header info), intvl must be set.
    loc : string, optional
        The loc of data ('superior') to be read in.  If filename does
        not contain loc metadata (header info), loc must be set.
    set : string, optional
        The set name to be assigned to the data.  This is useful when
        we have multiple data sets with the same metadata but different 
        values; e.g. multiple forecasts of the same variable for the
        same time period (ensemble members). If not provided, the default
        value will be assigned.

    Returns
    -------
    dataseries : object
        an object of class DataSeries with attributes dataKind,
        dataUnits, dataInterval, dataLocation, startDate,
        endDate, dataSet, and dataVals (the actual data)

    Notes
    -----
    If keywords kind, units, intvl, or loc are included in function call
    AND included in the metadata of filename, then they must match,
    otherwise an exception will be raised.

    Classic CGLRRM files were NOT required to include loc or kind in
    the metadata (header info) so these keyword args will likely
    be required when processing these files.
    '''


    #
    #  Get all of the non-comment content from the file as
    #  a list of strings
    #
    all_lines = __get_valid_lines(filename)

    #
    #  Create objects for each of the supplied metadata items
    #
    try:
        mkind  = None
        munits = None
        mintvl = None
        mloc   = None
        if kind:
            mkind = databank.getPrimaryName(meta='kind', name=kind)
        if units:
            munits = databank.getPrimaryName(meta='units', name=units)
        if intvl:
            mintvl = databank.getPrimaryName(meta='interval', name=intvl)
        if loc:
            mloc = databank.getPrimaryName(meta='location', name=loc)
    except:
        raise Exception('Invalid metadata passed to read_file()')

    #
    #  Assign correct value to mset variable
    #
    if set:
        mset = set
    else:
        mset = 'na'      # default value
        
    #
    #  Parse the data file lines to get any metadata
    #  specified within the file.
    #
    try:
        fkind, funits, fintvl, floc = __get_meta_data(all_lines)
    except:
        raise Exception('Unable to process metadata from ' + filename)

    #
    #  If any of the metadata items from the file are undefined,
    #  set that metadata object to None
    #
    if fkind  == 'na': fkind  = None
    if funits == 'na': funits = None
    if fintvl == 'na': fintvl = None
    if floc   == 'na': floc   = None

    #
    #  Do we have sufficient metadata to continue?
    #
    if not mkind  and not fkind:  raise Exception('Missing metadata (kind) for'
                                                  ' processing ' + filename)
    if not munits and not funits: raise Exception('Missing metadata (units) for'
                                                  ' processing ' + filename)
    if not mintvl and not fintvl: raise Exception('Missing metadata (interval)'
                                                  ' for processing ' + filename)
    if not mloc   and not floc:   raise Exception('Missing metadata (location)'
                                                  ' for processing ' + filename)

    #
    #  If we have metadata from both caller and file, do they match?
    #
    ok = True
    if mkind and fkind:
        if mkind != fkind:
            ok = False
    if munits and funits:
        if munits != funits:
            ok = False
    if mintvl and fintvl:
        if mintvl != fintvl:
            ok = False
    if mloc and floc:
        if mloc != floc:
            ok = False
    if not ok:
        raise Exception('Inconsistent metadata for processing ' + filename)

    #
    #  If needed, assign metadata from file to our primary variables
    #
    if not mkind:  mkind  = fkind
    if not munits: munits = funits
    if not mintvl: mintvl = fintvl
    if not mloc:   mloc   = floc

    #
    #  Process the file appropriately.  Note that the
    #  DataSeries object that is returned may have some missing
    #  metadata items if they were left blank in the file, but
    #  specified by the caller of this procedure.


    if mintvl == 'dy' or mintvl == 'wk' or mintvl == 'qm' or mintvl == 'mn':
        # beginning of what procData did (as far as I can tell)
        start = util.MISSING_DATE       # "missing" date value
        end = util.MISSING_DATE       # "missing" date value
        data_lines = []
        for line in all_lines:
            try:
                if line.find(':') > 0:
                    items = line.split(':')
                    if len(items) != 2:
                        raise Exception('Invalid header content or misplaced colon :')
                else:
                    data_lines.append(line)
            except:
                raise Exception('Error separating header and data lines.')

        format_name = __detect_format(filename, data_lines, mintvl)
        # end of what procData did
    else:
        raise Exception('Invalid interval')

    if format_name == 'unknown':
        raise Exception('Format of ' + filename + ' could not be recognized')
    try:
        start, end, datavals = __parse_data(filename, format_name, mintvl,
                                          data_lines, missing_value)
    except:
        raise Exception('Error calling parse_data for '+filename+';  format:'
                        +format_name+'; interval:'+mintvl)
    if not datavals:
        raise Exception('No data read in for '+filename+': format:'+format_name
                        +'; interval:'+mintvl)
    if start > end:
        raise Exception('endDate precedes startDate (likely that no dates were'
                        'recognized) for '+filename+': format:'+format_name
                        +'; interval:'+mintvl)

    ds = databank.DataSeries(kind=mkind, units=munits, intvl=mintvl, loc=mloc,
                             first=start, last=end, set=mset, values=datavals)
    return ds


#--------------------------------------------------------------------
def __parse_data(filename, format_name, intvl, data_lines, missing_value=None):
    '''
    --------------------------------------------------------------------
    This function returns a complete DataSeries object.
    --------------------------------------------------------------------
    Psuedo code:
            parse_data(args)
                 selection for interval
                    selection for format
                        loop thru to get start, end
                        loop thru to get data values
                return data start,end, data
                
    If missing_value is specified, it will be used to evaluate each data value
    read from the file. Values <= missing_value will be treated as missing.
    If the value is not specified, then the default value of -999999 will be
    used as the criterion value.
    '''

    # initialize return values
    start = dt.date(2999, 12, 31)    # ridiculously far into future
    end = dt.date(1000, 1, 1)    # ridiculously far into past
    datavals = None
    
    # set up for handling missing values
    if missing_value:
        missval = missing_value
    else:
        missval = -999999

    #
    if intvl == 'dy':
        #      ACCEPTED DAILY FORMATS
        #  CGLRRM: YYYY MM QQ VAL1 VAL2 ... VAL7 (VAL8)
        #  TABLE: YYYY-MM, VAL1, VAL2, ... VAL31(,)
        #  COLUMN: YYYY-MM-DD, VAL1(,)

        # first loop through to get the first/last date
        if format_name == 'cglrrm':
            for line in data_lines:
                items = [s.strip() for s in line.split() if s]
                yy = int(items[0])
                mm = int(items[1])
                qq = int(items[2])
                sd, ed = util.getQtrMonthStartEnd(year=yy, month=mm, qtr=qq)
                d1 = dt.date(yy, mm, sd)
                d2 = dt.date(yy, mm, ed)
                start = min(start, d1)
                end = max(end, d2)

        if format_name == 'table':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                ndays = util.days_in_month(year=yy, month=mm)
                d1 = dt.date(yy, mm, 1)
                d2 = dt.date(yy, mm, ndays)
                start = min(start, d1)
                end = max(end, d2)

        if format_name == 'column':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                dd = int(items[0].split('-')[2])
                d1 = dt.date(yy, mm, dd)
                d2 = d1
                start = min(start, d1)
                end = max(end, d2)


        ndays = (end - start).days + 1
        if ndays < 1:
            return start, end, datavals

        #  Now parse the lines to get data values and assign.
        datavals = [util.MISSING_REAL] * ndays
        if format_name == 'cglrrm':
            for line in data_lines:
                items = [s.strip() for s in line.split() if s]
                yy = int(items[0])
                mm = int(items[1])
                qq = int(items[2])
                sd, ed = util.getQtrMonthStartEnd(year=yy, month=mm, qtr=qq)
                d1 = dt.date(yy, mm, sd)
                for i in range(3, len(items)):
                    try:
                        val = float(items[i])
                    except:
                        val = util.MISSING_REAL
                    ndx = (d1 - start).days + i - 3
                    datavals[ndx] = val

        if format_name == 'table':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                d1 = dt.date(yy, mm, 1)
                ndays = util.days_in_month(year=yy, month=mm) + 1
                for i in range(1, ndays):
                    try:
                        val = float(items[i])
                    except:
                        val = util.MISSING_REAL
                    ndx = (d1 - start).days + i - 1
                    datavals[ndx] = val

        if format_name == 'column':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                dd = int(items[0].split('-')[2])
                d1 = dt.date(yy, mm, dd)
                #ndx = (d1 - start).days - 1
                ndx = (d1 - start).days 
                try:
                    datavals[ndx] = float(items[1])
                except:
                    datavals[ndx] = util.MISSING_REAL



    if intvl == 'wk':
        #      ACCEPTED WEEKLY FORMATS
        #  CGLRRM: YYYY MM DD VAL1
        #  COLUMN: YYYY-MM-DD, VAL1(,)
        # note: for weekly, all formats have same formatting (just different delims)
        # lets replace the delimiters so we can process them all the same way
        # This might not be very pythonic/good style but eliminates redundant code

        if format_name == 'cglrrm':
            for i, line in enumerate(data_lines):
                # replace YYYY MM DD with YYYY-MM-DD
                line = re.sub(r"(^[0-9]{4})\s+([0-9]{1,2})\s+([0-9]{1,2})", r"\1-\2-\3", line)
                # replace all white spaces with ", "
                line = re.sub(r"\s+", ", ", line)
                data_lines[i] = line

        # loop to get start/end
        for line in data_lines:
            items = [s.strip() for s in line.split(',') if s]
            yy = int(items[0].split('-')[0])
            mm = int(items[0].split('-')[1])
            dd = int(items[0].split('-')[2])
            d1 = util.getFridayDate(year=yy, month=mm, day=dd)
            start = min(start, d1)
            end = max(end, d1)


        nweeks = int((end - start).days/7) + 1
        if nweeks < 1:
            return start, end, datavals
        datavals = [util.MISSING_REAL] * nweeks

        # loop to get data values
        for line in data_lines:
            items = [s.strip() for s in line.split(',') if s]
            yy = int(items[0].split('-')[0])
            mm = int(items[0].split('-')[1])
            dd = int(items[0].split('-')[2])
            d1 = util.getFridayDate(year=yy, month=mm, day=dd)
            ndx = int((d1 - start).days/7)
            try:
                datavals[ndx] = float(items[1])
            except:
                datavals[ndx] = util.MISSING_REAL

    if intvl == 'qm':
       #      ACCEPTED QM FORMATS
       #  CGLRRM: YYYY QQ VAL1 VAL2 ... VAL12
       #  TABLE: YYYY-MM, VAL1, VAL2, VAL3, VAL4(,)
       #  COLUMN: YYYY-MM-QQ, VAL1(,)

    # first loop through to get start/end dates
        if format_name == 'cglrrm':
            start_qtr = 1  # this format can NOT have partial years
            for line in data_lines:
                items = [s.strip() for s in line.split() if s]
                yy = int(items[0].split()[0])
                sd = util.getQtrMonthStartEnd(year=yy, month=1, qtr=1)[0]
                ed = util.getQtrMonthStartEnd(year=yy, month=12, qtr=4)[1]
                d1 = dt.date(yy, 1, sd)
                d2 = dt.date(yy, 12, ed)
                start = min(start, d1)
                end = max(end, d2)
            nqtrs = 48*(end.year - start.year + 1)
            datavals = [util.MISSING_REAL] * nqtrs


        if format_name == 'table':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                sd = util.getQtrMonthStartEnd(year=yy, month=mm, qtr=1)[0]
                ed = util.getQtrMonthStartEnd(year=yy, month=mm, qtr=4)[1]
                d1 = dt.date(yy, mm, sd)
                d2 = dt.date(yy, mm, ed)
                if d1 <= start:
                    start = d1
                    start_mon = mm
                if d2 >= end:
                    end = d2
                    end_mon = mm
            nqtrs = 48*(end.year - start.year) + 4*(end_mon - start_mon + 1)
            datavals = [util.MISSING_REAL] * nqtrs


        if format_name == 'column':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                qq = int(items[0].split('-')[2])
                sd, ed = util.getQtrMonthStartEnd(year=yy, month=mm, qtr=qq)
                d1 = dt.date(yy, mm, sd)
                d2 = dt.date(yy, mm, ed)
                if d1 <= start:
                    start = d1
                    start_mon = mm
                    start_qtr = qq
                if d2 >= end:
                    end = d2
                    end_mon = mm
                    end_qtr = qq
            nqtrs = 48*(end.year - start.year) + 4*(end_mon - start_mon) + (end_qtr - start_qtr) + 1
            datavals = [util.MISSING_REAL] * nqtrs

        # now loop through to get datavals
        if format_name == 'cglrrm':
            for line in data_lines:
                items = [s.strip() for s in line.split() if s]
                yy = int(items[0])
                qq = int(items[1])

                yy_delta = yy - start.year
                qq_delta = qq - start_qtr

                for i in range(2, len(items)):
                    ndx = 48*yy_delta + qq_delta + 4*(i - 2)
                    try:
                        datavals[ndx] = float(items[i])
                    except:
                        datavals[ndx] = util.MISSING_REAL


        if format_name == 'table':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])

                yy_delta = yy - start.year
                mm_delta = mm - start_mon

                for i in range(1,len(items)):
                    ndx = 48*yy_delta + 4*mm_delta + (i - 1)
                    try:
                        datavals[ndx] = float(items[i])
                    except:
                        datavals[ndx] = util.MISSING_REAL

        if format_name == 'column':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                qq = int(items[0].split('-')[2])

                yy_delta = yy - start.year
                mm_delta = mm - start_mon
                qq_delta = qq - start_qtr

                ndx = 48*yy_delta + 4*mm_delta + qq_delta
                try:
                    datavals[ndx] = float(items[i])
                except:
                    datavals[ndx] = util.MISSING_REAL



    if intvl == 'mn':
       #      ACCEPTED MN FORMATS
       #  CGLRRM: YYYY VAL1 VAL2 ... VAL12
       #  TABLE: YYYY, VAL1, VAL2, ..., VAL12
       #  COLUMN: YYYY-MM, VAL1(,)

    # first loop through to get start/end dates
        if format_name == 'column':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])
                d1 = dt.date(yy, mm, 1)
                ndays = util.days_in_month(yy, mm)
                d2 = dt.date(yy, mm, ndays)
                start = min(start, d1)
                end = max(end, d2)
                nmonths = 12*(end.year - start.year) + (end.month - start.month + 1)

        # change cglrrm data_lines to be comma delimited
        if format_name == 'cglrrm':
            for i, line in enumerate(data_lines):
                line = re.sub(r"\s+", ", ", line)
                data_lines[i] = line

        if format_name == 'cglrrm' or format_name == 'table':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0])
                d1 = dt.date(yy, 1, 1)
                d2 = dt.date(yy, 12, 31)
                start = min(start, d1)
                end = max(end, d2)
                nmonths = 12*(end.year - start.year + 1)

        datavals = [util.MISSING_REAL] * nmonths



        # now loop through to get datavals
        if format_name == 'column':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0].split('-')[0])
                mm = int(items[0].split('-')[1])

                yy_delta = yy - start.year
                mm_delta = mm - start.month

                ndx = 12*yy_delta + mm_delta
                try:
                    datavals[ndx] = float(items[1])
                except:
                    datavals[ndx] = util.MISSING_REAL


        if format_name == 'cglrrm' or format_name == 'table':
            for line in data_lines:
                items = [s.strip() for s in line.split(',') if s]
                yy = int(items[0])
                yy_delta = yy - start.year
                for i in range(1, len(items)):
                    ndx = 12*yy_delta + i - 1
                    try:
                        datavals[ndx] = float(items[i])
                    except:
                        datavals[ndx] = util.MISSING_REAL


    #  Once all data values have been read and the datavals[] list 
    #  has been populated, check all values for meeting the criteria
    #  as "missing".  Replace those values with MISSING_REAL.
    for i in range(len(datavals)):
        if datavals[i] <= missval:
            datavals[i] = util.MISSING_REAL
                    
    #
    return start, end, datavals



# =================================================================
# ==============   BEGIN WRITE SECTION ============================
# =================================================================

#--------------------------------------------------------------------
def write_file(filename, file_format, series, overwrite=False, 
               width=9, prec=2, missing_value=None):
    '''

    Write dataseries (metadata and datavals) to filename

    Parameters
    ----------
    filename: string
        The desired filename to write to.
    file_format: string
        The format of the output file ("table" or "column")
        note that weekly data only can be written as column format
    series: object
        The dataseries object created by databank or read_file
        to be written to filename.
    overwrite: boolean, optional
        Flag to indicate whether or not to overwrite existing file
        named "filename". default is False (i.e don't overwrite files)
    width: string, optional
        Used to specify the formatting (column width) of values  written to filename. If
        set, should be used in conjunction with prec. MINIMUM VALUE = 6 to properly print
        the fill value of -999.9
    prec: string, optional
        Used to specify the formatting (decimal precision) of values written to filename.
        if set, it should be used in conjunction with width. MINIMUM VALUE = 1  to properly
        print the fill value of -999.9
    missing_value: numeric, optional
        Used to specify the numeric value to output when the data value is missing.
        Note that the value will be formatted according to width and prec.
    '''



    # check input arguments for validity
    if file_format not in ['table', 'column']:
        raise Exception('invalid file_format: use \"table\" or \"column\"')

    if series.dataInterval == 'wk' and file_format == 'table':
        raise Exception('weekly data must use \"column\" formatting')
        
    try:
        __write_metadata(filename, series, overwrite)
    except:
        raise Exception('Unable to write metadata to ' + filename)

    try:
        __write_datavals(filename, file_format, series, width, prec, missing_value)
    except:
        raise Exception('Unable to write datavals to ' + filename)

        
#--------------------------------------------------------------------
def __missing_value_string(width, prec, value=None):
    ''' constructs a missing data value appropriate for the specified width/precision 
    width: integer in range 3 to 20
        Used to specify the formatting (column width) of missing values for output.
        Minimum value = 3
        Maximum value = 20 (a little arbitrary, but reasonable)
    prec: integer in range 0 to 8
        Used to specify the formatting (decimal precision) of missing values for output.
        Minimum value = 0, which results in an integer value
        Maximum value = 8 (a little arbitrary, but reasonable)
        Values from 1 to 8 result in an output with that many digits after the decimal
        prec must be less than or equal to (width - 3)
    value: optional numeric argument. By convention, it should be something like
        -99999 or -9.99e5 or some other string that looks like a minus 99999.
        The value will be formatted according to width & prec.  
    
    e.g. (3,0)  -> -99 
         (4,1)  -> -9.9
         (9,4)  -> -999.9999
         (10,2) -> -999999.99
         (11,0) -> -9999999999
    '''
    if (width < 3):  raise Exception('Invalid width specifier < 3')
    if (width > 20): raise Exception('Invalid width specifier > 20')
    if (prec < 0):   raise Exception('Invalid precision specifier < 0')
    if (prec > 8):   raise Exception('Invalid precision specifier > 8')
    if (prec > (width-3)):  raise Exception('Invalid precision specifier > (width-3)')

    #
    #  No value specified
    #
    if value==None:
        s = '-'
        l = width - prec - 1
        if prec > 0: 
            l = l - 1
        for i in range(l):
            s = s + '9'
        if prec > 0:
            s = s + '.'
            for i in range(prec):
                s = s + '9'
        return s

    #
    f = '%' + str(width) + '.' + str(prec) + 'f'
    s = f % value
    return s
        
        
        
#--------------------------------------------------------------------
def __write_metadata(filename, dataseries, overwrite):
    ''' writes metadata to filename '''

    # check if file exists already
    if not overwrite and os.path.isfile(filename): raise Exception(filename + ' already exists')

    # create new file for writing
    try:
        f = open(filename, 'w')
    except:
        raise Exception('could not open ' + filename + ' for writing')

    # write metadata
    line = ':'.join(['KIND', dataseries.dataKind])
    f.write(line + '\n')
    line = ':'.join(['UNITS', dataseries.dataUnits])
    f.write(line + '\n')
    line = ':'.join(['INTERVAL', dataseries.dataInterval])
    f.write(line + '\n')
    line = ':'.join(['LOCATION', dataseries.dataLocation])
    f.write(line + '\n')


    # write comment line with date of generation
    now_str = dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = ''.join(['# file created: ', now_str])
    f.write(line + '\n')

    f.close()


#--------------------------------------------------------------------
def __write_datavals(filename, file_format, dataseries, width, prec, missing_value=None):
    ''' Write the dataseries.dataVals to filename
    optional parameter missing_value specifies the value to output as a missing value
    '''
    intvl = dataseries.dataInterval
    vals = dataseries.dataVals.copy()    # make an actual copy of the list
    start = dataseries.startDate
    end = dataseries.endDate
    fmt = '{0:' + str(width) + '.' + str(prec) + 'f}'

    # define fill value
    if missing_value:
        fill = missing_value
    else:
        fill = __missing_value_string(width-1, prec)
    
    # go through the temporary copy of the dataVals list, replacing all missing data
    # with the fill value.  By doing this now, we can process things efficiently
    # in all of the subsequent code.
    if (prec == 0):
        mval = int(fill)
    else:
        mval = float(fill)
    l = len(vals)
    for i in range(l):
        if vals[i] < util.MISSING_TEST:
            vals[i] = mval

    #
    f = open(filename, 'a')

    # daily
    if intvl == 'dy':
        if file_format == 'table':
            # prepare and write header (YYYY-MM, D1, D2, ...D31)
            days = range(1, 32)
            hdr = ', '.join(['{:>{wid}}'.format('D'+str(d), wid=width) for d in days])
            hdr_line = ','.join(['#YYYY-MM', hdr])
            f.write(hdr_line + '\n')

            # initialize date for loop
            date = start
            dom = util.days_in_month(date.year, date.month)
            first = 0          # vals index for first day of the month
            last = first + dom # vals index for last day of the month
            while date <= end:
                # construct the line to write to filename
                data_line = [fill] * 31
                data_line[0:dom] = vals[first:last]
                data_line = ', '.join([fmt.format(d) for d in data_line])
                full_line = ', '.join(['{:%Y-%m}'.format(date), data_line])
                f.write(full_line + '\n')

                # update date and indices for next iteration
                date = date + dt.timedelta(days=dom)
                dom = util.days_in_month(date.year, date.month)
                first = last
                last = first + dom

        if file_format == 'column':
            # prepare and write header (YYYY-MM-DD, VAL)
            hdr = '{:>{wid}}'.format('VAL', wid=width)
            hdr_line = ','.join(['#YYYY-MM-DD', hdr])
            f.write(hdr_line + '\n')

            # initialize date for loop
            date = start
            i = 0
            while date <= end:
                # construct a line to write to file
                data_line = fmt.format(vals[i])
                full_line = ', '.join(['{:%Y-%m-%d}'.format(date), data_line])
                f.write(full_line + '\n')

                # update
                i = i + 1
                date = date + dt.timedelta(days=1)


    # weekly
    if intvl == 'wk':
        # prepare and write header (YYYY-MM-DD, VAL)
        hdr = '{:>{wid}}'.format('VAL', wid=width)
        hdr_line = ','.join(['#YYYY-MM-DD', hdr])
        f.write(hdr_line + '\n')

        # init date for loop
        date = start
        i = 0
        while date <= end:
            # construct a line to write to file
            data_line = fmt.format(vals[i])
            full_line = ', '.join(['{:%Y-%m-%d}'.format(date), data_line])
            f.write(full_line + '\n')

            # update
            i = i + 1
            date = date + dt.timedelta(days=7)


    # quarter-monthly
    if intvl == 'qm':
        # this is not actually the format that we discussed, but I think its better
        if file_format == 'table':
            # prepare and write header (YYYY-MM, Q1, Q2, Q3, Q4)
            qtrs = range(1, 5)
            hdr = ', '.join(['{:>{wid}}'.format('Q'+str(q), wid=width) for q in qtrs])
            hdr_line = ','.join(['#YYYY-MM', hdr])
            f.write(hdr_line + '\n')

            date = start
            i = 0
            while date <= end:
                # construct line and write
                data_line = vals[i:i+4]
                data_line = ', '.join([fmt.format(d) for d in data_line])
                full_line = ', '.join(['{:%Y-%m}'.format(date), data_line])
                f.write(full_line + '\n')

                # update
                i = i + 4
                dom = util.days_in_month(date.year, date.month)
                date = date + dt.timedelta(days=dom)

        if file_format == 'column':
            # prepare and write header (YYYY-MM-QQ, VAL)
            hdr = '{:>{wid}}'.format('VAL', wid=width)
            hdr_line = ','.join(['#YYYY-MM-QQ', hdr])
            f.write(hdr_line + '\n')

            # init date for loop
            date = start
            i = 0
            q = 1  # quarter
            while date <= end:
                # construct line to write
                data_line = fmt.format(vals[i])
                date_str = '-'.join(['{:%Y-%m}'.format(date), '{:02d}'.format(q)])
                full_line = ', '.join([date_str, data_line])
                f.write(full_line + '\n')

                # update
                i = i + 1
                q = q + 1
                if i % 4 is 0:
                    q = 1
                    dom = util.days_in_month(date.year, date.month)
                    date = date + dt.timedelta(days=dom)


    if intvl == 'mn':
        if file_format == 'table':
            # prepare and write header (YYYY, M1, M2, M3, ..., M12)
            mons = range(1, 13)
            hdr = ',       '.join([month_names[m-1] for m in mons])
            hdr_line = ',      '.join(['#YYYY', hdr])
            f.write(hdr_line + '\n')

            # init date for loop
            date = int('{:%Y}'.format(start))
            end = int('{:%Y}'.format(end))
            i = 0
            while date <= end:
                # construct line to write
                data_line = vals[i:i+12]
                data_line = ', '.join([fmt.format(d) for d in data_line])
                full_line = ', '.join(['{:4d}'.format(date), data_line])
                f.write(full_line + '\n')

                # update
                i = i + 12
                date = date + 1


        if file_format == 'column':
            # prepare and write header (YYYY-MM, VAL)
            hdr = '{:>{wid}}'.format('VAL', wid=width)
            hdr_line = ','.join(['#YYYY-MM', hdr])
            f.write(hdr_line + '\n')

            # init date for loop
            date = start
            i = 0
            
            while date <= end:
                # construct line to write
                data_line = fmt.format(vals[i])
                full_line = ', '.join(['{:%Y-%m}'.format(date), data_line])
                f.write(full_line + '\n')

                # update
                i = i + 1
                dom = util.days_in_month(date.year, date.month)
                date = date + dt.timedelta(days=dom)

    f.close()

