"""All-in-one Spatial Activity Processor

By Joseph Broach <jbroach@pdx.edu>

An environment for processing, viewing, and analyzing spatial data. It was 
written with travel data in mind but may be useful in other fields of study. 
I created ASAP as part of my PhD work at Portland State University.

orcycle version only has Deployment class and limited psql output functionality.

"""

# Aug 8, 2013 - created

import csv
import time
import calendar
import math
import psycopg2 as ppg2
from datetime import datetime
from pprint import pprint
import pg_tools_orcycle as pg
import us_time_zones as ustz
#import am_tools as amt


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class Deployment():
    """A data collection deployment"""

    def __init__(self, time_start=None, time_stop=None, num=-9, 
                 gps=(None, None), 
                 am=(None, None), unit_id=None, group_id=None, phase_id=None,
                 activities=None, local_date_start=None, local_time_start=None,
                 local_time_stop=None, stats=None, **kwargs):
        """Initialize a Deployment instance.

        time_start : UTC time start [float]
        time_stop : UTC time stop [float]
        num : deployment number [int]
        ntd : number of travel days [int]
        gps : (gps_file, format)
        am : (am_file, format)
        activities : list of activities

        """
        self.time_start = time_start
        self.time_stop = time_stop
        self.num = num
        self.unit_id = unit_id
        self.group_id = group_id
        self.phase_id=phase_id
        self.activities = activities
        self.local_date_start = local_date_start
        self.local_time_start = local_time_start
        self.local_time_stop = local_time_stop
        if stats:
            self.stats = stats
        else:
            self.stats = {}
        if (type(gps[0]) == dict and 
              type(gps[0][gps[0].keys()[0]]) == pg.GPSCollection):
            # GPS data already processed
            self.gps = gps[0]
        elif gps[0]:
            gpsd = GPSData(gps[1])
            if gpsd.data_type == 'csv':
                self.gps = gpsd.load_csv(gps[0], 
                                         local_timezone='Pacific', 
                                         time_start=time_start, 
                                         time_stop=time_stop, 
                                         astable=kwargs.get('astable_gps', False))
            elif gpsd.data_type == 'psql':
                # gpsd.test(gpsd.data_type, db_info=gps[0], 
                          # unit_id=self.unit_id, deployment=self.num)
                # exit(0)
                self.gps = gpsd.load_psql(gps[0],
                                          local_timezone='Pacific', 
                                          time_start=time_start, 
                                          time_stop=time_stop,
                                          unit_id=self.unit_id,
                                          phase_id=self.phase_id,
                                          deployment=self.num)
        else:
            self.gps = pg.GPSCollection()
        if am[0]:
            self.am, header = amt.read_am(am[0], 
                                          time_start=time_start, 
                                          time_stop=time_stop)
            # self.attach_counts()
        else:
            #self.am = amt.AMCollection()
            self.am = None
        if activities is None:
            self.activities = []
        else:
            self.activities = activities
        self.kwargs = kwargs
        
    def get_stats(self):
        """Calculate statistics peculiar to deployment."""
        onedayseconds = 86400
        onedayminutes = 1440
        self.stats['gpstrips'] = {}
        self.stats['amhours'] = {}
        for day in range(self.kwargs['n_days']):
            self.stats['gpstrips'][day + 1] = 0
            self.stats['amhours'][day + 1] = 0.0
        self.stats['totalgpstrips'] = 0
        if self.gps and self.gps.trips:
            for trip in self.gps.trips:
                day = math.floor((trip.time_start - self.time_start) / 
                                 onedayseconds) + 1
                self.stats['gpstrips'][day] += 1
                self.stats['totalgpstrips'] += 1
        if self.am and self.am.counts:
            stop = 0
            for day in sorted(self.stats['amhours'].keys()):
                start = stop
                stop += 1440
                a = amt.AMCollection(counts=self.am.counts[start: stop])
                a.get_stats()
                self.stats['amhours'][day] = a.nvalid / 60.0
        
    def write_psql(self, db_conn, schema, GPSpoints=None, GPStrips=None, 
                   GPStripstages=None,
                   Activities=None, MOA=None, GPSroutes=None, GPSlinks=None,
                   Persons=None):
        """Output deployment to standard set of psql tables."""
        if GPSroutes and self.gps:
            query = ("INSERT INTO {}.{} (gpstripid,stageid,"
                     "aps,oddlinks,distancelink,distancepoint,"
                     "geom) VALUES ".format(schema, GPSroutes))
            data = []
            for trip_id, gpsc in self.gps.iteritems():
                for trip in gpsc.trips:
                    if trip.route:
                        r = trip.route
                        link_dst_mi = float(r.stats['link_dst']) * 0.000621371
                        query += ('(%s,%s,%s,%s,%s,%s,%s),')
                        data.extend([trip_id, trip.num,
                                     r.stats['aps'], r.stats['odd_links'],
                                     link_dst_mi, trip.stats['dist_GPS_mi'], 
                                     r.geom])
            query = query[:-1]
            if data:
                db_conn.run(query, data, commit=True)

        if GPSlinks and self.gps:
            query = ("INSERT INTO {}.{} (gpstripid,stageid,"
                     "featureorder,featureid,direction,"
                     "startpoint,endpoint,mps,odd,startpos,endpos) VALUES "
                     .format(schema, GPSlinks))
            data = []
            for trip_id, gpsc in self.gps.iteritems():
                for trip in gpsc.trips:
                    tk = trip.kwargs
                    if trip.route:
                        r = trip.route
                        j = 1
                        for i in r.links:
                            query += ('(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s),')
                            startpos = i.__dict__.get('pos0', 0.0)
                            # if 'pos0' in i.kwargs:
                                # startpos = i.kwargs['pos0']
                            # else:
                                # startpos = 0.0
                            endpos = i.__dict__.get('pos1', 1.0)
                            # if 'pos1' in i.kwargs:
                                # endpos = i.kwargs['pos1']
                            # else:
                                # endpos = 1.0
                            if i.id < 0:
                                direction = -1
                                i.id = -1 * i.id
                            else:
                                direction = 1
                            data.extend([tk['trip_id'], trip.num, 
                                         j, i.id, 
                                         direction, i.point_start, i.point_stop,
                                         i.mps, i.odd, startpos, endpos])
                            j += 1
            query = query[:-1]
            if data:
                db_conn.run(query, data, commit=True)

        if GPStripstages and self.gps:
            query = ("INSERT INTO {}.{} (gpstripid,"
                     "stageid,"
                     "starttime,endtime,startepoch,endepoch,startpoint,"
                     "endpoint,duration,npoints,pointdistance,airdistance,"
                     "avgspeed,avgpointspeed,maxpointspeed,startlong,"
                     "startlat,endlong,endlat,"
                     "indirectness,"
                     "geom) VALUES ".format(schema, GPStripstages))
            data = []
            for trip_id, gpsc in self.gps.iteritems():
                for trip in gpsc.trips:
                    tk = trip.kwargs
                    query += ('(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,'
                              '%s,%s,%s,%s,%s,%s),')
                    s = trip.stats
                    airdistance = float(s['dist_net_m']) * 0.000621371
                    data.extend([tk['trip_id'], trip.num, 
                                 trip.local_time_start,
                                 trip.local_time_end, trip.time_start, 
                                 trip.time_stop, trip.point_start, trip.point_stop,
                                 s['time_GPS_min'], len(trip.gps.points), 
                                 s['dist_GPS_mi'], airdistance, 
                                 s['spd_avg_seg_mph'], s['spd_avg_GPS_mph'],
                                 s['spd_max_GPS_mph'], 
                                 trip.lng_start, trip.lat_start, trip.lng_stop, 
                                 trip.lat_stop,
                                 s['indirectness'], trip.gps.geom]) # 21
            query = query[:-1]
            if data:
                db_conn.run(query, data, commit=True)
                
                
class GPSData(object):
    """Base class for raw GPS data formats.
    
    GPSData(format)
    
    Parameters
    ----
    format : a format module
    
    Methods
    ----
    test : Attempts to scan GPS data using kwargs and return a data snippet and
        summary statistics for troubleshooting.
    load_csv : Attempts to read GPS data using kwargs and return a GPSCollection
        object.
    
    """

    def __init__(self, format, **kwargs):
        def option(attribute):
            try:
                return format.__dict__[attribute]
            except KeyError:
                return None

        self.data_type = format.data_type
        self.srid = format.srid
        self.col_x = format.col_x
        self.col_y = format.col_y
        self.col_time = format.col_time
        self.col_speed = format.col_speed
        self.col_altitude = format.col_altitude
        self.timezone = ustz.zones[format.timezone]
        self.col_date = option('col_date')
        if self.col_date:
            self.time_fmt = format.date_fmt + ' ' + format.time_fmt
        else:
            self.time_fmt = format.time_fmt
        self.integer_datetime = option('integer_datetime')
        if self.integer_datetime:
            self.time_fmt = '%y-%m-%d %H:%M:%S'
        self.fieldnames = option('fieldnames')
        self.delimiter = option('delimiter')
        self.interval = option('interval')
        self.day_hour_start = option('day_hour_start')
        self.col_id_point = option('col_id_point')
        self.col_x_dir = option('col_x_dir')
        self.col_lat = option('col_lat')
        self.col_lng = option('col_lng')
        self.col_unit_id = option('col_unit_id')
        self.col_phase = option('col_phase')
        self.col_deployment = option('col_deployment')
        self.col_geom = option('col_geom')
        self.col_unit_id = option('col_unit_id')
        self.col_deployment = option('col_deployment')
        self.col_accel = option('col_accel')
        self.col_phase_id = option('col_phase_id')
        self.col_heading = option('col_heading')
        self.kwargs = kwargs

    def load_csv(self, path, local_timezone = 'UTC', test = False, 
                 time_start = None, time_stop = None, astable = False):
        """Load GPS data from text file and return a GPSCollection.
        
        path : full path to text file
        local_timezone :  Pacific, Mountain, Central, Eastern, UTC [UTC]
            This is still a work in progress. Currently, only four US time zones
            plus UTC are supported. See important notes in us_time_zones.py 
            about anomalies around DST changeover dates.
        test : whether to run in test mode [False]
        astable : False or kwarg to divide points by id 
        kwargs : unit_id, group_id, trip_id
        
        """

        def convert_integer_time(int_time):
            """Convert 24h integer time to hh:mm:ss string time
            
            int_time : [hh][mm][ss]
            
            Examples
            ----
            3558 = 00:35:58
            153201 = 15:32:01
            
            TODO what happens for exact midnight (blank) time?
            
            """
            h = int(int_time / 10000)
            m = int(int_time / 100) - int(int_time / 10000) * 100
            s = int(int_time) - int(int_time / 100) * 100
            time_str = '{h}:{m}:{s}'.format(h=h, m=m, s=s)
            return time_str

        def convert_integer_date(int_date):
            """Convert integer date to yyyy-mm-dd string date
            
             int_date : d[mm][yy]
            
            Examples
            ----
            10511 = 2011-05-01
            300411 = 2011-04-30
            """
            d = int(int_date / 10000)
            m = int(int_date / 100) - int(int_date / 10000) * 100
            y = int(int_date) - int(int_date / 100) * 100
            date_str = '{y}-{m}-{d}'.format(y=y, m=m, d=d)
            return date_str

        def _time_conv(gps_datetime, local_timezone):
            """Accept GPS datetime and return utc timestamp & local datetime."""
            utc_timestamp = calendar.timegm(gps_datetime.utctimetuple())
            local_datetime = datetime.fromtimestamp(utc_timestamp, local_timezone)
            return (utc_timestamp, local_datetime)

        try:
            f = open(path)
        except IOError:
            err = '{f} is in use or does not exist'.format(f=path)
            raise IOError(err)

        if self.fieldnames:
            data = csv.DictReader(f, delimiter=self.delimiter, fieldnames=self.fieldnames)
        else:
            data = csv.DictReader(f, delimiter=self.delimiter)
        if test:
            return data
        if astable:
            gtable = {}
        if time_start:
            deploy_start = time_start
        else:
            deploy_start = float('-inf')
        if time_stop:
            deploy_stop = time_stop
        else:
            deploy_stop = float('+inf')
        id_point_auto = 1
        floor = math.floor
        for row in data:
            if astable and row[astable] not in gtable:
                g = pg.GPSCollection(srid=self.srid, srid_gps=self.srid, interval=self.interval, deploy_start=time_start)
                gtable[row[astable]] = g
            try:
                x = float(row[self.col_x])
            except KeyError:
                raise KeyError('X column {col_x} not found'.format(col_x=self.col_x))

            try:
                y = float(row[self.col_y])
            except KeyError:
                raise KeyError('Y column {col_y} not found'.format(col_y=self.col_y))

            d = None
            if self.col_date:
                try:
                    d = row[self.col_date]
                    if self.integer_datetime:
                        d = convert_integer_date(int(d))
                except KeyError:
                    raise KeyError('Date column {col_date} not found'.format(col_date=self.col_date))

            try:
                t = row[self.col_time]
                if self.integer_datetime:
                    t = convert_integer_time(int(t))
                if d:
                    t = d + ' ' + t
                dt = datetime.strptime(t, self.time_fmt).replace(tzinfo=self.timezone)
            except KeyError:
                raise KeyError('Time column {col_time} not found'.format(col_time=self.col_time))

            if self.col_x_dir:
                try:
                    x_dir = row[self.col_x_dir]
                    if x_dir == 'W' or x_dir == 'w':
                        x = -x
                except KeyError:
                    raise KeyError('X direction column {col_date} not found'.format(col_x_dir=self.col_x_dir))

            if self.col_id_point:
                id_point = row[self.col_id_point]
            else:
                id_point = id_point_auto
                id_point_auto += 1
            try:
                col_speed, units = self.col_speed
                if units == 'mi/hr':
                    speed = float(row[self.col_speed[0]]) * 0.44704
                elif units == 'km/hr':
                    speed = float(row[self.col_speed[0]]) * 0.27778
                elif units == 'm/s':
                    speed = float(row[self.col_speed[0]])
            except TypeError:
                speed = None

            try:
                altitude = row[self.col_altitude]
            except KeyError:
                altitude = None

            heading = row.get('heading', None)
            if heading:
                heading = int(heading)
            utc_timestamp, local_datetime = _time_conv(dt, ustz.zones[local_timezone])
            if utc_timestamp >= deploy_start and utc_timestamp <= deploy_stop:
                try:
                    travelday = int((utc_timestamp - deploy_start) / 86400.0) + 1
                except OverflowError:
                    travelday = 0

                p = pg.GPSPoint(x, y, utc_timestamp, id_point=id_point, s=speed, local_datetime=local_datetime, z=altitude, lng=x, lat=y, travelday=travelday, heading=heading)
                g.points.append(p)

        return gtable

    def test(self, data_type, db_info = None, unit_id = None, deployment = None):
        """Scan GPS data, and return a data snippet and summary stats."""
        if data_type == 'text':
            data = self.load_csv(path, test=True)
            n = 0
            print '{v} columns found'.format(v=len(data.fieldnames))
            print ''
            try:
                pprint(data.next())
                n += 1
            except StopIteration:
                pass

            for row in data:
                n += 1
                if any((val in (None, '') for val in row.itervalues())):
                    print row
                    err = 'Problem found in data row {r} shown above'.format(n)
                    raise Error(err)

            print 'Scanned {n} rows of data successfully.'.format(n=n)
        elif data_type == 'psql':
            data = self.load_psql(db_info, unit_id=unit_id, deployment=deployment, test=True)
            n = 0
            print '{v} columns found'.format(v=len(data[0]))
            print ''
            while n < len(data):
                pprint(data[n])
                n += 1
                raw_input('Next?')

        else:
            print 'Data type {} not supported.'.format(data_type)          