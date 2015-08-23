"""GPS raw data format module.

    GPSData(srid, col_x, col_y, col_time, time_fmt, timezone='UTC', 
            col_speed=None, col_altitude=None, col_date=None, date_fmt=None,
            integer_datetime = False, **kwargs)

    Parameters
    ----
    srid : Spatial Reference System Identifier numeric code
    col_x,col_y : Column names for x & y coordinates
    col_date,col_time : Column names for GPS time and optional 
        date fields
    time_fmt : Python format string for time
    timezone : Pacific, Mountain, Central, Eastern, UTC ['UTC']  
    dst_adjust : whether GPS recorded time already adjusts for dst [False] 
    date_fmt : Python format string for date (optional)
    integer_datetime : A flag for integer time and (optional) date format
       If True, time_fmt and date_fmt user settings will be ignored.
       integer time (e.g. 3558 = 00:35:58, 153201 = 15:32:01)
       integer date (e.g. 10511 = 2011-05-01, 300411 = 2011-04-30)
       Available as a setting on Globalsat DG-100. Anything else? Are there
           other variants out there?
    col_speed : (Column name for speed, units) where units can be 'mi/hr', 
        'km/hr', or 'm/s'
    kwargs : {delimiter, fieldnames, pg_login, interval, day_hour_start, 
              col_x_dir, col_id_point}
      
"""

device = 'orcycle'
format = 'orcycle'
data_type = 'csv'
delimiter = ','
srid = 4326
timezone = 'Pacific'  ## Pacific, Mountain, Central, Eastern, UTC
time_fmt = '%Y-%m-%d %H:%M:%S'  ## Python standard time format OR integer
#date_fmt = 'integer'  ## Python standard data format OR integer
#col_date = 'date'
col_time = 'time'
col_x = 'long'  
col_y = 'lat'
#col_x_dir = 'x_dir'  ## Optional
col_speed = ('speed', 'm/s')
col_altitude = 'alt'
#integer_datetime = True
col_id_point = 'id'
### optional parameters ###
#fieldnames = ['v1', 'y', 'v2', 'x', 'x_dir', 'time', 'date', 'speed', 
#              'heading', 'altitude', 'v3', 'v4']
interval = 1             
day_hour_start = 3
col_lng = 'long'
col_lat = 'lat'



    