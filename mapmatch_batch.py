"""Test mapmatch."""

import csv
import os
import time
import copy
import errno
from pprint import pprint
import asap_orcycle as asap
import stats
import pg_tools_orcycle as pg
import orcycle_gps
import web_mapping as wm
import json


def gmap(gps, route, gmap_path, case, trip_id, srid=4326, nextmap=None):
    """Write a google map html file with points and route."""
    # create gmap output dir if not exists
    try:
        os.makedirs(gmap_path + str(trip_id))
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    gmap_file = gmap_path + trip_id + '/' + str(case) + '.html'
    if nextmap:
        nextmap = str(nextmap) + '.html'
    extent = gps.get_extent()
    gmap = wm.Gmap((extent[1] + extent[3]) / 2.0,
                    (extent[0] + extent[2]) / 2.0, 17)
    points = gps.points
    start = True
    for p in points:
        if p.kwargs.get('trip', -999) == case:
            if start:
                gmap.addpoint(p.lat, p.lng, color='#009933', radius=3)
                start = False
            else:
                gmap.addpoint(p.lat, p.lng)
        elif p.kwargs.get('trip'):
            gmap.addpoint(p.lat, p.lng, color='#808080')
        else:
            gmap.addpoint(p.lat, p.lng, color='#FFFFFF')

    if route:
        lines = pg.LineCollection(srid=srid, lines=route.links[:])
        lines.get_geom(db_conn, vertices=True)
        lines.transform(4326, db_conn)
        for line in lines.lines:
            path = [(y,x) for x, y in line.vertices]
            gmap.addpath(path)
        gmap.draw(gmap_file, extent, title=case, next=nextmap)
##    print ('gmap output written to {}'
##                   .format(gmap_file, tripnum))

#### USER SETTINGS ####
with open('/Users/bblanc/OneDrive/_BikeAppProject/ORcycle_Analysis_Tool_Suite/mapMatch_Broach/userSettings.json') as data_file:
    userSettings = json.load(data_file)

print(userSettings)
working_directory = userSettings["working_directory"]
db_pass = userSettings["db_credentials"]["db_pass"]
db_user = userSettings["db_credentials"]["db_user"]
db_name = userSettings["db_credentials"]["db_name"]
db_schema = userSettings["db_credentials"]["db_schema"]

gps_path = 'data/coordsForMatch.csv'  # relative path to GPS points csv file
## Network settings
net_table = 'wbnet2012'  # network table
net_sql = 'bik_rest != 1 '  # (optional) body of SQL where clause to select network links
                            # Ex) 'bik_rest != 1' or '' to omit
## Map match parameters
max_candidates = 20  # [40] number of unique route candidates to track
                     # lower increases speed at expense of accuracy
                     # 20-40 seems to be a reasonable range

max_radius = 1000  # [1000] maximum radius to search for initial link
                   # candidates. Little downside to larger search
                   # radius and can help if trip start is noisy.

initial_candidates = 25  # [25] target number of initial candidates.
                         # Radius will expand until either this number
                         # or max_radius reached

max_angle_flags = 999  # [2] Number of angle flags (heading perpendicular to
                     # link) before link abandoned. Set to 999 to ignore
                     # headings. Angle flags can be problematic if signal
                     # is noisy. Worth experimenting.

max_dist_ratio = 1.05  # [1.05] Max GPS distance/link distance before link
                       # abandoned. Has only small effect in our tests.
                       # Reasonable range is 1.0-1.25
## Output settings
gmap_output = True  # Output html google maps in specified directory
                    # Files will be named by caseid unless path ends
                    # with .html, in which case file will be overwritten
                    # each time (good for manual mod debugging).

gmap_path = ('results/google_gmaps/')  # relative path to write gmap files
                                # if path ends in .html,
                                # file will overwrite each
                                # time

psql_output = True  # Append route and link results to existing psql tables

db_conn = pg.db_conn({'user': db_user, 'pass': db_pass,
                      'db': db_name})

os.chdir(working_directory)

# check to be sure data accessible
try:
    open(gps_path)
except IOError:
    print 'data location {} not accessible'.format(working_directory + gps_path)
    raise IOError

if gmap_output:
    # create gmap output dir if not exists
    try:
        os.makedirs(gmap_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    print 'gmap output is ON'
    print 'gmaps will go to', gmap_path
else:
    print 'gmap output is OFF'
if psql_output:
    print 'psql output is ON'
else:
    print 'psql output is OFF'
if net_sql:
    print 'network restrictions in effect: {}'.format(net_sql)

print 'Loading network...'
network = pg.Network(db_conn, db_schema + '.' + net_table,
                            length_attr='length_m', id_attr='psuid',
                            sql=net_sql, geom_attr='geom')
print 'Loading data...'
deployments = asap.Deployment(gps=(gps_path, orcycle_gps),
                              astable_gps='trip_id')

index = 1
trip_id = None
t1 = time.time()
for trip_id in sorted(deployments.gps):
    print ('{}, {} of {}, {:.1f} minutes so far'
           .format(trip_id, index, len(deployments.gps),
                   (time.time() - t1) / 60.0))
    gps = deployments.gps[trip_id]
    #print '{} points loaded'.format(len(gps.points))
    # create dummy deployment for just this trip
    d = asap.Deployment(gps=({trip_id: deployments.gps[trip_id]}, orcycle_gps))
    #print 'Reprojecting...'
    gps.get_geom2(db_conn)
    gps.transform2(2838, db_conn)
    #print 'Screening...'
    # For now, this is in same units as GPS file
    gps.screen(max_spd=40, max_accel=10, min_altitude=-100,
               max_altitude=9999, min_spd=0.45, debug=False)

    #print 'Masking...'
    gps.mask(db_conn)
    gps.get_accel()
    gps.split_trips(deploy_start=None, min_trip_pts=10,
                    max_time_gap=300,
                    num_near_points=30, search_radius=15,
                    min_bundle_density=10.0, min_bundle_ratio=2.0/3.0,
                    min_bundle_time=120, sequence_length=15,
                    max_time_zero_spd=120, zero_spd=0.01,
                    attribute_points=True, debug=False)
    print '{} trip(s) detected'.format(len(gps.trips))

    tripnum = 1
    for trip in gps.trips:
        trip.get_stats()
        # Hack to keep track of original and split trip nums
        trip.kwargs['trip_id'] = trip_id
        trip.num = tripnum
        trip.gps.get_geom2(db_conn)
        #print 'Map matching {} of {}...'.format(tripnum, len(gps.trips))
        t2 = time.time()
        m = pg.MapMatch(db_conn, trip.gps, network, check=False)
        route = m.solve(max_radius=max_radius, max_dist_ratio=max_dist_ratio,
                        initial_candidates=25, max_candidates=20,
                        u_turn=False, max_angle_flags=999, oneway_wt=0.0,
                        attribute_points=False, min_match_links=2,
                        max_route_aps=75, max_link_mps=75, max_odd_links=3,
                        fallback='gps', pos_tol=0.001)
        if route:
            print ('possible route found: aps={:.1f}, odd={}, {:.1f} ms/pt'
                   .format(route.stats['aps'], route.stats['odd_links'],
                           (time.time() - t2) / len(trip.gps.points) * 1000))
            #print ('avg point score: {:.1f} '
            #       '(PSU <25 acceptable match, <75 Schuessler & Axhausen)'
            #       .format(route.stats['aps']))
            #print ('odd links: {} '
            #       '(3 is PSU and Schuessler & Axhausen max)'
            #       .format(route.stats['odd_links']))
            #print ('{:.1f} ms/pt (expected performance is 20-50 ms/pt)'
            #       .format((time.time() - t2) / len(trip.gps.points) * 1000))
            # sttach to dummy stage
            trip.route = route
            trip.route.get_geom(db_conn, srid=2838)
        else:
            print 'no route found with given parameters'
            trip.route = None
        if gmap_output:
            gps4326 = pg.GPSCollection(srid=2838,
                               points=[copy.copy(p) for p in d.gps[trip_id].points])
            gps4326.get_geom2(db_conn)
            gps4326.transform2(4326, db_conn)
            if tripnum == len(gps.trips):
                nextmap = None
            else:
                nextmap = tripnum + 1
            gmap(gps4326, route, gmap_path, tripnum, trip_id, srid=2838, nextmap=nextmap)
        tripnum += 1
    index += 1
    trip_id = None
    if psql_output:
        d.write_psql(db_conn, db_schema,
                     GPSroutes='orc_gpsroutes',
                     GPSlinks='orc_gpslinks',
                     GPStripstages='orc_gpstrips')

        print ''
