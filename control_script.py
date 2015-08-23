"""Load GPS points, map match, and ouptut to psql and gmaps.

2013-12-02 v01 Initial version sent to LCOG
2014-02-04 v02 Link count was showing from nodes instead of links (display only).
               Fixed gmap output bug when single trip to match.
               Fixed gmap bug not respecting SRID.
               Added db_conn.close() at end.
               
A lightweight version of the ASAP suite, this script loads GPS points
from a PostGIS enabled PostgreSQL (psql) database and matches to a
travel network using an algorithm based on Schuessler & Axhausen's (2009)
Multiple Hypothesis Technique
(https://edit.ethz.ch/ivt/vpl/publications/reports/ab568.pdf).

Basic Requirements

PostgreSQL server 9.1 and up
PostGIS 2.0 and up
Python 2.6.4 - 2.7.x (not 3.x)
psycopg2 python module (http://initd.org/psycopg/)
pg_tools_lite (provided with this package)

Other Requirements

pygmap python module (https://code.google.com/p/pygmaps/)
web_mapping (provided with this package)

Data Requirements

Travel network database with
1) a PostGIS geometry column (gist spatial index will help performance)
2) SRID should be planar geometry, preferably in meters
3) fields defining from-node and to-node connectivity
4) IMPORTANT (and weird): unique, symmetric link ids define a directed
   network. For example, if link 1 (fnode=0, tnode=1) defines one link,
   there should be a link -1 (fnode=1, tnode=0) with
   ST_Reverse(link1.geometry) unless link 1 is truly only available in
   the from-to direction. This may be relaxed in future versions.
5) length field in the same units as the geometry column

GPS table with
1) a PostGIS geometry column with same srid as travel network
2) fields defining x and y coordinates in same units as geometry
3) field with unique id for each point
4) field with unique id for each case (index to improve performance)
   Future versions may include tools for splitting trips and
   screening points, but for now, trips must be split and screened.

Basic Usage

Set variables defining database, travel network, and GPS point table.
Running this script will attempt to map match points for each case
and output results as routes and links to two PostgreSQL tables. You
can expect runtimes of 30-50ms per GPS point matched.
   
"""

import time
import pg_tools_lite as pg
import web_mapping as wm

def gmap(gpt, route, gmap_path, case, index, srid):
    """Write a google map html file with points and route."""
    if gmap_path.endswith('.html'):
        gmap_file = gmap_path
    else:
        gmap_file = gmap_path + str(case) + '.html'
    try:
        nextmap = str(gpt.cases[index + 1]) + '.html'
    except IndexError:
        nextmap = ''
    gpt.gps.get_geom2(db_conn)
    gpt.gps.transform2(4326, db_conn)
    extent = gpt.gps.get_extent()
    gmap = wm.Gmap((extent[1] + extent[3]) / 2.0,
                    (extent[0] + extent[2]) / 2.0, 17)
    points = gpt.gps.points
    gmap.addpoint(points[0].y, points[0].x, color='#009933', radius=3)
    gmap.addpoint(points[-1].y, points[-1].x, color='#330099', radius=3)
    for p in points[1:-1]:
        gmap.addpoint(p.y, p.x)
    if route:
        #route.get_geom(db_conn)
        lines = pg.LineCollection(srid=srid, lines=route.links[:])
        lines.get_geom(db_conn, vertices=True)
        lines.transform(4326, db_conn)
        for line in lines.lines:
            path = [(y,x) for x, y in line.vertices]
            gmap.addpath(path)
        gmap.draw(gmap_file, extent, title=case, next=nextmap)

#####################################
# User Settings

db_conn = pg.db_conn({'user': 'postgres', 'pass': 'asimov',
                      'db': 'fas', 'host': 'localhost'})  # psql connection details

# Information about psql point table
#
# Each sequence of points to be mapmatched must be uniquely identified
# by the col_caseid field. Geometry must be planar, not lat/long.
# All measurement variables will be assumed to match geometry units.
# Defaults are set for meters. ID column must indicate GPS track order.

table = 'orcycle_points'  # name of points table
schema = 'data'  # schema of data
col_id = 'gid'  # point id indicates order of GPS track
col_caseid = 'trip_id'  # unique caseid grouping points into route
                       # default max length is 50 characters
col_geom = 'geom'  # name of the geometry column
col_x = 'point_x'  # name of x (lng) value in same units as geom
col_y = 'point_y'  # name of y (lat) value in same unites as geom
srid = 2838  # SRID of data; must be planar
sql = 'gid >= 1849'  # (optional) body of SQL where clause to select points
          # Ex) 'phase=1' or '' to omit

# These settings control the mapmatch function. The defaults come
# closest to replicating Schuessler & Axhausen's Multiple Hypothesis
# Technique

net_table = 'wbnet2012'  # network table, same schema and srid as GPS table
net_sql = 'bik_rest !=1 '  # (optional) body of SQL where clause to select network links
                         # Ex) 'bik_rest != 1' or '' to omit
net_geom = 'geom'  # name of network geom column
net_id = 'gid'  # unique, directional link_id (positive FT, negative TF)
net_length = 'length_m'  # length attribute must be same units as SRID
net_fnode = 'fnode'  # from node attribute defining connectivity
net_tnode = 'tnode'  # to node attribute defining connectivity

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

psql_output = False  # Append route and link results to psql tables
                    # Tables will be created if not existing
tables_like = 'match'  # prefix for output tables

gmap_output = True  # Output html google maps in specified directory
                    # Files will be named by caseid unless path ends
                    # with .html, in which case file will be overwritten
                    # each time (good for manual mod debugging).
gmap_path = 'c:/data/gmap_test/'  # path to write gmap files
                                                    # if path ends in .html,
                                                    # file will overwrite each
                                                    # time


#####################################
# Program script
# comment raw_input prompt at end to run through all cases

if psql_output:
    # Create route table if not existing
    db_conn.create_table(schema, tables_like + '_routes',
                         [('caseid', 'varchar', 50),
                          ('startpos', 'double precision', 0),
                          ('endpos', 'double precision', 0),
                          ('aps', 'double precision', 0),
                          ('oddlinks', 'integer', 0),
                          ('distancelink', 'double precision', 0),
                          ('distancepoint', 'double precision', 0),
                          ('npoints', 'integer', 0)], srid=srid,
                          col_geom='geom')

    # Create link table if not existing
    db_conn.create_table(schema, tables_like + '_links',
                         [('caseid', 'varchar', 50),
                          ('featureorder', 'integer', 0),
                          ('featureid', 'integer', 0),
                          ('seltype', 'integer', 0),
                          ('direction', 'integer', 0),
                          ('startpoint', 'integer', 0),
                          ('endpoint', 'integer', 0),
                          ('startpos', 'double precision', 0),
                          ('endpos', 'double precision', 0),
                          ('mps', 'double precision', 0),
                          ('odd', 'integer', 0)])
    
print 'Loading Network...'
network = pg.Network(db_conn, schema + '.' + net_table, geom_attr=net_geom,
                     id_attr=net_id, length_attr=net_length,
                     fnode_attr=net_fnode, tnode_attr=net_tnode, sql=net_sql)
print '{} links loaded'.format(len(network.lines.lines))
gpt = pg.GPSPointTable(db_conn, table, srid, schema=schema, col_id=col_id,
                       col_caseid=col_caseid, col_geom=col_geom, col_x=col_x,
                       col_y=col_y, sql=sql)
print '{} cases found like {}'.format(len(gpt.cases), gpt.cases[0])
index = -1
for case in sorted(gpt.cases):
    index += 1
    print 'case:', case
    gpt.load_gps(db_conn, case)
    n = len(gpt.gps.points)
    print '{} points loaded'.format(n)
    t1 = time.time()
    m = pg.MapMatch(db_conn, gpt.gps, network)
    route = m.solve(max_candidates=max_candidates,
                    max_radius=max_radius,
                    initial_candidates=initial_candidates,
                    max_angle_flags=max_angle_flags,
                    max_dist_ratio=max_dist_ratio)
    if route:
        print 'possible route found'
        print ('avg point score: {:.1f} '
               '(PSU <25 acceptable match, <75 Schuessler & Axhausen)'
               .format(route.stats['aps']))
        print ('odd links: {} '
               '(3 is PSU and Schuessler & Axhausen max)'
               .format(route.stats['odd_links']))
        route.get_geom(db_conn, srid=srid)
        print 'links:', [link.id for link in route.links]
        if psql_output:
            route.write_psql(db_conn, schema, tables_like, case)
            print ('psql output written to {s}.{t}_routes and {s}.{t}_links'
                   .format(s=schema, t=tables_like))
        else:
            print 'psql output is turned off'
        if gmap_output:
            gmap(gpt, route, gmap_path, case, index, srid)
            print ('gmap output written to {}{}.html'
                   .format(gmap_path, case))
        else:
            print 'gmap output is turned off'
    else:
        print 'no route found with given parameters'
        
    print ('{:.1f} ms/pt (expected performance is 20-50 ms/pt)'
               .format((time.time() - t1) / n * 1000))
    #raw_input('Press enter to continue->')
db_conn.close()
