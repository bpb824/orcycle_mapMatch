"""PostGIS functions and classes for Python"""

# orcycle - limited version that splits trips and mapmatches only 

import math
import csv
import numbers
import types
import datetime
import os
import StringIO
from copy import deepcopy, copy
from pprint import pprint
import psycopg2 as pg
import psycopg2.extras as pg_extras
import stats


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class db_conn(object):
    """A PostgreSQL wrapper around psycopg2"""

    def __init__(self, login):
        """Initialize a new pg_conn instance

        login - {'user'='', 'pass'='', 'db':'', 'host':'loclahost'}

        """
        d = login
        if 'host' not in d:
            d['host'] = 'localhost'
        self.conn = pg.connect(host=d['host'], database=d['db'], user=d['user'], 
                               password=d['pass']) 
        self.cur = self.conn.cursor()
        self.dict_cur = self.conn.cursor(cursor_factory=
                                         pg_extras.RealDictCursor)
    def update(self, schema, table, col, set_val, where, where_val, 
               debug=False, commit=False):
        """Update column in table based on condition.
        
        col : name of column
        set_val : single value to set to if condition is True
        where : body of where clause with placeholders
          Ex) 'pid=%s AND phase=%s'
        where_val : ordered list of values for where clause
          Ex) ['101-1', 1]
        commit : if True, updates are committed immediately
        
        """
        
        query = ('UPDATE {}.{} SET {}=%s WHERE {};'
                 .format(schema, table, col, where))
        data = [set_val]
        data.extend(where_val)
        if not debug:
            self.run(query, data, commit=commit)
        else:
            print self.debug(query, data)
                 
        
    def create_table(self, schema, table, attributes, srid,
                     col_geom='the_geom', dim_geom=2):
        """Create a table with geometry column and gid primary key field
        
        attributes : list, excluding geometry field and gid primary key
          [('varname1', 'integer', 0), ('varname2', 'varchar', 5)...]
        
        """
        valid_types = set(['integer', 'bigint', 'char', 'varchar', 
                           'double precision', 'numeric'])
        # Check whether table exists already
        query = ('SELECT EXISTS(SELECT * from information_schema.tables where '
                 'table_schema=\'{}\' AND table_name=\'{}\');'
                 .format(schema, table))
        result = self.run(query)[0][0]
        if not result:
            query = 'CREATE TABLE {}.{} (gid serial primary key,'.format(schema, table) 
            for name, type, length in attributes:
                if type == 'char' or type == 'varchar':
                    query += '{} {}({}),'.format(name, type, length)
                elif type in valid_types:
                    query += '{} {},'.format(name, type, length)
                else:
                    print '{} type is not yet implemented'.format(type)
            if query[-1] == ',':
                query = query[:-1]
            query += ');'
            self.run(query)
            query = ("SELECT AddGeometryColumn('{}','{}','{}',{},'POINT',{});"
                     .format(schema, table, col_geom, srid, dim_geom))
            self.run(query, commit=True)
        else:
            print 'Warning: Table already exists! Create skipped.'
            
    def run(self, query, data=None, commit=False, dict_cursor=False):
        """Execute query and return result as fetchall list of tuples.

        query - string statement with %s insertion points for variables
        data - oredered tuple of variables for query
        commit : commit changes to database
        dictcursor : get result with DictCursor
          Ex) row = {'col1': value, 'col2': value,...}

        """
        if not dict_cursor:
            self.cur.execute(query, data)
            try:
                result = self.cur.fetchall()
            except pg.ProgrammingError:
                result = None
        else:
            self.dict_cur.execute(query, data)
            try:
                result = self.dict_cur.fetchall()
            except pg.ProgrammingError:
                result = None
        if commit:
            self.conn.commit()
        return result
        
    def debug(self, query, data=None):
        """Return cursor query as it would be sent to database."""
        return self.cur.mogrify(query, data)
        
    def close(self, commit=True):
        """Close the connection."""
        self.conn.commit()
        self.cur.close()
        self.conn.close()
    
        
class Shape(object):
    """Base class for spatial geometry objects."""
    pass


class ShapeCollection(object):
    """Base class for shape collections."""
    pass


class Line(Shape):
    """A line feature.

    Line(wkt=None, geom=None, length=None, **kwargs)

    wkt : well-known text representation of linestring
    geom : psql geometry (inlcudes srid)
    length : geometric length of line segment
    kwargs : {}

    Either wkt or geom must be supplied. Other will be calculated if missing.
    """

    def __init__(self, id, wkt=None, geom=None, length=None, **kwargs):
        self.id = id
        self.wkt = wkt
        self.geom = geom
        self.length = length
        self.kwargs = kwargs

    def substring(self, pos0, pos1, db_conn):
        """Return geom, wkt, length of substring from pos0 to pos1.
        
        pos0 : fraction along line to use as starting point
        pos1 : fraction along line to use as ending point
        

        """
        # query = ('SELECT ST_Line_Substring(%s,%s,%s), '
                 # 'ST_AsText(ST_Line_Substring(%s,%s,%s));')
        query = ("WITH stub as (select ST_Line_Substring(geometry(%s),%s,%s)) "
                 "SELECT (SELECT * from stub),ST_AsText((SELECT * from stub)),"
                 "ST_Length((SELECT * from stub));")
        data = [self.geom, pos0, pos1]
        geom, wkt, length = db_conn.run(query, data)[0]
        return geom, wkt, length
        
class Link(Line):
    """MapMatch candidate link.

    id : unique link id
      Links with id set to None are used as off-network links by MapMatch
    fnode : from node id
    tnode : to node id
    points : related point features
    kwargs : {transit_lines, pos0, pos1}
      pos0 & pos1 here refer to final positions along link
    #profile = 10% of mapmatch
    
    """
    def __init__(self, id, fnode, tnode, wkt=None, geom=None, length=None, 
                 points=None, **kwargs):
        super(Link, self).__init__(id=id, wkt=wkt, geom=geom, length=length,
                                   **kwargs)
        self.fnode = fnode
        self.tnode = tnode
        if points is None:
            self.points = []
        else:
            self.points = points
        self.odd = 0  ## set by Route.get_stats()
        self.flg = 0  ## set by MapMatch.check_route()
        self.pos = 0.0  ## current position along link

class LineCollection(ShapeCollection): 
    """A collection of line features with a common spatial reference."""
    def __init__(self, srid=4326, lines=None, **kwargs): 
        """Initialize a LineCollection. 

        srid - Spatial Reference System Identifier [integer] 
        lines - optional list of Line objects [list] 
        kwargs - {}

        """
        self.srid = srid 
        if lines is None: 
            self.lines = []
        else: 
            self.lines = lines 
        self.kwargs = kwargs
        
    def get_geom(self, db_conn, vertices=False, collect_geom=False): 
        """Create geometry attributes for xy data.

        Either geom or wkt attribute must exist. Missing attribute will be
        calculated.
        
        vertices : If True, make a list of vertices [(x, y),...(X, Y)]

        """

        # todo : Port faster get_geom2 from PointCollection
        for i in self.lines:
            if i.wkt:
                query = ("SELECT ST_GeomFromText(%s, %s);")
                data = [i.wkt, self.srid]
                i.geom = db_conn.run(query, data)[0][0]
            elif i.geom:
                query = ("SELECT ST_AsText(geometry(%s));")
                data = [i.geom]
                i.wkt = db_conn.run(query, data)[0][0]
            else:
                err = 'Neither geom nor wkt found for {}'.format(i.id)
                raise Error(err)
        if vertices:
            for i in self.lines:
                v = [(float(x), float(y)) for x,y in [x.split(' ') for \
                  x in i.wkt[len('LINESTRING') + 1:-1].split(',')]]
                i.vertices = v
        else:
            i.vertices = None
        if collect_geom:
            linearray = []
            query = "SELECT ST_Collect(%s),ST_AsText(ST_Collect(%s));"
            for i in self.lines:
                linearray.append(i.geom)
            data = [linearray, linearray]
            # print db_conn.debug(query, data)
            self.geom, self.wkt = db_conn.run(query, data)[0]
            # return self.geom
                
    def transform(self, new_srid, db_conn):
        """Transform (re-project) geometry and wkt coordinates of a
        LineCollection

        srid - SRID to transform collection to [int]
        pg_login - database connection info [dict]

        """
        for i in self.lines:
            # todo : test whether single query faster
            query = "SELECT ST_Transform(geometry(%s), %s);"
            data = [i.geom, new_srid]
            # print db_conn.cur.mogrify(query, data)
            i.geom = db_conn.run(query, data)[0][0]
            query = ("SELECT ST_AsText(geometry(%s));")
            data = [i.geom]
            i.wkt = db_conn.run(query, data)[0][0]
            if i.vertices:
                v = [(float(x), float(y)) for x,y in [x.split(' ') for \
                  x in i.wkt[len('LINESTRING') + 1:-1].split(',')]]
                i.vertices = v
        self.srid = new_srid

    def write_shp(self, name, data_dir, attributes=None):
        """Use PyShp module to write LineCollection to a shapefile.

        attributes : [(name, type, fmt),...] for shapefile.Writer.field

        """

        def write_prj(srid, location):
            """Generate prj file based on SRID"""
            # TODO: lookup projections from file or database on the fly
            f = open(location + '.prj', 'w')
            if self.srid == 4326:
                prj = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
            elif self.srid == 2838:
                prj = 'PROJCS["NAD83(HARN) / Oregon North",GEOGCS["NAD83(HARN)",DATUM["D_North_American_1983_HARN",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",46],PARAMETER["standard_parallel_2",44.33333333333334],PARAMETER["latitude_of_origin",43.66666666666666],PARAMETER["central_meridian",-120.5],PARAMETER["false_easting",2500000],PARAMETER["false_northing",0],UNIT["Meter",1]]'
            elif self.srid == 2913:
                prj = 'PROJCS["NAD83(HARN) / Oregon North (ft)",GEOGCS["NAD83(HARN)",DATUM["D_North_American_1983_HARN",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",46],PARAMETER["standard_parallel_2",44.33333333333334],PARAMETER["latitude_of_origin",43.66666666666666],PARAMETER["central_meridian",-120.5],PARAMETER["false_easting",8202099.738],PARAMETER["false_northing",0],UNIT["Foot",0.3048]]'
            f.write(prj)
            f.close()
        if attributes is None:
            attributes = []
        else:
            attributes = attributes
        write_prj(self.srid, data_dir + '/' + name)
        w = shapefile.Writer(shapeType=3)
        w.autoBalance = 1  # Forces every shape to have a record
        w.field('id', 'N', decimal=0)
        for attr, type, fmt in attributes:
            if type == 'N':
                w.field(attr, type, decimal=fmt)
            else:
                w.field(attr, type, fmt)
        for i in self.lines:
            vertices = [x.split(' ') for x in i.wkt[11:-1].split(',')]
            vertices = [[float(x), float(y)] for x, y in vertices]
            w.line(parts=[vertices])
            values = []
            for a in attributes:
                try:
                    values.append(i.__dict__[a[0]])
                except KeyError:
                    try:
                        values.append(i.kwargs[a[0]])
                    except KeyError:
                        values.append(-9)
            w.record(i.id, *values)
        w.save(data_dir + '/' + name)                

class Point(Shape): 
    """A point feature."""
    def __init__(self, x, y, z=None, id_point=None, geom=None, lat=None,
                 lng=None, **kwargs): 
        """Initialize a Point object. 

        x - x coordinate of point [str] 
        y - y coordinate of point [str] 
    	z - z coordinate of point
        geom - optional PostGIS geometry with embedded spatial reference [str]
        lat : latitude of point in srid=4326 (WGS84)
        lng : longitude of point in srid=4326 (WGS84)
        kwargs : {}

        """
        self.id = id_point 
        self.x = x
        self.y = y
        if z is None:
            self.z = z
        else:
            self.z = float(z)
        self.geom = geom
        self.lat = lat
        self.lng = lng
        self.kwargs = kwargs


### Mapmatch classes ###
class Network(object):
    """Travel network object

    A database representation of a network consisting of links and nodes.
    Node connectivity is assumed based on matching id's.

    Network(link_table, db_conn, geom_attr='the_geom', 
            id_attr='gid', fnode_attr='fnode', tnode_attr='tnode', 
            length_attr='length', **kwargs)

    link_table : schema.name of table storing network links
    db_conn : a db_conn connection object
    geom_attr : geometry column
    id_attr : unique link id column
    fnode_attr : name of from node column
    tnode_attr : name of to node column
    length_attr : name of link length column
    sql : optional sql WHERE clause (e.g. 'gid > 0')
    name : name for network
    kwargs : {transit_lines_attr, srid}

    todo : Enforce link uniqueness between node pairs.
    """
    # todo : re-write as LineCollection child?
    def __init__(self, db_conn, link_table, geom_attr='the_geom', 
                 id_attr='gid', fnode_attr='fnode', tnode_attr='tnode', 
                 length_attr='length', sql=None, name=None, **kwargs):
        self.table = link_table
        self.geom_attr = geom_attr
        self.id_attr = id_attr
        self.fnode_attr = fnode_attr
        self.tnode_attr = tnode_attr
        self.length_attr = length_attr
        self.kwargs = kwargs
        self.load(db_conn, sql, **kwargs)
        self.name = name

    def load(self, db_conn, sql, **kwargs):
        """Load Network from database into dict."""
        if not sql:
            sql = 'True'
        query = ("SELECT ST_SRID({}) FROM {} LIMIT 1;"
                 .format(self.geom_attr, self.table))
        srid = db_conn.run(query)[0][0]
        self.lines = LineCollection(srid=srid)
        query = ("SELECT {},{},{},{},{},ST_AsText({}) FROM {} WHERE {};"
                 .format(self.geom_attr, self.id_attr, self.fnode_attr, 
                         self.tnode_attr, self.length_attr, self.geom_attr,
                         self.table, sql))
        # print query
        # todo : Re-cast links as link collection (line collection?).
        result = db_conn.run(query)
        self.links = {}
        for geom, link_id, fnode, tnode, length, wkt in result:
        # todo : fix these links
            # if fnode == tnode:
                # print 'Warning: link {} fnode=tnode!'.format(link_id)
            if fnode not in self.links:
                self.links[fnode] = []
            link = Link(id=link_id, geom=geom, wkt=wkt, fnode=fnode, 
                        tnode=tnode, length=length)
            self.links[fnode].append(link)
            self.lines.lines.append(link)


def line_locate_point(db_conn, point, line):
    """Get position of Point along Line."""
    # todo : Handle missing geometry.
    query = "SELECT ST_Line_Locate_Point(%s,%s)"
    data = [line.geom, point.geom]
    pos = db_conn.run(query, data)[0][0]
    return pos

def line_interpolate_point(db_conn, line, pos=0.5):
    """Return a Point given line and position along.
    
    db_conn : PostgreSQL database connection
    line : Line object
    pos : position along line (0-1, default=0.5)
    
    """
    # todo : try ST_AsText instead of separate ST_X & ST_Y
    query = ("SELECT ST_X(ST_Line_Interpolate_Point(geometry(%s), %s)),"
             "ST_Y(ST_Line_Interpolate_Point(geometry(%s), %s));")
    data = [line.geom, pos, line.geom, pos]
    x, y = db_conn.run(query, data)[0]
    point = Point(x=x, y=y)
    return point
    
def link_angle_point(db_conn, p_prev, p, link, pos, ppos):
    """Return angle of point relative to current link"""
    def get_ab_link_angle(a, b, link):
        query = ('SELECT ST_Azimuth('
                 '(SELECT ST_Line_Interpolate_Point(%s, %s)),'
                 '(SELECT ST_Line_Interpolate_Point(%s, %s)))')
        data = [link.geom, a, link.geom, b]
        angle_link = (db_conn.run(query, data)[0][0] / 
                      (2 * math.pi) * 360)
        # link.angle = angle_link
        return angle_link
    query = 'SELECT ST_Azimuth(geometry(%s), geometry(%s))'
    data = [p_prev.geom, p.geom]
    try:
        angle_gps = (db_conn.run(query, data)[0][0] / 
                     (2 * math.pi) * 360)
        p.angle = angle_gps
    except TypeError:
        # points have identical geometry
        # raw_input('Warning!')
        return 0
    # Measure angle of link using gps track as "shadow segment" unless
    # results in zero-length segment
    # todo : consider behavior around turns; should check be restricted to
    #   points matched to same link?
    if pos > ppos:
        a = ppos
        b = pos
    elif pos < ppos:
        # against link movement should be > 90deg angle 
        a = pos
        b = ppos
    else:
        if pos > 0.01:
                a = pos - 0.01
                b = pos
        else:
            a = pos
            b = pos + 0.01
    angle_link = get_ab_link_angle(a, b, link)
    angle = abs(angle_gps - angle_link)
    if angle > 180:
        angle = 360 - angle
    return angle


class Route(object):
    """MapMatch candidate route (hypothesis)"""

    def __init__(self, links=None, linkset=None, score=0.0, point_dist=0.0,
                 points=None, initial=False, **kwargs):
        """Initialize a Route object

        d - dict containing 

        """
        # todo : replace d with keyword arguments
        if links is None:
            self.links = []
        else:
            self.links = links
        self.linkset = linkset
        self.score = score
        self.point_dist = point_dist
        if points is None:
            self.points = []
        else:
            self.points = points
        self.end = False
        self.angle_end = False
        self.angle_flags = 0
        self.kwargs = kwargs
        self.initial = initial  ## initial link conditions in effect

    def extend(self, db_conn, route, stitch_links=None):
        """Extend an existing route with a second route and 
        stitch links, if necessary.
        
        """
        if stitch_links:
            route.links.extend(stitch_links)
        self.links.extend(route.links[:])
        self.score += route.score
        self.points.extend(route.points[:])
        if 'pos0' not in self.__dict__:
            self.pos0 = route.pos0
        self.pos1 = route.pos1
        #print len(self.points), 'route points'
        self.get_stats(db_conn, adjust_geom=False)
        # pprint(self.stats)

    def score_point(self, db_conn, p):
        """Calculate score for given GPS point and match link."""
        link = self.links[-1]
        query = "SELECT ST_Distance(geometry(%s), geometry(%s))"
        data = [p.geom, link.geom]
        match_dst = db_conn.run(query, data)[0][0]
        self.score += match_dst
        link.points.append((p, match_dst))
        p.link = link.id

    def update(self, db_conn, p_prev, p, point_dist, max_angle_flags, 
               pos_tol=0.001, max_angle=85, max_dist_ratio=1.05):
        """Update Route with new GPS point.

        point_dist - distance from previous point along GPS trace
        pos_tol - tolerance value for position along link (default=0.001)
        """
        def check_for_end(p_prev, p):
            """Check to see if end of current match link has been reached."""
            link = self.links[-1]
            pos = line_locate_point(db_conn, p, link)
            if (link_angle_point(db_conn, p_prev, p, link, pos, 
                                 link.pos) > max_angle):
                self.angle_flags += 1
            else:
                self.angle_flags = 0
            # TODO: Formalize pos_tol.
            # pos_tol adjusts for bug in postgis ST_Line_Locate_Point
            # if pos > 0.0 + pos_tol:  # allows for "catch up" after early jump
            self.point_dist += point_dist
            if (pos < 1.0 - pos_tol and
                self.point_dist <= max_dist_ratio * link.length and
                self.angle_flags <= max_angle_flags):
                pass
            else:
                self.end = True
                if self.angle_flags > max_angle_flags:
                    self.angle_end = True
            link.pos = pos
            
        check_for_end(p_prev, p)

    def next_links(self, network, link, u_turn):
        """Return links connected to current end_node. 

        Route is valid if candidate link has not been used in route. If 
        candidate link results in u-turn, u-turns must be permitted. Note
        that u-turn check assumes unique link between any two nodes, a
        common network assumption enforced by Network class by default. A 
        u-turn is never allowed on the initial link of a segment.
        

        """

        candidates = []
        if link.tnode in network.links:
            for candidate in network.links[link.tnode]:
                if (candidate.id not in self.linkset and 
                   (candidate.tnode != link.fnode)):
                    candidates.append(candidate)
        # if link.tnode in network.links:
            # for candidate in network.links[link.tnode]:
                # if (candidate.tnode != link.fnode or 
                    # (u_turn and not self.initial)):
                    # candidates.append(candidate)
        else:
            pass
        return candidates
         
    def advance(self, db_conn, network, point, link, u_turn):
        """Spawn new candidate routes on forward links.
        
        #profile = 18% of mapmatch
        
        """
        candidate_links = self.next_links(network, link, u_turn)
        linkset = self.linkset
        routes_update = []
        for c_link in candidate_links:
            links_new = [Link(id=i.id, geom=i.geom, wkt=i.wkt, fnode=i.fnode,
                          tnode=i.tnode, length=i.length, 
                          points = [(p, s) for p, s in i.points],
                          **i.kwargs) for i in self.links]
            j = c_link
            i = Link(id=j.id, geom=j.geom, wkt=j.wkt, 
                     fnode=j.fnode, tnode=j.tnode, 
                     length=j.length, pos=0.0, seltype=1)
            links_new.append(i)
            linkset_new = linkset.copy()
            linkset_new.add(links_new[-1].id)
            # Create new Route object
            routes_update.append(Route(links=links_new, linkset=linkset_new,
                                       score=self.score))
            routes_update[-1].score_point(db_conn, point)
        if not candidate_links:
            # If dead end, freeze link at dead end
            routes_update.append(self)
            routes_update[-1].score_point(db_conn, point)
        return routes_update
    
    def get_geom(self, db_conn, srid=2838):
        """Set route geometry as MultiLinestring."""
        lines = LineCollection(srid, lines=self.links)
        lines.get_geom(db_conn, collect_geom=True)
        self.geom = lines.geom
        self.wkt = lines.wkt
        
    def get_stats(self, db_conn, attribute_points=True, 
                  max_link_mps=75, network=None, adjust_geom=True,
                  pos_tol=0.001):
        """Calculate route stats."""
        self.stats = {}
        npoints = sum([len(i.points) for i in self.links])
        try:
            self.stats['network'] = network.name
        except AttributeError:
            self.stats['network'] = None
        self.stats['aps'] = float(self.score) / npoints
        self.stats['matched_points'] = npoints  ## Reduced as needed by check_route. 
        mls = 0  ## Maximum link_mps; not part of Schuessler & Axhausen
        odd_links = 0
        if adjust_geom:
            # Store match positions of starting and ending links
            link0, link1 = (self.links[0], self.links[-1])
            # links.points list of tuples like [(p, match_dst),]
            self.pos0 = line_locate_point(db_conn, link0.points[0][0], link0)
            self.pos1 = line_locate_point(db_conn, link1.points[-1][0], link1)
            # We don't want to include last link if not reached.
            # todo : Consider shifting link1.points to new link1
            # print [i.id for i in self.links]
            # print self.pos0, self.pos1
            if link0 != link1 and self.pos1 <= 0.0 + pos_tol:
                self.pos1 = 1.0
                points = link1.points
                self.links = self.links[:-1]
                link1 = self.links[-1]
                link1.points.extend(points)
            # print [i.id for i in self.links]
            # print self.pos0, self.pos1
            # link0 length already adjusted during mapmatch
            # Adjust geometries of start & end lengths to reflect partial use.
            if link0 != link1:
                link0.geom, link0.wkt, link0.length = \
                  link0.substring(self.pos0, 1.0, db_conn)
                link1.geom, link1.wkt, link1.length = \
                  link1.substring(0.0, self.pos1, db_conn)
                # link1.length = self.pos1 * link1.length
            elif self.pos1 >= self.pos0:
                link0.geom, link0.wkt, link0.length = \
                  link0.substring(self.pos0, self.pos1, db_conn)
                # link0 already reflects pos0
                # try:
                    # link0.length = ((self.pos1 - self.pos0) * 
                                    # (link0.length / (1.0 - self.pos0)))
                # except ZeroDivisionError:
                    # link0.length = 0
            else:
                # GPS points are running opposite link; flip it!
                try:
                    candidate_links = network.links[link0.tnode]
                except TypeError:
                    candidate_links = []
                except KeyError:
                    candidate_links = []
                match = False
                for clink in candidate_links:
                    if clink.tnode == link0.fnode:
                        match = True
                        break
                if match:
                    length = link0.length / (1.0 - self.pos0)
                    link0 = Link(id=clink.id, fnode=clink.fnode, 
                                 tnode=clink.tnode,
                                 points=link0.points, geom=clink.geom, 
                                 wkt=clink.wkt, 
                                 length=link0.length / (1.0 - self.pos0),
                                 seltype=1)
                    p0 = self.pos0
                    p1 = self.pos1
                    self.pos0 = 1.0 - p0
                    self.pos1 = 1.0 - p1
                    
                    link0.geom, link0.wkt, link0.length = \
                      link0.substring(self.pos0, self.pos1, db_conn)
                    # link0.length = (self.pos1 - self.pos0) * link0.length
                    self.links[0] = link0
                else:
                    print 'Warning: Could not flip link {}!'.format(link0.id)
                # print [i.id for i in self.links]
                # print self.pos0, self.pos1
        link_dst = 0.0
        self.links[0].kwargs['pos0'] = self.pos0
        self.links[-1].kwargs['pos1'] = self.pos1
        for link in self.links:
            # link.raps = self.aps
            try:
                link.mps = min([s for (p, s) in link.points])
                if link.mps > mls:
                    mls = link.mps
                if link.mps > max_link_mps:
                    link.odd = 1
                    odd_links += 1
                link.point_start = link.points[0][0].id
                link.point_stop = link.points[-1][0].id
            except ValueError:
                # stitch link, no points assigned
                link.mps = None
                link.odd = None
                link.point_start = None
                link.point_stop = None
            link_dst += link.length
            
        self.stats['link_dst'] = link_dst
        self.stats['odd_links'] = odd_links
        self.stats['mls'] = mls
        # link = self.links[-1]
        # link.length_adj = link.pos1 * link.length
        link = self.links[0]
        link.length_adj = link.length
        if attribute_points is True:
            self.points = []
            # index = 0
            # self.pt_dst.reverse()
            for link in self.links:
                mps = link.mps
                for p, s in link.points:
                    self.points.append(p)
                    p.link = link.id
                    # todo : Grab this in MapMatch.solve()?
                    # p.pos = line_locate_point(db_conn, p, link)
                    p.score = s
                    # p.link_mps = mps
                    # p.rte_scor = self.stats['aps']
                    # p.pt_dst = self.pt_dst.pop()

    def write_csv(self, path='STRING', fieldnames_r=None, fieldnames_l=None, 
                  user_dict=None, none_val=-999):
        """Write Route details as csv.

        path : path to write csv (links file will be appended with _links)
          string - return text (e.g. to delayed write)
        fieldnames_r : optional ordered list of str, num, or datetime.datetime 
          fieldnames to write to route file
        fieldnames_l : optional ordered list of str, num, or datetime.datetime 
          fieldnames to write to link file
        user_dict : user-supplied dict to include in fieldnames search

        Note that for now top-level attribute dict keys named in fieldnames 
        must be unambiguous across all attribute dicts, or an error will be 
        raised.

        """
        def _parse(val):
            """Return True if attribute is a valid type."""
            if isinstance(val, types.BooleanType):
                return int(val)
            elif (isinstance(val, types.StringTypes) or
                isinstance(val, numbers.Number)):
                 return val
            #todo : handle individual time and date objects
            elif isinstance(val, datetime.datetime):
                return val.strftime('%Y-%m-%d %H:%M:%S')
            elif isinstance(val, types.NoneType):
                return none_val
            else:
                return 'INVALID'
                
        if fieldnames_r is None:
            fieldnames_r = ['hid', 'pid', 'phase', 'trip', 'stage',
                            'segment', 'pos0', 'pos1', 'npoints', 
                            'aps', 'odd_links', 'link_dst', 
                            'pt_dst']
        if fieldnames_l is None:
            fieldnames_l = ['pid', 'phase', 'trip', 'stage',
                            'segment', 'id', 'mps', 'point_start', 'point_stop',
                            'odd']
        if user_dict is None:
            user_dict = {}
        fieldmap = {}
        for attr in fieldnames_r:
            if attr in user_dict:
                val = _parse(user_dict[attr])
                if val != 'INVALID':
                    fieldmap[attr] = val
                else:
                    err = ('{0} is invalid for row {1}'
                           .format(attr, self.num))
                    raise TypeError(err)
            elif attr in self.__dict__:
                val = _parse(self.__dict__[attr])
                if val != 'INVALID':
                    fieldmap[attr] = val
                else:
                    err = ('{0} is invalid for row {1}'
                           .format(attr, self.num))
                    raise TypeError(err)
            else:
                for k, v in self.__dict__.iteritems():
                    if isinstance(v, types.DictType):
                        if attr in v:
                            if attr in fieldmap:
                                err = ('{0} is ambiguous for row {1}'
                                       .format(attr, self.num))
                                raise Error(err)
                            else:
                                val = _parse(v[attr])
                                if v != 'INVALID':
                                    fieldmap[attr] = v[attr]
                                else:
                                    err = ('{0} is invalid type for row {1}'
                                           .format(attr, self.num))
                                    raise TypeError(err)
            if attr not in fieldmap:
                err = ('{0} is not a valid attribute or sub-attribute'
                       .format(attr))
                raise AttributeError(err)
        row = {}
        for field in fieldnames_r:
            row[field] = fieldmap[field]
        if path == 'STRING':
            output_r = StringIO.StringIO()
            w = csv.DictWriter(output_r, fieldnames=fieldnames_r)
        else:
            try:
                f = open(path + '.csv')
                f.close()
                f = open(path + '.csv', 'ab')
                w = csv.DictWriter(f, fieldnames=fieldnames_r)
            except IOError:
                f = open(path + '.csv', 'ab')
                w = csv.DictWriter(f, fieldnames=fieldnames_r)
                w.writeheader()
        w.writerow(row)
        
        if path != 'STRING':
            f.close()
            try:
                f = open(path + '_links.csv')
                f.close()
                f = open(path + '_links.csv', 'ab')
                w = csv.DictWriter(f, fieldnames=fieldnames_l)
            except IOError:
                f = open(path + '_links.csv', 'ab')
                w = csv.DictWriter(f, fieldnames=fieldnames_l)
                w.writeheader()
        else:
            output_l = StringIO.StringIO()
            w = csv.DictWriter(output_l, fieldnames=fieldnames_l)
        fieldmap = {}
        rows = []
        for link in self.links:
            for attr in fieldnames_l:
                if attr in user_dict:
                    val = _parse(user_dict[attr])
                    if val != 'INVALID':
                        fieldmap[attr] = val
                    else:
                        err = ('{0} is invalid for row {1}'
                               .format(attr, user_dict['segment']))
                        raise TypeError(err)
                elif attr in link.__dict__:
                    val = _parse(link.__dict__[attr])
                    if val != 'INVALID':
                        fieldmap[attr] = val
                    else:
                        err = ('{0} is invalid for row {1}'
                               .format(attr, user_dict['segment']))
                        raise TypeError(err)
                else:
                    for k, v in link.__dict__.iteritems():
                        if isinstance(v, types.DictType):
                            if attr in v:
                                if attr in fieldmap:
                                    err = ('{0} is ambiguous for row {1}'
                                           .format(attr, user_dict['segment']))
                                    raise Error(err)
                                else:
                                    val = _parse(v[attr])
                                    if v != 'INVALID':
                                        fieldmap[attr] = v[attr]
                                    else:
                                        err = ('{0} is invalid type for row {1}'
                                               .format(attr, user_dict['segment']))
                                        raise TypeError(err)
                if attr not in fieldmap:
                    err = ('{0} is not a valid attribute or sub-attribute'
                    .format(attr))
                    raise AttributeError(err)
            row = {}
            for field in fieldnames_l:
                row[field] = fieldmap[field]
            rows.append(row)
        w.writerows(rows)
        if path != 'STRING':
            f.close()
        else:
            return output_r.getvalue(), output_l.getvalue()

        
class MapMatch():
    """A map matching problem and functions

    MapMatch(db_conn, segment, network, check=False)

    db_conn : db_conn connection object
    segment : Segment instance
    network : Network instance
    check : Whether to check if route is valid
    
    Based on the Multiple Hypthesis Technique, Schuessler & Axhausen (xxxx)

    """
    def __init__(self, db_conn, segment, network, check=False):
        # self.segment = segment
        self.points = segment.points
        # self.start = segment.point_start
        # self.end = segment.point_stop
        ##        self.max_candidates = max_candidates
        ##        self.u_turn = u_turn
        ##        self.max_angle_flags = max_angle_flags
        ##        self.oneway_wt = oneway_wt
        self.network = network
        self.routes = []
        self.db_conn = db_conn
        self.check = check

    def start_search(self, radius, max_radius, initial_candidates, 
                     pos_tol=0.001):
        """Find initial set of candidate links within radius of segment 
        start point.
        
        Initial candidate links must be within search radius and the relevant
        perpendicular distance to the start point must not be the link end.

        pos_tol : tolerance value for position along link (default=0.001)
        
        """
        # todo : re-write using distance-ordered list of start links?
        start_point = self.points[0]
        net = self.network
        o_candidates = set()
        while len(o_candidates) < initial_candidates and radius <= max_radius:
            query = ("SELECT {g},ST_AsText({g}),{i},{f},"
                     "{t},{e}, ST_Line_Locate_Point({g},%s) " 
                     "FROM {a} WHERE ST_DWithin(geometry(%s),"
                     "geometry({g}),%s) and ST_Line_Locate_Point({g},%s) < "
                     "1 - %s and bik_rest=0;"
                     .format(g=net.geom_attr, i=net.id_attr, f=net.fnode_attr, 
                             t=net.tnode_attr, e=net.length_attr, a=net.table))
            data = [start_point.geom, start_point.geom, radius, 
                    start_point.geom, pos_tol]
            result = self.db_conn.run(query, data)
            for geom, geom_text, gid, fnode, tnode, length, pos in result:
                # Don't add link if endponit nearest point
                if pos < 1.0 - pos_tol:
                    # Adjust length of initial link by start position
                    o_candidates.add((geom, geom_text, gid, fnode, tnode, 
                                      (1.0 - pos) * length))
            radius += 100
        return o_candidates

    def prune_candidates(self, max_candidates):
        """Remove candidate routes with worst (highest) scores."""
        cr = self.routes[:]
        cr.sort(key=lambda cr: cr.score)
        while len(cr) > max_candidates:
            cr.pop()
        self.routes = cr[:]

    def setup(self, db_conn, radius, max_radius, 
              initial_candidates, max_candidates):
        """Setup a map matching problem."""
        o_candidates = self.start_search(radius, max_radius, 
                                         initial_candidates)
        
        for geom_link, geom_text, gid_link, fnode, tnode, length in \
            o_candidates:
            # Note that length is adjusted in start_search() based on point
            # position.
            i = Link(id=gid_link, geom=geom_link, wkt=geom_text, fnode=fnode, 
                     tnode=tnode, length=length, pos=0.0, seltype=1)
            self.routes.append(Route(links=[i], linkset=set([gid_link]),
                                     initial=True))
            self.routes[-1].score_point(db_conn, self.points[0])
        # Pruning here may give initial point too much weight.
        # self.routes.prune_candidates(max_candidates)    

    def solve(self, radius=100, max_radius=1000, 
              initial_candidates=25, max_candidates=20, 
              u_turn=False, max_angle_flags=2, oneway_wt=0.0, 
              attribute_points=False, min_match_links=2, max_route_aps=75,
              max_link_mps=75, max_odd_links=3, fallback='gps', pos_tol=0.001,
              max_dist_ratio=1.05):
        """Solve a map matching problem
        start_point : GPSPoint instance
        radius : radius to search for initial link from start point, in 
          units of point data
        max_radius : maximum radius to search for initial link; overrides
          initial_candidates
        initial_candidates : requested minimum number of initial candidates 
          before search is stopped
        min_links - if candidate links < min_links, radius increased by 
          100m until candidates >= min_links or radius >= max_radius - 100
        max_candidates : maximum number of route candidates in the queue
        # moved to route.update() for some reason
        max_dist_ratio : maximum ratio of point dist to assigned link distance
                         before forced to join to following link
        u_turn : True if u-turns allowed
        max_angle_flags : maximum headings perpendicular to link before link is
                          abandoned
        oneway_wt : penalty, in measurement units of spatial reference system,
                    given to links traversed in restricted direction
        attribute_points : ??
        min_match_links : minimum number of links for a valid route
        max_route_aps : maximum average point distance for a valid route
        max_link_mps : minimum point score for a conforming link
        max_odd_links : maximum number of links not meeting max_link_mps in
          valid route
        fallback : what to return if no valid route found
          'gps' - Use GPS trace and attribute point distance. Set links=None.
          
        """
        def check_route(db_conn, route, fallback, kill_route=True):
            """Validate best candidate route, and try to complete bad routes.

            kill_route : If True, invalid route is removed entirely. If
              False, route is flagged and odd links are replaced using
              fallback method from MapMatch.solve().

            """
            route.get_stats(db_conn, max_link_mps, 
                            attribute_points, network=self.network)

            if (len(route.links) < min_match_links or
                route.aps > max_route_aps or self.odd_links > max_odd_links):
                route.flg = 1
                if kill_route:
                    self.segment.route = None
                else:
                    # Find first problem spot, and try to re-match from there
                    for link in route.links[1:]:
                        # todo : better handle initial link problems
                        if link.odd:
                            for p in link.points:
                                p.link = None
                            route.stats['matched_points'] -= len(link.points)
            else:
                pass

        db_conn = self.db_conn
        seg_points = self.points[:]
        seg_points.reverse()
        p = seg_points.pop()  ## get first GPS point
        self.setup(db_conn, radius, max_radius, initial_candidates, 
                   max_candidates) 
        # get_distance = dist_pythag  ## localize function 
        # routes_update = []

        # advance through points
        debug_id = 26374
        debug = False
        while seg_points:
            if not self.routes:
                # route list is empty
                break
            # if p.id >= id_debug:
                # debug=True
            p_prev = p
            p = seg_points.pop()
            point_dist = ((p.x - p_prev.x) ** 2 + (p.y - p_prev.y) ** 2) ** 0.5
            # print self.routes[0].point_dist_sum
            # print '+', point_dist
            routes_update = []
            if debug and p.id >= debug_id:
                print '**', p.id
            for route in self.routes:
                # for link in route.links:
                    # if link.kwargs['uturn']:
                        # print link.id, link.kwargs
                        # raw_input()
                if debug and p.id >= debug_id:
                    print route.score
                    print route.point_dist
                    print route.links[-1].length
                    print [i.id for i in route.links]
                # try:
                    # print [(i.id, i.pos, i.pos0, i.pos1, i.uturn) for i in route.links]
                # except AttributeError:
                    # print [(i.id, i.pos, i.uturn) for i in route.links]
                link = route.links[-1]
                route.update(db_conn, p_prev, p, point_dist,
                             max_angle_flags, max_dist_ratio=max_dist_ratio)
                if not route.end:
                    routes_update.append(route)
                    routes_update[-1].score_point(db_conn, p)
                else:
                    routes_update.extend(route.advance(db_conn, self.network, 
                                                       p, link, u_turn))
            self.routes = routes_update[:]
            self.prune_candidates(max_candidates)
            # for route in self.routes:
                # print route.links[-1].id
            if debug and p.id >= debug_id:
                raw_input('?')
        if len(self.routes) > 0:
            route = self.routes[0] 
            if self.check:
                check_route(db_conn, route, fallback)
            else:
                route.get_stats(db_conn,
                                max_link_mps=max_link_mps, 
                                attribute_points=attribute_points,
                                network=self.network)
                # Catch routes matched only to single link vertex
                if len(route.links) == 1 and route.pos0 == route.pos1:
                    route = None
        else:
            # print 'Warning: No route found.'
            route = None
        return route
        
        
class MapMatch2():
    """A map matching problem and functions, shortest path based.
    
    buffer : Only network links completely inside buffer around gps
      will be considered in solve.
      
    Uses the Network.solve2 algorithm modified to restruct search
    to a subset of links.
    """
    def __init__(self, db_conn, segment, network, buffer=20, check=False):
        if not segment.gps.geom:
            segment.gps.get_geom2
        # Create a subset of network links within buffer.
        query = ("WITH buffer as (SELECT ST_Buffer(%s, %s)) "
                 "SELECT {} from {} WHERE ST_Within({},(SELECT * FROM buffer)"
                 .format(network.id_attr, network.table, network.geom_attr))
        data = [segment.gps.geom, buffer]
        result = db_conn.run(query, data)
        self.blinks = set(result[0])
        self.segment = segment
        self.points = segment.gps.points
        self.start = segment.point_start
        self.end = segment.point_stop
        self.network = network
        self.routes = []
        self.db_conn = db_conn
        self.check = check
    
    def solve(self, db_conn, origin, destinat, cost_params, pg_table=None, 
              phi=3.14):
        """Solve a buffered, shortest-path mapmatch problem,
        
        Based on Dalumpines & Scott (2011) but uses link selection instead
        of route barriers.
        
        """
        def nearest_link(db_conn, point, max_candidates=100):
            """Return nearest link and its doppelganger to given point and 
            point locations along links.
            
            """
            query = ("with index_query as (select st_distance({}, %s) as dist,"
                     "{},{},{},{} from {} order by {} <#> %s "
                     "limit %s) "
                     "select dist,the_geom,ST_Reverse(the_geom),gid,fnode,"
                     "tnode from index_query where dist order by dist limit 1;"
                     .format(self.geom_attr, self.geom_attr, self.id_attr, 
                             self.fnode_attr, self.tnode_attr, self.table, 
                             self.geom_attr))
            data = [point.geom, point.geom, max_candidates]
            # print db_conn.debug(query, data)
            result = db_conn.run(query, data)[0]
            match_dst, geom, geomd, id, fnode, tnode = result
            link = Link(id, fnode, tnode, geom=geom)
            linkd = Link(-1 * id, tnode, fnode, geom=geomd)
            link.pos = line_locate_point(db_conn, point, link)
            linkd.pos = 1.0 - link.pos
            return match_dst, {link.id: link, linkd.id: linkd}
            
        links = self.links
        ## containers for queues and previous link table
        queue = dict()      ## {link_id:util}
        queue_sort = list() ## [(link_id,util)]
        bucket = list()
        prev = dict()       ## {link_id:previous link}
        ## Find nearest network link to origin and destination points
        match_dst, olinks = nearest_link(db_conn, origin)
        match_dst, dlinks = nearest_link(db_conn, destinat)
        # Adjust origin/destination links to represent only relevant section
        for id, link in olinks.iteritems():
            link.geom, link.wkt, link.length = link.substring(link.pos, 1.0, 
                                                              db_conn)
        for id, link in dlinks.iteritems():
            link.geom, link.wkt, link.length = link.substring(0.0, link.pos, 
                                                              db_conn)
            # Always include dlinks in buffered area.
            self.blinks.add(id)
        ## initialize the queue with origin links
        for id, link in olinks.iteritems():
            cost = calc_cost(link.__dict__)
            queue[link] = cost
            bucket.append((cost, link))
            prev[link] = None  ## no preceding link for origin links
        
        # (3) SELECT LINK WITH LOWEST KNOWN COST FROM ORIGIN TO SCAN NEXT
        iter_num = 0
        topsort = 0
        labeled = dict() ## holds true costs
        ## main shortest path loop
        while True:
            iter_num+= 1
            ####### priority queue (queue_sort) #######
            if queue_sort and topsort == 0:
                pass
            elif topsort == 1:

                queue_sort.sort(reverse = True)
                topsort = 0
            elif bucket:

                bucket.sort(reverse = True)
                queue_size = len(bucket)
                p_queue_size = int(math.ceil(queue_size/phi))
                queue_sort = bucket[queue_size-p_queue_size:]
                bucket = bucket[:queue_size-p_queue_size]
                max_queue = queue_sort[0][0]    # max cost in queue_sort
            else:
                ## minimum spanning tree complete, but no route found
                print 'No route found from {o} to {d}.'.format(o=origin,
                                                               d=destinat)
                break
            ###############################
            ## grab the next link from the bottom of the (reversed)
            ## priority queue
            cost_i, link_i = queue_sort.pop()
            if link_i not in labeled:
                ## this link is now permanently labeled
                labeled[link_i] = cost_i

                ## TERMINATION CHECK
                if link_i.id in dlinks:
                    break

                # (4) UPDATE ("RELAX") COSTS ON LINKS CONNECTING TO CURRENT LINK
                # ??????if link_i in links:
                fnode_j = link_i.tnode
                for link_j in links[fnode_j]:
                    ## if link is destination link, use modified link
                    if link_j.id in dlinks:
                        link_j = dlinks[link_j.id]
                    ## If link already labeled or outside buffered area, 
                    ## ignore it.
                    if link_j.id not in self.blinks or link_j in labeled:
                        pass
                    ## if link in queue, update if lower cost path just found
                    elif link_j in queue:
                        cost_j = cost_i + calc_cost(link_j.__dict__)
                        if cost_j < queue[link_j]:
                            queue[link_j] = cost_j
                            prev[link_j] = link_i
                            if cost_j < max_queue:
                                queue_sort.append((cost_j, link_j))
                                topsort = 1
                            else:
                                bucket.append((cost_j, link_j))
                    ## if newly discovered link, add it to the appropriate queue
                    else:
                        cost_j = cost_i + calc_cost(link_j.__dict__)
                        queue[link_j] = cost_j
                        prev[link_j] = link_i
                        if cost_j < max_queue:
                            queue_sort.append((cost_j, link_j))
                            topsort = 1
                        else:
                            bucket.append((cost_j, link_j))
                            
        # RECONSTRUCT LEAST COST PATH
        path = []
        link_j = link_i
        while link_j:
            path.append(link_j)
            link_j = prev[link_j]
        path.reverse()
        lines = LineCollection(srid=self.kwargs['srid'], lines=path,
                               cost=cost_i)
        return lines
        
       

### End Mapmatch classes ###        


class GPSPoint(Point): 
    """A Point feature with time attribute(s)"""
    def __init__(self, x, y, utc, z=None, id_point=None, geom=None, s=None, 
	             a=None, valid=1, local_datetime=None, lat=None, lng=None,
                 mode=None, **kwargs): 
        """Initialize a GPSPoint. 

        s - speed as dst/time, None if no preceding point [double] 
        a - average acceleration as change speed/change time, None if 
            no speed for preceding point [double] 
        valid : 1 if point location considered trustworthy
        local_datetime - Python datetime using local time zone 

        Note that lat & lng are not touched by geometric methods like
        transform. They always apply to the point as initially attributed.

        """
        super(GPSPoint, self).__init__(x, y, z=z, id_point=id_point, geom=geom, 
                                       lat=lat, lng=lng, **kwargs)
        self.t = utc 
        self.s = s 
        self.a = a 
        self.valid = valid 
        self.local_datetime = local_datetime
        if local_datetime:
            self.local_time = local_datetime.strftime('%m/%d/%Y %I:%M:%S %p')
        else:
            self.local_time=None
        # self.time_local = local_datetime.strftime('%I:%M:%S %p')
        self.mode = mode

class PointCollection(object): 
    """A collection of point features with a common spatial reference."""
    def __init__(self, srid=4326, points=None, geom=None, index=False, 
                 **kwargs): 
        """Initialize a PointCollection. 

        srid - Spatial Reference System Identifier [integer] 
        points - optional list of Point objects [list] 
        geom : geometry of collection
        kwargs - {}

        """
        self.srid = srid 
        if points is None: 
            self.points = [] 
        else: 
            self.points = points
        self.geom = geom 
        if index:
            self.index_points()
    
    def index_points(self):
        """Attribute points with position index for lookups."""
        index = 0
        for p in self.points:
            p.index = index
            index += 1
        
    def get_centroid(self, db_conn, latlng=False):
        """Set centroid for PointCollection.
        
        latlng : get lat and lng along with srid xy and geom
        
        """
        if not self.geom:
            self.get_geom2(db_conn)
        data = [self.geom]
        if latlng:
            query = ("with centroid as (SELECT ST_Centroid(%s) geom) "
                     "SELECT (SELECT geom FROM centroid),"
                     "ST_X((SELECT geom FROM centroid)),"
                     "ST_Y((SELECT geom FROM centroid)),"
                     "ST_X((SELECT ST_Transform((SELECT geom FROM centroid),4326))),"
                     "ST_Y((SELECT ST_Transform((SELECT geom FROM centroid),4326)));")
            # print db_conn.debug(query, data)
            # raw_input()
            geom, x, y, lng, lat = db_conn.run(query, data)[0]
            self.centroid = Point(geom=geom, x=x, y=y, lng=lng, lat=lat)
        else:
            query = ("with centroid as (SELECT ST_Centroid(%s) geom) "
                     "SELECT (SELECT geom FROM centroid),"
                     "ST_X((SELECT geom FROM centroid)),"
                     "ST_Y((SELECT geom FROM centroid));")
            # print db_conn.debug(query, data)
            # raw_input()
            geom, x, y = db_conn.run(query, data)[0]
            self.centroid = Point(geom=geom, x=x, y=y)
        # data = [self.geom]
        # geom = db_conn.run(query, data)[0][0]
        # query = "SELECT ST_X(%s), ST_Y(%s);"
        # data = [geom, geom]
        # x, y = db_conn.run(query, data)[0]
        
    
    def get_extent(self):
        """Return extent for use in web_mapping module."""
        extent = [float('inf'), float('inf'), float('-inf'), float('-inf')]
        for p in self.points:
            x = p.x
            y = p.y
            if x < extent[0]:
                extent[0] = x
            if y < extent[1]:
                extent[1] = y
            if x > extent[2]:
                extent[2] = x
            if y > extent[3]:
                extent[3] = y
        return extent
    
    def write_psql(self, db_conn, schema, table, attributes=None, 
                   col_geom='the_geom', dim_geom=2, create=False, 
                   user_dict=None, missing=None):
        """Write PointCollection to new or existing PostgreSQL table.
        
        For now, attribute names must match PostgreSQL column names.
        
        table : schema.table
        attributes : list, excluding geometry field
          [('varname1', 'integer'), ('varname2', 'numeric')...]
        col_geom : name of geometry column
        dim_geom : dimensions of geometry
        create : if true, create table if not exists
        user_dict : loook here if attribute not found in __dict__ or kwargs
        
        """
        if self.points:
            query = 'INSERT INTO {}.{} ('.format(schema, table)
            query += ', '.join([a for a in attributes])
            query += ', {}) VALUES '.format(col_geom)
            data = []
            query += ', '.join(['(' + ','.join(['%s' for i in 
                               range(len(attributes) + 1)]) + 
                               ')' for p in self.points])
            query += ';'
            
            if not user_dict:
                user_dict = {}
            for p in self.points:
                values = []
                for name in attributes:
                    try:
                        values.append(p.__dict__[name])
                    except KeyError:
                        try:
                            values.append(p.kwargs[name])
                        except KeyError:
                            try:
                                values.append(user_dict[name])
                            except KeyError:
                                values.append(missing)
                values.append(p.geom)
                data.extend(values)
            i = 0
            while False:
                print db_conn.debug(query, data)[i: i + 1000]
                i += 1000
                raw_input()
            db_conn.run(query, data, commit=True)
            # raw_input('?')
            
    def write_shp(self, name, data_dir, attributes=None):
        """Use PyShp module to write PointCollection to a shapefile.

        attributes : [(name, type, fmt),...] for shapefile.Writer.field

        """

        def write_prj(srid, location):
            """Generate prj file based on SRID"""
            # TODO: lookup projections from file or database on the fly
            f = open(location + '.prj', 'w')
            if self.srid == 4326:
                prj = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
            elif self.srid == 2838:
                prj = 'PROJCS["NAD83(HARN) / Oregon North",GEOGCS["NAD83(HARN)",DATUM["D_North_American_1983_HARN",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",46],PARAMETER["standard_parallel_2",44.33333333333334],PARAMETER["latitude_of_origin",43.66666666666666],PARAMETER["central_meridian",-120.5],PARAMETER["false_easting",2500000],PARAMETER["false_northing",0],UNIT["Meter",1]]'
            elif self.srid == 2913:
                prj = 'PROJCS["NAD83(HARN) / Oregon North (ft)",GEOGCS["NAD83(HARN)",DATUM["D_North_American_1983_HARN",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",46],PARAMETER["standard_parallel_2",44.33333333333334],PARAMETER["latitude_of_origin",43.66666666666666],PARAMETER["central_meridian",-120.5],PARAMETER["false_easting",8202099.738],PARAMETER["false_northing",0],UNIT["Foot",0.3048]]'
            f.write(prj)
            f.close()
        if attributes is None:
            attributes = []
        else:
            attributes = attributes
        write_prj(self.srid, data_dir + '/' + name)
        w = shapefile.Writer(shapeType=1)
        w.autoBalance = 1  # Forces every shape to have a record
        w.field('id', 'N', decimal=0)
        for attr, type, fmt in attributes:
            if type == 'N':
                w.field(attr, type, decimal=fmt)
            else:
                w.field(attr, type, fmt)
        for p in self.points:
            w.point(p.x, p.y)
            values = []
            for a in attributes:
                try:
                    values.append(p.kwargs[a[0]])
                except KeyError:
                    values.append(-9)
            w.record(p.id, *values)
        w.save(data_dir + '/' + name)            

    def get_geom(self, db_conn): 
        """Create geometry attributes for xy data."""

        for p in self.points: 
            query = ("SELECT ST_GeomFromText('POINT(%s %s)', %s);")
            data = [p.x, p.y, self.srid]
            p.geom = db_conn.run(query, data)[0][0]
        # query = "SELECT ST_Collect(Array%s);"
        # data = [','.join([g for g in self.geom])]
        # self.geom = db_conn.run(query, data)[0][0]

    def get_geom2(self, db_conn, recalc_points=True): 
        """Create geometry attributes for xy data."""
        if not self.points:
            self.geom = None
            print 'Warning: get_geom2 called on empty GPSCollection'
            return False
        self.geom_wkt = 'MULTIPOINT('
        for p in self.points: 
            self.geom_wkt += '{x} {y},'.format(x=p.x, y=p.y)
        # if self.geom_wkt.endswith(','):
        self.geom_wkt = self.geom_wkt[:-1]
        self.geom_wkt += ')'
        query = "SELECT ST_GeomFromText(%s, %s);"
        data = [self.geom_wkt, self.srid]
        self.geom = db_conn.run(query, data)[0][0]
        if recalc_points:
            query = "SELECT ST_Dump(%s);"
            data = [self.geom]
            point_geoms = db_conn.run(query, data)
            i = 0
            # todo : parsing the result set seems fragile
            # ex row: ('({1},0101000020E610000093C6681D55A95EC092B3B0A71DC64640)',)
            # try (ST_Dump(geom)).geom as cleaner solution
            for row in point_geoms:
                geom = row[0].split(',')[1][:-1]
                self.points[i].geom = geom
                i += 1

    def transform(self, new_srid, db_conn, xy_update=True):
        """Transform (re-project) geometry and xy coordinates of a
        GPSCollection

        srid - SRID to transform collection to [int]
        pg_login - database connection info [dict]
        xy_update - update XY coordinates as well as geometry [boolean] 

        """
        for p in self.points:
            if p.geom is not None:
                query = "SELECT ST_Transform(geometry(%s), %s);"
                data = [p.geom, new_srid]
                # print db_conn.cur.mogrify(query, data)
                p.geom = db_conn.run(query, data)[0][0]
            else:
                query = ("SELECT ST_Transform("
                         "(SELECT ST_GeomFromText('POINT(%s %s)', %s))"
                         ",%s);")
                data = [p.x, p.y, self.srid, new_srid]
                try:
                    p.geom = db_conn.run(query, data)[0][0]
                except Exception, e:
                    print e.pgerror
                    raw_input('any key to continue')
            if xy_update is True:
                query = "SELECT ST_x(%s), ST_y(%s);"
                data = [p.geom, p.geom]
                p.x, p.y = db_conn.run(query, data)[0]
        self.srid = new_srid
        data = [[p.geom for p in self.points]]
        query = "SELECT ST_Collect(%s);"
        self.geom = db_conn.run(query, data)[0][0]

    def transform2(self, new_srid, db_conn):
        """Transform (re-project) geometry and xy coordinates of a
        GPSCollection

        srid - SRID to transform collection to [int]
        pg_login - database connection info [dict]
        xy_update - update XY coordinates as well as geometry [boolean] 

        """
        if not self.geom:
            print 'Warning: transform2 called on empty PointCollection.'
            return False
        data = [self.geom, new_srid]
        query = "SELECT ST_Transform(geometry(%s), %s);"
        self.geom = db_conn.run(query, data)[0][0]
        self.srid = new_srid
        data = [self.geom]
        query = "SELECT ST_Dump(%s);"
        point_geoms = db_conn.run(query, data)
        data = [self.geom]
        query = "SELECT ST_AsText(%s);"
        self.geom_wkt = db_conn.run(query, data)[0][0]
        xy_text = self.geom_wkt[11:-1].split(',')
        point_xy = [xy.split(' ') for xy in xy_text]
        i = 0
        for row in point_geoms:
            geom = row[0].split(',')[1][:-1]
            self.points[i].geom = geom
            self.points[i].x = float(point_xy[i][0])
            self.points[i].y = float(point_xy[i][1])
            i += 1

class GPSCollection(PointCollection):
    """A collection of sequential GPS points with a common spatial reference.

    GPSCollection(srid=4326, points=None, unit_id=None, group_id=None, 
	              phase_id=None, activities=None, trips=None, **kwargs)

	Parameters
    ----
    srid : Spatial Reference System Identifier numeric code
    points : Sequential list of GPSPoint objects
	unit_id : Unique individual (e.g. participant) identifier
	group_id : Group (e.g. household) identifier
    phase_id : Phase identifier, for multi-period data
    activities : non-travel periods
    trips : travel sequences between activities
    kwargs : {interval}

    """
    def __init__(self, srid=4326, points=None, unit_id=None, group_id=None, 
	             phase_id=None, activities=None, trips=None, geom=None,
                 **kwargs): 
        # todo : move activities to Deployment attribute
        super(GPSCollection, self).__init__(srid=srid, points=points, 
                                            geom=geom, **kwargs)
        #self.unit_id = unit_id
        #self.group_id = group_id		
        self.phase_id = phase_id 
        if trips is None:
            self.trips = []
        else:
            self.trips = trips
        if activities is None:
            self.activities = []
        else:
            self.activities = activities
        self.kwargs = kwargs
        self.stats = {}
        if self.points:
            self.time_start = self.points[0].t
            self.time_stop = self.points[-1].t
        else:
            self.time_start = None
            self.time_stop = None

    def get_stats(self):
        """Calculate statistics for collection."""
        points = self.points
        self.stats['time_GPS_min'] = (points[-1].t - points[0].t) / 60.0
        meters = 0.0
        index = 1
        # todo : seems costly to loop through points again here
        for p in points[1:]:
            meters += (((points[index].x - points[index - 1].x) ** 2 + 
                        (points[index].y - points[index - 1].y) ** 2) ** 0.5)
            index += 1
        self.stats['dist_GPS_m'] = meters    
        self.stats['dist_GPS_mi'] = meters * 0.000621371
        self.stats['dist_net_m'] = (((points[-1].x - points[0].x) ** 2 + 
                                     (points[-1].y - points[0].y) ** 2) ** 0.5)
        try:
            self.stats['spd_avg_GPS_mph'] = stats.mean([p.s for p in points]) * 2.23694
        except TypeError:
            self.stats['spd_avg_GPS_mph'] = None
        try:
            self.stats['spd_avg_seg_mph'] = ((self.stats['dist_GPS_mi'] /
                                                self.stats['time_GPS_min']) * 60.0)
            self.stats['spd_avg_seg_ms'] = (self.stats['dist_GPS_m'] / 
                                       (self.stats['time_GPS_min'] * 60.0))
        except ZeroDivisionError:
            self.stats['spd_avg_seg_mph'] = 0
            self.stats['spd_avg_seg_ms'] = 0
        try:
            self.stats['spd_max_GPS_mph'] = max([p.s for p in points]) * 2.23694
            self.stats['spd_avg_GPS_ms'] = stats.mean([p.s for p in points])
            self.stats['spd_med_GPS_ms'] = stats.median([p.s for p in points])
            self.stats['spd_med_GPS_ms_sq'] = self.stats['spd_med_GPS_ms'] ** 2
            self.stats['spd_nf_GPS_ms'] = stats.percentile([p.s for p in points], 0.95)
            self.stats['spd_nf_GPS_ms_sq'] = self.stats['spd_nf_GPS_ms'] ** 2
            self.stats['spd_var_GPS_ms'] = stats.variance([p.s for p in points])
        except TypeError:
            self.stats['spd_max_GPS_mph'] = None
            self.stats['spd_avg_GPS_ms'] = None
            self.stats['spd_med_GPS_ms'] = None
            self.stats['spd_med_GPS_ms_sq'] = None
            self.stats['spd_nf_GPS_ms'] = None
            self.stats['spd_nf_GPS_ms_sq'] = None
            self.stats['spd_var_GPS_ms'] = None
        try:
            self.stats['acc_nf_GPS_mss'] = stats.percentile([abs(p.a) for p in points], 0.95)
        except TypeError:
            pass
        
        try:
            self.stats['spd_cv_GPS_ms'] = ((self.stats['spd_var_GPS_ms'] ** 0.5) /
                                            self.stats['spd_avg_GPS_ms'])

        except ZeroDivisionError:
            # zero distance
            self.stats['spd_cv_GPS_ms'] = 0
        except TypeError:
            # no speed
            self.stats['spd_cv_GPS_ms'] = None
        try:
            self.stats['indirectness'] = (self.stats['dist_GPS_m'] / 
                                          self.stats['dist_net_m'])
        except ZeroDivisionError:
            # zero distance segment; set to 0
            self.stats['indirectness'] = 0
            
    def get_density(self, num_near_points, search_radius):

        """Calculate number of points within given radius of each GPSPoint.

        num_near_points limits number of sequential points to search in each 
        direction.

        """
        # todo : Consider calculating with PostGIS DWithin()
        
        points = self.points
        npoints = len(points)
        i = 0
        for p in points:
            n_within = 0
            for pp in points[max(0, i - num_near_points): i + num_near_points + 1]:
                d = ((p.x - pp.x) ** 2 + (p.y - pp.y) ** 2) ** 0.5
                if d <= search_radius:
                    n_within += 1
            i += 1
            p.kwargs['density'] = n_within - 1  ## exclude self

    def get_accel(self):
        """Calculate average acceleration values between points."""
        if not self.points:
            print 'Warning: get_accel called on empty GPSCollection.'
            return False
        s0 = None
        pp = self.points[0]
        pp.a = 0.0
        for p in self.points[1:]:
            try:
                p.a = (p.s - pp.s) / (p.t - pp.t)
            except ZeroDivisionError:
                p.a = 0.0
            pp = p

    def load_trips(self, tf, unit_id, attribute_points=False):
        """Load trip profiles from a Trip.write_csv file.

        Note that points must already be loaded with geom.
        This will erase any existing trips from memory.

        tf : TripFile object
        """
        # todo : handle assigning AM counts during loading
        self.trips = []
        if unit_id not in tf.segments:
            print 'Warning: no trips loaded for unit_id {0}'.format(unit_id)
            return False
        n_points = len(self.points)
        index = 0
        for trip_num, d in sorted(tf.segments[unit_id].iteritems()):
            start = d['time_start']
            end = d['time_stop']
            mode = d['mode']
            match = False
            points = []
            # Use times instead of point ids for version safety
            # todo : Test behavior with different point masks applied.
            while index < n_points:
                #print self.points[index].t, start
                if self.points[index].t == start:
                    match = True
                    try:
                        while index < n_points and self.points[index].t <= end:
                            points.append(self.points[index])
                            index += 1
                    except IndexError:
                        print num, start, end, index
                        raise IndexError
                    g = GPSCollection(srid=self.srid, points=points[:], 
                                      unit_id=unit_id)
                    self.trips.append(Trip(trip_num, 
                                           gps=g, mode=mode))
                    break
                else:
                    index += 1
            if not match:
                print self.points[0].t, start, end, self.points[-1].t
                err = ('No matching points found for {0}-{1}'
                       .format(unit_id, trip_num))
                raise Error(err)
        if attribute_points:
            for trip in self.trips:
                for p in trip.gps.points:
                    p.kwargs['trip'] = trip.num
                    p.kwargs['trip_mode'] = trip.mode
        return True

    def split_trips(self, deploy_start=None, min_trip_pts=10, 
                    max_time_gap=300,
                    num_near_points=30, search_radius=15, 
                    min_bundle_density=10.0, min_bundle_ratio=2.0/3.0,
                    min_bundle_time=120, sequence_length=15,
                    max_time_zero_spd=120, zero_spd=0.01,
                    attribute_points=False, debug=False):
        """Divide GPS time into discrete trips and activity participation 
        periods. Based on method described in Schuessler & Axhausen (xxxx).
        
        deploy_start : deployment start time (epoch, UTC)
        min_trip_pts : minimum number of points for valid trip
        max_time_gap : max time gap (seconds) before trip splits
        num_near_points : number of preceding and succeeding points included in
            density calculation and succeeding points in point bundle search
        search_radius : max distance for point to be included in density count 
        min_bundle_density : min point density defining a bundle
        min_bundle_time : elapsed time with at least min_bundle_density
            before trip splits
        min_bundle_ratio : min p.kwargs['density'] >= min_point_density / 
            n_points defining a bundle
        sequence_length : incremental number of points in search_bundle()
        attribute_points : add trip attribute to point.kwargs

        """

        def validate(max_indirectness=3.0, min_net_distance=100,
                     max_slope=0.4, max_speed=40, max_time_gap_ratio=0.3,
                     max_alt_gap_ratio=0.7):
            """Flag trips based on data quality and point behavior.



            """
            for trip in self.trips:
                points = trip.gps.points
                count_time_gaps = 0.0
                count_alt_gaps = 0.0
                count_spd_gaps = 0.0 
                distance_point = 0.0
                indirect_flag = 0
                net_dst_flag = False
                time_gap_flag = False
                alt_gap_flag = False
                kill_flag = False
                # todo : autodetect planar geometry
                distance_direct = (((points[-1].x - points[0].x) ** 2 + 
                                    (points[-1].y - points[0].y) ** 2) ** 0.5)
                time_elapsed = trip.time_stop - trip.time_start
                points_n = len(trip.gps.points)
                i = 1
                for p in trip.gps.points[1:]:
                    k = p.kwargs
                    pp = trip.gps.points[i - 1]
                    kp = pp.kwargs
                    if (p.t - pp.t) > self.kwargs['interval']:
                        count_time_gaps += 1
                    # todo: handle non-planar geometries
                    d = ((p.x - pp.x) ** 2 + (p.y - pp.y) ** 2) ** 0.5
                    distance_point += d
                    try:
                        if abs((p.z - pp.z) / d) > max_slope:
                            count_alt_gaps += 1
                    except ZeroDivisionError:
                        pass
                    # todo: find out what could cause this [mask set]
                    try:
                        if d / (p.t - pp.t) > max_speed:
                            count_spd_gaps += 1
                            if debug:
                                print 'spd_gap:', pp.id, p.id, d, p.t - pp.t
                    except ZeroDivisionError:
                        pass
                    i += 1
                if count_time_gaps / points_n > max_time_gap_ratio:
                    time_gap_flag = True
                if count_alt_gaps / points_n > max_alt_gap_ratio:
                    alt_gap_flag = True
                try:
                    if distance_point / distance_direct > max_indirectness:
                        indirect_flag = 1
                except ZeroDivisionError:
                    if distance_point > 0:
                        indirect_flag = 1
                if distance_direct < min_net_distance:
                    net_dst_flag = True

                if time_gap_flag and alt_gap_flag:
                    kill_flag = True
                elif net_dst_flag and (time_gap_flag or alt_gap_flag):
                    kill_flag = True
                elif indirect_flag:
                    loop_flag = True
                trip.kwargs['kill'] = kill_flag

                if 1:
                    trip.stats['spd_flg'] = count_spd_gaps
                    trip.stats['indirect_flg'] = indirect_flag
                    trip.stats['dist_net_flg'] = net_dst_flag
                    # todo : handle other units
                    trip.stats['spd_avg_seg_mph'] = ((distance_point / time_elapsed) * 
                                               2.23694)  # m/s to mi/hr
                    trip.stats['spd_avg_GPS_mph'] = \
                      ((sum([p.s for p in trip.gps.points]) / points_n) *
                       2.23694)
                    trip.stats['spd_max_GPS_mph'] = \
                      max([p.s for p in trip.gps.points]) * 2.23694
                    trip.stats['dist_GPS_mi'] = distance_point * 0.000621371
                    trip.stats['dist_net_m'] = distance_direct
                    ##########################
                    try:
                        trip.stats['indirectness'] = distance_point / distance_direct
                    except ZeroDivisionError:
                        trip.stats['indirectness'] = -999
                    trip.stats['time_gaps'] = count_time_gaps
                    trip.stats['alt_gaps'] = count_alt_gaps
                    trip.stats['spd_gaps'] = count_spd_gaps
                    trip.stats['dst_point'] = distance_point
                    trip.stats['npoints'] = points_n
                    trip.stats['time_GPS_min'] = time_elapsed / 60.0
                    ##########################    
            if not debug:
                valid_trips = []
                inspect_trips = []
                killed_trips = []
                v = 1
                i = 1
                k = 1
                j = 1
                travelday = 0
                # Separate trips into bins, and renumber.
                for trip in self.trips:
                    if trip.kwargs['kill'] is False:
                        if trip.kwargs['travelday'] > travelday:
                            travelday = trip.kwargs['travelday']
                            j = 1
                        trip.num = v
                        trip.kwargs['travday_num'] = j
                        valid_trips.append(trip)
                        v += 1
                        j += 1
                    elif trip.stats['dist_net_m'] >= 100:
                        trip.num = i
                        inspect_trips.append(trip)
                        i += 1
                    else:
                        trip.num = k
                        killed_trips.append(trip)
                        k += 1
                self.trips = valid_trips[:]
                self.trips_inspect = inspect_trips[:]
                self.trips_killed = killed_trips[:]

        def get_trip_end(i, trip_start, max_gap=5.0, N=2):
            """Find trip end based on last valid trip end sequence.
            This trims noisy points at the end of trips.

            i - proposed trip end/activity start
            trip_start - start point of this trip
            gap - max time gap (seconds) for points in trip start sequence
            N - number of points less one in trip end sequence (# gaps)

            """
            n = 0
            for h in range(i, trip_start, -1):
                if (self.points[h].t - self.points[h - 1].t <= max_gap
                      and self.points[h - 1].kwargs['density'] < min_bundle_density):
                    n += 1
                else:
                    n = 0
                if n == N:
                    return h + 1
            return trip_start

        def get_activity_end(i, min_spd=0.75, max_gap=5.0, N=3):
            """Find activity end based on next valid trip start sequence.

            i - proposed activity end/trip start
            spd - minimum speed for trip start sequence
            gap - max time gap (seconds) for points in trip start sequence
            n - number of points in trip start sequence

            """
            h = i
            n = 0
            for p in self.points[i + 1:]:
                if p.s >= min_spd and p.t - self.points[h].t <= max_gap:
                    n += 1
                else:
                    n = 0
                h += 1
                if n == N:
                    return h - N
            return h

        def get_activity_gaps():
            """Find potential activity starts and stops based on
            time gaps between GPS points and speed close to zero.

            """
            activities_prov = ap = []
            i = 1
            time_zero_spd = 0.0
            start_zero_spd = None
            act_end = 0
            for p in self.points[i:]:
                if p.s <= zero_spd:
                    time_zero_spd += p.t - self.points[i - 1].t
                    if start_zero_spd is None:
                        start_zero_spd = i - 1
                    if p.t - self.points[i - 1].t > max_time_gap:
                        # 1
                        if time_zero_spd > max_time_zero_spd:
                            act_start = get_trip_end(start_zero_spd, act_end)
                            act_end = get_activity_end(i)
                            ap.append((act_start, 0))
                            ap.append((act_end, 1))
                            # print '1',act_start,act_end
                        #2
                        else:
                            act_start = get_trip_end(i - 1, act_end)
                            act_end = get_activity_end(i)
                            ap.append((act_start, 0))
                            ap.append((act_end, 1))
                            # print '2',act_start,act_end
                        time_zero_spd = 0.0
                        start_zero_spd = None
                elif p.s > zero_spd:
                    #3
                    if time_zero_spd > max_time_zero_spd:
                        act_start = get_trip_end(start_zero_spd, act_end)                        
                        act_end = get_activity_end(i - 1)
                        ap.append((act_start, 0))
                        ap.append((act_end, 1))
                        # print '3',act_start,act_end
                    #4
                    elif p.t - self.points[i - 1].t > max_time_gap:
                        act_start = get_trip_end(i - 1, act_end)
                        act_end = get_activity_end(i)
                        ap.append((act_start, 0))
                        ap.append((act_end, 1))
                        # print '4',act_start,act_end
                    time_zero_spd = 0.0
                    start_zero_spd = None
                i += 1
            return activities_prov

        def search_bundles(min_density, min_bundle_ratio, sequence_length):
            """Return bundles of high-density points.

            Bundles always end on dense point.
            Minimum bundle time is min_bundle_time.

            """
            bundles = []
            index = 0
            # todo : could speed up by skipping all leading dense points
            #   in failed bundles when starting next sequence search
            end = False
            while index < len(self.points):
                if self.points[index].kwargs['density'] >= min_density:
                    start = index
                    search_index = index
                    last_dense_point = index
                    sequence_points = 0.0
                    sequence_dense_points = 0.0
                    bundle = True
                    while bundle:
                        for i in range(sequence_length):
                            try: 
                                sequence_points += 1
                                if (self.points[search_index].kwargs['density'] 
                                    >= min_density):
                                    sequence_dense_points += 1
                                    bundle_ratio = (sequence_dense_points / 
                                      sequence_points)
                                    if bundle_ratio >= min_bundle_ratio:
                                        last_dense_point = search_index                                            
                                search_index += 1         
                            except IndexError:
                                end = True
                                break
                        if ((sequence_dense_points / sequence_points) >=
                            min_bundle_ratio and not end):
                            sequence_points = 0.0
                            sequence_dense_points = 0.0
                            index = last_dense_point  # outer while index                                
                        elif (self.points[last_dense_point].t - 
                              self.points[start].t >= min_bundle_time):
                            bundles.extend([(start, 0), (last_dense_point, 1)])
                            index = last_dense_point  # outer while index
                            bundle = False
                        else:
                            bundle = False
                index += 1
            return bundles
        
        def resolve_trips(activities_prov):
            "Finalize provisional trip and activity sequence."""                
            activities = []
            trips = []
            act_start = []
            act_end = []
            act_num = 1
            act_travday_num = 1
            act_index = 0
            trip_start = 0
            trip_end = None
            trip_num = 1
            trip_num_travday = 1
            max_index = len(self.points) - 1
            travelday = 0
            deploy_start = self.kwargs['deploy_start']
            for index, code in sorted(activities_prov):
                if code == 0:
                    act_start.append(index)
                else:
                    act_end.append(index)
                if len(act_end) == len(act_start):
                    start = trip_end = act_start[0]
                    end = act_end[-1]
                    if trip_end - trip_start >= min_trip_pts:
                        if self.points[start].kwargs['travelday'] > travelday:
                            travelday = self.points[start].kwargs['travelday']
                            trip_num_travday = 1
                        activities.append((start, end))
                        # self.activities.append(Activity(act_num, self.points[start].geom,
                                       # self.points[end].geom,
                                       # self.points[start].t,
                                       # self.points[end].t,
                                       # gps=GPSCollection(self.srid,
                                                     # points=[pt for pt in self.points[start: end]])))
                        act_index += 1
                        # TODO: Double check trip_end vs trip_end + 1
                        self.trips.append(Trip(trip_num, 
                                               gps=GPSCollection(self.srid,
                                                                 points=[pt for pt in self.points[trip_start: trip_end + 1]]),
                                               travelday=travelday, 
                                               travday_num=trip_num_travday))
                        trip_num += 1
                        trip_num_travday += 1
                    elif act_index > 0:
                        activities[act_index - 1] = ((activities[act_index - 1][0],
                                                      end))
                    else:
                        activities.append((start, end))
                    trip_start = min(end + 1, max_index) # Next trip starts where current activity ends
                    act_start = []
                    act_end = []
            trip_end = len(self.points) - 1
            #if ((not self.trips or trip_end - trip_start >= min_trip_pts) and
            if trip_end - trip_start >= min_trip_pts:
                if self.points[start].kwargs['travelday'] > travelday:
                    travelday = self.points[start].kwargs['travelday']
                    trip_num_travday = 1
                self.trips.append(Trip(trip_num, 
                                       gps=GPSCollection(self.srid,
                                                         points=[pt for pt in self.points[trip_start: trip_end + 1]]),
                                       travelday=travelday, 
                                       travday_num=trip_num_travday))
            # print [(t.num, t.kwargs['travelday'], t.kwargs['travday_num']) for t in self.trips]
           
        # START MAIN FUNCTION
        #print num_near_points, search_radius
        self.get_density(num_near_points, search_radius)  # Calculate and store point densities.
        # ap stores provisional activity periods
        activities_prov = ap = get_activity_gaps()
        bundles = search_bundles(min_bundle_density, min_bundle_ratio,
                                 sequence_length)
        #print bundles
        ap.extend(bundles)
        #pprint(sorted(ap))
        # print len(bundles), len(ap)
        if ap:
            resolve_trips(ap)
            validate()
        elif self.points:
            # entire sequence is a trip?
            trip = Trip(1, 
                        gps=GPSCollection(self.srid,
                          points=[pt for pt in self.points]),
                        travelday=self.points[0].kwargs['travelday'], 
                        travday_num=1)
            self.trips.append(trip)
        else:
            # GPSCollection is empty
            pass
        if attribute_points:
            for trip in self.trips:
                for p in trip.gps.points:
                    p.kwargs['trip'] = trip.num
                    p.kwargs['trip_travday'] = trip.kwargs['travday_num']
        ##################################
        if debug:
            print ('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13}'
                   .format('num', 'start', 'end', 'n_points', 'indirect_flag',
                           'indirectness', 'loop_flag', 'net_dst', 'dst_point',
                           'time_gaps', 'alt_gaps', 'spd_gaps', 't1', 't2'))
            for trip in self.trips:
                k = trip.stats
                print ('{0} {1} {2} {3} {4} {5:.2f} {6} {7:.0f} {8:.0f} {9} {10} '
                       '{11} {12} {13}'
                       .format(trip.num, trip.point_start,  
                               trip.point_stop,
                               k['npoints'],
                               k['indirect_flg'], k['indirectness'],
                               k['dist_net_flg'], k['dist_net_m'], k['dist_GPS_mi'],
                               k['time_gaps'],
                               k['alt_gaps'], k['spd_gaps'], trip.time_start,
                               trip.time_stop))
        ###################################


    def write_shp(self, name, data_dir, attributes=None, missing_val=-9):
        """Use PyShp module to write PointCollection to a shapefile.

        attributes : [(name, type, fmt),...] for shapefile.Writer.field

        """

        def write_prj(srid, location):
            """Generate prj file based on SRID"""
            # TODO: lookup projections from file or database on the fly
            f = open(location + '.prj', 'w')
            if self.srid == 4326:
                prj = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
            elif self.srid == 2838:
                prj = 'PROJCS["NAD83(HARN) / Oregon North",GEOGCS["NAD83(HARN)",DATUM["D_North_American_1983_HARN",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",46],PARAMETER["standard_parallel_2",44.33333333333334],PARAMETER["latitude_of_origin",43.66666666666666],PARAMETER["central_meridian",-120.5],PARAMETER["false_easting",2500000],PARAMETER["false_northing",0],UNIT["Meter",1]]'
            elif self.srid == 2913:
                prj = 'PROJCS["NAD83(HARN) / Oregon North (ft)",GEOGCS["NAD83(HARN)",DATUM["D_North_American_1983_HARN",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",46],PARAMETER["standard_parallel_2",44.33333333333334],PARAMETER["latitude_of_origin",43.66666666666666],PARAMETER["central_meridian",-120.5],PARAMETER["false_easting",8202099.738],PARAMETER["false_northing",0],UNIT["Foot",0.3048]]'
            f.write(prj)
            f.close()
        if attributes is None:
            attributes = []
        else:
            attributes = attributes
        write_prj(self.srid, data_dir + '/' + name)
        w = shapefile.Writer(shapeType=1)
        w.autoBalance = 1  # Forces every shape to have a record
        w.field('id', 'N', decimal=0)
        w.field('utc_time', 'N', decimal=0)
        w.field('local_datetime', 'C', 19)
        w.field('speed', 'N', decimal=2)
        w.field('accel', 'N', decimal=2) 
        w.field('alt', 'N', decimal=1)
        w.field('valid', 'N', decimal=0)
        for attr, type, fmt in attributes:
            if type == 'N':
                w.field(attr, type, decimal=fmt)
            else:
                w.field(attr, type, fmt)
        for p in self.points:
			# todo : automate srid detection, and transform if needed
            w.point(p.x, p.y)
            # w.point(p.kwargs['x_gps'], p.kwargs['y_gps'])
            values = []
            for a in attributes:
                try:
                    values.append(p.__dict__[a[0]])
                except KeyError:
                    try:
                        values.append(p.kwargs[a[0]])
                    except KeyError:
                        values.append(missing_val)
            w.record(p.id, p.t, p.local_datetime, p.s, p.a, p.z, int(p.valid),
			         *values)
        w.save(data_dir + '/' + name)         

    def screen(self, max_spd=40, max_accel=10, min_altitude=-100, 
               max_altitude=9999, min_spd=1.0, debug=False):
        """Screen points based on time, speed, and altitude."""
        n_masked = 0
        i = 1
        j = 1  # distance back to last valid point in time
        time_warn = True
        # todo : the units thing again
        for p in self.points[1:]:
            pp = self.points[i - j]
            spd = (((p.x - pp.x) ** 2 + (p.y - pp.y) ** 2) ** 0.5) / (p.t - pp.t)
            p.kwargs['spd'] = spd
            if p.t == pp.t:
                p.valid = 0
                n_masked += 1
                j = 1
            elif p.t < pp.t:
                # Try to find next point in future
                p.valid = 0
                n_masked += 1
                j += 1
                if time_warn:
                    print ('Warning: Time runs backward after point '
                           '{0}!'.format(pp.id))
                    time_warn = False
                # err = ('Time runs backward between points {0} and {1}!'
                       # .format(pp.id, p.id))
                # raise Error(err)
            elif p.z > max_altitude or p.z < min_altitude:
                p.valid=0
                n_masked += 1
                j = 1
            elif p.a and p.a > max_accel:
                p.valid = 0
                n_masked += 1
                j = 1
            elif spd > max_spd or spd < min_spd:
                p.valid = 0
                n_masked += 1
                j = 1
            else:
                p.valid = 1
                j = 1
            i += 1
        if debug and n_masked > 0:
            print 'masked =', n_masked

    def mask(self, db_conn):
        """Remove masked points from point stream."""
        # todo : consider adding a renumber function
        valid_points = []
        self.points_masked = []
        self.points_all = self.points[:]  ## Store full list of points
        for p in self.points:
            if p.valid:
                valid_points.append(p)
            else:
                p.index = None
                self.points_masked.append(p)
        self.points = valid_points[:]
        # todo : Is there a way around this? Seems costly.
        self.get_geom2(db_conn)
        self.index_points()

    def smooth(self, db_conn, bandwidth=10, max_points=20, recalc=False):
        """Smooth x,y,[z] values using a Gaussian filter.

           The kernel density function is taken from Schuessler & Axhausen 
           (xxxx). max_points was not included in their specification. Original
           coordinate values are stored in new x_raw and y_raw attributes.

           gps - GPSCollection object
           bandwidth - Sigma in denominator of smoothing function,
               seconds; note that large values may give weird results.
           max_points - Max points before or after to include in
               filtered position calculation. This is used to keep computation
               time reasonable without substantially altering the result.
           recalc - If True, recalculate speed and acceleration using
                    smoothed points.

        """
        t0 = None
        i = 0  # counter
        n_points = len(self.points)
        for p in self.points:
            p.x_raw = p.x
            p.y_raw = p.y
            t = p.t
            x_num_sum = 0.0
            y_num_sum = 0.0
            den_sum = 0.0
            # look only max_points points before and after
            # make endpoint smoothing symmetric to left & right
            nearest_tail = min(i - 0, n_points - i)
            lbound = max(i - max_points, i - nearest_tail)
            rbound = min(i + nearest_tail, i + max_points)
            for j in xrange(lbound, rbound):
                t_j = self.points[j].t
                # use raw values for already smoothed points
                if j < i:
                    x_j = float(self.points[j].x_raw)
                    y_j = float(self.points[j].y_raw)
                else:
                    x_j = float(self.points[j].x)
                    y_j = float(self.points[j].y)
                w_tj = math.exp(-((t - t_j)**2) / (2 * bandwidth**2))
                try:
                    x_num_sum += w_tj * x_j
                except TypeError:
                    j, w_tj, x_j
                y_num_sum += w_tj * y_j
                den_sum += w_tj
            try:
                x = x_num_sum / den_sum
                y = y_num_sum / den_sum
            except ZeroDivisionError:
                x = p.x_raw
                y = p.y_raw
            p.x = x
            p.y = y
            x_prev = x
            y_prev = y
            t0 = t
            i += 1
        self.get_geom2(db_conn)
        if recalc is True:
            self.get_speed()
            self.get_accel()


class Segment(object):
    """Most basic subset of a GPSCollection.

    Segment(num, geom_start, geom_stop, time_start, time_stop,
            gps=None, modes=None, mode=None, am=None, **kwargs)
    num : trip number, unique within Person Deployment [int]
    gps : GPSCollection defining this segment
    modes : MNL or fuzzy model of modes
    mode : mode of segment
    am : AMCollection of AMCount objects contained by this segment
    stats : dict of segment statistics
    kwargs : {}

    """
    def __init__(self, num, gps=None, modes=None, mode=None, am=None, 
                 trip=None, stage=None, stats=None, route=None, geom=None,
                 **kwargs):
        self.num = int(num)
        if modes is None:
            self.modes = {}
        else:
            self.modes = modes
        self.mode = mode
        if stats is None:
            self.stats = {}
        else:
            self.stats = stats
        self.gps = gps
        self.geom = geom
        if self.gps:
            start = self.gps.points[0]
            end = self.gps.points[-1]
            self.geom_start = start.geom
            self.geom_stop = end.geom
            self.time_start = float(start.t)
            self.time_stop = float(end.t)
            self.point_start = start.id
            self.point_stop = end.id
            self.local_datetime_start = start.local_datetime
            self.local_datetime_stop = end.local_datetime
            try:
                self.local_time_start = \
                  start.local_datetime.strftime('%m/%d/%Y %I:%M:%S %p')
                # self.local_time_start = start.local_datetime.strftime('%I:%M:%S %p')
                self.local_time_end = \
                  end.local_datetime.strftime('%m/%d/%Y %I:%M:%S %p')
                # self.local_time_end = end.local_datetime.strftime('%I:%M:%S %p')
            except AttributeError:
                self.local_time_start = None
                self.locl_time_stop = None
            self.lat_start = start.lat
            self.lng_start = start.lng
            self.lat_stop = end.lat
            self.lng_stop = end.lng
            self.npoints = len(self.gps.points)
        else:
            self.point_start = None
            self.point_stop = None
            self.local_datetime_start = None
            self.local_datetime_stop = None
            self.local_datetime_start_str = None
            # self.local_time_start = None
            self.local_datetime_end_str = None
            # self.local_time_end = None
            self.lat_start = None
            self.lng_start = None
            self.lat_stop = None
            self.lng_stop = None
            # todo : Should this be zero or None?
            self.npoints = None
            self.gps.stats = None
        self.am = am
        self.route = route
        self.kwargs = kwargs
        # self.get_stats()

    def get_stats(self):
        """Calculate and collect gps, am, and route stats."""
        if self.gps:
            if not self.gps.stats:
                self.gps.get_stats()
            self.stats.update(self.gps.stats)
        if self.am:
            if not self.am.stats:
                self.am.get_stats()
            self.stats.update(self.am.stats)
        if self.route:
            if not self.route.stats:
                self.route.get_stats()
            self.stats.update(self.route.stats)
        else:
            # todo : generalize this
            self.stats['taps'] = 9999

    def mapmatch(self, db_conn, network, attribute_points=False, u_turn=False, 
                 transit=False, initial_candidates=25, max_candidates=40, 
                 check=False, max_radius=1000, max_angle_flags=2):
        """Match GPS points to a Network."""
        if not transit:
            m = MapMatch(db_conn, self.gps, network, check)
        else:
            m = TransitMapMatch(db_conn, self.gps, network, check)
        self.route = m.solve(attribute_points=attribute_points,
                             u_turn=u_turn, initial_candidates=initial_candidates, 
                             max_candidates=max_candidates, max_radius=max_radius,
                             max_angle_flags=max_angle_flags)

    def write_csv(self, path, fieldnames=None, user_dict=None, none_val=-999):
        """Append Segment as row to csv file.

        path : path to write csv
        fieldnames : optional ordered list of str, num, or datetime.datetime 
          fieldnames to write
        user_dict : user-supplied dict to include in fieldnames search

        Note that for now top-level attribute dict keys named in fieldnames 
        must be unambiguous across all attribute dicts, or an error will be 
        raised.

        """
        def _parse(val):
            """Return True if attribute is a valid type."""
            if (isinstance(val, types.StringTypes) or
                isinstance(val, numbers.Number)):
                 return val
            #todo : handle individual time and date objects
            elif isinstance(val, datetime.datetime):
                return val.strftime('%Y-%m-%d %H:%M:%S')
            elif isinstance(val, types.NoneType):
                return none_val
            else:
                return False

        if fieldnames is None:
            fieldnames = ['hid', 'pid', 'phase_id', 'trip_num', 'stage_num', 
                          'num', 'point_start', 'point_stop', 
                          'local_datetime_start', 'local_time_stop', 
                          'time_start', 'time_stop', 'lat_start', 'lng_start',
                          'lat_stop', 'lng_stop']
        if user_dict is None:
            user_dict = {}
        fieldmap = {}
        for attr in fieldnames:
            if attr in user_dict:
                val = _parse(user_dict[attr])
                if val:
                    fieldmap[attr] = val
                else:
                    err = ('{0} is invalid for row {1}'
                           .format(attr, self.num))
                    raise TypeError(err)
            elif attr in self.__dict__:
                val = _parse(self.__dict__[attr])
                if val is not False:
                    fieldmap[attr] = val
                else:
                    err = ('{0} is invalid for row {1}'
                           .format(attr, self.num))
                    raise TypeError(err)
            else:
                for k, v in self.__dict__.iteritems():
                    if isinstance(v, types.DictType):
                        if attr in v:
                            if attr in fieldmap:
                                err = ('{0} is ambiguous for row {1}'
                                       .format(attr, self.num))
                                raise Error(err)
                            else:
                                val = _parse(v[attr])
                                if v is not False:
                                    fieldmap[attr] = v[attr]
                                else:
                                    err = ('{0} is invalid type for row {1}'
                                           .format(attr, self.num))
                                    raise TypeError(err)
            # if attr not in fieldmap:
                # err = ('{0} is not a valid attribute or sub-attribute'
                       # .format(attr))
                # raise AttributeError(err)
        row = {}
        for field in fieldnames:
            try:
                row[field] = fieldmap[field]
            except KeyError:
                row[field] = none_val
        try:
            f = open(path)
            f.close()
        except IOError:
            f = open(path, 'ab')
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
        try:
            w.writerow(row)
        except NameError:
            f = open(path, 'ab')
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writerow(row)
        f.close()
    
    def write_psql(self, db_conn, schema, table, create=False):
        """Write Segment to PostgreSQL table
        
        db_conn : psql database connection
        schema : schema holding table
        table : name of table to use
        create : crate table if not exists
        
        """
        pass

class MNL(object):
    """Multinomial logit model"""

    def __init__(self, data, alts):
        self.data = data  # dict of utility function attribute values by alt
        self.choice = {} 
        self.probs = {}  ## dict of choice prob by alt
        self.alts = alts

    def choose(self):
        """Return highest utility candidate and choice probabilities."""
        
        max_util = float(-inf)
        data = self.data
        alts = self.alts
        max_v = float('-inf')
        for alt in alts:
            alt.v = alt.eval(data)
            if alt.v > max_v:
                max_v = alt.v
                self.choice = alt.name
        sum_exp_v = sum([exp(alt.v) for alt in alts])
        for alt in alts:
            self.probs[alt.name] = alt.v / sum_exp_v
            
class Stage(Segment):
    """A collection of GPSPoints defining a single-mode stage.

    Next level up is a Trip. Next level down is a Segment.

    """
    def __init__(self, num, gps=None, modes=None, mode=None, segments=None,
                 am=None, stats=None, **kwargs):
        """Initialize a Stage instance.

        segments : list of Segment objects contained by this Stage
        kwargs : {}

        """
        super(Stage, self).__init__(num, gps=gps, modes=modes, mode=mode, 
                                    am=am, stast=stats, **kwargs)
        if segments is None:
            self.segments = []
        else:
            self.segments = segments
        if self.gps:
            self.gps.get_stats()
        # print 'stage {} from {} to {}'.format(num, self.gps.points[0].id, self.gps.points[-1].id)
        # raw_input()
    
    def get_stats(self):
        """Collect stats from GPS, AM, route, and mode data."""
        super(Stage, self).get_stats()
        # todo : re-write to remove hardcoded limit
        if self.stats['indirectness'] > 3.0:
            self.stats['indirect_flg'] = 1
        else:
            self.stats['indirect_flg'] = 0
                
    def mode_mnl(self, model, attribute_points=False):
        """Assign travel mode probabilities to Stage.

        model : an MNL mode model module
        
        """
        # todo : handle missing data
        alts = model.alts
        m = dc_tools.MNL(self.stats, alts)
        m.choose()
        self.kwargs['model'] = model.name
        self.modes = m.probs
        self.mode = m.choice.num
        # todo : generalize modes
        self.kwargs['pw'] = self.modes['1']
        self.kwargs['pb'] = self.modes['2']
        self.kwargs['pa'] = self.modes['3']
        if '4' in self.modes:
            self.kwargs['pt'] = self.modes['4']
            self.kwargs['tmodes'] = ';'.join(list(self.route.tmodes))
        if 'pt' not in self.kwargs:
            self.kwargs['pt'] = -999
            self.kwargs['tmodes'] = -999
        if attribute_points:
            for p in self.gps.points:
                p.mode = self.mode

    def mode_fuzzy(self, modes, fuzzy_vars):
        """Use fuzzy logic to detect travel mode fuzzy membership for each stage."""
        fv = fuzzy_vars
        # TODO: Clean up handling of missing data
        missing_data = False
        for mode in modes:
            max_rule_score = 0.0
            for rule in mode.rules:
                min_score = 1.0
                for var_nm, val_nm in rule:
                    if var_nm in self.gps.stats:
                        x = self.gps.stats[var_nm]
                    else:
                        x = self.am.stats[var_nm]
                    # Missing values cause variable to be ignored.
                    if x is None:
                        pass
                    else:
                        y = fv[var_nm].evaluate(x, val_nm)
                        if y < min_score:
                            min_score = y
                if missing_data is True:
                    break
                elif min_score > max_rule_score:
                    max_rule_score = min_score
            if missing_data is True:
                break
            else:
                self.modes[mode.name] = max_rule_score

    def stitch_routes(self, db_conn, network, transit=False, 
                      stitch_type='path', 
                      min_match=0.75, missing_aps=50):
        """Stitch segment routes together.

        stitch_type : how to handle gpas between segments
          'path' - use shortest path to connect 
          'gap' - leave gaps in route (e.g. when skipping unmatched segments)
        min_match : minimum portion of valid route that is on-network
        missing_aps : score assigned for off-network segments
        Placeholder for method that will stitch segments together to form
        stage.route. Need to decide whether to use shortest path imputation
        or treat inter-segment travel as off-network        

        """
        def make_dummy():
            """Return a dummy route for off-network travel."""
            pass
            
        if len(self.segments) == 1:
            # no stitching needed
            self.route = self.segments[0].route
        elif transit:
            self.route = self.segments[0].route
            for segment in self.segments[1:]:
                if self.route and segment.route:
                    self.route.extend(db_conn, segment.route, segment.npoints)
                    # except:
                        # print 'Warning: Problem stitching route.'
                else:
                    self.route = None
            if not self.route:
                self.stats['taps'] = missing_aps
        # elif 1:
            # # For now, ignore multi-segment routes
            # self.route = None
        else:
            self.route = None
            index = 0
            segments = self.segments
            for segment in segments:
                stitch_links = []
                if index > 0 and stitch_type == 'path':
                    # Stitch route between segments
                    path = network.solve2(db_conn, 
                                          segments[index - 1].gps.points[-1],
                                          segment.gps.points[0],
                                          {'length': 1.0})
                    if path:
                        for i in path.lines:
                            i.kwargs['seltype'] = 2
                        stitch_links = path.lines
                    else:
                        # no shortest path between segments
                        self.route = None
                        break
                elif index > 0 and stitch_type == 'gap':
                    # Leave gaps between segments to represent 
                    # unmatched travel.
                    pass
                if segment.route:
                    if not self.route:
                        self.route = Route()
                    # try to connect segments with shortest path
                    # todo : rewind to non-odd match link
                    # todo : take cost attribute as variable
                    self.route.extend(db_conn, segment.route, 
                                      stitch_links=stitch_links)
                    # below all handled in Route.extend
                    # dst = sum([i.length for i in stitch_links])
                    # self.route.stats['link_dst'] += dst
                    # self.route.links.extend(stitch_links)
                    # self.route.points.extend(segment.route.points)
                    # self.route.score += segment.route.score
                    # if index == 0:
                        # self.route.pos0 = segment.route.pos0
                    # self.route.pos1 = segment.route.pos1
                else:
                    # If any segment has no route, stage has no route
                    self.route = None
                    break
                index += 1
        if not transit and self.route:
            if (self.route.stats['aps'] < 25 and 
                self.route.stats['odd_links'] < 4):
                self.kwargs['mapmatch'] = 1
                self.route.get_geom(db_conn)
            elif not self.route:
                self.kwargs['mapmatch'] = 3
                self.geom = None
            else:
                self.kwargs['mapmatch'] = 2
                self.route.get_geom(db_conn)
                
    def split_segments(self, time_gap=120, dist_gap=500, 
                       attribute_points=True):
        """Divide Stage into segments based on GPS gaps.

        Segments represent uninterrupted slices of GPS points. Segments are
        split where gaps or potential changes of direction are detected. These
        are meant to capture the effects of GPS signal problems (e.g. tunnels)
        and multi-loop trips (e.g. home -> park -> home -> elsewhere). The
        multi-hypothesis map matching cannot handle repeated traversals 
        of the same link in the same direction.
          
        
        time_gap : maximum time (s) between GPS points before segment splits
        dist_gap : maximum distance (m) between GPS points before segment 
          splits
          
        """
             
        # todo : check if meaningful at split_stages defaults
        # todo : consider attributing points with cumulative distance
        self.segments = []
        # points = self.gps.points[:]
        points = self.gps.points
        # points.reverse()
        # pp = points.pop()  ## previous point
        pp = points[0]
        pp.angle = None
        queue = [pp]
        num = 1
        i = 1
        while i < len(points):
            p = points[i]
            # print p.angle
            # print p.kwargs['dheading']
            # raw_input()
            # todo : make dheading limit parameter
            # if (p.t - pp.t <= time_gap and 
                # ((p.x - pp.x) ** 2 + (p.y - pp.y) ** 2) ** 0.5 <= dist_gap and
                # (dheading(pp, p) < 175)):
            if (p.t - pp.t <= time_gap and 
                ((p.x - pp.x) ** 2 + (p.y - pp.y) ** 2) ** 0.5 <= dist_gap):
                queue.append(p)
            else:
                g = GPSCollection(srid=self.gps.srid, points=queue)
                # try: 
                s = Segment(num, gps=g)
                # except IndexError:
                    # print queue
                    # raw_input()
                self.segments.append(s)
                queue = [p]
                num += 1
            i += 1
            pp = p
        if queue:
            g = GPSCollection(srid=self.gps.srid, points=queue)
            s = Segment(num, gps=g)
            self.segments.append(s)
        if attribute_points:
            for seg in self.segments:
                for p in seg.gps.points:
                    p.kwargs['segment'] = seg.num

    def write_csv(self, path, fieldnames=None, user_dict=None, none_val=-999):
        """Append Stage as row to csv file.

        path : path to write csv
        fieldnames : optional ordered list of str, num, or datetime.datetime 
          fieldnames to write
        user_dict : user-supplied dict to include in fieldnames search

        Note that for now top-level attribute dict keys named in fieldnames 
        must be unambiguous across all attribute dicts, or an error will be 
        raised.

        """
        if fieldnames is None:
            fieldnames = ['hid', 'pid', 'phase_id', 'deploy_id', 'trip_num', 
                          'num', 'point_start', 'point_stop', 'npoints',
                          'local_datetime_start', 'local_time_stop', 
                          'time_start', 'time_stop', 'time_GPS_min', 
                          'lat_start', 'lng_start',
                          'lat_stop', 'lng_stop', 'mode', 'model',
                          'pw', 'pb', 'pa', 'pt', 'dist_GPS_mi', 
                          'dist_net_m', 'spd_avg_GPS_mph', 
                          'spd_avg_rte_mph', 'spd_max_GPS_mph',
                          'spd_med_GPS_ms', 'spd_var_GPS_ms', 'acc_nf_GPS_mss',
                          'taps', 'tmodes', 'med_vc', 'med_pc', 'med_hc',
                          'med_st', 'ncounts', 'valid_counts', 
                          'notworn_counts']
        super(Stage, self).write_csv(path, fieldnames=fieldnames, 
                                     user_dict=user_dict, none_val=-999)


class SegmentFile(object):
    """Segment as generated by Segment.write_csv()."""
    def __init__(self, path, col_segment_num='segment_id', 
                 col_time_start='time_start', 
                 col_time_stop='time_stop', col_unit_id='pid', 
                 col_trip_num='trip_id', col_stage_num='stage_id',
                 col_mode='mode'):
        self.segments = {}
        reader = csv.DictReader(open(path))
        for row in reader:
            unit_id = row[col_unit_id]
            time_start = float(row[col_time_start])
            time_stop = float(row[col_time_stop])
            if col_trip_num in row:
                trip_num = int(row[col_trip_num])
            else:
                trip_num = None
            if col_stage_num in row:
                stage_num = int(row[col_stage_num])
            else:
                stage_num = None
            if col_segment_num in row:
                segment_num = int(row[col_segment_num])
            else:
                segment_num = None
            if col_mode in row:
                mode = row[col_mode]
            else:
                mode = None
            self.update(unit_id, time_start, time_stop, trip_num, stage_num,
                        segment_num, mode)

    def update(self, unit_id, time_start, time_stop, trip_num, stage_num,
               segment_num, mode):
        """Update the segment profile."""

        if unit_id not in self.segments:
                self.segments[unit_id] = {}
        if trip_num not in self.segments[unit_id]:
            self.segments[unit_id][trip_num] = {}
        if stage_num not in self.segments[unit_id][trip_num]:
            self.segments[unit_id][trip_num][stage_num] = {}
        self.segments[unit_id][trip_num][stage_num][segment_num] = \
          {'time_start': time_start, 'time_stop': time_stop, 'mode': mode}


class RouteFile(object):
    """Map matched Route file"""

    def __init__(self, path):
        pass

    def update(self):
        """Update the Route profile."""
        pass
        
class Trip(Segment):
    """A collection of GPSPoints defining a trip."""
    def __init__(self, num, gps=None, am=None, stages=None, **kwargs):
        """Initialize a Trip instance.

        num - trip number, unique within Person Deployment [int]
        gps - GPSCollection
        stages - collection of Stage objects defining Trip [list]
        kwargs - {}

        """
        super(Trip, self).__init__(num, gps=gps, 
                                   am=am, **kwargs)

        self.minutes = (self.time_stop - self.time_start) / 60.0
        if stages is None:
            self.stages = []
        else:
            self.stages = stages
    
    def get_stats(self):
        """Collect stats from GPS, AM, route, and mode data."""
        super(Trip, self).get_stats()
        

class Activity(Segment):
    """Time and location period defining segments between Trips.

    num : activity number, unique within unit and phase [int]
    geom_start : PostGIS geometry at activity start [str]
    geom_start : PostGIS geometry at activity end [str]
    time_start : time at activity start [float]
    time_stop : time at activity end [float]
    gps : GPSCollection of points during activity
    
    Note that a single point defines trip start/activity end 
    and likewise for trip end/activity start so that a point 
    may be assigned both a trip and an activity.
    
    """
    def __init__(self, num, time_start, time_stop, gps=None, am=None, 
                 stats=None, modes=None, mode=None, trip=None, stage=None,
                 geom=None, db_conn=None, route=None, **kwargs):
        super(Activity, self).__init__(num, gps=gps, modes=modes, mode=mode, 
                                       am=am, trip=trip, stage=stage, 
                                       stats=stats, route=route, geom=geom,
                                       **kwargs)
        self.time_start = float(time_start)
        self.time_stop = float(time_stop)
        if 'local_time_start' in kwargs:
            self.local_time_start = kwargs['local_time_start']
        if 'local_time_stop' in kwargs:
            self.local_time_end = kwargs['local_time_stop']
        if self.gps:
            if not self.geom:
                self.gps.get_geom2(db_conn)
            self.gps.get_centroid(db_conn) 

def distance(g1, g2, db_conn):
    """Returns 2d cartesian distance between any two features with a geom
    attribute in projected units.

    """
    query = "SELECT ST_Distance(%s, %s);"
    data = [g1.geom, g2.geom]
    result = db_conn.run(query, data)
    return result[0][0]

def dist_pythag(g1, g2):
    """Returns planar distance between any two Points with xy coordinates
    in projected units

    """
    d = ((g2.x - g1.x) ** 2 + (g2.y - g1.y) ** 2) ** 0.5
    return d

def centroid(g, db_conn):
    """Returns centroid for multipoint object (maybe others)

    g - multipoint object with geom attribute

    """
    query = "SELECT ST_Centroid(%s);"
    data = [g.geom]
    geom = db_conn.run(query, data)[0][0]
    query = "SELECT ST_X(%s), ST_Y(%s);"
    data = [geom, geom]
    x, y = db_conn.run(query, data)[0]
    c = Point(geom=geom, x=x, y=y)
    return c

def collect(g_list, db_conn):
    """Returns multi geometry or geometry collection

    g_list - list of objects with geom attributes

    """
    g_list = [g.geom for g in g_list]
    query = "SELECT ST_Collect(%s);"
    data = [g_list]
    result = db_conn.run(query, data)
    #todo check this result
    return result[0][0]

def distance_matrix(points, directional=True):
    """Calculate an OD distance matrix using pythagorean distance.

    points - a list of points with unique ids
    directional - If False, then only upper diagonal will be calculated

    """
    # todo test speed improvement of directional=False and consider implementing
    d = {}
    for p1 in points:
        if p1.id in d:
            print 'Warning: multiple points with id {0}'.format(p1.id)
        d[p1.id] = {}
        for p2 in points:
            d[p1.id][p2.id] = dist_pythag(p1, p2)
    return d

def lines_interpolate_points(db_conn, lines, interval, output, 
                             endpoint='ignore'):
    """Sample points at regular intervals along line.

    db_conn : psql database connection
    lines : LineCollection
    interval : sampling interval in projected units
    # todo : Implement psql & csv output
    output : tuple of output type and location
      type - shp, psql, or csv
      location - (data_dir, name), schema.table, or file path
    endpoint : how to handle remainder of link near endpoint
      'ignore' - Ignore the remainder sections in sampling. This means the
        endpoint will likely not be included in sample. Lines shorter than
        <interval> will have only the startpoint sampled.
      # todo : Implement other endpoint options.
      'include' - Ignore interval, and include endpoint as sample point. This
        means samples will not always be <interval> distance apart.
      'average_directions' - Sample in both directions and average. Ensures
        aggregate stats symmetrical i n opposite directions

    """

    if output[0] == 'shp':
        os.listdir(output[1][0])
    # elif output[0] == 'csv':
       # f = open(output[1], 'wb')
       # writer = csv.DictWriter(f, fieldnames=['row', 'line_id', 'point_id'])
    else:
        err = 'Output type {} is not supported.'.format(output[0])
        raise Error(err)

    nlines = len(lines.lines)  ## for progress update 
    interval = float(interval)
    PC = PointCollection(srid=lines.srid) 
    n = 1
    nl = 1
    for line in lines.lines:
        if nl % 5000 == 0:
            print '{} of {}'.format(nl, nlines)
        geom = line.geom
        line_id = line.id
        np = 1
        try:
            r = interval / line.length  ## sample interval proportion
        except ZeroDivisionError:
            print 'Warning: zero length link {}'.format(line.id)
            r = 1.0
        pos = 0.0
        while pos <= 1:
            query = ("SELECT ST_X(ST_Line_Interpolate_Point(geometry(%s), %s)),"
                     "ST_Y(ST_Line_Interpolate_Point(geometry(%s), %s));")
            data = [geom, pos, geom, pos]
            x, y = db_conn.run(query, data)[0]
            PC.points.append(Point(x=x, y=y, id_point=n, line_id=line_id, 
                                   sample_n=np))
            pos += r
            n += 1
            np += 1
        nl += 1
    print 'Writing shapefile...'
    PC.write_shp(data_dir=output[1][0], name=output[1][1], 
                 attributes=[('line_id', 'N', 0), ('sample_n', 'N', 0)])
    print 'Wrote {} records to {}'.format(len(PC.points), output[1])

def nearest_neighbor(db_conn, point, neighbor_table, col_geom, 
                     col_attr=None, within=None):
    """Return nearest neighbor from psql table to given Point

    db_conn : psql database connection
    point : Point instance
    neighbor_table : schema.name of neighbor point table (needs index)
    col_geom : name of geometry column on neighbor table
    col_attr : optional second attribute to return
    within : max distance to assign neighbor

    This function uses the PostGIS <-> operator to make use of a spatial index
    on the neighbor table.

    """ 
    if not col_attr:
        col_attr = ''
    if within:
        query = ("SELECT ST_X({}),ST_Y({}){} FROM {} "
                 "WHERE ST_DWithin(geometry({}),geometry(%s),%s) "
                 "ORDER BY {} <-> geometry(%s) LIMIT 1;"
                 .format(col_geom, col_geom, ',' + col_attr, neighbor_table, 
                         col_geom, col_geom))
        data = [point.geom, within, point.geom]
    else:
        query = ("SELECT ST_X({}),ST_Y({}){} FROM {} "
                 "ORDER BY {} <-> geometry(%s) LIMIT 1;"
                 .format(col_geom, col_geom, ',' + col_attr, neighbor_table, 
                         col_geom))
        [point.geom, point.geom]

    result = db_conn.run(query, data)
    if len(result) > 0:
        return result[0]
    else:
        return None

def nearest_neighbor2(db_conn, geom, neighbor_table, col_geom, 
                      col_attr=None, max_candidates=100):
    """Return nearest neighbor distance and attribute from psql 
    table to given feature.

    db_conn : psql database connection
    geom : geometry of feature to neighbor
    neighbor_table : schema.name of neighbor point table (needs index)
    col_geom : name of geometry column on neighbor table
    col_attr : optional second attribute to return
    max_candidates : max number of approximate nearest neighbors to search

    This function uses the PostGIS <#> operator to make use of a spatial index
    on the neighbor table. <#> uses bounding box distances, so it isn't exact
    for anything but point to point neighbors. To overcome this, candidates are
    returned and then the nearest is returned from candidates.
    
    http://boundlessgeo.com/2011/09/indexed-nearest-neighbour-search-in-postgis/

    """ 
    
    query = ("with index_query as (select st_distance({}, %s) as dist"
             .format(col_geom))
    if col_attr:
        query += ",{}".format(col_attr)
    query += (" from {} order by {} <#> %s limit %s) "
              "select * from index_query order by dist limit 1;"
             .format(neighbor_table, col_geom))
    data = [geom, geom, max_candidates]
    result = db_conn.run(query, data)[0]
    return result
        
def point_slope(db_conn, p1, p2, attribute_points=True, skip=False):
    """Calculate slope between consecutive points.

    skip : not yet implemented

    """
    miss = 0
    d = dist_pythag(p1, p2)
    try:
        dz = p2.z - p1.z
    except TypeError:
        dz = None
    try:
        return dz / d
    except ZeroDivisionError:
        return None
    except TypeError:
        return None

def point_distance(points):
    """Return cumulative distance in for point sequence.
    
    points : sequence of Point-like objects with x and y values in a planar
      geometry
    
    """
    if len(points) > 1:
        index = 1
        d = 0
        for p in points[: 1]:
            pp = points[index - 1]
            d += ((p.x - pp.x) ** 2 + (p.y - pp.y) ** 2) ** 0.5
            index += 1
        return d
    elif len(points) == 1:
        return 0
    else:
        return None