"""Classes and functions for generating web maps.

By Joseph Broach <jbroach@pdx.edu>

"""

import pygmaps
    
class Gmap(pygmaps.maps):
    """Extends base class
    
    Define map extent to cover point features.
    Draw points as circles instead of pushpin markers.
    
    
    """
    def drawmap(self, f, extent):
        f.write('\t\tvar centerlatlng = new google.maps.LatLng(%f, %f);\n' % (self.center[0],self.center[1]))
        f.write('\t\tvar bounds = new google.maps.LatLngBounds();\n')
        f.write('\t\tvar extent = [new google.maps.LatLng(%s, %s), new google.maps.LatLng(%s, %s)];\n' % (extent[1], extent[0], extent[3], extent[2]))
        f.write(('\t\tfor (var i = 0; i < 2; i++) {var coords = extent[i];\n'
                 'bounds.extend (coords);}\n'))
        f.write('\t\tvar myOptions = {\n')
        f.write('\t\t\tzoom: %d,\n' % (self.zoom))
        f.write('\t\t\tcenter: centerlatlng,\n')
        f.write('\t\t\tmapTypeId: google.maps.MapTypeId.ROADMAP\n')
        f.write('\t\t};\n')
        f.write('\t\tvar map = new google.maps.Map(document.getElementById("map_canvas"), myOptions);\n')
        f.write('\t\tmap.fitBounds( bounds );\n')
        f.write('\n')
    
    def drawpoint(self,f,lat,lon,color,radius):
        f.write('\t\tvar latlng = new google.maps.LatLng(%f, %f);\n'%(lat,lon))
        f.write(('\t\tvar point = new google.maps.Circle('
        '{{strokeColor: "{color}", strokeOpacity: .8, strokeWeight: {radius},'
        'fillColor: "{color}", fillOpacity: .35, map: map,' 
        'center: latlng, radius: {radius}}});\n'.format(color=color, 
                                                        radius=radius)))
    
    def addpoint(self, lat, lng, color='#FF0000', radius=1):
        self.points.append((lat,lng,color,radius))
        
    def drawpoints(self,f):
        for point in self.points:
            self.drawpoint(f,point[0],point[1],point[2],point[3])
            
    def draw(self, htmlfile, extent, title='Untitled Map', next=''):
        f = open(htmlfile,'w')
        f.write('<html>\n')
        f.write('<head>\n')
        f.write('<meta name="viewport" content="initial-scale=1.0, user-scalable=no" />\n')
        f.write('<meta http-equiv="content-type" content="text/html; charset=UTF-8"/>\n')
        f.write('<title>{}</title>\n'.format(title))
        f.write('<script type="text/javascript" src="http://maps.google.com/maps/api/js?sensor=false"></script>\n')
        f.write('<script type="text/javascript">\n')
        f.write('\tfunction initialize() {\n')
        self.drawmap(f, extent)
        self.drawgrids(f)
        self.drawpoints(f)
        self.drawradpoints(f)
        self.drawpaths(f,self.paths)
        f.write('\t}\n')
        f.write('</script>\n')
        f.write('</head>\n')
        f.write('<body style="margin:0px; padding:0px;" onload="initialize()">\n')
        if next:
            f.write('\t<div id="navigate" style="width: 100%; text-align: center;">')
            f.write('<a href="{}">Next ></a></div>'.format(next))
        f.write('\t<div id="map_canvas" style="width: 100%; height: 100%;"></div>\n')
        f.write('</body>\n')
        f.write('</html>\n')		
        f.close()
