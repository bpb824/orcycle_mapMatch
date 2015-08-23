--drop and re-create tables
drop table if exists data.orc_GPSroutes;

create table data.orc_GPSroutes (
gid serial primary key,
gpstripid	Integer,
stageid	Integer,
aps	double precision,
oddlinks	Integer,
distancelink	double precision,
distancepoint	double precision,
geom 	geometry(MultiLinestring, 2838)
);

drop table if exists data.orc_GPSlinks;

create table data.orc_GPSlinks (
gid serial primary key,
gpstripid	Integer,
stageid	Integer,
featureorder	Integer,
featureid	integer,
direction	Integer,
startpoint	Integer,
endpoint	Integer,
mps	Double precision,
odd	Integer,
startpos	Double precision,
endpos	Double precision
);

drop table if exists data.orc_GPStrips;

create table data.orc_GPStrips (
gid serial primary key,
gpstripid	Integer,
stageid		Integer,
starttime       timestamp,
endtime         timestamp,
startepoch	integer,
endepoch	integer,
startpoint	integer,
endpoint	integer,
duration	double precision,
npoints		integer,
pointdistance	double precision,
airdistance	double precision,
avgspeed	double precision,
avgpointspeed	double precision,
maxpointspeed	double precision,
startlong	double precision,
startlat	double precision,
endlong	double precision,
endlat	double precision,
indirectness double precision,
geom geometry(MultiPoint, 2838)
)
