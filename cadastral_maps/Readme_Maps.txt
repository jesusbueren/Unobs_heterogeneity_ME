The 4 zipfiles contain:

1.	village_plots_kmz.zip has one kmz file for each map. you should be able to open in googleearth by doubleclick.

2.	tbl.zip has one csv file for each map, with following attributes:

		a.	plot – plot number from map, multiples where the plot has been subdivided
		b.	sub or subplot – usually not populated, but in some cases is subplot id
		c.	flag – notes from the digitizing
		d.	pid – unique (within map) id, links to neighbor matris
		e.	area_ac – area of plot in acres
		f.	bnd_km – distance to map edge (outer polygon boundary), zero for edge plots. you should be able to use this to identify your ‘center’ plots

3.	nbs.zip has one file for each map containing the neighbor matrix. rowname and colname are the plot id (pid), populated with 1 for neighbors, else 0

4.	shp.zip – the polygon shapefiles of plots

Codes for villages number (map_village variable in _match_map database):

Ayyavaripalli			1
Dharmapur			2
Itukalapalli			3
Jajapur				4
M. Venkata Puram		5
Manesamudram			6
Muddireddy Palli		7
Pamireddypalli			8
Ramachandrapuram		9
Reddipalli			10
Siddarampuram			11
Suddakuntapalli			12
Thipparasipalli			13
Y.B. Halli			14
