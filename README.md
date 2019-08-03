
# corder

The corder ( correlation and order analyzor package ) was made, on July 17 2017, with the aim of analysis structures and its properties from the trajectory of the VASP's first-princiles molecular dynamics ( FPMD ) calculation.

The required files are follows:

| file | role |
|:-----|:-----|
| XDATCAR | trajectory of the VASP's FPMD calculation |
| param.in | parameters that user-controllable for analyze XDATCAR |


### _dependencies_ ###

The corder performs with depends on these packages:

| file | role |
|:-----|:-----|
| includes | some auxiliary functions |
| linalg | linear algebra functions that overwritting its operators |
| voro++ | package for calculate Voronoi diagram |
| qhull | package for create convex fulls by Qhull algorithm |

Where, the packages voro++ was cloned from [GitHub](https://github.com/spatialfruitsalad/pomelo/tree/master/lib/voro%2B%2B), and qhull was download on its [site](http://www.qhull.org).

*Note that*, some lines were added into the voro++ source code for perform the corder and draw Voronoi diagrams which is formatted to the POV-Ray style.


### _installation_ ###

Download and install:
```
tar -zxvf corder-1.0.0.tar.gz
cd ./corder
./INSTALL.sh
```
*Note that*, the install prefix is set to current directory. If you specify the install location, change this line:
```
export PREFIX=/path/to/install/location
```
