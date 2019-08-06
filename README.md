# corder #

---

The corder ( correlation and order analyzor package ) was made, on July 17 2017, with the aim of analysis structures and its properties from the trajectory of the VASP's first-princiles molecular dynamics ( FPMD ) calculation.

The required files are follows:

| file | role |
|:-----|:-----|
| XDATCAR | trajectory of the VASP's FPMD calculation |
| param.in | parameters that user-controllable for analyze XDATCAR |


### _dependencies_ ###

-

The corder performs with depends on these packages:

| file | role |
|:-----|:-----|
| includes | some auxiliary functions |
| linalg | linear algebra functions that overwritting its operators |
| [voro++](https://github.com/spatialfruitsalad/pomelo/tree/master/lib/voro%2B%2B) | package for calculate Voronoi diagram |
| [qhull](http://www.qhull.org) | package for create convex fulls by Qhull algorithm |

*Note that*, some lines were added into the voro++ source code for perform the corder and draw Voronoi diagrams which is formatted to the POV-Ray style.


### _installation_ ###

-

Download and install:
```
git clone git@github.com:Cetus-K/corder.git
cd ./corder
./INSTALL.sh
```
*Note that*, the install prefix is set to current directory. If you specify the install location, change this line:
```
export PREFIX=/path/to/install/location
```


# Supporting methods #

---

The corder supports several geometric algorithms as follows:

| function | calculates |
|:-----|:-----|
| `gofr` | pair correlation function |
| `sofq` | structure factor |
| `dfc`  | mean-squared displacement for diffusion coefficient, and velocity correlation |
| `lboo` | local bond-orientational order parameter |
| `bacf` | bond-angle correlation function |
| `pofqw` | bond-orientational probability distribution |
| `bofree` | bond-orientational ( Landau ) free energy |
| `csform` | distrubution of cluster's sharing formations w. r. t. its sharing type |
| `povoro` | make trajectory animtion of voronoi cell |

### _pair correlation function_ ###

-

Correlations between each atom related to positions are directly calculated by

<div style="text-align: center;">
<img src="https://latex.codecogs.com/gif.latex?g(r),\&space;g_{\alpha\beta}(r),\&space;g_{AB}(r)">.
</div>

The <img src="https://latex.codecogs.com/gif.latex?r"> indicates radius between two atoms with unit of <img src="https://latex.codecogs.com/gif.latex?\AA">, suffixes <img src="https://latex.codecogs.com/gif.latex?\alpha,\&space;\beta,\&space;A,\&space;B"> mean atomic specie with respected to. Each <img src="https://latex.codecogs.com/gif.latex?A"> and <img src="https://latex.codecogs.com/gif.latex?B"> is regarded as the identical specified, by users.

### _<font color="OrangeRed">This document is in the middle way of my writting..., partially completed.</font>_ ###


### _structure factor_ ###

-


### _mean-squared displacement and velocity correlation_ ###

-


### _local bond-orientational order parameter_ ###

-

### _bond-angle correlation function_ ###

-

### _bond-orientational probability distribution_ ###

-

### _bond-orientational ( Landau ) free energy_ ###

-

### _distrubution of cluster's sharing formations_ ###

-
