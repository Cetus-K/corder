# corder

***

The corder ( correlation and order analyzor package ) was made, on July 17 2017, with the aim of analysis structures and its properties from the trajectory of the VASP's first-princiles molecular dynamics ( FPMD ) calculation.

The required files are follows:

| file | role |
|:-----|:-----|
| XDATCAR | trajectory of the VASP's FPMD calculation |
| param.in | parameters that user-controllable for analyze XDATCAR |


### _dependencies_

--

The corder performs with depends on these packages:

| file | role |
|:-----|:-----|
| includes | some auxiliary functions |
| linalg | linear algebra functions that overwritting its operators |
| [voro++](https://github.com/spatialfruitsalad/pomelo/tree/master/lib/voro%2B%2B) | package for calculate Voronoi diagram |
| [qhull](http://www.qhull.org) | package for create convex fulls by Qhull algorithm |

*Note that*, some lines were added into the voro++ source code for perform the corder and draw Voronoi diagrams which is formatted to the POV-Ray style.


### _installation_

--

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


# Supporting methods

***

The corder supports several geometric algorithms as follows:

| function | calculates |
|:-----|:-----|
| `gofr` | Pair correlation function |
| `sofq` | Static structure factor |
| `dfc`  | Mean-squared displacement for diffusion coefficient, and averaged velocity |
| `lboo` | Local bond-orientational order parameter |
| `bacf` | Bond-angle correlation function |
| `pofqw` | Bond-orientational probability distribution |
| `bofree` | Bond-orientational ( Landau ) free energy coefficient |
| `csform` | Distrubution of cluster's sharing formations w. r. t. its sharing type |
| `povoro` | Make trajectory animtion of Voronoi cell |


### _For related to the atomic specie_

--

Each functions are represented by the total, partial and "specified"-partial component. For example, put a function <img src="https://latex.codecogs.com/gif.latex?f(x)"/>, the each component are expressed as follows:

| symbol | component |
|:-----|:-----|
| <img src="https://latex.codecogs.com/gif.latex?f(x)"/> | total |
| <img src="https://latex.codecogs.com/gif.latex?f^{\alpha\beta}(x)"/> | partial by <img src="https://latex.codecogs.com/gif.latex?\alpha"/>-<img src="https://latex.codecogs.com/gif.latex?\beta"/> |
| <img src="https://latex.codecogs.com/gif.latex?f^{AB}(x)"/> | partial by <img src="https://latex.codecogs.com/gif.latex?A"/>-<img src="https://latex.codecogs.com/gif.latex?B"/>, as user defined; <img src="https://latex.codecogs.com/gif.latex?e.g.\&space;A=\alpha_{1}\beta_{2},\&space;B=\alpha_{2}\beta_{1}\beta_{3}"/> |

Since this term, I abbreviate the symbol of each component.

### _Pair correlation function_

--

Correlations between each two-atoms-pair related to positions are directly calculated by

<img src="https://latex.codecogs.com/gif.latex?g(r)=\frac{n(r)}{4\pi&space;r^{2}\rho\,{\rm&space;d}r}"/>.

The <img src="https://latex.codecogs.com/gif.latex?r"/> indicates radius between two atoms with unit of <img src="https://latex.codecogs.com/gif.latex?\AA"/>, and <img src="https://latex.codecogs.com/gif.latex?n(r)"/> means the number of other atoms at the distance <img src="https://latex.codecogs.com/gif.latex?r"/> from the center atom.

### _Static structure factor_

--

The spectrums of the interference between neutrons in the structure are calculated by

<img src="https://latex.codecogs.com/gif.latex?s(q)=1+4\pi\int{\rm&space;d}r\,\{\rho(r)-\rho_{0}\}r^{2}\frac{\sin(qr)}{qr}"/>,

This function is derived by integrate, equal to the Fourier transform of, the pair correlation function as liquid state. Where the <img src="https://latex.codecogs.com/gif.latex?q"/> means the norm of the wave vector with unit of <img src="https://latex.codecogs.com/gif.latex?\AA^{-1}"/>. The <img src="https://latex.codecogs.com/gif.latex?\rho(r)=\langle\,\rho_{0}(t)\,g(r,t)\,\rangle"/>, and the average atomic density of the structure <img src="https://latex.codecogs.com/gif.latex?\rho_{0}"/>. The angular bracket indicates the time averaging related to the instantaneous variable <img src="https://latex.codecogs.com/gif.latex?t"/>.


### _Mean-squared displacement and averaged velocity_

--

The mean-squared displacement <img src="https://latex.codecogs.com/gif.latex?\delta^{2}"/>, defined by its each instantaneous location from <img src="https://latex.codecogs.com/gif.latex?\vec{r}(t)"/> to <img src="https://latex.codecogs.com/gif.latex?\vec{r}(t+T)"/>, derives its diffusion coefficient <img src="https://latex.codecogs.com/gif.latex?D"/> using the Einstein equation:

<img src="https://latex.codecogs.com/gif.latex?\delta^{2}=6tD"/>.


### _Local bond-orientational order parameter_

--

The ligand forms a geometric shape for an atom, which can be represented by the identical values <img src="https://latex.codecogs.com/gif.latex?q_{\ell},\&space;w_{\ell}"/> that is likes the spectrum as follows:

<img src="https://latex.codecogs.com/gif.latex?q_{\ell}=\sqrt{\frac{4\pi}{2\ell+1}\sum_{|m|\leq\ell}|q_{\ell&space;m}|^{2}},\&space;w_{\ell}=\sum_{m_{1}+m_{2}+m_{3}=1}\begin{pmatrix}\&space;\ell&\ell&\ell\&space;\\&space;\&space;m_{1}&m_{2}&m_{3}\&space;\end{pmatrix}\frac{q_{\ell&space;m_{1}}q_{\ell&space;m_{2}}q_{\ell&space;m_{3}}}{\left(\sum_{|m|\leq\ell}|q_{\ell&space;m}|^{2}\right)^{3/2}}"/>,

where,

<img src="https://latex.codecogs.com/gif.latex?q_{\ell&space;m}=\sum_{f\in\mathcal{F}}\frac{S_{f}}{S}Y_{\ell&space;m}^{*}({\it\hat{\Omega}}_{\!f}),\&space;S=\sum_{f\in\mathcal{F}}S_{f}"/>.

The coefficients in the summation of <img src="https://latex.codecogs.com/gif.latex?w_{\ell}"/> is the Wigner <img src="https://latex.codecogs.com/gif.latex?3j"/> symbols, and the <img src="https://latex.codecogs.com/gif.latex?{\it\hat{\Omega}}_{\!f}"/> is the normal vector of the surface <img src="https://latex.codecogs.com/gif.latex?f"/>.
For robust calculating of the cluster symmetries, there is the method that takes summation with Voronoi facet area <img src="https://latex.codecogs.com/gif.latex?{S}_{f}"/>, and the facets <img src="https://latex.codecogs.com/gif.latex?\mathcal{F}"/>, as a weighting parameter. 
These values <img src="https://latex.codecogs.com/gif.latex?q_{\ell}"/> and <img src="https://latex.codecogs.com/gif.latex?w_{\ell}"/> are rotationally invariant so it give us the interpretation for the distinguish cluster symmetry in any coordinates and also perspectives.


### _Bond-angle correlation function_

--

The cluster symmetries are correlated to the other symmetries, and the sustainability directly reflects its order as the correlation length. Also, the length can statistically estimate that the clusters are connected by being tied in a row. The definition is follows:

<img src="https://latex.codecogs.com/gif.latex?G_{\ell}(r)=\frac{4\pi}{2\ell+1}\sum_{|m|\leq\ell}\frac{\langle\,Y_{\ell&space;m}^{*}({\it\hat{\Omega}}_{\!f})Y_{\ell&space;m}({\it\hat{\Omega}}_{\!f'})\,\rangle(r)}{G_{0}(r)}"/>,

where,

<img src="https://latex.codecogs.com/gif.latex?r=\lVert\,\vec{r}_{f}-\vec{r}_{f'}\,\rVert,\&space;G_{0}(r)=4\pi\langle\,Y_{00}^{*}({\it\hat{\Omega}}_{\!f})Y_{00}({\it\hat{\Omega}}_{f'})\,\rangle(r),"/>

<img src="https://latex.codecogs.com/gif.latex?\langle\,Y_{\ell&space;m}^{*}({\it\hat{\Omega}}_{\!f})Y_{\ell&space;m}({\it\hat{\Omega}}_{f'})\,\rangle(r)=\sum_{f,f'}\frac{S_{f}S_{f'}}{S_{ff'}}Y_{\ell&space;m}^{*}({\it\hat{\Omega}}_{\!f})Y_{\ell&space;m}({\it\hat{\Omega}}_{f'})"/>.

This angular bracket takes ensamble average, and swapping the order <img src="https://latex.codecogs.com/gif.latex?m"/>-summation derives follows by generalized addition theorem, this package actually evaluate it:

<img src="https://latex.codecogs.com/gif.latex?G_{\ell}(r)=\frac{\langle\,P_{\ell}(z_{ff'})\,\rangle(r)}{\langle\,P_{0}(z_{ff'})\,\rangle(r)}"/>,

where <img src="https://latex.codecogs.com/gif.latex?z_{ff'}=\cos(\omega_{ff'})"/> indicates the inner product between bonds which parallel with each normal vector. This returns the value of the <img src="https://latex.codecogs.com/gif.latex?\ell"/>-th order correlation related to angular <img src="https://latex.codecogs.com/gif.latex?\omega_{ff'}"/>.


### _Bond-orientational ( Landau ) free energy coefficient_

--

Structures have a "geometric" free energy caused by cluster symmetries, and bond-orientational order describes that characteristics as the order parameter of the Landau theory. This function calculates internal and primary external free energy:

<img src="https://latex.codecogs.com/gif.latex?f_{\ell}=I_{\ell}+\sum_{\ell'\neq\ell}J_{\ell'\ell},"/>

where,

<img src="https://latex.codecogs.com/gif.latex?I_{\ell}=\sum_{m_{1}+m_{2}+m_{3}=1}\begin{pmatrix}\&space;\ell&\ell&\ell\&space;\\&space;\&space;m_{1}&m_{2}&m_{3}\&space;\end{pmatrix}q_{\ell&space;m_{1}}q_{\ell&space;m_{2}}q_{\ell&space;m_{3}},"/>

<img src="https://latex.codecogs.com/gif.latex?J_{\ell'\ell}=\sum_{m_{1}+m_{2}+m_{3}=1}\begin{pmatrix}\&space;\ell'&\ell'&\ell\&space;\\&space;\&space;m_{1}&m_{2}&m_{3}\&space;\end{pmatrix}q_{\ell'&space;m_{1}}q_{\ell'&space;m_{2}}q_{\ell&space;m_{3}}."/>

The full Landau free energy is constructed from two more additional terms:

<img src="https://latex.codecogs.com/gif.latex?K_{\ell'\ell''\ell}=\sum_{m_{1}+m_{2}+m_{3}=1}\begin{pmatrix}\&space;\ell'&\ell''&\ell\&space;\\&space;\&space;m_{1}&m_{2}&m_{3}\&space;\end{pmatrix}q_{\ell'&space;m_{1}}q_{\ell''&space;m_{2}}q_{\ell&space;m_{3}},"/>

<img src="https://latex.codecogs.com/gif.latex?K_{\ell'\ell''\ell}^{p}=\sum_{m_{1}+m_{2}+m_{3}=1}\begin{pmatrix}\&space;\ell'&\ell''&\ell\&space;\\&space;\&space;m_{1}&m_{2}&m_{3}\&space;\end{pmatrix}q_{\ell'&space;m_{1}}q_{\ell''&space;m_{2}}q_{\ell&space;m_{3}}^{*}."/>

The first and second of additional term are secondary external free energy and periodic external free energy. I assume that the first term is usually small value and second term is always zero due to the <img src="https://latex.codecogs.com/gif.latex?{\it\Gamma}"/>-point sampling during the FPMD calculation.

### _Distrubution of cluster's sharing formations_

--

Clusters are linked to each other in the system as various shape, and then there are some sharing part on the clusters. These sharing forms are categorized as follows:

| form | sharing part |
|:-----|:-----|
| isolated | nothing |
| corner sharing | a vertex |
| edge sharing | an edge |
| surface sharing | a surface |
| bicap sharing | some part of volume |

This function returns the distribution of these amount.
