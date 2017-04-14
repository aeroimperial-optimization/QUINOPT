# QUINOPT (QUadratic INtegral OPTimisation)
An open-source add-on for YALMIP to solve optimisation problems with polynomial quadratic integral inequality constraints. Below is a quick guide to QUINOPT, but details, examples, and much more can be found in the [full online documentation](http://quinopt.readthedocs.io/).

**Latest release:** 1.4  
**Release date:** 02 March 2017  
**Known bugs in version 1.4:** No known bugs (yet!)

## Contents
- [System requirements](#Requirements)
- [Installation](#Install)
- [Getting started](#GettingStarted)
- [How to cite](#Cite)
- [Copyright](#Copyright)

## System requirements<a name="Requirements"></a>

In order to use QUINOPT, you will need:

1. A working version of [YALMIP](https://yalmip.github.io/), the MATLAB optimization modelling software by J. L&ouml;fberg
2. A suitable SDP solver. Choices include [SeDuMi](https://github.com/sqlp/sedumi), [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), [SDPA](http://sdpa.sourceforge.net/), [Mosek](https://www.mosek.com/) (free for
    users in academia).

QUINOPT has been succesfully tested on MATLAB 7.10  (R2010a) and higher. If you have a different version of MATLAB, use at your own risk!

## Installation<a name="Install"></a>

To install QUINOPT:

1. [Download](https://yalmip.github.io/download/) and [install](https://yalmip.github.io/tutorial/installation/) YALMIP 
2. Install a semidefinite programming (SDP) solver compatible with YALMIP. [Click here for a complete list of YALMIP-compatible SDP solvers](https://yalmip.github.io/allsolvers/).  
3. Install QUINOPT by running the MATLAB installer:

```Matlab
>> installQUINOPT
```

The installer should add the required folders to the MATLAB path and run some test problems to make sure everything is working.
Please report any installation problems to Giovanni Fantuzzi (gf910[at]ic.ac.uk).

_**NOTE:** QUINOPT has been tested with [SeDuMi](https://github.com/sqlp/sedumi), 
  [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), 
  [SDPA](http://sdpa.sourceforge.net/) and 
  [Mosek](https://www.mosek.com/) (free for users in academia). 
  Any other YALMIP-compatible SDP solver should work, but use at your own risk!_
  
## Getting started<a name="GettingStarted"></a>

To get started with QUINOPT, please look at the sample scripts provided the folder "./examples/". A description of the problems being solved can be found in the document "examples.pdf".

For more information on the main functions in QUINOPT, type

```Matlab
>> help quinopt
>> help indvar
>> help depvar
>> help parameters
```

at the MATLAB command prompt.


## How to cite<a name="Cite"></a>
  
If you find QUINOPT useful, or have used it in your own work, please reference
it by citing the following papers:

* G. Fantuzzi, A. Wynn, _Semidefinite relaxation of a class of quadratic
 integral inequalities_, [55th IEEE Conference on Decision and Control, 2016]
 (http://dx.doi.org/10.1109/CDC.2016.7799221).
 
 ```
@inproceedings{FW2016CDC,
    address = {Las Vegas, USA},
    author = {Fantuzzi, G. and Wynn, A.},
    booktitle = {Proc. 55th IEEE Conf. Decis. Control},
    doi = {10.1109/CDC.2016.7799221},
    pages = {6192--6197},
    publisher = {IEEE},
    title = {{Semidefinite relaxation of a class of quadratic integral inequalities}},
    url = {http://dx.doi.org/10.1109/CDC.2016.7799221},
    volume = {2},
    year = {2016}
	}
 ```
 
* G. Fantuzzi, A. Wynn, P. Goulart, A. Papachristodoulou, _Optimization 
with affine homogeneous quadratic integral inequality constraints_,
[arXiv:1607.04210 [math.OC]](https://arxiv.org/abs/1607.04210#).

 ```
 @article{FWGP2016,
	archivePrefix = {arXiv},
	eprint = {1607.04210},
	primaryClass = "math-OC",
	author = {Fantuzzi, Giovanni and Wynn, Andrew and Goulart, Paul and Papachristodoulou, Antonis},
	title = {{Optimization with affine homogeneous quadratic integral inequality constraints}}
	}
 ```

A selection of BibTex styles that support arXiv preprints can be found [here](http://arxiv.org/hypertex/bibstyles/).
Should you wish to cite the code directly, please use the following BibTeX entry:

```
@misc{QUINOPTv1.4,
    author       = {Fantuzzi, Giovanni and Wynn, Andrew and Goulart, Paul and Papachristodoulou, Antonis},
    title        = {{QUINOPT}, version 1.4},
    howpublished = {\url{https://github.com/aeroimperial-optimization/QUINOPT}},
    month        = Mar,
    year         = 2017
    }
```

## Copyright<a name="Copyright"></a>
- Giovanni Fantuzzi (Department of Aeronautics, Imperial College London, UK. Email: gf910[at]ic.ac.uk)  
- Andrew Wynn (Department of Aeronautics, Imperial College London, UK. Email: a.wynn[at]imperial.ac.uk)
- Paul Goulart (Department of Engineering Science, University of Oxford, UK. Email: paul.goulart[at]eng.ox.ac.uk)
- Antonis Papachristodoulou (Department of Engineering Science, University of Oxford, UK. Email: antonis[at]eng.ox.ac.uk)
