# QUINOPT (QUadratic INtegral OPTimisation)
An open-source add-on for YALMIP to solve optimisation problems with polynomial quadratic integral inequality constraints.

- [Latest release](#LatestRelease)
- [Copyright](#Copyright)
- [System requirements](#Requirements)
- [Installation](#Install)
- [Getting started](#GettingStarted)
- [How to cite](#Cite)

## Latest release<a name="LatestRelease"></a>
Version 1.0  
19 July 2016  

## Copyright<a name="Copyright"></a>
- Giovanni Fantuzzi (Department of Aeronautics, Imperial College London, UK. Email: gf910[at]ic.ac.uk)  
- Andrew Wynn (Department of Aeronautics, Imperial College London, UK. Email: a.wynn[at]imperial.ac.uk)
- Paul Goulart (Department of Engineering Science, University of Oxford, UK. Email: paul.goulart[at]eng.ox.ac.uk)
- Antonis Papachristodoulou (Department of Engineering Science, University of Oxford, UK. Email: antonis[at]eng.ox.ac.uk)

## System requirements<a name="Requirements"></a>

In order to use QUINOPT, you will need:

1. A working version of [YALMIP](http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.WhatIsYALMIP), the MATLAB optimization modelling software by J. L&ouml;fberg
2. A suitable SDP solver. Choices include [SeDuMi](https://github.com/sqlp/sedumi), [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), [SDPA](http://sdpa.sourceforge.net/), [Mosek](https://www.mosek.com/) (free for
    users in academia).

QUINOPT has been succesfully tested on MATLAB 7.10  (R2010a) and higher. If you have a different version of MATLAB, use at your own risk!

## Installation<a name="Install"></a>

To install QUINOPT:

1. Install YALMIP (download from [here](http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.Download) 
   and follow the [installation instructions](http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.Installation))
2. Install a semidefinite programming (SDP) solver compatible with YALMIP. A complete list of YALMIP-compatible SDP solvers can be found [here](http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.Solvers).  
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

To get started with QUINOPT, please look at the sample scripts provided the folder "examples/". A description of the problems being solved can be found in the document "examples.pdf".

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
if by citing the following paper:

G. Fantuzzi, A. Wynn, P. Goulart, A. Papachristodoulou, _Optimization 
with affine homogeneous quadratic integral inequality constraints_,
[arXiv:1607.04210v1 [math.OC]](http://arxiv.org/abs/1607.04210v1).
