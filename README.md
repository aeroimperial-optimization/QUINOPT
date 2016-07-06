# QUINOPT (QUadratic INtegral OPTimisation)
An open-source add-on for YALMIP to solve optimisation problems with polynomial quadratic integral inequality constraints.

## Release
Version 1.0  
16 May 2016  

## Copyright
- Giovanni Fantuzzi (Department of Aeronautics, Imperial College London, UK. Email: gf910[at]ic.ac.uk)  
- Andrew Wynn (Department of Aeronautics, Imperial College London, UK. Email: a.wynn[at]imperial.ac.uk)
- Paul Goulart (Department of Engineering Science, University of Oxford, UK. Email: paul.goulart[at]eng.ox.ac.uk)
- Antonis Papachristodoulou (Department of Engineering Science, University of Oxford, UK. Email: antonis[at]eng.ox.ac.uk)

## System requirements

In order to use QUINOPT, you will need:

1. A working version of [YALMIP](http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.WhatIsYALMIP), the MATLAB optimization modelling software by J. L&ouml;fberg
2. A suitable SDP solver. Choices include [SeDuMi](https://github.com/sqlp/sedumi), [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), [SDPA](http://sdpa.sourceforge.net/), [Mosek](https://www.mosek.com/) (free for
    users in academia).

QUINOPT has been succesfully tested on MATLAB 7.6  (R2008a) and higher. If you have a different version of MATLAB, use at your own risk!

## Installation

A typical installation of QUINOPT requires the following steps:

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
