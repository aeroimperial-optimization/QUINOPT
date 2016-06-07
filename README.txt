% ======================================================================= %
                          WELCOME TO QUINOPT!
% ======================================================================= % 

A free add-on for YALMIP to solve optimisation problems with polynomial 
quadratic integral inequality constraints.



% ======================================================================= %
                           RELEASE DETAILS
% ======================================================================= % 

QUINOPT Version 1.0
Release date: 18/05/2016

Copyright:    
- Giovanni Fantuzzi (Department of Aeronautics, Imperial College London)  
- Andrew Wynn (Department of Aeronautics, Imperial College London)
- Paul Goulart (Department of Engineering Science, University of Oxford)
- Antonis Papachristodoulou (Department of Engineering Science, University 
of Oxford)



% ======================================================================= %
                              CONTENTS
% ======================================================================= %

(1) Introduction
(2) System requirements
(3) Licence
(4) Installation
(5) How to use QUINOPT
(6) References
(7) User support



% ======================================================================= %
                         (1) INTRODUCTION
% ======================================================================= %

Welcome to QUINOPT, a free add-on for YALMIP to solve optimisation problems
with polynomial quadratic integral inequality constraints.



% ======================================================================= %
                      (2) SYSTEM REQUIREMENTS
% ======================================================================= %

In order to use QUINOPT, you will need:

1. A working version of YALMIP, the MATLAB optimization modelling software by
   J. L�fberg (http://users.isy.liu.se/johanl/yalmip/)

2. A YALMIP.-compatible SDP solver.

QUINOPT has been succesfully tested on MATLAB 7.10 (R2010a) and higher with
the following common SDP solvers:

	 * SeDuMi v1.3   (free) 			   http://sedumi.ie.lehigh.edu	
	 * SDPT3  v4.0   (free) 			   http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
	 * SDPA   v7.3.8 (free) 			   http://sdpa.sourceforge.net
	 * Mosek  v7.0   (free for academia)   https://www.mosek.com	 

If you have a different version of MATLAB or a different SDP solver, use at
your own risk!


% ======================================================================= %
                            (3) LICENCE
% ======================================================================= %
QUINTOPT is distrubuted under the terms of the GNU Lesser General Public 
Licence (LGPL) v3.0. You should have received a file named "LICENCE.txt." 
containing the licence terms; if not, you can find the licence at

http://www.gnu.org/licenses/lgpl-3.0.en.html

or

https://opensource.org/licenses/LGPL-3.0



% ======================================================================= %
                         (4) INSTALLATION
% ======================================================================= %
A typical installation of QUINOPT requires the following steps:

1. Install YALMIP. You can download YALMIP from 

	http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.Download
	
	and follow the installation instructions at 
	
	http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.Installation
	

2. Install a semidefinite programming (SDP) solver compatible with YALMIP. 
   A complete list of YALMIP-compatible SDP solvers can be found at 
   
   http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.Solvers
   
   For a list of solvers that have been tested with QUINOPT, please see 
   Section (2) of this README file.


3. Install QUINOPT by running the MATLAB installer:

	>> installQUINOPT


The installer should add the required folders to the MATLAB path and run 
some test problems to make sure everything is working. Please report any
installation problems to Giovanni Fantuzzi (gf910[at]ic.ac.uk).



% ======================================================================= %
                        (5) HOW TO USE QUINOPT
% ======================================================================= %

To get started with QUINOPT, please look at the sample scripts provided
the folder "resources". A description of the problems solved by the scripts
is given in the reference paper and in "AdditionalExamples.pdf".

For more information on the main functions in QUINOPT, type

>> help quinopt
>> help indvar
>> help depvar
>> help parameters

at the MATLAB command prompt.



% ======================================================================= %
                           (6) REFERENCES
% ======================================================================= %

If you find QUINOPT useful, or have used it in your own work, please reference
if by citing the following paper:

    G. Fantuzzi, A. Wynn, P. Goulart, A. Papachristodoulou, "Solution of 
    polynomial integral inequalities using semidefinite programming", in 
    preparation.


Since QUINOPT is an add-on for YALMIP, we also recommend that you consult 
the following additional references:

[2] J. L�fberg., "YALMIP: A Toolbox for Modeling and Optimization in MATLAB",
    Proceedings of the CACSD Conference, Taipei, Taiwan, 2004.

[3] J. L�fberg, "Pre- and post-processing sum-of-squares programs in practice",
    IEEE Transactions on Automatic Control, 54(5):1007-1011, 2009.



% ======================================================================= %
                           (7) USER SUPPORT
% ======================================================================= %

QUINOPT is distributed as free software in the hope that it will be useful.
The software is provided "as is", and although we have tested it, we cannot
guarantee its functionality. If you need specific help with QUINOPT or find 
a bug, please contact:

gf910@ic.ac.uk



% ======================================================================= %
                    Thanks for downloading QUINOPT!
% ======================================================================= %