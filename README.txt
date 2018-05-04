% ======================================================================= %
                          WELCOME TO QUINOPT!
% ======================================================================= %


Version 2.2
04 May 2018

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
(6) How to cite
(7) References
(8) User support



% ======================================================================= %
                         (1) INTRODUCTION
% ======================================================================= %

QUINOPT is an open-source add-on for YALMIP to solve optimisation problems
with polynomial quadratic integral inequality constraints.



% ======================================================================= %
                      (2) SYSTEM REQUIREMENTS
% ======================================================================= %

In order to use QUINOPT, you will need:

1. A working version of YALMIP, the MATLAB optimization software by J.
   Löfberg (http://users.isy.liu.se/johanl/yalmip/)

2. A YALMIP-compatible SDP solver.

Instructions on how to obtain YALMIP and an SDP solver are given in Section (4)
of this README file.

QUINOPT has been succesfully tested on MATLAB 7.10 (R2010a) and higher with
the following common SDP solvers:

  * SeDuMi v1.3        http://sedumi.ie.lehigh.edu
  * SDPT3  v4.0        http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
  * SDPA   v7.3.8      http://sdpa.sourceforge.net
  * Mosek  v7.0        https://www.mosek.com

If you have a different version of MATLAB or a different SDP solver, use at
your own risk!



% ======================================================================= %
                            (3) LICENCE
% ======================================================================= %

QUINOPT is distributed under the following licence:

Copyright 2016, G. Fantuzzi, A. Wynn, P. Goulart, A. Papachristodoulou

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.



% ======================================================================= %
                         (4) INSTALLATION
% ======================================================================= %

To install QUINOPT:

1. Install YALMIP. You can download YALMIP from

	http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.Download

   and follow the installation instructions at

	http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.Installation


2. Install a solver for semidefinite programs (SDPs) compatible with YALMIP.
   A complete list of YALMIP-compatible SDP solvers can be found at

   	http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.Solvers

   For a list of solvers that have been tested with QUINOPT, please see
   Section (2) of this README file.


3. Install QUINOPT by running the MATLAB installer:

	>> installQUINOPT


The installer should compile the required files, add the required folders to the
MATLAB path, and run some test problems to make sure everything is working.
Please report any installation issues to Giovanni Fantuzzi (gf910[at]ic.ac.uk).

IMPORTANT NOTES
---------------
1. During installation, you may receive the following warning:

        "Warning: Compilation of mex files by installQUINOPT failed.
         QUINOPT will still work without compiled mex files, but
         it will be slower. To resolve the issue, make sure that
         a supported compiler is installed and re-run the installer."

   QUINOPT should still work, but you may wish to resolve the issue with the
   mex file compilation. You can find a list of supported compilers for
   MATLAB's latest version at

   https://uk.mathworks.com/support/compilers.html

   For all other versions of MATLAB please consult

   https://uk.mathworks.com/support/sysreq/previous_releases.html

2. QUINOPT has been tested with the following SDP solvers: SeDuMi, SDPT3, SDPA,
   and Mosek (free for users in academia). Any other YALMIP-compatible SDP
   solver should work, but use them at your own risk!



% ======================================================================= %
                        (5) HOW TO USE QUINOPT
% ======================================================================= %

To get started with QUINOPT, please look at the sample scripts provided
the folder "examples". A step-by-step description of the examples can be found
in the online documentation, available at:

http://quinopt.readthedocs.io/04_examples/index.html

For more information on the main functions in QUINOPT, check the online
documentation or type

>> help quinopt
>> help indvar
>> help depvar
>> help parameters

at MATLAB's command prompt.



% ======================================================================= %
                           (6) HOW TO CITE
% ======================================================================= %

If you find QUINOPT useful, or have used it in your own work, please reference
it by citing the following papers:

* G. Fantuzzi, A. Wynn (2016). Semidefinite relaxation of a class of quadratic
 integral inequalities. In: Proceedings of the 55th IEEE Conference on Decision
 and Control. DOI: 10.1109/CDC.2016.7799221.

* G. Fantuzzi, A. Wynn, P. Goulart, A. Papachristodoulou (2017). Optimization
 with affine homogeneous quadratic integral inequality constraints, IEEE
 Transactions on Automatic Control 62(12), 6221-6236.
 DOI: 10.1109/TAC.2017.2703927.

Should you wish to cite the code directly, please use the following BibTeX
entry, replacing X.X.X with the appropriate version of QUINOPT:

@misc{QUINOPTvX.X.X,
    author = {Fantuzzi, Giovanni and Wynn, Andrew and Goulart, Paul and Papachristodoulou, Antonis},
    title = {{QUINOPT}, version X.X.X},
    howpublished = {\url{https://github.com/aeroimperial-optimization/QUINOPT}},
    year = 2017
    }



% ======================================================================= %
                           (7) REFERENCES
% ======================================================================= %

Since QUINOPT is an add-on for YALMIP, we also recommend that you consult
the following additional references:

[2] J. Löfberg., "YALMIP: A Toolbox for Modeling and Optimization in MATLAB",
    Proceedings of the CACSD Conference, Taipei, Taiwan, 2004.

[3] J. Löfberg, "Pre- and post-processing sum-of-squares programs in practice",
    IEEE Transactions on Automatic Control, 54(5):1007-1011, 2009.



% ======================================================================= %
                           (8) USER SUPPORT
% ======================================================================= %

QUINOPT is distributed as free software in the hope that it will be useful.
The software is provided "as is", and although we have tested it, we cannot
guarantee its functionality. If you need specific help with QUINOPT or find
a bug, please contact:

giovanni.fantuzzi10[at]ic.ac.uk



% ======================================================================= %
                    Thanks for downloading QUINOPT!
% ======================================================================= %
