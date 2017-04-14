Install QUINOPT
===============

QUINOPT is easily installed by running the installer ``installQUINOPT.m`` in MATLAB. Below is a detailed installation guide.

-----------------------
Step 1: Install YALMIP
-----------------------

QUINOPT is an add-on for YALMIP, the optimization modeling software by J. Löfberg. If you already have YALMIP installed, you cak skip this step. Otherwise, `download YALMIP's latest version <https://yalmip.github.io/download/>`_ and install it by adding the following folders to MATLAB's path:

.. code:: matlab

   YALMIP-master
   YALMIP-master/extras
   YALMIP-master/solvers
   YALMIP-master/modules
   YALMIP-master/modules/parametric
   YALMIP-master/modules/moment
   YALMIP-master/modules/global
   YALMIP-master/modules/sos
   YALMIP-master/operators

You can test your YALMIP installation by running

.. code:: matlab

   >> yalmiptest
   
at the MATLAB command prompt. More details on how to install or update YALMIP be found `on YALMIP's website <https://yalmip.github.io/tutorial/installation/>`_.

     
.. note::

   The folder names above are the default when YALMIP is downloaded from GitHub. Should you wish to use a different folder name, simply replace ``<YALMIP-master>`` with the appropriate path.
   
  
------------------------------
Step 2: Install an SDP solver
------------------------------

To be able to use QUINOPT, you need to install a semidefinite programming (SDP) solver compatible with YALMIP. A complete list of YALMIP-compatible SDP solvers can be found `here <https://yalmip.github.io/allsolvers/>`_.

.. warning::

   QUINOPT has been tested with `SeDuMI`_, `SDPT3`_, `SDPA`_, and `Mosek`_ 
   (free for users in academia). Any other YALMIP-compatible SDP solver should 
   work, but use at your own risk!
   
   .. _SeDuMi: https://github.com/sqlp/sedumi
   .. _SDPT3: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
   .. _SDPA: http://sdpa.sourceforge.net/
   .. _Mosek: https://www.mosek.com/
     
	 
-----------------------
Step 3: Install QUINOPT
-----------------------
	 
If you have successfully installed YALMIP and a compatible SDP solver, you are ready to install QUINOPT. First, download QUINOPT's latest version (you can `find in on GitHub <https://github.com/aeroimperial-optimization/QUINOPT/releases>`_, or simply `click here to download <https://github.com/aeroimperial-optimization/QUINOPT/archive/v1.4.zip>`_). 
After unzipping the downloaded folder, navigate to it in MATLAB and simply run the installer:

.. code:: matlab

   >> installQUINOPT

The installer should add the required folders to the MATLAB path and run 
some test problems to make sure everything is working. If you experience any installation problems, please `contact us`_ or file an issue issues via the `Github issue tracker <https://github.com/aeroimperial-optimization/QUINOPT/issues>`_.

.. _contact us: mailto:giovanni.fantuzzi10@imperial.ac.uk?Subject=QUINOPT%20installation%20issue



`Back to Table of Contents <http://quinopt.readthedocs.io/>`_