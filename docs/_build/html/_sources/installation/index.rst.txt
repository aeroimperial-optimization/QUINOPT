Install QUINOPT
===============

To install QUINOPT:

1. `Download <https://yalmip.github.io/download/>`_ and `install <https://yalmip.github.io/tutorial/installation/>`_ YALMIP 
2. Install a semidefinite programming (SDP) solver compatible with YALMIP. Click `here <https://yalmip.github.io/allsolvers/>`_ for a complete list of YALMIP-compatible SDP solvers.

.. note::

   QUINOPT has been tested with `SeDuMI`_, `SDPT3`_, `SDPA`_, and `Mosek`_ 
   (free for users in academia). Any other YALMIP-compatible SDP solver should 
   work, but use at your own risk!
   
   .. _SeDuMi: https://github.com/sqlp/sedumi
   .. _SDPT3: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
   .. _SDPA: http://sdpa.sourceforge.net/
   .. _Mosek: https://www.mosek.com/
     
3. Install QUINOPT by running the MATLAB installer:

.. code:: matlab

   >> installQUINOPT

The installer should add the required folders to the MATLAB path and run 
some test problems to make sure everything is working. If you experience any installation problems, please `contact us`_.

.. _contact us: mailto:giovanni.fantuzzi10@imperial.ac.uk?Subject=QUINOPT%20installation%20issue

