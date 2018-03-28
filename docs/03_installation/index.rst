.. include:: ../substitutions.txt


Install QUINOPT
===============

QUINOPT is easily installed by running the installer ``installQUINOPT.m`` in MATLAB. Below is a detailed installation guide.

-----------------------
Step 1: Install YALMIP
-----------------------

QUINOPT is an add-on for YALMIP, the optimization modeling software by J. Löfberg. If you already have YALMIP installed, you cak skip this step. Otherwise, `download YALMIP <https://yalmip.github.io/download/>`_ and install it by adding the following folders to MATLAB's path:

.. code::

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

.. code-block:: matlabsession

   >> yalmiptest

at the MATLAB command prompt. More details on how to install or update YALMIP be found `on YALMIP's website <https://yalmip.github.io/tutorial/installation/>`_.


.. note::

   The folder names above are the default when YALMIP is downloaded from GitHub. Should you wish to use a different folder name, simply replace ``<YALMIP-master>`` with the appropriate path.

.. warning::

   YALMIP is regularly updated, and changes in YALMIP may sometimes affect the functionality of QUINOPT. If you have downloaded YALMIP's latest version and are experiencing installation problems, please `contact us`_ or file an issue via the `GitHub issue tracker <https://github.com/aeroimperial-optimization/QUINOPT/issues>`_. 

------------------------------
Step 2: Install an SDP solver
------------------------------

To be able to use QUINOPT, you need to install a semidefinite programming (SDP) solver compatible with YALMIP. A complete list of YALMIP-compatible SDP solvers can be found `here <https://yalmip.github.io/allsolvers/>`_.

.. warning::

   QUINOPT has been tested with `SeDuMi <https://github.com/sqlp/sedumi>`_, `SDPT3 <http://www.math.nus.edu.sg/~mattohkc/sdpt3.html>`_, `SDPA <http://sdpa.sourceforge.net/>`_, and `Mosek <https://www.mosek.com/>`_ (free for users in academia). Other suitable YALMIP-compatible SDP solver should work, but use them at your own risk!


-----------------------
Step 3: Install QUINOPT
-----------------------

If you have successfully installed YALMIP and a compatible SDP solver, you are ready to install QUINOPT. First, `download QUINOPT's latest stable version`_
(for the "developer" version or previous versions, visit the `Download`_ page).

After unzipping the downloaded folder, navigate to it in MATLAB and simply run the installer:



.. code-block:: matlabsession

   >> installQUINOPT

The installer should compile the required files, add the required folders to the MATLAB path, and run some test problems to make sure everything is working. If you experience any installation problems, please `contact us`_ or file an issue via the `GitHub issue tracker <https://github.com/aeroimperial-optimization/QUINOPT/issues>`_.

.. warning::

    During installation, you may receive the following warning:

    .. code-block:: matlabsession

        Warning: Compilation of mex files by installQUINOPT failed.
        QUINOPT will still work without compiled mex files, but
        it will be slower. To resolve the issue, make sure that
        a supported compiler is installed and re-run the installer.

    QUINOPT should still work, but you may wish to resolve the issue with the mex file compilation. You can find a list of supported compilers for MATLAB's latest version `on this webpage <https://uk.mathworks.com/support/compilers.html>`_; for all other versions of MATLAB please `look at this webpage <https://uk.mathworks.com/support/sysreq/previous_releases.html>`_.

.. Define links (ignore syntax highlighting in Atom?):

.. _contact us: mailto:giovanni.fantuzzi10@imperial.ac.uk?Subject=QUINOPT%20installation%20issue

.. _download QUINOPT's latest stable version: https://github.com/aeroimperial-optimization/QUINOPT/archive/master.zip

.. _Download: ../01_download/index.html

----------------------

* :doc:`Back to Table of Contents <../index>`
