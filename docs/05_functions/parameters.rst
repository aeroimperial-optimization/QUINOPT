.. include:: ../substitutions.txt

``parameters()``
================================================

Create symbolic optimization variables for QUINOPT.

---------------------

:Syntax: ``P = parameters(m,n)``
:Description: reates an :math:`m \times n` matrix of parameters ``P`` to be used
        as optimization variables with QUINOPT. The parameters are
        YALMIP variables (class ``sdpvar``).

---------------------

:Syntax: ``parameters p1 p2 p3 p4``
:Description: creates multiple scalar optimization parameters, ``p1``, ..., ``p4``.

------------

* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
