.. include:: ../substitutions.txt

``depvar()``
================================================

Define dependent variables to define integral inequalities in QUINOPT.

--------------

:Syntax: ``U = depvar(x)``
:Description: sets up a symbolic variable modelling a generic function :math:`U(x)`, where ``x`` is a valid independent variable with domain :math:`[a,b]` created with the command indvar. ``U`` behaves like a function handle, and is used with the syntax

    .. code-block:: matlab

            U(POINT)

    where ``POINT`` is either the independent variable ``x``, the lower extremum :math:`a` of the domain of ``U``, or the upper extremum :math:`b` of the domain of ``U``. Moreover, derivatives of ``U`` can be created/accessed using the syntax

    .. code-block:: matlab

        U(POINT,DERIVATIVE)

    where ``POINT`` is ``x``, :math:`a` or :math:`b` and ``DERIVATIVE`` is the desired derivative order. Note that ``U(POINT,0)`` is equivalent to ``U(POINT)``.

------------

:Syntax: ``[U1,U2,...Uq] = depvar(x)``
:Description: sets up multiple dependent variables, ``U1``, ..., ``Uq``. Each
        dependent variable depends on the independent variable ``x``.

------------
        
* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
