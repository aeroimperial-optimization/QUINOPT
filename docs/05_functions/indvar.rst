.. include:: ../substitutions.txt

``indvar()``
================================================

Create an independent variable of integration to define integral inequalities in QUINOPT.

---------------------

:Syntax: ``x = indvar(a,b)``
:Description: creates an independent variable of integration with
        domain :math:`[a,b]` used to set up an integral inequality constraint with the toolbox QUINOPT.

.. warning::
    An integral inequality constraint defined with QUINOPT can only have one independent variable - multivariable integrals are not allowed. However, one independent variable can be used to define multiple integral inequalities, and different inequalities can have different independent variables (although this is not necessary).

------------

* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
