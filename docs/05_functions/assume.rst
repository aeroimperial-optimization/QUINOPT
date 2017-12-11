.. include:: ../substitutions.txt

``@depvar/assume()``
================================================

Add assumption on dependent variables in QUINOPT.

---------------------

:Syntax: ``assume(U,STR)``
:Description: adds the assumption specified by the character string ``STR`` on the dependent variable ``U`` (class ``depvar``). Currently, allowed values for ``STR`` are:

    * ``'even'``: assume that ``U`` is symmetric with respect to the midpoint of the domain :math:`[a,b]`  in which the dependent variable ``U`` is defined.

    * ``'odd'``: assume that ``U`` is anty-symmetric with respect to the midpoint of the domain :math:`[a,b]`  in which the dependent variable ``U`` is defined.

    * ``'none'``: remove all previous assumptions on the dependent variable ``U``.

---------------------

    :Syntax: ``assume(U1,STR1,U2,STR2,...)``
    :Description: adds the assumptions specified by the character strings ``STR1``, ``STR2``, ..., on the dependent variables ``U1``, ``U2``, ..., as if set by the commands ``assume(U1,STR1)``, ``assume(U2,STR2)``, and so on.

------------

* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
