.. include:: ../substitutions.txt

``@legpoly/legpolyval()``
================================================

 Evaluate polynomial in Legendre basis (class ``legpoly``).

---------------------

:Syntax: ``F = legpolyval(p,x)``
:Description: evaluates the polynomial ``p`` at the points specified by the vector (or matrix) ``x``. The points in ``x`` must be in the interval :math:`[a,b]` where ``p`` is defined (this can be recovered by calling ``getDomain(p)``).

------------

* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
