.. include:: ../substitutions.txt

``@legpoly/int()``
================================================

Integrate a polynomial in the Legendre basis (class ``legpoly``).

---------------------

:Syntax: ``P = INT(p)`` or ``P = INT(p,x)``
:Description: integrates the polynomial ``p`` in Legendre basis  (class ``legpoly``" with respect to its independent variable ``x``. Indefinite integration is performed such that :math:`P(0)=0`. If a different behaviour is required, please use the function ``legpolyint()``.

------------

:Syntax: ``P = INT(p,x,a,b)``
:Description: computes the integral of ``p`` from ``a`` to ``b``. The integration limits ``a`` and ``b`` must be contained within the domain of definition of the polynomial ``p``, as returned by calling ``getDomain(p)``.

------------

* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
