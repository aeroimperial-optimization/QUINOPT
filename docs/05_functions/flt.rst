.. include:: ../substitutions.txt

``flt()``
================================================

Fast Legendre transform according to the algorithm presented in `Iserles, Numer. Math. 117, 529-553 (2010) <http://link.springer.com/article/10.1007%2Fs00211-010-0352-1>`_.

---------------------

:Syntax: ``C = flt(FUN,N,DOMAIN)``
:Description: computes the first ``N`` Legendre coefficients of the expansion of the function ``FUN``, defined over the domain ``DOMAIN``, and returns them in the vector ``C``. The input ``FUN`` must be a handle to the function to be projected onto the first ``N`` Legendre polynomials, ``N`` must be a non-negative integer, and ``DOMAIN`` must be a  vector with 2 elements, i.e. ``DOMAIN = [a,b]`` with ``a<b``.

------------

* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
