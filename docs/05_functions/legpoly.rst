.. include:: ../substitutions.txt

``legpoly()``
================================================

Create polynomial in the Legendre basis to use with QUINOPT.

---------------------

:Syntax: ``P = legpoly(x,DEG)``
:Description: creates a polynomial ``P`` in the independent variable ``x`` of degree ``DEG``, expressed in Legendre basis. That is, the polynomial :math:`P(x)` is expressed as

        .. math::

            P(x) = C_1\,\mathcal{L}_0[z(x)] + C_2\,\mathcal{L}_1[z(x)] + \cdots + C_{{\rm DEG}+1}\,\mathcal{L}_{\rm DEG}[z(x)]

        where :math:`\mathcal{L}_n(z)` is the Legendre polynomial of degree :math:`n`. Since Legendre polynomials are defined over the standard domain :math:`[-1,1]`, the original independent variable :math:`x` with domain :math:`[a,b]` is rescaled to

            .. math::

                z(x) = \frac{2\,x-b-a}{b-a}

        The input ``x`` must be a valid independent variable (class ``indvar``), and ``DEG`` should be a non-negative integer. The coefficients of the polynomial are YALMIP variables (class ``sdpvar``) and can be recovered with the command ``C = coefficients(P)``. Finally, ``P`` can be displayed symbolically in the standard monomial basis using the command ``sdisplay(P)``.

------------

:Syntax: ``[P,C] = legpoly(x,DEG)``
:Description: also returns the Legendre coefficients of the polynomials in the vector ``C``. These are YALMIP variables (class ``sdpvar``). The coefficients in ``C`` are listed in order of increasing degree of the corresponding Legendre polynomial (see above).

---------------

:Syntax: ``P = legpoly(x,DEG,COEF)``
:Description: creates a polynomial ``P`` expressed in Legendre basis whose coefficients are specified by ``COEF``. ``COEF`` can be a numeric/sdpvar vector, or an :math:`M\times N` cell array whose entries are numeric/sdpvar vectors. When ``COEF`` is a cell array, an :math:`M\times N` matrix of Legendre polynomials is created such that the coefficients of the entry ``P(i,j)`` are given by ``COEF{i,j}``.

---------------

:Syntax: ``[P,C] = legpoly(x,DEG,M,N)``
:Description: creates an :math:`M\times N` matrix of
        Legendre polynomials of degree ``DEG``. The output ``C``, containing the coefficients of each entry of ``P``, is optional.

--------------

* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
