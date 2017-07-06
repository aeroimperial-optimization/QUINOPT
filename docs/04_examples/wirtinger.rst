Wirtinger's inequality
======================================================

Wirtinger's inequality states that for a :math:`2\pi`-periodic function :math:`v` that satisfies

.. math::

	\int_0^{2\pi} v(x) {\rm d}x = 0

there exists a constant :math:`C` such that

.. math::

	\int_0^{2\pi} \left[
    C \vert v'(x) \vert^2 - \vert v(x) \vert^2
    \right]{\rm d}x \geq 0.

In this example, we verify that the smallest possible constant is :math:`C = 1`. We do so by computing a sequence of convergent upper and lower bounds on the smallest :math:`C` with QUINOPT. The objective of this example is to demonstrate how to override the default options in QUINOPT.

:download:`Download the MATLAB file for this example <./downloads/example04.m>`

--------------------------------
1. Reformulate the problem
--------------------------------
Before inputting the integral inequality into QUINOPT, we need to remove the integral constraint on :math:`v(x)`. This can be done by defining

.. math::

    u(x) = \int_0^x v(t) \, {\rm d}t

so the zero-integral and periodicity conditions on :math:`v` become

.. math::

    u'(0)-u'(2\pi)=0,\quad u(2\pi) = 0.

Moreover, by definition of :math:`u(x)` we have the additional boundary condition :math:`u(0)=0`. With this change of variables, Wirtinger's inequality becomes

.. math::

	\int_0^{2\pi} \left[
    C \vert u''(x) \vert^2 - \vert u'(x) \vert^2
    \right]{\rm d}x \geq 0.

--------------------------
2. Set up the problem
--------------------------

As usual, we start by clearing the model and creating the variables:

.. code-block:: matlabsession

    >> clear;
    >> yalmip clear;
    >> quinopt clear;
    >> x = indvar(0,2*pi);
    >> u = depvar(x);
    >> parameters C;

To set up the integral inequality constraint, we specify the integrand and the vector of boundary conditions:

.. code-block:: matlabsession

	>> expr = C*u(x,2)^2 - u(x,1)^2;
	>> bc = [u(0); u(2*pi); u(0,1)-u(2*pi,1)];


--------------------------
3. Solve the problem
--------------------------

We seek upper and lower bounds on the lowest possible :math:`C`. QUINOPT's default behaviour is to compute upper bounds on the optimal objective value by formulating and inner approximation of the feasible set of the integral inequality constraints. A lower bound is found by overriding the default method with an ``options`` structure:

.. code-block:: matlabsession

    >> options.method = 'outer';

This makes QUINOPT formulate an outer approximation of the feasible set of the integral inequality constraints. We can also tell QUINOPT to run quietly by setting YALMIP's "verbose" option to 0, and cache YALMIP's available solvers to improve YALMIP's performance:

.. code-block:: matlabsession

    >> options.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);

To compute a lower bound on the optimal :math:`C` we can then simply run

.. code-block:: matlabsession

    >> quinopt(expr,bc,C,options);
    >> LB = value(C);                      % extract the lower bound on the optimal C

To compute an upper bound, we need to reset QUINOPT's default behaviour:

.. code-block:: matlabsession

    >> options.method = 'inner';           % reset the default behaviour: inner approximation
    >> quinopt(expr,bc,C,options);
    >> UB = value(C);                      % extract the upper bound on the optimal C

.. note::

    The commands above return an upper bound``UB = 1.000618`` , but a lower bound ``LB = NaN``. This is because QUINOPT's default outer approximation is always feasible, and so the optimization problem that is solved has an unbounded objective value. This issue is resolved in the next section.

--------------------------
4. Improve the results
--------------------------

As we have seen, the lower bound obtained with QUINOPT's default outer approximation is not good. This issue can be resolved by refining the approximation that QUINOPT builds. Roughly speaking, QUINOPT builds such approximations by expanding the dependent variables as polynomials of degree :math:`N`. By default, QUINOPT determines :math:`N` based on the problem (`see our paper for details <https://arxiv.org/pdf/1607.04210.pdf>`_): for Wirtinger's inequality, the default value is :math:`N=2`. Fortunately, we can tell QUINOPT to use a larger value by specifying the option ``options.N``:

.. code-block:: matlabsession

    >> options.N = 3;                     % use polynomial expansions of degree N=3
    >> options.method = 'outer';          % use outer approximation
    >> quinopt(expr,bc,C,options);
    >> LB = value(C);                     % extract the improved lower bound on the optimal C
    >> options.method = 'inner';          % use inner approximation
    >> quinopt(expr,bc,C,options);
    >> UB = value(C);                     % extract the improved upper bound on the optimal C

The lower bound obtained with these settings is ``LB = 0.657974``, and the upper bound improves to ``UB=1.000034``. Increasing :math:`N` further, we see that the two values converge to 1:

=========== ============ ============ =============
:math:`N`   Lower bound  Upper bound  Difference
=========== ============ ============ =============
2 (default) NaN          1.000618     NaN
3           0.657974     1.000034     3.42e-01
4           0.939960     1.000001     6.00e-02
5           0.992796     1.000000     7.20e-03
6           0.999413     1.000000     5.87e-04
7           0.999966     1.000000     3.35e-05
8           0.999999     1.000000     1.17e-06
9           1.000000     1.000000     1.97e-08
=========== ============ ============ =============

.. note::

    If ``options.N`` is lower than the minimum value (again, `see our paper for details <https://arxiv.org/pdf/1607.04210.pdf>`_), QUINOPT issues a warning and uses the minimum value of :math:`N` instead.


-----------------------
5. Summary
-----------------------

In summary, the optimal constant for Wirtinger's inequality can be determined with the following simple lines of code:


.. code-block:: matlabsession

    >> % Clean up
    >> clear;
    >> yalmip clear;
    >> quinopt clear;
    >> % Set up the problem
    >> x = indvar(0,2*pi);
    >> u = depvar(x);
    >> parameters C;
    >> expr = C*u(x,2)^2 - u(x,1)^2;
    >> bc = [u(0); u(2*pi); u(0,1)-u(2*pi,1)];
    >> % Set options for YALMIP
    >> options.YALMIP = sdpsettings('verbose',0,'cachesolvers',1);
    >> % Compute a lower bound on the optimal C
    >> options.method = 'outer';
    >> quinopt(expr,bc,C,options);
    >> LB = value(C);
    >> % Compute an upper bound on the optimal C
    >> options.method = 'inner';
    >> quinopt(expr,bc,C,options);
    >> UB = value(C);
    >> % Improve the solution by setting options.N
    >> options.N = 9;
    >> options.method = 'outer';
    >> quinopt(expr,bc,C,options);
    >> LB = value(C);
    >> options.method = 'inner';
    >> quinopt(expr,bc,C,options);
    >> UB = value(C);


----------------------

* :doc:`Back to Table of Contents <../index>`
