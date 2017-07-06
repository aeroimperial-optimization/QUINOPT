.. include:: ../substitutions.txt

Feasibility of an integral inequality
================================================

One of the most basic applications of QUINOPT is to find one value of an optimization parameter such that a homogeneous quadratic integral functional is positive. Here, we demonstrate how to use QUINOPT to find a value :math:`\gamma` such that

.. math::

	\int_0^1 \left[
	\vert u'(x) \vert^2 +\gamma\,u'(x)\,u(x) + \vert u(x) \vert^2
	\right] {\rm d}x
	\geq 0

for all differentiable functions :math:`u(x)`. Clearly, possible choices are :math:`\gamma=0`, and :math:`\gamma=-2` or :math:`\gamma=2` (when the integrand is a perfect square).

:download:`Download the MATLAB file for this example <./downloads/example01.m>`

--------------------------
1. Create the variables
--------------------------

The first step to use QUINOPT is to set up the problem variables. These are the integration variable :math:`x\in[0,1]` (the *independent variable*), the unknown function :math:`u(x)` (the *dependent variable*), and the optimization parameter :math:`\gamma`.

First, we create the independent variable :math:`x\in[0,1]` using the command ``indvar()``, as

.. code-block:: matlabsession

	>> x = indvar(0,1);			% Create the independent variable with domain [0,1]

Then, we set up the dependent variable :math:`u(x)` using the command ``depvar()``:

.. code-block:: matlabsession

	>> u = depvar(x);			% Create the dependent variable u(x)


Finally, we set up the optimization parameter :math:`\gamma` using the command ``parameters``

.. code-block:: matlabsession

	>> parameters gamma;			% Create the optimization variable gamma



.. note::

		The commands ``indvar()`` and ``depvar()`` return MATLAB objects of class ``@indvar`` and ``@depvar``, respectively. While the ``@indvar`` class behaves like a usual YALMIP variable, the ``@depvar`` class is specific to QUINOPT and does **not** behave like a YALMIP variable. Instead, it is intended to be used only as shown in the following.


-----------------------------
2. Set up the inequality
-----------------------------

Once the variables have been set up, we can set up the inequality. This is done in QUINOPT by constructing the integrand expression.

.. code-block:: matlabsession

	>> EXPR = u(x,1)^2 + gamma*u(x,1)*u(x) + u(x)^2;	% Create the integrand

In the expression above, the syntax ``u(x,DER)`` is used to specify the derivative of :math:`u(x)` of order ``DER``. In other words, ``u(x,1)`` is the first derivative of :math:`u(x)`.


.. note::

		The integration interval has already been specified when defining the independent variable.

-----------------------------------
3. Solve the problem with QUINOPT
-----------------------------------

Once the variables and the integrand of the inequality have been set up, a value of :math:`\gamma` for which the integral functional is positive semidefinite can be found using the command ``quinopt()``, together with YALMIP's command ``value()``

.. code-block:: matlabsession

	>> quinopt(EXPR);	% Solve the problem
	>> value(gamma) 	% Extract the value of gamma


-----------------------------------
4. Summary
-----------------------------------
In summary, a feasible value :math:`\gamma` such that the integral inequality at the top of the page holds can be found using the following simple commands:

.. code-block:: matlabsession

	>> x = indvar(0,1); 					% Create the independent variable with domain [0,1]
	>> u = depvar(x);					% Create the dependent variable u(x)
	>> parameters gamma;					% Create the optimization variable gamma
	>> EXPR = u(x,1)^2 + gamma*u(x,1)*u(x) + u(x)^2;       % Create the integrand
	>> quinopt(EXPR);					% Solve the problem
	>> value(gamma) 					% Extract the value of gamma



----------------------

* :doc:`Back to Table of Contents <../index>`
