.. include:: ../substitutions.txt

Feasibility of an integral inequality with boundary conditions
==============================================================

We now revisit the :doc:`previous example <./feas_nobc>` to show how to specify boundary conditions on the dependent variables. Specifically, we demonstrate how to use QUINOPT to find the minimum value :math:`\gamma` such that

.. math::

	\int_0^1 \left[
	\vert u'(x) \vert^2 +\gamma\,u'(x)\,u(x) + \vert u(x) \vert^2
	\right] {\rm d}x
	\geq 0

for all differentiable functions :math:`u(x)` that satisfy the boundary conditions

.. math::

	u'(0) = 0, \quad \text{and} \quad u(1)=0.

:download:`Download the MATLAB file for this example <./example02.m>`

---------------------------
1. Clearing the workspace
---------------------------

Before we solve a new example, it is good practice to clear the variables from the workspace using MATLAB's ``clear`` command. Moreover, to prevent the build-up of unused internal variables in QUINOPT, it is useful to clear QUINOPT's internal variables using the command ``quinopt clear``:

.. code-block:: matlabsession

	>> clear
	>> quinopt clear


---------------------------
2. Creating the variables
---------------------------

As in the :doc:`previous example <./feas_nobc>`, we create the independent variable :math:`x\in[0,1]`, the dependent variable :math:`u(x)`, and the optimization parameter :math:`\gamma`:

.. code-block:: matlabsession

	>> x = indvar(0,1);
	>> u = depvar(x);
	>> parameters gamma;


------------------------------
3. Setting up the inequality
------------------------------

As in the :doc:`previous example <./feas_nobc>`, we begin by constructing the integrand expression.

.. code-block:: matlabsession

	>> EXPR = u(x,1)^2 + gamma*u(x,1)*u(x) + u(x)^2;	% Create the integrand

The boundary condition can be specified through a vector ``BC``, which is interpreted internally as the element-wise condition ``BC=0``:

.. code-block:: matlabsession

	>> BC(1) = [u(0,1)];		% Create the boundary condition u'(0)=0
	>> BC(2) = [u(1)]; 		% Create the boundary condition u(0)-u(1)=0


.. note::

		The syntax ``u(POINT,DER)`` is used to specify the derivative of :math:`u(x)` of order ``DER``, evaluated at the point ``POINT``. Possible values of ``POINT`` are the independent variable of integration ``x``, or the extrema of the domain of integration, in this case 0 and 1.



-----------------------------------
4. Solving the problem with QUINOPT
-----------------------------------

Once the variables and the integrand of the inequality have been set up, a value of :math:`\gamma` for which the integral functional is positive semidefinite can be found using the command ``quinopt()`` with two inputs:

.. code-block:: matlabsession

	>> quinopt(EXPR,BC);	% Solve the problem
	>> value(gamma) 	% Extract the value of gamma


-----------------------------------
5. Summary
-----------------------------------
In summary, a feasible value :math:`\gamma` such that the integral inequality at the top of the page holds can be found using the following simple commands:

.. code-block:: matlabsession

	>> clear						% Clear workspace
	>> quinopt clear					% Clear QUINOPT internals
	>> x = indvar(0,1); 					% Create the independent variable with domain [0,1]
	>> u = depvar(x);					% Create the dependent variable u(x)
	>> parameters gamma;					% Create the optimization variable gamma
	>> EXPR = u(x,1)^2 + gamma*u(x,1)*u(x) + u(x)^2;       % Create the integrand
	>> BC(1) = [u(0,1)];					% Create the boundary condition u'(0)=0
	>> BC(2) = [u(1)]; 					% Create the boundary condition u(1)=0
	>> quinopt(EXPR,BC);					% Solve the problem
	>> value(gamma) 					% Extract the value of gamma


`Back to Table of Contents <http://quinopt.readthedocs.io/>`_
