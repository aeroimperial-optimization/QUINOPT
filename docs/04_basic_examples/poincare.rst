Poincaré's inequality with Dirichlet boundary conditions
=========================================================

Poincaré's inequality for functions :math:`u:[0,1]\to\mathbb{R}` that satify the Dirichlet boundary conditions

.. math::

	u(0)=0=u(1)

states that

.. math::

	\int_0^1 \vert u'(x) \vert^2 {\rm d}x \geq \pi^2 \int_0^1 \vert u(x) \vert^2 {\rm d}x.

In this example, we verify that the constant :math:`\pi^2` on the right-hand side is optimal, in the sense that it solve the optimization problem

.. math::

	\begin{aligned}
	\max_{\gamma} \quad &\gamma\\
	\text{subject to} \quad
	&\int_0^1 \left[
	\vert u'(x) \vert^2 -\gamma \vert u(x) \vert^2
	\right] {\rm d}x \geq 0,
	\quad u\in C^1([0,1],\mathbb{R}),\quad u(0)=0=u(1).
	\end{aligned}

:download:`Download the MATLAB file for this example <./example03.m>`

--------------------------
1. Creating the variables
--------------------------

As usual, we start by clearing the model and creating the variables:

.. code-block:: matlabsession

	>> clear
	>> quinopt clear
	>> x = indvar(0,1);
	>> u = depvar(x);
	>> parameters gamma;


------------------------------
2. Setting up the inequality
------------------------------

To set up Poincaré's inequality constraint, first we specify the integrand:

.. code-block:: matlabsession

	>> EXPR = u(x,1)^2 - gamma*u(x)^2;

and then we set up the vector of boundary conditions (this can be a row vector as in the previous example, or a column vector):

.. code-block:: matlabsession

	>> BC = [u(0); u(1)];


--------------------------
3. Solving the problem
--------------------------

To solve the problem and maximize :math:`\gamma`, we use once again the command ``quinopt()``, this time with three arguments: ``EXPR``, ``BC`` and the objective function.

.. note::
	When calling ``quinopt(EXPR,BC,OBJECTIVE)``, QUINOPT **minimizes** the specified objective function.

Since QUINOPT minimizes the specified objective function, instead of maximizing :math:`\gamma` we minimize :math:`-\gamma`:

.. code-block:: matlabsession

	>> quinopt(EXPR,BC,-gamma);		% Maximize gamma (by minimizing -gamma)
	>> value(gamma)/pi^2			% Get the optimal value (in units of pi^2)

With the default parameters in QUINOPT, we obtain :math:`\gamma_{\rm opt} = 0.9994 \,\pi^2`, i.e. the optimal solution returned by QUINOPT is within 99.9% of true optimum :math:`\gamma_{\rm opt}=\pi^2`.


-----------------------
4. Summary
-----------------------

In summary, the optimal constant for Poincaré's inequality can be determined with the following simple lines of code:


.. code-block:: matlabsession

	>> clear
	>> quinopt clear
	>> x = indvar(0,1);
	>> u = depvar(x);
	>> parameters gamma;
	>> EXPR = u(x,1)^2 - gamma*u(x)^2;
	>> BC = [u(0); u(1)];
	>> quinopt(EXPR,BC,-gamma);
	>> value(gamma)/pi^2


`Back to Table of Contents <http://quinopt.readthedocs.io/>`_
