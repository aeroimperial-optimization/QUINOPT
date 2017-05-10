Poincaré's inequality for odd and periodic functions
=========================================================

Poincaré's inequality for *odd* functions :math:`u:[-1,1]\to\mathbb{R}` that satify the periodicity condition

.. math::

	u(-1)=u(1)

states that

.. math::

	\int_{-1}^1 \vert u'(x) \vert^2 {\rm d}x \geq \pi^2 \int_0^1 \vert u(x) \vert^2 {\rm d}x.

In this example, we verify that the constant :math:`\pi^2` on the right-hand side is optimal, in the sense that it is the optimal value for the optimization problem

.. math::

	\begin{aligned}
	\max_{\nu} \quad &\nu\\
	\text{subject to} \quad
	&\int_{-1}^1 \left[
	\vert u'(x) \vert^2 -\nu \vert u(x) \vert^2
	\right] {\rm d}x \geq 0,
	\quad u\in C^1([0,1],\mathbb{R}),\,
    \begin{cases}
    u(-1)=u(1), \\ u(-x)=-u(x).
    \end{cases}
	\end{aligned}

The aim of this example is to show how QUINOPT can handle symmetry constraints on the dependent variables.

:download:`Download the MATLAB file for this example <./downloads/example05.m>`

--------------------------
1. Create the variables
--------------------------

As usual, we start by clearing the workspace, YALMIP's and QUINOPT's internal variables, and creating the variables for the problem:

.. code-block:: matlabsession

    >> clear
    >> yalmip clear
    >> quinopt clear
    >> x = indvar(-1,1);     % the integration variable with domain [-1,1]
    >> u = depvar(x);        % the dependent variable u(x)
    >> parameters nu;        % the optimization variable nu


------------------------------
2. Set up the inequality
------------------------------

To set up Poincaré's inequality constraint, first we specify the integrand:

.. code-block:: matlabsession

    >> EXPR = u(x,1)^2 - nu*u(x)^2;

Then, we set the boundary and symmetry conditions on :math:`u(x)`. The periodic boundary conditions is enforced as :math:`u(-1)-u(1)=0`, while the symmetry condition can be enforced using the command ``assume()``:

.. code-block:: matlabsession

    >> BC = [ u(-1)-u(1) ];
    >> assume(u,'odd')

.. note::

	Other valid assumptions are ``assume(u,'even')`` to assume that :math:`u(x)` is even, and ``assume(u,'none')`` to remove any previous assumption. Moreover, when the domain of the independent variable used to construct the dependent variable ``u`` is a generic interval :math:`[a,b]` rather than a symmetric interval :math:`[-a,a]`, the symmetry condition set with the command ``assume()`` is relative to the midpoint of the domain, :math:`(a+b)/2`.


--------------------------
3. Solve the problem
--------------------------

To solve the problem and maximize :math:`\nu`, we use the command ``quinopt()`` with three arguments: ``EXPR``, ``BC`` and the objective function. Since QUINOPT minimizes the specified objective function, instead of maximizing :math:`\nu` we minimize :math:`-\nu`.

.. code-block:: matlabsession

    >> quinopt(EXPR,BC,-nu);       % Maximize nu (by minimizing -nu)
    >> value(nu)/pi^2              % Get the optimal value (in units of pi^2)

With the default parameters in QUINOPT, we obtain :math:`\nu_{\rm opt} = 0.8184 \,\pi^2`, i.e. the optimal solution returned by QUINOPT is within 81.8% of true optimum :math:`\nu_{\rm exact}=\pi^2`. In fact, we can refine the approximation of the integral inequalities by increasing the number of Legendre coefficients used by QUINOPT to expand them. We do this by setting the option ``options.N``, as in :doc:`the previous example <./wirtinger>`:

.. code-block:: matlabsession

    >> options.N = 5;              % Use N=5 expansion coefficients
    >> quinopt(EXPR,BC,-nu);       % Maximize nu (by minimizing -nu)
    >> value(nu)/pi^2              % Get the optimal value (in units of pi^2)

The optimal value of :math:`\nu` returned by QUINOPT in this case is :math:`\nu_{\rm opt} = 0.999965 \,\pi^2`, meaning that the numerical optimum is essentially indistinguishable from the true optimum :math:`\nu_{\rm exact}=\pi^2`. Setting ``options.N`` to larger values further improves the numerical optimum (note that roundoff errors might result in a numerical optimum that is slightly larger than the exact solution).


-----------------------
4. Summary
-----------------------

In summary, the optimal constant for Poincaré's inequality for odd, periodic functions can be determined with the following simple lines of code:


.. code-block:: matlabsession

    >> % Set up the variables
    >> clear
    >> quinopt clear
    >> x = indvar(-1,1);
    >> u = depvar(x);
    >> parameters nu;
    >> % Build the inequality
    >> EXPR = u(x,1)^2 - nu*u(x)^2;
    >> BC = [ u(-1)-u(1) ];
    >> assume(u,'odd')
    >> % Solve with the default parameters
    >> quinopt(EXPR,BC,-nu);
    >> value(nu)/pi^2
    >> % Refine the solution: solve with N=5 expansion coefficients
    >> options.N = 5;
    >> quinopt(EXPR,BC,-nu,options);
    >> value(nu)/pi^2


----------------------

* :doc:`Back to Table of Contents <../index>`
