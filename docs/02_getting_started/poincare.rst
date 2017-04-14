Poincaré inequality with Dirichlet boundary conditions
======================================================

The Poincaré inequality for functions :math:`u:[0,1]\to\mathbb{R}` that satify the Dirichlet boundary conditions

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


-----------------------
Creating the variables
-----------------------


--------------------------
Setting up the inequality
--------------------------


-----------------------
Defining the objective
-----------------------


-----------------------
Solving the problem
-----------------------




`Back to Table of Contents <http://quinopt.readthedocs.io/>`_