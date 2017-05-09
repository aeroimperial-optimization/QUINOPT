.. include:: ../substitutions.txt

What is QUINOPT?
================

QUINOPT (QUadratic INtegral OPTimisation) is an open-source add-on for `YALMIP <https://yalmip.github.io/>`_ to compute rigorous upper and lower bounds on the optimal value of optimization problems with infinite-dimensional polynomial quadratic integral inequality constraints. Such problems commonly arise from stability analysis of linear PDEs using the so-called *energy method* (also known as :math:`\mathcal{L}^2` stability), and in bounding time-average properties of turbulent fluid flows using the so-called *background method*.

In the simplest form, given a bounded interval :math:`[a,b] \subset \mathbb{R}`, a :math:`k`-times continuously differentiable function :math:`u \in C^k([a,b],\mathbb{R})`, and a vector of optimization variables :math:`\gamma \in \mathbb{R}^s`, QUINOPT computes upper and/or lower bounds on the optimal value of the optimization problem

.. math::

	\begin{aligned}
	\min_{\gamma} \quad &c^T \gamma\\
	\text{subject to} \quad &\int_a^b Q_{\gamma}(x,u(x),u'(x),...,u^k(x))
	\,{\rm d}x \geq 0 \quad \forall u(x) \in \mathcal{H}
	\end{aligned}

by constructing SDP-representable inner and outer approximations of its feasible set.
In the problem above, :math:`Q_{\gamma}(x,u(x),u'(x),...,u^k(x))` is

* a quadratic polynomial in :math:`u(x),u'(x),...,u^k(x)`;
* a polynomial in in :math:`x`;
* an affine function of the optimization variable :math:`\gamma`.

Moreover, :math:`\mathcal{H}` is the subspace of functions that satisfy :math:`m` homogeneous boundary conditions,  i.e.

.. math::

	\mathcal{H} := \left\{ u \in C^k([a,b],\mathbb{R})
	\quad
	a_1 u(a) + a_2 u(b) + a_3 u'(a) + \cdots + a_{2k} u^k(b) = 0\right\},

where :math:`a_0,\,\ldots,\,a_{2k} \in \mathbb{R}^m` are known vectors.

.. note::

	Inhomogeneous boundary conditions can be "lifted" by changing variables according to :math:`u(x)=v(x)+p(x)`, where :math:`p(x)` is a polynomial of sufficiently high degree satisfying the inhomogeneous boundary conditions.

A particularly simple example of an optimization problem with an integral inequality is to determine the best Poincaré constant, i.e. the largest :math:`\gamma > 0` such that

.. math::

   \int_0^1 \left[
   \vert u'(x)\vert^2 - \gamma \vert u(x)\vert^2 \right] {\rm d}x \geq 0
   \quad \forall u \in C^2([0,1],\mathbb{R}),\quad u(0)=0=u(1).

Upper and lower bounds on the largest :math:`\gamma` are found by QUINOPT upon solving two SDPs.

.. note::
   QUINOPT can also handle problems with more dependent variables, i.e. :math:`u:[a,b]\to\mathbb{R}^q`, and problems in which the boundary values of the dependent variables and their derivatives appear explicitly in the integrand of the inequality constraint.




----------------------

* :doc:`Back to Table of Contents <../index>`
