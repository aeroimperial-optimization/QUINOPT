What is QUINOPT?
================

QUINOPT (QUadratic INtegral OPTimisation) is an open-source add-on for `YALMIP <https://yalmip.github.io/>`_ to solve optimisation problems with polynomial quadratic integral inequality constraints. In the simplest form, these have the form

.. math::

	\begin{aligned}
	\min_{\gamma} \quad &c^T \gamma\\
	\text{subject to} \quad &\int_a^b Q_{\gamma}(x,u(x),u'(x),...,u^k(x)) 
	\,{\rm d}x \geq 0 \quad \forall u(x) \in \mathcal{H}
	\end{aligned}
	
where :math:`Q_{\gamma}(x,u(x),u'(x),...,u^k(x))` is 

1. polynomial in :math:`x \in [a,b]`
2. homogeneous quadratic in :math:`u(x),u'(x),...,u^k(x)`
3. affine the optimization (vector) variable :math:`\gamma \in \mathbb{R}^s`, and

.. math::

	\mathcal{H} := \left\{ u \in C^k([a,b],\mathbb{R})
	\quad
	a_1 u(a) + a_2 u(b) + a_3 u'(a) + \cdots + a_{2k} u^k(b) = 0\right\}
	
is the space of functions that satisfy the :math:`m` homogeneous boundary conditions specified through the vectors :math:`a_0,\,\ldots,\,a_{2k} \in \mathbb{R}^m`.


-----------
Quick Links
-----------
* :ref:`mastertoc`