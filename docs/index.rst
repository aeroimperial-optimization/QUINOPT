.. QUINOPT documentation master file, created by
   sphinx-quickstart on Thu Apr 13 18:14:55 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QUINOPT's documentation!
===================================

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


-------
License
-------
QUINOPT is distributed under the `Apache 2.0 License <https://www.apache.org/licenses/LICENSE-2.0>`_

------------------------
Bug reports and support
------------------------

Please report any issues via the `Github issue tracker <https://github.com/aeroimperial-optimization/QUINOPT/issues>`_, or `contact us <mailto:giovanni.fantuzzi10@imperial.ac.uk?Subject=QUINOPT%20issue>`_. All types of issues are welcome, including bug reports, documentation typos, feature requests and so on.


.. toctree::
   :titlesonly:
   :maxdepth: 1
   :caption: Contents

   index
   requirements/index
   installation/index
   examples/index
