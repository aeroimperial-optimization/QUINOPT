.. include:: ../substitutions.txt

``quinopt()``
================================================

Minimize a cost function subject to quadratic integral inequality constraints, or reset QUINOPT.

------------

:Syntax: ``quinopt clear`` or ``quinopt('clear')``
:Description: clears QUINOPT's internal variables, to be used in combination with MATLAB's ``clear`` command.

------------

:Syntax: ``quinopt(EXPR)``
:Description: tests whether a quadratic integral inequality with integrand specified by ``EXPR`` is feasible using a finite dimensional relaxation based on semidefinite programming. If multiple integral inequalities must be tested simultaneously, ``EXPR`` can be a vector such that the i-th entry specifies the integrand of the i-th integral inequality. Each entry of ``EXPR`` must be a polynomial of the integration variable returned by the command ``indvar()``, and a quadratic polynomial of the dependent variables returned by the function ``depvar()``.

------------

:Syntax: ``quinopt(EXPR,BC)``
:Description: determines whether the integral inequalities, with integrand specified by ``EXPR``, are feasible for all dependent variables satisfying the homogeneous boundary conditions specified by the vector ``BC``. Specifically, ``BC`` is interpreted as the list of boundary conditions ``BC(1)=0``, ..., ``BC(end)=0``. Like ``EXPR``, ``BC`` must be created using the variables returned  by the commands ``indvar()`` and ``depvar()``.

------------

:Syntax: ``quinopt(EXPR,BC,OBJ)``
:Description: **minimizes** the objective function ``OBJ`` constrained by the integral inequalities specified by ``EXPR`` and ``BC``.

------------

:Syntax: ``quinopt(EXPR,BC,OBJ,OPTIONS)``
:Description: overrides the default options. ``OPTIONS`` is a structure containing any of the following fields:

    - ``OPTIONS.YALMIP``: a substructure containing the options for YALMIP, set with YALMIP's command ``sdpsettings()``.

    - ``OPTIONS.N``: an integer specifying the number of Legendre coefficients to use in the expansion of the dependent variable to obtain an SDP-representable relaxation of the quadratic integral inequality.

    - ``OPTIONS.method``: if set to ``'inner'`` (default), QUINOPT generates an inner approximation of the feasible set of the integral inequalities specified by ``EXPR``, i.e. the quadratic integral inequality is strenghtened. If set to ``'outer'``, an outer approximation is constructed, i.e. the integral inequality is relaxed into a weaker constraint.

    - ``OPTIONS.BCprojectorBasis``: string specifying which basis to use for the projection on the boundary conditions. If set to ``'rref'`` (default), use a "rational" basis. If set to ``'orth'``, use an orthonormal basis. The orthonormal basis may be preferable numerically, but it may destroy sparsity of the data.

    - ``OPTIONS.sosdeg``: the degree of the sum-of-squares polynomials used in the S-procedure to localize SOS constraints from the integral inequality to the integration domain. Default value: 6.

    - ``OPTIONS.solve``: if set to ``'true'`` (default), QUINOPT calls the solver specified by the YALMIP options (or YALMIP's default solver). If set to ``'false'``, QUINOPT does not call the solver, but simply sets up the YALMIP problem structure. In this case, additional outputs to QUINOPT are required (see below).

------------

:Syntax: ``quinopt(EXPR,BC,OBJ,OPTIONS,CNSTR)`` or ``quinopt(EXPR,BC,OBJ,OPTIONS,CNSTR,PARAMETERS)``
:Description: minimizes the objective function ``OBJ`` subjet to the integral inequalities specified by ``EXPR`` and ``BC``, and the additional constraints given by ``CNSTR``. ``CNSTR`` is a constraint object built with YALMIP. If ``CNSTR`` contains sum-of-square constraints, then the variable parameters in the polynomial expressions **must** be specified in the input vector ``PARAMETERS``. See YALMIP's function ``sos()`` and ``solvesos()`` for more details on specifying sum-of-squares constraints with YALMIP.

------------

:Syntax: ``[SOL,CNSTR,DATA] = quinopt(...)``
:Description: returns solution informations in the structure ``SOL``, the set of YALMIP constraint ``CNSTR`` used to solve the optimization problem, and a structure ``DATA`` containing all variables used to set up the constraints in ``CNSTR``. The solution structure ``SOL`` contains the following fields:

    - ``SOL.setupTime``: the time taken to set up the problem

    - ``SOL.solutionTime``: the time taken by YALMIP to solve the problem

    - ``SOL.problem``: code of problem encountered during setup. Values are:

        * ``0``: no problem
        * ``1``: ill-posed inequality
        * ``2``: infeasible relaxation

    - ``SOL.FeasCode``: code for the feasibility of the solution returned by YALMIP. Run ``quinoptFeasCode`` at the MATLAB command line to obtain a complete list.

    - ``SOL.YALMIP``: the solution structure returned by YALMIP. See YALMIP's functions ``optimize()`` and ``solvesos()`` for more details.

------------


* :doc:`Back to the list of main functions <./index>`
* :doc:`Back to Table of Contents <../index>`
