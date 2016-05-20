function quinoptFeasCode
% 
% QUINOPTFEASCODE.m Guide to output codes returned by QUINOPT.
% 
% The codes returned in SOL.FeasCode by QUINOPT are the same as the problem
% codes in YALMIP. Specifically:
% 
%      -9 Specified solver name not recognized
%      -8 Problem does not satisfy geometric programming rules
%      -7 Solver does not return error codes
%      -6 Search space not bounded (bound all variables)
%      -5 License problems in solver
%      -4 Solver not applicable
%      -3 Solver not found in MATLAB path
%      -2 Successfully solved
%      -1 Unknown error
%       0 Successfully solved
%       1 Infeasible problem
%       2 Unbounded objective function
%       3 Maximum #iterations or time-limit exceeded
%       4 Numerical problems
%       5 Lack of progress
%       6 Initial solution infeasible
%       7 YALMIP sent incorrect input to solver
%       8 Feasibility cannot be determined
%       9 Unknown problem in solver
%      10 bigM failed (obsolete)  
%      11 Other identified error (use savesolveroutput and refer to solver
%      documentation)
%      12 Infeasible or unbounded
%      13 YALMIP cannot determine status in solver
%      14 Model creation failed
%      15 Problem either infeasible or unbounded
%      16 User terminated
%      17 Presolve recovery failed
%      18 Missing non-negativity bounds in GP formulation
%      19 Convexity requirements not met
%
% See also QUINOPT, YALMIPERROR

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

help quinoptFeasCode