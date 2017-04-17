function L = setLinearForm(coeff,monom,ivars,allvars)

% SETLINEARFORM.m Set up linear form from polynomial
%
% Set up vector that represents the linear part of a quadratic polynomial. 
% The order of the variables is:
% 1) First, all dependent variables and derivatives
% 2) Then, all boundary variables with the same order
%
% e.g. with 2 dependent variables with 1 derivative each, Q multiplies vector 
% [u u_x v v_x u(-1) u(1) u_x(-1) u_x(1) v(-1) v(1) v_x(-1) v_x(1)]

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    17/04/2017
% Last Modified:    17/04/2017
% ----------------------------------------------------------------------- %

nterms = length(coeff);
nvars = length(allvars);
L = zeros(nvars,1);

% UGLY FIX FOR LEGPOLY OBJECTS
% Otherwise assignments does not work
if isa(coeff,'legpoly')
    L = legpoly(L);
end

% LOOP TO SET QUADRATIC TERM
for i = 1:nterms
    
    vars = ivars(monom(i,:)~=0);    % variables in the monomial
    ind = find(allvars==vars);      % position in global indexing
    
    % Assign linear form vector
    if isa(coeff(i),'sdpvar') && ~isa(L,'sdpvar')
        % Due to behaviour of YALMIP variables, assignment of sdpvar to
        % double results in NaN, so need to add full matrices.
        L = L + sparse(ind(1),1,coeff(i),nvars,1);
    else
        L(ind) = coeff(i);
    end
    
    
end
