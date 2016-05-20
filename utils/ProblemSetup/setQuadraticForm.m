function Q = setQuadraticForm(coeff,monom,ivars,allvars)

% SETQUADRATICFORM.m Set up quadratic form from polynomial
%
% Set up matrix that represents a quadratic polynomial. Indices are sorted
% so Q is always upper triangular. Order of variables is:
% 1) First, all dependent variables and derivatives
% 2) Then, all boundary variables with the same order
%
% e.g. with 2 dependent variable with 1 derivative each, Q multiplies vector 
% [u u_x v v_x u(-1) u(1) u_x(-1) u_x(1) v(-1) v(1) v_x(-1) v_x(1)]

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    11/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

nterms = length(coeff);
nvars = length(allvars);
Q = zeros(nvars);

% UGLY FIX FOR LEGPOLY OBJECTS
% Otherwise assignments does not work
if isa(coeff,'legpoly')
    Q = legpoly(Q);
end

% LOOP TO SET QUADRATIC TERM
for i = 1:nterms
    
    vars = ivars(monom(i,:)~=0);    % variables in the monomial
    if length(vars)==1
        % square monomial
        ind = [find(allvars==vars), find(allvars==vars)];
    else
        % quadratic monomial
        ind = sort( [find(allvars==vars(1)), find(allvars==vars(2))] );
    end
    
    % Assign quadratic form
    if isa(coeff(i),'sdpvar') && ~isa(Q,'sdpvar')
        % Due to behaviour of YALMIP variables, assignment of sdpvar to
        % double results in NaN, so need to add full matrices.
        Q = Q + sparse(ind(1),ind(2),coeff(i),nvars,nvars);
    else
        Q(ind(1),ind(2)) = coeff(i);
    end
    
    
end
