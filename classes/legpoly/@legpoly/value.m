function [p,warn] = value(p,silent)

% @legpoly/VALUE.m Overloaded sdpvar/value on legpolys
%
% q = value(p), where p is a polynomial in the Legendre basis (class <a href="matlab:help('legpoly')">legpoly</a>),
%               returns a legendre polynomial q whose coefficients are the
%               evaluations of the coefficients of p when these are variables of
%               class <a href="matlab:help('sdpvar')">sdpvar</a>. For example:
%
%               x = indvar(-1,1);
%               p = legpoly(x,4);
%
%               Defines a degree-4 polynomial p with variable coefficients, to
%               be chosen so as to minimize an objective function. Then, after
%               the optimization procedure has been carried out, the optimal
%               choice of p is given by
%
%               q = value(p);
%
% See also LEGPOLY, INDVAR, SDPVAR, @SDPVAR/VALUE
%

if nargin<2; silent = 0; end

warn = 0;
for i = 1:numel(p)
    p.coef = value(p.coef);
    if any(isnan(p.coef))
        warn=1;
    end
end

if warn & ~silent
   warning('legpoly:value',['Some coefficients evaluated to NaN.\n',...
       'If there are unused variables in your polynomial, please consider removing them.']) 
end
