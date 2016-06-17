function prd = legProd(c1,c2,DOMAIN)

%% LEGPROD coefficients of Legendre polynomial after multiplication by x
%
% prd = LEGPROD(c1,c2,DOMAIN) returns the coefficients of the product P1*P2,
% where P1, P2 are the polynomials whose Legendre coefficients are given by c1
% and c2. The domain of the independent variable does not matter at all.
%
% NOTE: The multiplication uses the recursion relationship for Legendre
%       polynomials in the form
%
%       xP_i(x) = ((i + 1)*P_{i + 1}(x) + i*P_{i - 1}(x))/(2i + 1)
%
% Adapted from the Python package <a href
% ="https://github.com/numpy/numpy/blob/v1.10.0/numpy/polynomial/legendre.py">numpy</a>.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    17/06/2016
% Last Modified:    17/06/2016
% ----------------------------------------------------------------------- %


%% Remove trailing zeros from coefficients
c1 = removeTrailingZeros(c1(:));
c2 = removeTrailingZeros(c2(:));
degc1 = length(c1)-1;
degc2 = length(c2)-1;

%% Find poly of lower order to make shortest loop
if degc1 > degc2
    c = c2(:);
    xs = c1(:);
else
    c = c1(:);
    xs = c2(:);
end

%% Multiply
if length(c) == 1
    c0 = c(1)*xs;
    c1 = 0;
    dc0 = length(c0)-1;
    dc1 = length(c1)-1;
elseif length(c) == 2
    c0 = c(1)*xs;
    c1 = c(2)*xs;
    dc0 = length(c0)-1;
    dc1 = length(c1)-1;
else
    nd = length(c);
    c0 = c(end-1)*xs;
    c1 = c(end)*xs;
    dc0 = length(c0)-1;
    dc1 = length(c1)-1;
    for i=3:length(c)
        tmp = c0;
        nd = nd-1;
        dxs = length(xs)-1;
        c0 = [c(end-i+1)*xs;zeros(dc1-dxs,1)]-[c1.*((nd-1)/nd);zeros(dxs-dc1,1)];
        c1 = [tmp;zeros(dc1+1-dc0,1)]+[legMulx(c1).*((2*nd-1)/nd);zeros(dc0-1-dc1,1)];
        dc0 = length(c0)-1;
        dc1 = length(c1)-1;
    end
end
prd = [c0;zeros(dc1+1-dc0,1)] + [legMulx(c1); zeros(dc0-1-dc1,1)];
        