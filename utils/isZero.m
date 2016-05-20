function y = isZero(x)

%% ISZERO.m Determine if input is matrix of zeros

% y = ISZERO(x) returns true if all elements in x are zero, false otherwise. 
%     Only works for matrix inputs (up to 2 dimensional arrays).

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

if ndims(x)>2
    error('iszero only works for matrix inputs (at most 2D arrays).')
end

if isempty(x)
    y = 1;

elseif isnumeric(x)
    y = all(~any(x));
    
elseif isa(x,'sdpvar')
    y = all(~any(any(x)));       % any(x) gives sparsity pattern of x
    
elseif isa(x,'legpoly')
    x = coefficients(x);
    if iscell(x);
        % x was multidimensional legpoly
        x = cellfun(@any,x,'uniformoutput',0);
        x = cellfun(@any,x);
        y = all(~any(x));
    elseif isa(x,'sdpvar')
        % x was a scalar legpoly with sdpvar coefficients
        y = all(~any(any(x)));
    else
        % x was a scalar legpoly with numeric coefficients
        y = all(~any(x));
    end
    
end

% ----------------------------------------------------------------------- %
end