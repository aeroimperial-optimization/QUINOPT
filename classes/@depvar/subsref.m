function B = subsref(u,S)

% OVERLOADER: depvar/subsref
%
% Overload subsref so that the depvar object u behaves like a function
% handle. Desired behaviour: the command u(x,d) should return the
% derivative d of u evaluated at x.
% Inputs are:
% - u: the depvar object
% - S: a structure with two fields, type and subs, as per MATLAB subsref.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

if strcmpi(S.type,'()')
    
    % Extract inputs & check
    x = S.subs{1};
    if length(S.subs)==1
        d = 0;
    elseif length(S.subs)==2
        d = S.subs{2};
        % Check if valid derivative
        if ~isnumeric(d) || ~isscalar(d) || d<0
            error('When calling U(POINT,DERIVATIVE), DERIVATIVE must be a nonnegative integer.');
        end
    else
        error('Too many input arguments.')
    end
    B = u.h(x,d);
    
else
    error('Indexing of type %s not defined for depvar objects.',S.type)
end

