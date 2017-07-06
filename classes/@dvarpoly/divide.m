function p = divide(x,y)

%% OVERLOADED: dvarpoly/divide
% Element-wise division

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2017
% Last Modified:    25/04/2017
% ----------------------------------------------------------------------- %

% Check for size mismatch
[mx,nx] = size(x);
[my,ny] = size(y);
if (numel(x)~=1 && numel(y)~=1) && (mx~=my || nx~=ny)
    error('Size mismatch.')
end

% Re-use code
p = x.*(1./y);

% END FUNTION
end