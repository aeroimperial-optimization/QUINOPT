function X = blkdiag(varargin)

%% OVERLOADED: legpoly/blkdiag

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

if nargin < 1
    error('Not enough input arguments.')
    
elseif nargin > 20
    % avoid to many recursions
    A = blkdiag(varargin{1:20});
    X = blkdiag(A,varargin{21:end});
    
elseif nargin == 1
    X = varargin{1};
    
elseif nargin == 2
    [m1,n1] = size(varargin{1});
    [m2,n2] = size(varargin{2});
    X = [varargin{1}, zeros(m1,n2); zeros(m2,n1), varargin{2}];
    
else
    % recursion
    X = blkdiag(varargin{1},blkdiag(varargin{2:end})) ;
end
