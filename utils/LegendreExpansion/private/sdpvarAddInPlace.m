function Z = sdpvarAddInPlace(X,Y,varargin)

% SDPVARADDINPLACE Add entries to submatrix of SDPVAR without calling subsasn
% 
% Add Y to X(I) in place. I is a linear index.

% Get variables and basis
Y = Y(:);
Y_basis = getbase(Y);
Y_vars = getvariables(Y);
Y_length = length(Y);

% Get input indices
if nargin==3
    % Linear indices
    I = varargin{1};
elseif nargin==4
    % Matrix indices
    l1 = length(varargin{1});
    l2 = length(varargin{2});
    
    if Y_length==l1 &&  l1==l2 
        I = sub2ind(size(X),varargin{1}(:),varargin{2}(:));
    elseif Y_length==l1*l2
        [C,R] = meshgrid(varargin{2}(:),varargin{1}(:));
        I = sub2ind(size(X),R(:),C(:));
    end
end

% Index within limits?
if max(I)>numel(X) || any(I<0)
    error('Index exceeds matrix dimensions.')
end

% Expand Y
[nX,mX] = size(X);
[Iy,Jy,Vy] = find(Y_basis);
Y_basis =  sparse(I(Iy),Jy,Vy,numel(X),1+length(Y_vars));
Y = sdpvar(nX,mX,[],Y_vars,Y_basis);

% Set output
Z = X+Y;




