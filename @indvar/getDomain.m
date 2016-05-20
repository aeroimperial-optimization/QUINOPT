function [domain,varID] = getDomain(x)

% GETDOMAIN.m Recover domain of indvar object
%
% DOMAIN = GETDOMAIN(x) returns the domain of the independent variable x
%       used by QuadIntIneq. DOMAIN is a 2-by-1 row vector such as [a,b].
%
% See also INDVAR, @INDVAR/SETDOMAIN

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

model = qiimodel('query');
varID = depends(x);
i = ismember(model.INDVARMODEL.IVARID,varID);
domain = vertcat(model.INDVARMODEL.DOMAIN(i,:));