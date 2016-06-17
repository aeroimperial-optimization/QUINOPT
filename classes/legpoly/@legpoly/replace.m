function P = replace(P,variables,values)

% OVERLOADED: legpoly/replace
% Replace sdpvar object inside a legpoly with a value

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

%% Input checks
if ~isnumeric(values)
    error('Input ''values'' must be numeric.')
    
elseif ~isa(variables,'sdpvar')
    error('Input ''variables'' must be an sdpvar object.')
    
end

%% Replace
% Check if asked to replace the independent variable - not allowed!
ivar = depends(getivar(P)); 
domn = getDomain(P);
if ismember(ivar,depends(variables))
   error('Replacing independent variable of a legpoly is not allowed.')
end

% Replace - check to remove trailing zeros and constant entries in P to
% replace the independent variable
[m,n] = size(P);

coef = cellfun(@(x)replace(x,variables,values),{P.coef},'uniformoutput',0);  
coef = cellfun(@removeTrailingZeros,coef,'uniformoutput',0);
[P(1:m,1:n).coef] = deal(coef{:});

ivar =  cellfun(@(y)(length(y)>1),{P.coef},'uniformoutput',0);
domn =  cellfun(@(y)(length(y)>1)*domn,{P.coef},'uniformoutput',0);
ivar =  cellfun(@(x,y)x*y,ivar,{P.ivar},'uniformoutput',0);
[P(1:m,1:n).ivar] = deal(ivar{:});
[P(1:m,1:n).domn] = deal(domn{:});

        