function [Fi,Fb] = integrateByPartsFiDiagBlock(Fi,Fb,x)

% INTEGRATEBYPARTSFIDIAGBLOCK.m Integrate by parts diagonal block of Fi
%
% Integrate by parts a diagonal block (assume diagonal block is upper
% triangular): one dependent variable, infer derivative orders from size of
% the matrices. Create a matrix Fb of appropriate size if empty.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    12/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

n = size(Fi,2);

% Make sure Fb is an sdpvar - otherwise assignments produce NaN
if isempty(Fb)
    % Reinitialise to sdpvar with dummy entry
    clear Fb
    Fb(2*n,2*n) = x;
elseif isnumeric(Fb)
    % Add to an sdpvar with dummy entry to obtain sdpvar
    A(2*n,2*n) = x;
    Fb = A+Fb;
end

% Integrate by parts
for i = n:-1:2
    
    % Get column above the diagonal
    P = Fi(1:i-1,i);     
    if isZero(P) % Nothing to do!
        continue
    end
    
    % Update Fi
    J = jacobian(P,x);
    Fi(1:i-1,i-1) = Fi(1:i-1,i-1) - [0; P(1:end-1)] - J.*[ones(i-2,1); 0.5];
    if isa(Fi(1:i-1,i),'legpoly')
        Fi(1:i-1,i) = legpoly(0);
    else
        Fi(1:i-1,i) = 0;
    end
    
    % Find boundary values
    if isa(P,'legpoly')
        bv = getbv(P);
        BVatp1 = bv(:,1).*[ones(i-2,1); 0.5];
        BVatm1 = bv(:,2).*[ones(i-2,1); 0.5];
    else
        BVatp1 = replace(P,x,1).*[ones(i-2,1); 0.5];
        BVatm1 = replace(P,x,-1).*[ones(i-2,1); 0.5];
    end
    
    Fb(2:2:2*i-2,2*i-2) = Fb(2:2:2*i-2,2*i-2) + BVatp1;
    Fb(1:2:2*i-3,2*i-3) = Fb(1:2:2*i-3,2*i-3) - BVatm1;
    
end

% ----------------------------------------------------------------------- %
% NESTED FUNCTIONS
% ----------------------------------------------------------------------- %

    function qbv = getbv(q)
        % Only works for vectors - very specific to this application.
        % First BV at 1, then at -1
        qbv = [];
        for j = 1:numel(q)
            qc = coefficients(q(j));
            %qbv = [qbv; sum(qc), (1.^(0:length(qc)-1))*qc(:)];        % Wrong?
            qbv = [qbv; sum(qc), ( (-1).^(0:length(qc)-1) )*qc(:)];
        end
    end

%
% END
end