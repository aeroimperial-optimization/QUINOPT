function [Li,Lb] = integrateByPartsLi(Li,Lb,x)

% INTEGRATEBYPARTSLI.m Integrate by parts linear term matrix Li
%
% Integrate by parts vector of linear terms Li.
% This assumes that the matrix Li multiplies only one dependent variable
% and its derivatives on the right (Li in this function is a submatrix of
% the full mixed term matrix for a problem with multiple dependent
% variables). Inputs:
%
% - Li, Lb: matrices to integrate by parts
% - x: variable of integration
%
% Outputs:
%
% - Li, Lb: the updated matrices after the integration by parts

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    17/04/2017
% Last Modified:    17/04/2017
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% PRELIMINARIES
% ----------------------------------------------------------------------- %
[ind] = length(Li);

% Make sure Fb is an sdpvar - otherwise assignments produce NaN
if isempty(Lb)
    % Reinitialise to sdpvar with dummy entry
    clear Lb
    Lb(2*ind,1) = x;
elseif isnumeric(Lb)
    % Add to an sdpvar with dummy entry to obtain sdpvar
    A(2*ind,1) = x;
    Lb = A+Lb;
end

% ----------------------------------------------------------------------- %
% INTEGRATE BY PARTS ELEMENTS THAT NEED IT
% ----------------------------------------------------------------------- %
for i = ind:-1:2
    
    % Get elements to integrate
    P = Li(2:i);
    if ~isZero(P)
        
        % Update Li
%         Li(1:i-1) = Li(1:i-1) - jacobian(P,x);
%         if isa(Li(i),'legpoly')
%             % can only assign object of same class since "subsasgn" was not
%             % redefined for legpoly class
%             Li(i) = legpoly(0);
%         else
%             % subsasgn works for sdpvars!
%             Li(i) = 0;
%         end

        % First remove what I integrated
        if isa(Li(i),'legpoly')
            % can only assign object of same class since "subsasgn" was not
            % redefined for legpoly class
            Li(2:ind) = legpoly(0);
        else
            % subsasgn works for sdpvars!
            Li(2:ind) = 0;
        end
        
        % Then add jacobian
        Li(1:i-1) = Li(1:i-1) - jacobian(P,x);
        
        % Find boundary values
        if isa(P,'legpoly')
            bv = getbv(P);
            BVatp1 = bv(:,1);
            BVatm1 = bv(:,2);
        else
            % P is an sdpvar so use YALMIP commands
            BVatp1 = replace(P,x,1);
            BVatm1 = replace(P,x,-1);
        end
        
        Lb(2:2:2*i-2) = Lb(2:2:2*i-2) + BVatp1;
        Lb(1:2:2*i-3) = Lb(1:2:2*i-3) - BVatm1;
        
    end
    
end

% END CODE
end


% ----------------------------------------------------------------------- %
% NESTED FUNCTIONS
% ----------------------------------------------------------------------- %
function qbv = getbv(q)
    % Only works for vectors - very specific to this application.
    % First col gives, BV at 1, second col gives BV at -1
    qbv = [];
    for j = 1:numel(q)
        qc = coefficients(q(j));
        %             qbv = [qbv; sum(qc), (1.^(0:length(qc)-1))*qc(:)];        % Wrong?
        qbv = [qbv; sum(qc), ( (-1).^(0:length(qc)-1) )*qc(:)];
    end
end
% ----------------------------------------------------------------------- %

