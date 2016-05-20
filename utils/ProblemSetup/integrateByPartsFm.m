function [Fm,Fb] = integrateByPartsFm(Fm,Fb,x,rmcols)

% INTEGRATEBYPARTSFM.m Integrate by parts mixed term matrix Fm
%
% Integrate by parts matrix of mixed terms Fm to set required columns to 0.
% This assumes that the matrix Fm multiplies only one dependent variable
% and its derivatives on the right (Fm in this function is a submatrix of
% the full mixed term matrix for a problem with multiple dependent
% variables). Inputs:
%
% - Fm, Fb: matrices to integrate by parts
% - x: variable of integration
% - rmcols: first col that should be set to zero, that is cols rmcols:end
%           should become 0 at the end of the routine.
%
% Outputs:
%
% - Fm, Fb: the updated matrices after the integration by parts
% - equalities: a list of variables that should be set to 0.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    12/04/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% PRELIMINARIES
% ----------------------------------------------------------------------- %
[row,col] = size(Fm);

% Make sure Fb is an sdpvar - otherwise assignments produce NaN
if isempty(Fb)
    % Reinitialise to sdpvar with dummy entry
    clear Fb
    Fb(row,2*col) = x;
elseif isnumeric(Fb)
    % Add to an sdpvar with dummy entry to obtain sdpvar
    A(row,2*col) = x;
    Fb = A+Fb;
end

% ----------------------------------------------------------------------- %
% INTEGRATE BY PARTS COLUMNS THAT NEED IT
% ----------------------------------------------------------------------- %
for i = col:-1:rmcols
    
    % Get column without last row
    P = Fm(:,i);
    if ~isZero(P)
        
        % Update Fm
        Fm(:,i-1) = Fm(:,i-1) - jacobian(P,x);
        if isa(Fm(1:row-1,i),'legpoly')
            % can only assign object of same class since "subsasgn" was not
            % redefined for legpoly class
            Fm(:,i) = legpoly(0);
        else
            % subsasgn works for sdpvars!
            Fm(:,i) = 0;
        end
        
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
        
        Fb(:,2*i-2) = Fb(:,2*i-2) + BVatp1;
        Fb(:,2*i-3) = Fb(:,2*i-3) - BVatm1;
        
    end
    
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
            qbv = [qbv; sum(qc), (1.^(0:length(qc)-1))*qc(:)];
        end
    end
% ----------------------------------------------------------------------- %

% END CODE
end