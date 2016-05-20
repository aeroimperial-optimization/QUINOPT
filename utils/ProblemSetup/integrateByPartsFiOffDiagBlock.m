function [Fi,Fb,equalities,FLAG] = integrateByPartsFiOffDiagBlock(Fi,Fb,x,rmrows,rmcols)

% INTEGRATEBYPARTSFIOFFDIAGBLOCK.m Integrate by parts off-diagonal block of Fi
%
% Integrate by parts an off-diagonal block coupling two dependent variables.
% Infer derivative orders from size of the matrix. Inputs are
%
% - Fi, Fb: matrices describing the quadratic form to integrate by parts
% - x: variable of integration
% - rmrows, rmcols: first row/col that should be set to zero, that is rows
%                   rmrows:end and cols rmcols:end should become 0 at the
%                   end of the routine.
%
% Outputs:
%
% - Fi, Fb: the updated matrices after the integration by parts
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
FLAG = 0;
[row,col] = size(Fi);
equalities = [];

% Make sure Fb is an sdpvar - otherwise assignments produce NaN
if isempty(Fb)
    % Reinitialise to sdpvar with dummy entry
    clear Fb
    Fb(2*row,2*col) = x;
elseif isnumeric(Fb)
    % Add to an sdpvar with dummy entry to obtain sdpvar
    A(2*row,2*col) = x;
    Fb = A+Fb;
end

if rmrows==1 || rmcols==1
    % ----------------------------------------------------------------------- %
    % REMOVE ENTIRE MATRIX
    % ----------------------------------------------------------------------- %
    equalities = Fi(:);
    
else
    
    % ----------------------------------------------------------------------- %
    % INTEGRATE BY PARTS COLUMNS THAT NEED IT
    % ----------------------------------------------------------------------- %
    for i = col:-1:rmcols
        
        % Get column without last row
        P = Fi(1:row-1,i);
        if isZero(P) % Nothing to do!
            continue
        end
        
        % Update Fi
        J = jacobian(P,x);
        Fi(:,i-1) = Fi(:,i-1) - [0; P] - [J; 0];
        if isa(Fi(1:row-1,i),'legpoly')
            % can only assign object of same class since "subsasgn" was not
            % redefined for legpoly class
            Fi(1:row-1,i) = legpoly(0);
        else
            % subsasgn works for sdpvars!
            Fi(1:row-1,i) = 0;
        end
        
        % Find boundary values
        if isa(P,'legpoly')
            bv = getbv(P);
            BVatp1 = bv(:,1);
            BVatm1 = bv(:,2);
        else
            BVatp1 = replace(P,x,1);
            BVatm1 = replace(P,x,-1);
        end
        
        Fb(2:2:2*row-2,2*i-2) = Fb(2:2:2*row-2,2*i-2) + BVatp1;
        Fb(1:2:2*row-3,2*i-3) = Fb(1:2:2*row-3,2*i-3) - BVatm1;
        
    end
    
    % ----------------------------------------------------------------------- %
    % INTEGRATE BY PARTS ROWS THAT NEED IT
    % ----------------------------------------------------------------------- %
    % Integrate each row
    for i = row:-1:rmrows
        
        % Get row
        P = Fi(i,:);
        if isZero(P) % Nothing to do!
            continue
        end
        
        if ~isZero(P(rmcols-1:col))
            % Set last cols-rmcols+1 to zero since:
            % 1) cannot integrate by parts las column with existing derivatives
            % 2) if IP possible, would re-populate columns that have been zeroed earlier
            equalities = [equalities; P(rmcols-1:col)'];
            if isa(Fi,'legpoly')
                % See above for this distinction
                Fi(i,rmcols-1:col) = legpoly(0);
            else
                Fi(i,rmcols-1:col) = 0;
            end
        end
        
        % Integrate by parts the other columns (statement not entered if
        % rmcols==2, so no worries)
        if ~isZero(P(1:rmcols-2))
            
            % Update Fi
            J = jacobian(P(1:rmcols-2),x);
            Fi(i-1,1:rmcols-1) = Fi(i-1,1:rmcols-1) - [0 P(1:rmcols-2)] - [J 0];
            if isa(Fi,'legpoly')
                % Like above...
                Fi(i,1:rmcols-2) = legpoly(0);
            else
                % subsasgn works for sdpvars!
                Fi(i,1:rmcols-2) = 0;
            end
            
            % Find boundary values
            if isa(P,'legpoly')
                bv = getbv(P(1:rmcols-2));
                BVatp1 = bv(:,1)';
                BVatm1 = bv(:,2)';
            else
                BVatp1 = replace(P(1:rmcols-2),x,1);
                BVatm1 = replace(P(1:rmcols-2),x,-1);
            end
            
            Fb(2*i-2,2:2:2*rmcols-4) = Fb(2*i-2,2:2:2*rmcols-4) + BVatp1;
            Fb(2*i-3,1:2:2*rmcols-5) = Fb(2*i-3,1:2:2*rmcols-5) - BVatm1;
            
        end
    end
end

% ----------------------------------------------------------------------- %
% SET EQUALITY CONSTRAINTS FROM LIST IN "EQUALITIES"
% ----------------------------------------------------------------------- %
% Get coefficients of polynomial entries - "equalities" becomes an SDPVAR
if isempty(equalities)
    % All done
    return
    
elseif isa(equalities,'sdpvar')
    equalities = coefficients(equalities,x);
    
elseif isa(equalities,'legpoly')
    equalities = coefficients(equalities,x);
    if iscell(equalities)
        % need to convert (this is the case if have more than one equality constraint)
        equalities = vertcat(equalities{:});
    end
    
elseif ~isnumeric(equalities)
    error('Unknown class of variables for the coefficients.')
    
end

% Find the numeric entries & check for infeasible equality
numentries = ~cellfun('isclass',num2cell(equalities),'sdpvar');
if any(equalities(numentries))
    % Woops, we require nonzero number == 0
    FLAG = 1;
    return
    
else
    % Remove obvious equalities 0==0, keep nontrivial equalities
    equalities = equalities(~numentries);
    FLAG = 0;
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

% END
end