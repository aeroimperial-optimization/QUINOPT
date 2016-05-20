function p = mpower(x,d)

%% OVERLOADED : dvarpoly/mpower

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

if ~isnumeric(d) || numel(d)~=1 || rem(d,1)~=0 || d<0
    error('Power d in x^d must be a non-negative integer.')
end

% Operate - Not efficient!
[m,n] = size(x);
if d == 0
    % Easy!
    p = ones(m,n);
    return
    
elseif d==1
    % Just return input
    p = x;
    return
    
else
    % Loop over entries
    p = struct('coeff',[],'ivars',[],'monom',[]);
    for i=1:m
        for j=1:n
            
            if size(x(i,j).monom,1)==1
                % A monomial - easy
                p(i,j).coeff = x(i,j).coeff.^d;
                p(i,j).monom = x(i,j).monom.*d;
                p(i,j).ivars = x(i,j).ivars;
                
            else
                dfact = factor(d);
                t = x(i,j);
                temp = x(i,j);               % dummy
                for k1 = 1:length(dfact)
                    for k2 = 2:dfact(k1)
                        temp = temp*t;
                    end
                    t = temp;
                end
%                 p(i,j) = setstructfields(p(i,j),dvarpoly2struct(t));
                p(i,j) = dvarpoly2struct(t);
                
            end
        end
    end
    
    % Assign class
    p = dvarpoly(p);
    return
    
end