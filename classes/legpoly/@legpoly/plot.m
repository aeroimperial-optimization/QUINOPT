function plot(varargin)

% @legpoly/PLOT.m  Plot legendre polynomials
%
% PLOT(X,P) plots the Legendre polynomial P at the points specified by the
%         vector X. All points in X must be within the domain of definition
%         [a,b] of the Legendre polynomial
%
% PLOT(X,P,LINESPEC) and PLOT(P,LINESPEC) apply the formatting specified by
%         LINESPEC, as in the normal plot command.
%
% See also LEGPOLY, INDVAR, PLOT
%

if nargin<2; error('Not enough inputs.'); end
if nargin<3; varargin(3:4) = {[], []}; end

if isa(varargin{1},'legpoly') && isnumeric(varargin{2})
    plot(varargin{2},varargin{1},varargin{3:end});
    
elseif isnumeric(varargin{1}) && isa(varargin{2},'legpoly')
    
    % Good range?
    x = sort(varargin{1}(:));
    DOM = getDomain(varargin{2});
    if x(1)<DOM(1)-1e-15 || x(end)>DOM(2)+1e-15
        error('Range is outside the domain of the Legendre polynomial')
    end
    
    % plot
    plotlegpoly(x,varargin{2},varargin(3:end));
    
else
    error('The inputs should be a numeric vector and a legendre polynomial (class legpoly).')
end


% Nested function
    function plotlegpoly(x,q,format)
        
        hold_on = [];
        for i = 1:numel(q)
            [p,warn] = value(q(i),1);
            if ~warn
                % Get hold status
                if isempty(hold_on)
                    hold_on = ishold;
                end
                
                % Plot
                plot(x,legpolyval(p,x),format{:})
                if numel(q) > 1
                    hold on
                end
            else
                [r,c] = ind2sub(size(q),i);
                warning('legpoly:plot',...
                    'Could not plot entry (%i,%i) due to NaNs.',r,c)
            end
        end
        
        % Reset hold option
        if ~hold_on
            hold off
        end
    end

end