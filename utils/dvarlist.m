function varargout = dvarlist(varargin)

% DVARLIST Set up persistent list of dependent variable monomials

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    01/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Declare persistent & initialise if empty
persistent VLIST
if isempty(VLIST)
    VLIST = 0;
end

% Set outputs
switch nargin
    
    % ------------------------------------------------------------------- %
    case 0
        
        varargout{1} = VLIST;
        return
        
        % ------------------------------------------------------------------- %
    case 1
        
        if ischar(varargin{1}) && strcmpi(varargin{1},'clear')
            % Clear all depvarpoly variables and reset list
            W = evalin('caller','whos');
            for i = 1:size(W,1)
                if strcmp(W(i).class,'dvarpoly')
                    evalin('caller', ['clear ' W(i).name ';']);
                end
            end
            VLIST = 0;
            return
            
        elseif isnumeric(varargin{1})
            % Add new variables to VLIST if they does not exist
            newvars = varargin{1}(~ismember(varargin{1},VLIST));
            VLIST = [VLIST, newvars];
            
        end
        
        % ------------------------------------------------------------------- %
    case 2
        
        if isnumeric(varargin{1}) && strcmpi(varargin{2},'remove')
            % Remove dependent variable from list and clear all depvarpoly
            % that depend on it
            VLIST = VLIST(~ismember(VLIST,varargin{1}));
            W = evalin('caller','whos');
            for i = 1:size(W,1)
                if strcmp(W(i).class,'dvarpoly')
                    thevar = evalin('caller',W(i).name);
                    if ismember(thevar.vars,varargin{1})
                        evalin('caller', ['clear ' W(i).name ';']);
                    end
                end
            end
            return
        end
        
end
end