function symb = sdisplay(varargin)

%% @LEGPOLY/SDISPLAY.m Symbolic display of legpoly expression
%
% SDISPLAY(P) displays the polynomial P of class legpoly to screen. The
%       polynomial is displayed using the standard monomial basis for
%       convenience.
%
% S = SDIPSLAY(P) returns the output to the cell array S, where S{i,j}
%       is the output of SDISPLAY(P(i,j)).
%
% See also LEGPOLY, @YALMIP/yalmip/extras/SDISPLAY

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    07/10/2015
% Last Modified:    08/04/2016
% ----------------------------------------------------------------------- %

% Get sdpvar variable names in base workspace
W = evalin('caller','whos');
for i = 1:size(W,1)
    if strcmp(W(i).class,'sdpvar')
        
        % Get the sdpvar variable name
        thevar = evalin('caller',W(i).name);
        eval([W(i).name,'=sdpvar(thevar);']);
        
    end
end
clear thevar;
tmp = sdisplay(sdpvar(varargin{1}));
if nargout
    symb = tmp;
else
    disp(tmp)
end