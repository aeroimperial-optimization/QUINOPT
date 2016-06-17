function varargout = horzcat(varargin)

% OVERLOADED: depvar/horzcat
% Prevent concatenation of depvar objects

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    25/04/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

error('Concatenation of depvar objects is not allowed.')