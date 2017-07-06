function display(x)

%% OVERLOADED: depvarpoly/display

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    16/04/2015
% Last Modified:    16/04/2016
% ----------------------------------------------------------------------- %

[m,n] = size(x);
deg = degree(x);
try
    vars = unique(horzcat(x.ivars));
catch
    vars = unique(vertcat(x.ivars));
end
classification = sprintf('%ix%i dvarpoly (%i variables, degree %i)',m,n,length(vars),deg);
disp([classification])