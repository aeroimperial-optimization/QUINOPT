function runAll

%% RUNALL.m Run all QUINOPT examples
%
% Run all examples included with QUINOPT to test the installation.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    29/05/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% Get example names
d=dir('example*.m');
demotime = zeros(length(d),1);

% Loop over all examples
for k=1:length(d)
    time = tic;
    evalin('caller',['run(''',d(k).name(1:end-2),''')']);
    demotime(k) = toc(time);
end

% End
evalin('caller','close(gcf); clear;')
disp('All tests completed. QUINTOPT has been successfully installed.')

end