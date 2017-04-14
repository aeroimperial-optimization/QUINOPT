function installQUINOPT

%% INSTALLQUINOPT.m Install QUINOPT
%
% Install QUINOPT and removes previous versions from the MATLAB path.

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    29/05/2015
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

% CODE

currDir = pwd;

% ----------------------------------------------------------------------- %
% CHECK MATLAB VERSION
% ----------------------------------------------------------------------- %
if verLessThan('matlab','R2008a')
    warning('QUINTOPT has only been tested for compatibility for Matlab versions from 7.6 onwards.  Use at your own risk.');
end

% ----------------------------------------------------------------------- %
% CHECK IF YALMIP IS INSTALLED
% ----------------------------------------------------------------------- %
detected = which('yalmip.m','-all');
% Copied and edited from yalmiptest.m, YALMIP
if isa(detected,'cell') && ~isempty(detected)
    if length(detected)>1
        disp(['You seem to have multiple installations of YALMIP in your path. '...
            'Please correct this, then run the installation again.']);
        detected
        return
    end
    
    % ------------------------------------------------------------------- %
    % CHECK IF SDP SOLVER IS INSTALLED
    % ------------------------------------------------------------------- %
    % Check for general purpose solvers interfaced with YALMIP:
    % SeDuMi
    % SDPT3
    % Mosek
    % SDPA
    % CSDP (from OPTIToolbox)
    isavailable = 0;
    s = exist('sedumi.m','file'); s = (s~=0) & (s~=7); isavailable = isavailable | s;
    s = exist('sdpt3.m','file');  s = (s~=0) & (s~=7); isavailable = isavailable | s;
    s = exist('mosekopt','file'); s = (s~=0) & (s~=7); isavailable = isavailable | s;
    s = exist('sdpam.m','file');   s = (s~=0) & (s~=7); isavailable = isavailable | s;
    s = exist('csdp','file');   s = (s~=0) & (s~=7); isavailable = isavailable | s;
    if ~isavailable
        warn = sprintf(['\n(1) WARNING: none of the general purpose SDP solvers recommended for YALMIP could be found.\n',...
            '    You need an SDP solver to use QUINTOPT, so please install one.\n'...
            '    QUINTOPT has been tested using SeDuMi, SDPT3, Mosek and SDPA; click <a href="http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.Solvers">here</a>\n'...
            '    for a complete list of SDP solvers for YALMIP. If you have already\n',...
            '    installed a YALMIP-compatible SDP solver, please ignore this warning.\n']);
    else
        warn = [];
    end
    
    
    % ------------------------------------------------------------------- %
    % IF YALMIP INSTALLED, INSTALL
    % ------------------------------------------------------------------- %
    
    % Find all versions, cd to this versions folder and test if any other
    % versions on path
    allversions = which('installQUINOPT','-all');
    allversions = cellfun(@(x)fileparts(x),allversions,'uniformoutput',0);
    cd(allversions{1})
    onPath = ~cellfun(@isempty,allversions); % 1 if folder is on path
    
    % Remove old versions (if any) to avoid conflicts
    for i=2:length(allversions)
        if onPath(i)
            rmpath(genpath(allversions{i})); % remove current version from path
        end
    end
    
    % Add new version - platform dependent
    % Also try to remove .git folder by default (disable warnings first, then
    % enable them again to avoid long list of warnings when .git not on path)
    addpath(genpath(allversions{1}));
    rmpath(genpath([allversions{1},filesep,'examples']));
    rmpath(genpath([allversions{1},filesep,'docs']));
    rmpath(genpath([allversions{1},filesep,'papers']));
    wrn =  warning('query','all');
    warning('off','all');
    rmpath(genpath([allversions{1},filesep,'.git'])); 
    warning(wrn);
    savepath;
    
else
    error(['A working version of YALMIP is required. '...
        'Please correct this, then run the installation again.']);
end

% ----------------------------------------------------------------------- %
% RUN ALL DEMOS
% ----------------------------------------------------------------------- %
if isavailable
    cd('examples/')
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
else
    disp('QUINTOPT has been installed, but the following warning was generated:')
    disp(warn);
    disp('Please try to resolve this issue, then test QUINTOPT by running the demos in the folder "demos."')
end
cd(currDir)

%% END SCRIPT