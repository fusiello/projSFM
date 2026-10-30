function setup_projSFM()
%SETUP_PROJSFM Add the projSFM folders to the MATLAB path.
%
%   Run once per session (or add a call to your startup.m):
%       setup_projSFM
%   Then try the test scripts, e.g.  test_solvability

root = fileparts(mfilename('fullpath'));
dirs = {'basicCV', 'graphs', 'projectivity', 'simulations', 'solvability'};
for k = 1:numel(dirs)
    addpath(fullfile(root, dirs{k}));
end
addpath(root);
fprintf('projSFM added to the MATLAB path (%s)\n', root);
end
