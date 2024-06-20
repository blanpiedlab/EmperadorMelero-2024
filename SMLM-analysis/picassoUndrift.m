%% [outpath,yamlpath,driftpath] = picassoUndrift(locpath,varargin)
% Function picassoUndrift runs Undrift by RCC from Picasso via the command line.
% Requires Picasso installed and set correctly to be on the main cmd
% path variable. 
% REQUIRED INPUTS
% - locpath: full path to the Picasso hdf5 loc file you want to run
%    undrift on. ie 'Y:\Aaron\picassotest\dbscan weeks\New Folder\GluN2A_locs_filter_render_render.hdf5'
% NAME/VALUE OPTIONAL INPUTS
% - seg (num): set to 1000 by default, define the number of frames to bin
% - d (num): set to 0 by default, suppresses showing the drift (set to 1 to
%       show; doesn't work great in this code)
% - mode (char): choose from {'std','render','framepair'} to change how
%       drift correct works, leave unset for default
% - fromfile (char): input full file path to existing drift file to apply
%       that drift rather than calculate it
% - verbose: default 1, set to 0 to suppress Python cmd line outputs
%       and only print matlab outputs
% OUTPUTS
% - outpath = full path to '_undrift.hdf5' file
% - yamlpath = full path to '_undrift.yaml' file
% - drfitpath = full path to '_drift.txt' file

function [outpath,yamlpath,driftpath] = picassoUndrift(locpath,varargin)

    % Check for picasso command to appear in cmd line
    [status,~] = system('picasso undrift -h');
    if status ~= 0
        error('Something is not right in the system environment variable for cmd to find the path to picassopip env.')
    end

    p = inputParser;
    addRequired(p,'locpath',@ischar)
    addParameter(p,'seg',1000,@isnumeric)
    addParameter(p,'d',0,@isnumeric)
    opts = {'std','render','framepair'};
    addParameter(p,'mode','null',@(x)mustBeMember(x,opts))
    addParameter(p,'fromfile','null',@(x)mustBeFile(x))
    addParameter(p,'verbose',1,@isnumeric)
    parse(p,locpath,varargin{:})
    
    locpath = p.Results.locpath;
    
    s = p.Results.seg;
    d = p.Results.d;
    m = p.Results.mode;
    f = p.Results.fromfile;

    if d == 1
        d = '';
    else
        d = ' -d';
    end
     
    if strcmp(m,'null') && strcmp(f,'null')
        S = ['picasso undrift' d ' -s ' num2str(s) ' "' locpath '"'];
    elseif ~strcmp(m,'null') && strcmp(f,'null')
        S = ['picasso undrift' d ' -s ' num2str(s) ' -m ' m ' "' locpath '"'];
    elseif strcmp(m,'null') && ~strcmp(f,'null')
        S = ['picasso undrift' d ' -s ' num2str(s) ' -f "' f '" "' locpath '"'];
    elseif ~strcmp(m,'null') && ~strcmp(f,'null')
        S = ['picasso undrift' d ' -s ' num2str(s) ' -m ' m ' -f "' f '" "' locpath '"'];
    end
    
    switch p.Results.verbose
        case 1
            [status,~] = system(S,'-echo'); %RCC drift correct
        case 0
            [status,~] = system(S); %RCC drift correct suppress python output
    end

    if status == 0
        fprintf('Drift correct complete!\n')
    else
        fprintf('There was an error in the input, try again!\n')
    end

    % output files
    [folder,file,~] = fileparts(locpath);
    outpath = fullfile(folder,[file '_undrift.hdf5']);
    yamlpath = fullfile(folder,[file '_undrift.yaml']);
    driftpath = fullfile(folder,[file '_drift.txt']);

end