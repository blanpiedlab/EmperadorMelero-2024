%% function joinpath = picassoJoin(pathcell,varargin)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Function to join multiple picasso localization files using cmd
% Inputs
% Required
%   pathcell = cell array of absolute paths you want to join. note it will
%   merge in order of filenames in pathcell array so if you want it
%   numerical order you'll need to adjust manually before input
% Optional
%   keep = default 0; set to 1 to not reindex frame number
%   verbose = defaut 1; set to 0 to suppress any cmd line info. minimal for
%   this code
% Outputs
% joinpath = absolute path to result file

function joinpath = picassoJoin(pathcell,varargin)

    % Check for picasso command to appear in cmd line
    [status,~] = system('picasso join -h');
    if status ~= 0
        error('Something is not right in the system environment variable for cmd to find the path to picassopip env.')
    end

    % Parse some inputs
    p = inputParser;
  
    addRequired(p,'pathcell',@iscell)
    addParameter(p,'keep',0,@isnumeric)
    addParameter(p,'verbose',1,@isnumeric)
    parse(p,pathcell,varargin{:})
    pathcell = p.Results.pathcell;
    
    % Check input file is tif or raw
    [~,~,ext] = fileparts(pathcell);
    if sum(strcmp(ext,{'.hdf5'}))<size(ext,1)
        error('You did not feed this function hdf5 files try again.')
    end
    
    S = 'picasso join';
    for i = 1:size(pathcell,1)
        
        S = [S ' "' pathcell{i} '"'];
        
    end
   
    if p.Results.keep == 1
        
        S = [S ' -k'];
        
    end

    switch p.Results.verbose
        case 1
            [status,~] = system(S,'-echo'); % do localization and print python cmd
        case 0
            [status,~] = system(S); % do localization without print
    end
    
    % Output whether it worked
    if status == 0
        fprintf('Joining complete!\n')
    else
        fprintf('There was an error in joining, try again!\n')
    end
    
    [fold,fname,~] = fileparts(pathcell{1});
    joindir = dir(fullfile(fold,[fname '_join.hdf5']));
    joinpath = fullfile({joindir.folder},{joindir.name});
    joinpath = joinpath{:};
end
