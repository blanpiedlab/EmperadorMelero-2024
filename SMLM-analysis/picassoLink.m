%% [outpath,yamlpath] = picassoLink(locpath,radius,darkframes,varargin)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Function picassoLink runs the Link command from Picasso via the command line.
% Requires Picasso installed and set correctly to be on the main cmd
% path variable. 
% REQUIRED INPUTS
% - locpath: full path to the Picasso hdf5 loc file you want to link
%    ie 'Y:\Aaron\picassotest\dbscan weeks\New Folder\GluN2A_locs_filter_render_render.hdf5'
% - radius: radius in pixels to link within
% - darkframes: number of allowed dark frames for a linked localization to
%    still be linked
% NAME/VALUE OPTIONAL INPUTS
% - verbose: default 1, set to 0 to suppress Python cmd line outputs
%    and only print matlab outputs
% OUTPUTS
% - outpath = full path to '_link.hdf5' file
% - yamlpath = full path to '_link.yaml' file

function [outpath,yamlpath] = picassoLink(locpath,radius,darkframes,varargin)
    
    % Check for picasso command to appear in cmd line
    [status,~] = system('picasso link -h');
    if status ~= 0
        error('Something is not right in the system environment variable for cmd to find the path to picassopip env.')
    end
    
    p = inputParser;
    addRequired(p,'locpath',@ischar)
    addRequired(p,'radius',@isnumeric)
    addRequired(p,'darkframes',@isnumeric)
    addParameter(p,'verbose',1,@isnumeric)
    parse(p,locpath,radius,darkframes,varargin{:})
    
    % string to pass to cmd
    S = ['picasso link -d ' num2str(p.Results.radius) ' -t ' num2str(p.Results.darkframes) ' "' p.Results.locpath '"']; 
    
    out = 'Running Picasso Link with a %1.2f pixel radius and %d dark frames allowed.\n';
    fprintf(out,p.Results.radius,p.Results.darkframes);

    switch p.Results.verbose
        case 1
            [status,~] = system(S,'-echo'); % linking
        case 0
            [status,~] = system(S); % linking without cmd output
    end
    
    if status == 0
        fprintf('Linking is complete!\n')
    else
        fprintf('There was an error in the input, try again!\n')
    end

    %output file
    [folder,file,~] = fileparts(locpath);
    outpath = fullfile(folder,[file '_link.hdf5']);
    yamlpath = fullfile(folder,[file '_link.yaml']);

end
