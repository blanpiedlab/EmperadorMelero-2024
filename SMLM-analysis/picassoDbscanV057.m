%% [dbpath,dbyaml,clusterpath] = picassoDbscan(locpath,radius,minpts,psize,varargin)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Function picassoDbscan runs DBSCAN from Picasso via the command line.
% Requires Picasso installed and set correctly to be on the main cmd
% path variable. This version is required for picasso 0.5.7 and above
% including picasso 0.6.0
% REQUIRED INPUTS
% - locpath: full path to the Picasso hdf5 loc file you want to run
%    DBSCAN on. ie 'Y:\Aaron\picassotest\dbscan weeks\New Folder\GluN2A_locs_filter_render_render.hdf5'
% - radius: radius in pixels of DBSCAN, equivalent to epsilon vallue
% - minpts: minimum points in radius required for cluster for DBSCAN
% - psize: pixelsize in nm ie 160
% NAME/VALUE OPTIONAL INPUTS
% - verbose: default 1, set to 0 to suppress Python cmd line outputs
% and only print matlab outputs
% - spotcheck: default 0, set to 1 to plot a figure showing original and
% clustered locs.
% OUTPUTS
% - dbpath = full path to '_dbscan.hdf5' file
% - dbyaml = full path to '_dbscan.yaml' file
% - clusterpath = full path to '_dbclusters.hdf5' file
% NOTE this also works on v0.6.0

function [dbpath,dbyaml,clusterpath] = picassoDbscanV057(locpath,radius,minpts,psize,varargin)
    
    % Check for picasso command to appear in cmd line
    [status,~] = system('picasso dbscan -h');
    if status ~= 0
        error('Something is not right in the system environment variable for cmd to find the path to picassopip env.')
    end

    p = inputParser;
    addRequired(p,'locpath',@ischar)
    addRequired(p,'radius',@isnumeric)
    addRequired(p,'minpts',@isnumeric)
    addRequired(p,'psize',@isnumeric)
    addParameter(p,'verbose',1,@isnumeric)
    addParameter(p,'spotcheck',0,@isnumeric)
    parse(p,locpath,radius,minpts,psize,varargin{:})

    % String to pass to command line
    S = ['picasso dbscan "' p.Results.locpath '" ' num2str(p.Results.radius) ' ' num2str(p.Results.minpts) ' ' num2str(p.Results.psize)]; 
    
    out = 'Running Picasso DBSCAN with %1.2f pixel radius and %d min points.\n';
    fprintf(out,p.Results.radius,p.Results.minpts)

    switch p.Results.verbose
        case 1
            [status,~] = system(S,'-echo'); % Run the command
        case 0
            [status,~] = system(S); % Run the comand
    end

    if status == 0
        fprintf('DBSCAN is complete!\n')
    else
        warning('There was an error in the input, try again!\n')
    end

    %output files
    [folder,file,~] = fileparts(locpath);
    clusterpath = fullfile(folder,[file '_dbclusters.hdf5']);
    dbpath = fullfile(folder,[file '_dbscan.hdf5']);
    dbyaml = fullfile(folder,[file '_dbscan.yaml']);

    if p.Results.spotcheck

        [hs,loc] = Omniloader(locpath,'verbose',0);
        col = getColumns(hs);
        [dbhs,dbloc] = Omniloader(dbpath,'verbose',0);
        dbcol = getColumns(dbhs);
        figure
        hold on
        plotOneChannel(loc,col,'c','k')
        gscatter(dbloc(:,dbcol.x),dbloc(:,dbcol.y),dbloc(:,dbcol.groups))
        legend('off')
        sgtitle([file ' OG locs (black) and clusters (colors)'],'Interpreter','none')

    end

end
