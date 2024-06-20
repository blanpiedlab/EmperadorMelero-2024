%% [in, out, dbs_hs] = clusterFilter(dbscanfile,clusterfile,std_range,varargin)
% Function clusterFilter will take Picasso DBSCAN output and filter the
% clusters based on mean_frame and std_frame to remove clusters of spurious
% binding.
% REQUIRED INPUTS
% - dbscanfile: full path to the Picasso '_dbscan.hdf5' loc file you want to run
%     ie 'Y:\Aaron\picassotest\dbscan weeks\New Folder\GluN2A_locs_filter_render_render_dbscan.hdf5'
%       this file is just the loc file with added "group" column
% - clusterfile: full path to the Picasso '_dbclusters.hdf5'  file
%       associated with above
%      ie 'Y:\Aaron\picassotest\dbscan weeks\New Folder\GluN2A_locs_filter_render_render_dbclusters.hdf5'
%       this file has all the cluster statistics in it
% - std_range: a set of two values ie [5000 20000] representing the low and
%     high values to filter std_frame on
% NAME/VALUE OPTIONAL INPUTS
% - mean_range: default is 2, change to edit the multiplier on std of
%       mean_frame (ie, filtering is mean_range*std(mean_frame))
% - flag: default is 1, set to 0 to suppress figure showing cluster
% statistics
% - minn: default is 0, set to an integer to filter clusters based on a
%       minimum number of points in the cluster (ie, 'minn',50 will remove
%       clusters with <50 localizations
% - spotcheck: default is 0, set to 1 to show a plot of the original,
% rejected, and kept clusters.
% OUTPUTS
% - in = loc table of dbscanfile with only kept clusters
% - out = loc table of dbscanfile with only removed clusters
% - db_hs = headerstring from dbscan file
% updated 3/2023 to deal with picasso 0.6.0 which uses "frame" instead of
% "mean_frame" 

function [in, out, dbs_hs] = clusterFilterV060(dbscanfile,clusterfile,std_range,varargin)

    p = inputParser;
    addRequired(p,'dbscanfile',@isfile)
    addRequired(p,'clusterfile',@isfile)
    addRequired(p,'std_range',@isnumeric)
    addParameter(p,'mean_range',2,@isnumeric)
    addParameter(p,'flag',1,@isnumeric)
    addParameter(p,'minn',0,@isnumeric)
    addParameter(p,'maxn',10^10,@isnumeric)
    addParameter(p,'spotcheck',0,@isnumeric)
    parse(p,dbscanfile,clusterfile,std_range,varargin{:})
    
    % Load the files
    [dbs_hs, dbloc] = Omniloader(p.Results.dbscanfile,'verbose',0);
    dbcol = getColumns(dbs_hs);
    
    [clust_hs, clusters] = Omniloader(p.Results.clusterfile,'verbose',0);
    clustcol = getColumns(clust_hs);

    % Filter on mean frame value
    pdmeanframe = fitdist(clusters(:,clustcol.frame),'Normal');
    maxmeanframe = pdmeanframe.mu + p.Results.mean_range*pdmeanframe.sigma;
    minmeanframe = pdmeanframe.mu - p.Results.mean_range*pdmeanframe.sigma;
    mean_frame_keep = clusters(:,clustcol.frame)>minmeanframe & clusters(clustcol.frame)<maxmeanframe;

    % Filter on std frame range
    maxstdframe = p.Results.std_range(2);
    minstdframe = p.Results.std_range(1);
    std_frame_keep = clusters(:,clustcol.std_frame)>minstdframe & clusters(:,clustcol.std_frame)<maxstdframe;
    
    % Filter on max and min allowed n in a cluster. If not set in input, it
    % will filter from 0 to 10^10 so everything
    minn_keep = clusters(:,clustcol.n) >= p.Results.minn;
    max_keep = clusters(:,clustcol.n) <= p.Results.maxn;
    kept_cluster_tf = std_frame_keep & mean_frame_keep & minn_keep & max_keep;
    kept_clusters = clusters(kept_cluster_tf,clustcol.groups);

    in = dbloc(ismember(dbloc(:,dbcol.groups),kept_clusters),:);
    out = dbloc(~ismember(dbloc(:,dbcol.groups),kept_clusters),:);

    if p.Results.flag

        [~,file,~] = fileparts(dbscanfile);
        figure
        sgtitle(file,'Interpreter','none')
        subplot(1,2,1)
        hold on
        histogram(clusters(:,clustcol.frame),0:1000:max(dbloc(:,dbcol.frame)))
        xline(maxmeanframe,'-','max mean frame')
        xline(minmeanframe,'-','min mean frame')
        title('mean frame')
        xlabel('frame bin')
        ylabel('count')
        subplot(1,2,2)
        hold on
        histogram(clusters(:,clustcol.std_frame),0:1000:max(clusters(:,clustcol.std_frame)))
        xline(maxstdframe,'-','max std frame')
        xline(minstdframe,'-','min std frame')
        title('std frame')
        xlabel('frame bin')
        ylabel('count')

    end

    if p.Results.spotcheck

        [folder,file,~] = fileparts(dbscanfile);
        newfile = extractBetween(file,1,find(file=='_',1,'last')-1);
        [hs,loc] = Omniloader(fullfile(folder,[newfile{:} '.hdf5']),'verbose',0);
        col = getColumns(hs);
        figure
        hold on
        plotOneChannel(loc,col,'c','k')
        plotOneChannel(dbloc,dbcol,'c','r')
        plotOneChannel(in,dbcol,'c','g')
        legend({'Original locs','Rejected clusters','Kept clusters'})
        sgtitle([file 'Kept clusters (green) and rejected clusters (red)'],'Interpreter','none')

    end


end