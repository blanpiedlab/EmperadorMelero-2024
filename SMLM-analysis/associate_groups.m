%% associatedGroups = associate_groups(clusterLoc,testLoc,col,dbclust,clustcol,ringwidth,boxsize,psize)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Poorna Dharmasri, Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% This function takes a synaptic cluster and tests which groups of another
% protein have centroids within the dilated polyshape of that cluster
% Requires dilateCluster function
%
% REQUIRED INPUTS
% clusterLoc = entire loc array from Omniloader of this cluster to test (ie 1
% synapse)
% testLoc = entire loc array from Omniloader of test localizations for
% overlap
% col = col struct from getColumns for testLoc
% dbclust = the "loc" parameter of Omniloader for the testLoc-associated
% dbclusters file
% clustcol = the "col" parameter of the same
% ringwidth = distance you want to dilate polyshape of clusterLoc in nm; negative values contract
% boxsize = search distance around clusterLoc to look for testLoc groups
% psize = image pixelsize, ie 160
% OUTPUTS
% associatedGroups = group number(s) of testLoc that are within dilated
% polyshape of clusterLoc
% Written 4/2023 ADL; dilation code from PD

function associatedGroups = associate_groups(clusterLoc,testLoc,col,dbclust,clustcol,ringwidth,boxsize,psize)

 	p = inputParser;
    addRequired(p,'clusterLoc',@isnumeric)
    addRequired(p,'testLoc',@isnumeric)
    addRequired(p,'col',@isstruct)
    addRequired(p,'dbclust',@isnumeric)
    addRequired(p,'clustcol',@isstruct)
    addRequired(p,'ringwidth',@isnumeric)
    addRequired(p,'boxsize',@isnumeric)
    addRequired(p,'psize',@isnumeric)
    parse(p,clusterLoc,testLoc,col,dbclust,clustcol,ringwidth,boxsize,psize)
    clusterLoc = p.Results.clusterLoc;
    testLoc = p.Results.testLoc;
    col = p.Results.col;
    dbclust = p.Results.dbclust;
    clustcol = p.Results.clustcol;
    ringwidth = p.Results.ringwidth;
    boxsize = p.Results.boxsize;
    psize = p.Results.psize;
    
    % dilate the test cluster
    dilatedCluster = dilateCluster(clusterLoc(:,[col.x col.y]),ringwidth,psize);

    % find testloc groups within the box around this cluster
    ESboxwidth = boxsize/psize;
    ESbox_minx = min(clusterLoc(:,col.x)) - ESboxwidth;
    ESbox_maxx = max(clusterLoc(:,col.x)) + ESboxwidth;
    ESbox_miny = min(clusterLoc(:,col.y)) - ESboxwidth;
    ESbox_maxy = max(clusterLoc(:,col.y)) + ESboxwidth;

    EStestloc = testLoc(testLoc(:,col.x) < ESbox_maxx & testLoc(:,col.x) > ESbox_minx ...
        & testLoc(:,col.y) < ESbox_maxy & testLoc(:,col.y) > ESbox_miny,:); 
    inboundtest = unique(EStestloc(:,col.groups));
    
    % test whether the testloc group centers are within the expanded
    % cluster
    testcenters = nan(size(inboundtest,1),3);
    for o = 1:size(inboundtest,1)
        testcenters(o,:) = [inboundtest(o) dbclust(dbclust(:,clustcol.groups) == inboundtest(o),[clustcol.x clustcol.y])];
    end
    
    if ~isempty(testcenters)
        associatedGroups = testcenters(isinterior(dilatedCluster, testcenters(:,2:3)),1);
    elseif isempty(testcenters)
        associatedGroups = [];
    end
    
end
   
