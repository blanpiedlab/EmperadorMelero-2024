%% [filteredLoc, overlap, nonoverlap] = findPSDoverlap_nlocs(refloc,testloc,col,boxnm,psize,nlocs)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Function findPSDoverlap will take a reference loc table with cluster data
% and a variable number of test loc tables and filter the reference loc
% table to only clusters that have overlap with the test loc tables. For
% example, if refloc is PSD95 clusters and testloc is Bassoon and N2B
% localizations, this function will return filteredLoc, which is the PSD95
% clusters that had overlap with both Bassoon and N2B localizations. Note
% that overlap currently only requires 1 localization from testloc to
% overlap with the convex hull of the refloc cluster. Function will also
% return the cluster number of those refloc clusters, as well as the
% cluster numbers of the refloc clusters that did not overlap with each
% testloc. This function works by looping through each refloc cluster,
% making an alpha shape, testing whether testloc is inShape, and putting
% together a list of overlapped refloc clusters.
% REQUIRED INPUTS
% - refloc: loc table that has a col.groups from picasso dbscan or other.
% This is the reference set of clusters you are filtering.
% - testloc: cell array that has localization tables you want to test for
% overlap with refloc. Must be input as a cell array and can be any number
% of loc tables. ie, {loc}, or {bsnLoc, n2aLoc}.
% - col: output of getColumns for refloc
% - psize: camera pixel size in nm
% - nlocs: set the minimum number of locs in each testloc channel that need
% to be within the alpha shape of refloc groups to be kept
% NAME/VALUE OPTIONAL INPUTS
% - boxnm: default is 500, this defines how large a box you want to make
% around the refloc cluster to look for inShape. 
% - flag: default is 1, set to 0 to suppress figure
% OUTPUTS
% - filteredLoc = refloc with only kept clusters
% - overlap = m x 1 array of cluster numbers of refloc that did have
% overlap with refloc (filteredLoc is just refloc(ismember(refloc(:,col.groups),overlap),:))
% - nonoverlap = length(testloc) x 1 cell array of cluster numbers of
% refloc that did not overlap with refloc. ie if testloc = {bsnLoc, n2aLoc}
% then nonoverlap{1} is the refloc cluster numbers that did not have any
% overlap with bsnLoc, and nonoverlap{2} is the same for n2aLoc.
% - Will save filteredClusterNumbers.mat which holds overlap and nonoverlap 
% Written ADL 2/1/2022
% Requires subfunctions sortByField

function [filteredLoc, overlap, nonoverlap] = findPSDoverlap_nlocs(refloc,testloc,col,boxnm,psize,nlocs)
    
    % Input block
    p = inputParser;
    addRequired(p,'refloc',@isnumeric)
    addRequired(p,'testloc',@iscell)
    addRequired(p,'col',@isstruct)
    addRequired(p,'boxnm',@isnumeric)
    addRequired(p,'psize',@isnumeric)
    addRequired(p,'nlocs',@isnumeric)
    parse(p,refloc,testloc,col,boxnm,psize,nlocs)
    refloc = p.Results.refloc;
    testloc = p.Results.testloc;
    col = p.Results.col;
    boxnm = p.Results.boxnm;
    psize = p.Results.psize;
    nlocs = p.Results.nlocs;

    % initialize some useful things
    ntest = length(testloc);
    b = boxnm/psize;
    refloccell = sortByField(refloc,col.groups);
    nonoverlap = cell(ntest,1);
    g = nan(length(refloccell),1);

    % find overlaps between each refloccell and testlocs
    for p = 1:length(refloccell)
    
        % Get the bounding box
        thisp = refloccell{p};
        if isempty(thisp)
            continue;
        end
        box_minx = min(thisp(:,col.x)) - b;
        box_maxx = max(thisp(:,col.x)) + b;
        box_miny = min(thisp(:,col.y)) - b;
        box_maxy = max(thisp(:,col.y)) + b;

        % Find all locs in the box
        inrectlocs = cell(ntest,1);
        for l = 1:ntest
            thisloc = testloc{l};
            inrectlocs{l} = thisloc(thisloc(:,col.x) < box_maxx & thisloc(:,col.x) > box_minx & thisloc(:,col.y) < box_maxy & thisloc(:,col.y) > box_miny,:);
        end
        
%         % Diagnostic figure of locs in the box wile testing code
%         figure
%         hold on
%         plotOneChannel(inrectlocs{1},col,'c','k')
%         plotOneChannel(inrectlocs{2},col,'c','b')
%         plotOneChannel(inrectlocs{3},col,'c','r')
%         plotOneChannel(inrectlocs{4},col,'c','g')
%         rectangle('Position',[box_minx box_miny (box_maxx - box_minx) (box_maxy - box_miny)])
%         set(gca,'YDir','reverse')

        % Store the group number
        g(p) = thisp(1,col.groups);
  
        % Get alpha shape of reference loc and test which testloc overlap
        shp = alphaShape(thisp(:,col.x),thisp(:,col.y));
        for t = 1:ntest
            inshploc = sum(inShape(shp,inrectlocs{t}(:,col.x),inrectlocs{t}(:,col.y)));
            if sum(inshploc) < nlocs
                nonoverlap{t} = [nonoverlap{t}; g(p)];
            end
        end

    end

    % Find index of refloccells that had overlap with all of testloc
    c = cell2mat(nonoverlap);
    keptclusterstf = ~ismember(g,unique(c));
    overlap = g(keptclusterstf);

    % Generate loc table of filtered locs
    filteredLoc = cell2mat(refloccell(keptclusterstf));
    
    % Save the overlap variables for posterity
    % save('filteredClusterNumbers.mat','overlap','nonoverlap')


end
