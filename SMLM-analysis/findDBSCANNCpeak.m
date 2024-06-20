function peak = findDBSCANNCpeak(nclocs,psize,varargin)
%PAD for NDMAR project. This code borrows a strategy from our NC detection
%approach and finds a peak of highest density within a NC that has been
%detected using DBSCAN, which does not generate a point of highest density.
%It is advantageous for us to use the peak density point within a
%nanocluster to analyzing enrichments from, vs the center of the
%nanocluster, as the center is somewhat arbitrary whereas the peak density
%is likely where the protein of interest is most concentrated, even within
%our resolution limits. This code requires only the localizations belonging
%to a specific nanocluster to be input, though it takes optional inputs
%such as 'randomizations', to alter the number of random NCs we test below,
%and Tra, which is a factor that changes the search radius to look for
%nearest neighbor points in the real NC.
%PAD here to borrow this for LCv2KO analysis on 10.25.22. Added a section
%for this to handle 3D data.

p = inputParser;
addRequired(p,'nclocs',@isnumeric);
addRequired(p,'psize',@isnumeric);
addParameter(p,'randomizations','',@isnumeric);
addParameter(p,'Tra', '',@isnumeric);
addParameter(p,'is3D',0,@isnumeric);
parse(p,nclocs,psize,varargin{:});


randomizations = 20;
if ~isempty(p.Results.randomizations)
    randomizations = p.Results.randomizations;
end

Tra = 2.5;
if ~isempty(p.Results.Tra)
    Tra = p.Results.Tra;
end

%ncarea = alphaShape(nclocs(:,1), nclocs(:,2));
mmdout = NaN(randomizations,1);
for randnum = 1:randomizations
    if p.Results.is3D == 0
    randnc = get_cluster_randomized_ADL_PD2D_dbscanNC(p.Results.nclocs, 1, 0);
    locdnsrand = sort(pdist2(randnc,randnc));
    locdnsrand = locdnsrand(2:end,:);
    mmdrandnc = mean(locdnsrand(1,:));
    mmdout(randnum,1) = mmdrandnc;
    elseif p.Results.is3D == 1
    randnc = get_cluster_randomized_ADL_PD3D_dbscanNC(p.Results.nclocs, 1, p.Results.psize, 0);
    locdnsrand = sort(pdist2(randnc,randnc));
    locdnsrand = locdnsrand(2:end,:);
    mmdrandnc = mean(locdnsrand(1,:));
    mmdout(randnum,1) = mmdrandnc;
    end 
end


searchdistance = mean(mmdout)*Tra;

ncdns = sort(pdist2(p.Results.nclocs,p.Results.nclocs));
ncdns = ncdns(2:end,:);
ncdns = ncdns<searchdistance;
nearestneighbors = [];
for i = 1:size(ncdns,2)
    nearestneighbors = [nearestneighbors; sum(ncdns(:,i))];
end
potentialpeaks = p.Results.nclocs(nearestneighbors==max(nearestneighbors),:);
if p.Results.is3D == 0
peak = [mean(potentialpeaks(:,1)) mean(potentialpeaks(:,2))];
elseif p.Results.is3D == 1
    peak = [mean(potentialpeaks(:,1)) mean(potentialpeaks(:,2)) mean(potentialpeaks(:,3))];
end
end