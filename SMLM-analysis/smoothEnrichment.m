%% function smoothed = smoothEnrichment(data,cutoffs)
% This function will smooth enrichment curves by doing the following:
% 1) run prepEnrichmentForOutliers on data to first convert trailing 0 to
% nan and crop out first column
% 2) Run ROUT with 0.1% cutoff in PRISM
% 3) Get descriptive stats on cleaned data per column
% 4) run this function using the max values from the descriptive stats in
% the cutoffs variable
% smoothed is the data but with outliers values for each column replaced
% with the highest non-outlier value in each column.
% REQUIRED INPUTS
% data = m x n array of enrichment curves, where m is number of
% nanoclusters that have curves and n is the number of bins, usually 32 if
% the first bin has been removed
% cutoffs = 1 x n matrix where n is equal to the number of columns n in
% data; the cutoffs are the maximum non-outlier values from step 3 above
% OUTPUTS
% smoothed = same as data but with outlier values replaced with highest
% nonoutlier values. nan are not changed.

function smoothed = smoothEnrichment(data,cutoffs)

    p = inputParser;
    addRequired(p,'data',@isnumeric)
    addRequired(p,'cutoffs',@isnumeric)
    parse(p,data,cutoffs)
    
    if size(data,2) ~= size(cutoffs,2)
        error('You do not have matching cutoff and data sizes try again')
    end
    
    tf = data>cutoffs;
    
    [~,ncol] = size(data);
    smoothed = data;
    for thiscol = 1:ncol
        thistf = tf(:,thiscol);
        smoothed(thistf,thiscol) = cutoffs(thiscol);
    end


end