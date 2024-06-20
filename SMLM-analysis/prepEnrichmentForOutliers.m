%% function cleaned = prepEnrichmentForOutliers(data,crop)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Function to prepare enrichment data for outlier removal by converting
% trailing 0s to nan (presumably they come from being outside of real
% synapse data). Leading and bounded 0 are not converted. In addition, the
% first column is optionally deleted as this represents the self-enrichment bin (too
% small).
% REQUIRED INPUTS
% data = m x n array of enrichment data where m is number of clusters that
% have enrichments and n is the radii.
% crop = set to 0 to leave data intact, set to 1 to crop first column
% OUTPUTS
% cleaned = same m x n array as data, but with trailing 0 converted to nan
% and if crop == 1, with first column removed
function cleaned = prepEnrichmentForOutliers(data,crop)

    p = inputParser;
    addRequired(p,'data',@isnumeric)
    addRequired(p,'crop',@isnumeric)
    parse(p,data,crop)
    
    if crop~=0 && crop~=1
        error('Crop needs to be 0 or 1 try again')
    end
    
    if crop % crop first column if directed
        data = data(:,2:end);
    end
    
    checknan = sum(sum(isnan(data)));
    if checknan ~= 0
        error('Nans detected you need a different method')
    end

    [nrow,ncol] = size(data);
    cleaned = nan(size(data));
    for thisrow = 1:nrow

        thisenrich = data(thisrow,:);
        zerospots = thisenrich == 0;
        theseprops = regionprops(zerospots,'PixelIdxList');
        if size(theseprops,1) == 0 % if there are no zero regions...
            zeroblocks = []; % we'll return the original region
        elseif size(theseprops,1) == 1 % if there is one zeros region...
            if ismember(ncol,theseprops.PixelIdxList) % and it's the trailing zeros...
                zeroblocks = theseprops.PixelIdxList;   % then we need to convert to nan
            else % but if it's not trailing zeros...
                zeroblocks = []; % then we need to keep them the same
            end
        elseif size(theseprops,1) > 1 % if there is more than one groups of zeros...
            if ismember(ncol,theseprops(end).PixelIdxList) % check the last group to see if it includes the end
                zeroblocks = vertcat(theseprops(end).PixelIdxList); % then we need to convert to nan
            else % otherwise they're all bounded 0 and we don't convert to nan
                zeroblocks = [];
            end         
        end
        thisenrich(zeroblocks) = nan;
        cleaned(thisrow,:) = thisenrich;

    end

end % function
