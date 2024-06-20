%% function cleaned = cleanCrossCorr(data)
% Clean 2d cross correlation array. Leading nan and nan bounded by values
% are converted to 0 as these indicate non-overlap of the clusters;
% trailing nan are left as these indicate having gone paste the edge of the
% cluster.
% REQUIRED INPUTS
% data = m x n array of cross correlation data. Usually (m synapses) x (m
% bins) where (m bins) = 51.
% crop = 0 or 1, where 0 does not crop out the first bin from data, and 1
% does crop out the first bin (ie, usually the first bin is the 0 pixel
% size bin which we want to remove, generally with 51 bins you will only
% look at 50)
% OUTPUTS
% cleaned = m x n array same as data with leading and bounded nan converted
% to 0. If crop = 1, cleaned m = size(data,1)-1.
% Last updated ADL 5/24/23
function cleaned = cleanCrossCorr(data,crop)

    p = inputParser;
    addRequired(p,'data',@isnumeric)
    addRequired(p,'crop',@isnumeric)
    parse(p,data,crop);
   
    if crop ~=0 && crop ~= 1
        error('Set crop value to 0 or 1')
    end
    if crop
        data = data(:,2:end);
    end

    [nrow,ncol] = size(data);
    cleaned = nan(size(data));
    for thisrow = 1:nrow
        
        thisxc = data(thisrow,:);
        nanspots = isnan(thisxc);
        theseprops = regionprops(nanspots,'PixelIdxList');
        if size(theseprops,1) == 0 % if there are no nan regions...
            nanblocks = []; % we'll return the original region
        elseif size(theseprops,1) == 1 % if there is one nan region...
            if ismember(ncol,theseprops.PixelIdxList) % and it's the trailing nans...
                nanblocks = [];                         % then return original data
            else % but if it's not trailing nans...
                nanblocks = theseprops.PixelIdxList; % then we need to reset them to 0
            end
        elseif size(theseprops,1) > 1 % if there is more than one groups of nans...
            nanblocks = vertcat(theseprops(1:end-1).PixelIdxList); % we definitely want to reset all but the last block
            if ~ismember(ncol,theseprops(end).PixelIdxList) 
                % we also reset the last block if it is not a trailing nan 
                %(ie it doesn't include the last column)
                nanblock = [nanblocks; vertcat(theseprops(end).PixelIdxList)];
            end
        end
        thisxc(nanblocks) = 0;
        cleaned(thisrow,:) = thisxc;
    
    end

end