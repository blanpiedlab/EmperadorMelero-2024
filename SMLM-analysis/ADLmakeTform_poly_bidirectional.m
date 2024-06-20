%% [x,y,pairs,d] = ADLmakeTform_poly_bidirectional(filename,pixelsize,direction,varargin)
% Generates a tform(direction).mat file for correcting dual-view localizations
% This specific script uses the polynomial form of fitgeotrans function
% with polynomial of 2.
%
% In addition, this script takes an input "direction" which is definted as
% 'l2r' or 'r2l' and indicates whether the correction should be of the left
% side onto right or the right side onto left. It will save tforml2r.mat or
% tformr2l.mat, and should be used with ADLcorrectDUV_poly_bidirectional to
% correct in either direction.
%
% REQUIRED INPUTS
%     --filename: ex 'duv.csv' that is localization table for dual-view SMLM
%     acquisition of sub-diffraction beads.
%           ***set filename to 'null' and define 'loc' and 'headerstring'
%           name/value pairs to use already loaded data, rather than re-loading a file.
%     --pixelsize: nm per pixel of camera, ie 160
%     --direction: 'l2r' or 'r2l' indicates direction of correction, either
%     left onto right or right onto left.
%
% OPTIONAL INPUTS
%     as optional name/value pairs, you can change the split point, the photon
%     filter, the stretch filter, and the tolerance (below). If not
%     specified, they will be set to the defaults here.
%            name/value: 'SplitPoint', 256 (default) = splitpoint on the
%               camera in pixels
%            name/value: 'PhotonFilter', 500 (default)
%            name/value: 'ErrorFilter', 10 (deafult)
%            name/value: 'Tolerance', 3.6 (default) = distance away the
%               program looks for pairs
%            name/value: 'loc', [] (default). Set only if filename is null.
%               Loc must be in [px]
%            name/value: 'headerstring, '' (default). Set only if filename
%               is null. Make sure loc is in [px].
%            name/value: 'filterp',filterp (empty struct default). If using
%            a .txt matlab loc file, you can specify non-default filtering
%            parameters
%               filterp = (defaults below)
%         struct with fields:
%                   factor: 1
%                windowbot: 100
%                windowtop: 400
%               preccutoff: 10
%             photoncutoff: 400
%                 R2cutoff: 0.4500
%                ampcutoff: 60
%                  stretch: 0.5
%                   layer: 172 ****If you are just doing 2d cali, then set
%                   this to the number of images you took (ie, 50, or 100),
%                   otherwise it's the number of z-steps in the 3d version
%        Tolerance is the maximum pixel distance away from a loc on the left a pair
%        candidate on the right can be after approximate overlay to be considered
%        a pair. If there is more than one loc on the right <tolerance away, no
%        pair will be found. So as tolerance decreases you will only select for
%        paired locs that are already very close together. Thus the "accuracy"
%        will improve but you're really limiting the analysis to already well
%        paired beads.
%        --nKNNframes sets the number of frames to use for the nearest
%        neighbor search because its very slow with >100
%        --'sarah' allows you refilter the data based on DUV correction
%        residual error and redo the tform - ie, find bead pairs with error
%        less than sarah, and then recalculate tform using only those
%        pairs. not sure if this is helpful or kosher, but sarah wanted it.
%        --saveinplace set to 1 to save tform in same folder as filename,
%        rather than in whatever pwd is
% OUTPUTS:
%     --graphs of original bead overlay, best overlay using nearest neighbor minimization,
%     and corrected overlay after tform calculation, as well as histogram
%     of the residual error after DUVcali
%     --tform(direction).mat is saved to the active folder
% General strategy:
%     1) Import loc file and restructure columns to appropriate order (See below)
%     2) Filter for bright and non-stretched locs
%     3) Approximately overlay the locs by splitting the field in half at 256 px
%     4) Maximize the overlap by iterative nearest-neighbor search to minimize
%     distances between locs
%     5) Pair corresponding locs from left and right side of image that are in the same frame
%     6) Calculate tform using fitgeotrans on pairs in indicated direction
%     7) Check result by transforming in indicated direction; save tform.mat
% Modified and updated by Aaron Levy 3/26/20 during Coronocation 2K20
% Updated to be bidirectional ADL 1/6/22

function [x,y,pairs,d] = ADLmakeTform_poly_bidirectional(filename,pixelsize,direction,varargin)

% Parse inputs and optional name/value pairs
p = inputParser;
addRequired(p,'filename',@ischar);
addRequired(p,'pixelsize',@isfloat);
addRequired(p,'direction',@ischar);
addParameter(p,'SplitPoint',256);
addParameter(p,'PhotonFilter',500);
addParameter(p,'ErrorFilter',10);
addParameter(p,'Tolerance',3.6);
addParameter(p,'loc',[]);
addParameter(p,'headerstring','');
addParameter(p,'filterp',struct([]));
addParameter(p,'nKNNframes',101);
addParameter(p,'sarah',0)
addParameter(p,'saveinplace',0)
parse(p,filename,pixelsize,direction,varargin{:});

filename = p.Results.filename;
pixelsize = p.Results.pixelsize;
direction = p.Results.direction;
splitpoint = p.Results.SplitPoint;
photonfilter = p.Results.PhotonFilter;
errorfilter = p.Results.ErrorFilter;
tolerance = p.Results.Tolerance;
filterp = p.Results.filterp;
nKNNframes = p.Results.nKNNframes;
sarah = p.Results.sarah;
saveinplace = p.Results.saveinplace;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Load data from file or existing loc/headerstring in workspace
if strcmp(filename,'null') == 0
    
    [headerstring, loc] = Omniloader(filename,'convertToPixels',pixelsize,'verbose',0);
    col = getColumns(headerstring);
    
elseif strcmp(filename,'null') == 1
    
    loc = p.Results.loc;
    headerstring = p.Results.headerstring;
    
    if ~ischar(headerstring) || ~isfloat(loc) || isempty(loc) || isempty(headerstring) %throw an error if loading issue
        error(['You screwed up the loc/headerstring input loading'...
            'Either you are missing one or they are not formatted right. Try again.'])
    end
    
    % Convert columns to px if needed
    col = getColumns(headerstring);
    cols = [col.x col.y];
    headers = strsplit(headerstring,',');
    pmult = [1 1];    % This will be the values to multiply xyz by to convert to px or nm
    for i = 1:length(cols)
        thiscol = cols(i);
        if ~isempty(regexp(headers{thiscol},'(\[)nm(\])','match'))
            pmult(i) = 1/pixelsize;
            headers{thiscol} = replace(headers{thiscol},'nm','px');
        end
    end
    
    loc(:,col.x) = loc(:,col.x).*pmult(1);
    loc(:,col.y) = loc(:,col.y).*pmult(2);
    % Throw an error if needed
    if strcmp(headers{col.x},'x [nm]')
        error('Your positions are in nm, and need to be in pixels');
    end
    
    S = 'Using loc and headerstring that were provided in function input.\n Converted to pixels if needed.\n';
    fprintf(S)
    
end

mederror = sqrt(loc(:,col.lpx).^2 + loc(:,col.lpy).^2);

nlocs = size(loc,1);

if strcmp(filename,'all_localizations_found.txt') || strcmp(headerstring,['x [px],y [px],sigmax [px],amplitude [photons],offset [photons],'...
        'sigmay [px],ellipticity,frame,photconv,intensity [photons],bkgintensity [photons],bkgstd [photons],error [nm],R^2'])
    % Specifically filter matlab all_localizations_found based on
    % parameters in the 3dcali file.
    if isempty(filterp) == 0 % if filterp was specified use entered values
        loc = loc(loc(:,2) > filterp.windowbot & loc(:,2) < filterp.windowtop,:);
        loc = loc(loc(:,13) < filterp.preccutoff,:);
        loc = loc(loc(:,10) > filterp.photoncutoff,:);
        loc = loc(loc(:,14) > filterp.R2cutoff,:);
        loc = loc(loc(:,4) > filterp.ampcutoff,:);
        loc(:,2) = loc(:,2) * filterp.factor;
        loc(:,6) = loc(:,6) * filterp.factor;
        loc = loc(loc(:,8) > filterp.layer,:);
        loc = loc(abs(loc(:,3)-loc(:,6)) < filterp.stretch,:);
    else % if filterp was not specified run defaults
        S = ['You loaded a matlab .txt file but did not specify the filters. '...
        'Check help for this function to see appropriate filterp struct formatting.'];
        error(S);
    end
else
    % filter anything else for photons and error (unless error isn't a
    % column like for phasor).
    loc = loc(loc(:,col.phot) > photonfilter,:);
    if any(strcmp('error',fieldnames(col)))
        loc = loc(loc(:,col.error) > errorfilter,:);
    end
end

processloc = size(loc,1);
S = 'Started with %d locs and removed %d locs in filtering. \nWorking with %d locs. \n';
fprintf(S,nlocs,nlocs - processloc,processloc);

%~~~~CHANGE COLUMNS HERE AS NEEDED TO APPROPRIATE COLUMN IN INPUT LOC TABLE
loc = [loc(:,col.x) loc(:,col.y) loc(:,col.frame)];

%~~~~~~~~~~~Do approximate alignment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Roughly split into left/right and overlay points
left = loc(loc(:,1) < splitpoint-1,:); % just splits right down the middle to start
right = loc(loc(:,1) > splitpoint,:);
origxsplit = -splitpoint; % save for later
origysplit = 0;
off_x = origxsplit;
off_y = origysplit;

% Automatic search for best x/y re-alignment of beads
breaklimit = 100;
offsx = [off_x; zeros(breaklimit-1,1)];
offsy = [off_y; zeros(breaklimit-1,1)];
tf1 = 0;
tf2 = 0;
loopindex = 1;
while tf1 == 0 || tf2 == 0 % will continue to loop til x and y offset stabilize
    
    % find nearest neighbor and plot calculated shifted points 
    % UPDATE 6/8/21 using only the first nKNNframes frames to deal with
    % large numbers of locs
    left2 = [left(left(:,3)<nKNNframes,1) left(left(:,3)<nKNNframes,2)]; % crop left to just the xy coords
    right2 = [(right(right(:,3)<nKNNframes,1) + off_x) (right(right(:,3)<nKNNframes,2) + off_y)]; % crop y to just xy + offset
    IDX = knnsearch(right2, left2); % find nearest neighbor in right for every point in left
    diff_x = zeros(size(left2,1),1); % initialize some variables
    diff_y = zeros(size(left2,1),1);
    sizeleft = size(left2,1);
    % find the x and y distances from each point in y to its nearest neighbor
    % in x and append to diff_x/y
    for i = 1:sizeleft
        diff_x(i,1) = left2(i,1) - right2(IDX(i),1);
        diff_y(i,1) = left2(i,2) - right2(IDX(i),2);
    end
    % recalculate the offsets with the median differences added in
    off_x = off_x + median(diff_x);
    off_y = off_y + median(diff_y);
    % prepare variables to check whether offsets stabilized
    tf1 = off_x == offsx(loopindex);
    tf2 = off_y == offsy(loopindex);
    offsx(loopindex+1,1) = off_x;
    offsy(loopindex+1,1) = off_y;
    % prevent the loop from going on forever
    if loopindex == breaklimit
        break
    end
    
    loopindex = loopindex + 1;
    
end

% Plot the original split and the new alignment after auto calculation
h1 = figure;
subplot(2,2,1)
hold on
scatter(left(:,1), left(:,2), 10, 'green', 'filled')
scatter(right(:,1)+origxsplit, right(:,2)+origysplit, 10, 'red', 'filled')
axis square
title('Initial overlay of DUV')
legend('Left','Right')
hold off

subplot(2,2,2)
hold on
scatter(left(:,1), left(:,2), 10, 'green', 'filled')
scatter(right(:,1)+off_x, right(:,2)+off_y, 10, 'red', 'filled')
axis square
title(['Stable DUV alignment after NN search took ' num2str(loopindex) ' rounds'])
legend('Left','Right')
hold off

f = 0; % turn this to 1 to plot some diagnostic of what's happening to offsets
if f == 1
    h2 = figure;
    subplot(1,2,1)
    plot(offsx(offsx~=0),'blue');
    title('Rounds to offset_x stabilization')
    
    gca(h2);
    subplot(1,2,2)
    plot(offsy(offsy~=0),'red');
    title('Rounds to offset_y stabilization')
end


%~~~~~~~~~~~~Pair left and right side locs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Find paired locs. For each loc in left, finds whether there are any locs
% in the same frame in right, and then whether they are < some distance away. If
% only 1 loc fits this, it gets paired with the loc in left.
pairs = nan(length(left), 4);
table = tabulate(right(:,3));
a = table(:,2);
rightperframe = mat2cell(right(:,1:2),a);
zeroidx = table(1,1); % 1 if 1 idx, 0 if 0 idx; this is needed for dealing with the rightperframe formating and mismatch between rightperframe{1} having frame 0 not frame 1

for t = 1:length(left)

    x = left(t, 1);     % get left x position of loc t
    y = left(t, 2);     % get left y position of loc t
    pairs(t,1) = x;  % put the x and y values in pairs col 1/2
    pairs(t,2) = y;
    if zeroidx == 0
        rr = rightperframe{left(t,3)+1}; % get locs in same frame when frame is 0 indexed
    elseif zeroidx == 1
        rr = rightperframe{left(t,3)}; % same but when frame is 1 idx
    end
    if isempty(rr)
        continue
    end
    dist = sqrt((x - rr(:,1) - off_x).^2 + (y - rr(:,2) - off_y).^2); % calculate linear distance of each loc in rr to left loc t
    if min(dist) <= tolerance % if any loc is less than 3.6 px away
        f = dist < tolerance; % find the locs <3.6 px awy
        if sum(double(f)) > 1 % if there is more than one quit the loop
            continue
        end
        f = find(dist == min(dist)); % find the index of the loc that is closest
        pairs(t,3) = rr(f,1); % put the x/y values of the closest loc into pairs col 3 and 4
        pairs(t,4) = rr(f,2);
    end
end

%~~~~~~~~~~~~Filter the pairs and calculate tform~~~~~~~~~~~~~~~~~~~~~~~~%

% remove rows with nan values on the right (ie, there was no pair)
f = ~isnan(pairs(:, 3));
pairs = pairs(f, :);

% remove rows with duplicate values on the right side (ie, two different beads on the
% left matched to the same bead on the right)
[~,ia,~] = unique(pairs(:,3:4),'rows','stable');
pairs = pairs(ia,:);

S = ['%0.0f possible bead localization pairs processed. \n'...
    '%0.0f pairs removed for not being a pair. \n'...
    '%0.0f pairs removed for getting paired to more than one bead. \n'...
    '%0.0f bead pairs are going into the tform. \n'];
startpairs = length(left);
lostfornan = startpairs - size(pairs,1);
lostforunique = startpairs - lostfornan - length(ia);
totalpairs = startpairs - lostfornan - lostforunique;
fprintf(S,startpairs,lostfornan,lostforunique,totalpairs);

% Calculate the tform on the pairs
switch direction
    case 'r2l'
        fixedPoints = pairs(:,1:2);                     % left side are fixed
        movingPoints = pairs(:,3:4);                    % right side are moving
    case 'l2r'
        fixedPoints = pairs(:,3:4);                     % right side are fixed
        movingPoints = pairs(:,1:2);                    % left side are moving
    otherwise
        error('Need to input r2l or l2r as desired direction, try again.\n')
end
tform = fitgeotrans(fixedPoints(:,:),movingPoints(:,:),'polynomial',2); % calculate the tform
switch direction
    case 'r2l'
        if saveinplace == 0
            save tformr2l.mat tform;
        elseif saveinplace == 1
            [folder,~,~] = fileparts(filename);
            save(fullfile(folder,'tformr2l.mat'),'tform');
        end
    case 'l2r'
        if saveinplace == 0
            save tforml2r.mat tform;
        elseif saveinplace == 1
            [folder,~,~] = fileparts(filename);
            save(fullfile(folder,'tforml2r.mat'),'tform');
        end
end

% estimate transform error from median distance of shifted right to left side points
[x, y] = transformPointsInverse(tform, movingPoints(:,1), movingPoints(:,2));
dx = x - fixedPoints(:,1);
dy = y - fixedPoints(:,2);
d = sqrt(dx.^2 + dy.^2);
erroramt = median(d) * pixelsize;
S = 'Dual-view cali error is %0.1f nm \n';
fprintf(S,erroramt)

if sarah ~= 0
    S = 'You have chosen to use only filtered points at < %3f DUV error\n';
    fprintf(S,sarah*pixelsize)

    goodpts = d < sarah;
    fixedPoints = fixedPoints(goodpts,:);
    movingPoints = movingPoints(goodpts,:);
    tform = fitgeotrans(fixedPoints(:,:),movingPoints(:,:),'polynomial',2); % calculate the tform
    save tform.mat tform;

    % estimate transform error from median distance of shifted right to left side points
    [x, y] = transformPointsInverse(tform, movingPoints(:,1), movingPoints(:,2));
    dx = x - fixedPoints(:,1);
    dy = y - fixedPoints(:,2);
    d = sqrt(dx.^2 + dy.^2);
    erroramt = median(d) * pixelsize;
    S = 'Dual-view cali error is %0.1f nm \n';
    fprintf(S,erroramt)
end

% Display final figure
gca(h1);
subplot(2,2,3)
hold on;
scatter(fixedPoints(:,1), fixedPoints(:,2), 10,'green','filled');
scatter(x, y, 10, 'red','filled');
title('Aligned points from tform');
switch direction
    case 'r2l'
        legend('Left','Right')
    case 'l2r'
        legend('Right','Left')
end
axis equal

subplot(2,2,4)
hold on
histogram(d*pixelsize,0:1:50,'FaceColor','blue','normalization','probability');
try %needed because xline isn't in 2018a or earlier
    xline(median(mederror*pixelsize),'-',{'Med. loc error ',[ num2str(median(mederror*pixelsize)) ' nm']})
    xline(median(d*pixelsize),'-',{'Med. res. error ',[ num2str(median(d*pixelsize)) ' nm']})
catch ME
    if(strcmp(ME.identifier,'MATLAB:UndefinedFunction'))
        
        a = median(mederror*pixelsize);
        b = median(d*pixelsize);
        S = ['Printing residual errors due to old ML version.\n'...
        'the median loc error was %0.1f nm and the median residual after DUV was %0.1f nm\n'];       
        fprintf(S,a,b)
        text(0,1,'residual error printed in cmd line','Units','normalized')
    end
end 
xlabel('residual error in nm')
ylabel('count')
title('histogram of residual errors')

% Display final figure
% figure
% hold on;
% scatter(pairs(:,1), pairs(:,2), 10,'green','filled');
% scatter(x, y, 10, 'red','filled');
% title('Aligned points from tform');
% legend('Left','Right')
% axis equal


end