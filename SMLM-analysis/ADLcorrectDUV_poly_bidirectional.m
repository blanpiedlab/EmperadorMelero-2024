%% [headerstring, loc] = ADLcorrectDUV_poly_bidirectional(filename,pixelsize,varargin)
% “Distinct active zone protein machineries mediate Ca2+ channel clustering and vesicle priming at hippocampal synapses”
% Emperador-Melero et al, 2024
% authors: Aaron D Levy
% copyright: 2024 Blanpied Lab, Dept of Physiology, University of Maryland School of Medicine
%
% Uses tforml2r.mat or tformr2l.mat file generated from ADLmakeTform_poly_bidirectional
% to overlay left and right side of dual-view imaged SMLM.
%
% This assumes that tform(direction).mat and the file you want to process are in
% the active directory. 
%
% This bidirectional version of the function will load a tforml2r.mat or tformr2l.mat
% and correct the localization file in the appropriate direction based on the name.
% CHANNEL ONE IS ALWAYS THE LEFT SIDE OF THE FIELD
% 
% This version does not use the corner (aka, your DUVcali file and images 
% should have the same ROI) and does not zero the image. If you want to include 
% the corner you will need to modify this script to do that and uncomment out some lines
%
% Required Inputs: 
%     -filename: Filename of loc file to load, can be .txt, .csv, or .hdf5
%         Format should be char. Set to 'null' and define name/value pairs
%         'loc' and 'headerstring' to use data in the workspace rather than
%         loading in new locs.
%     -pixelsize: nm per pixel of camera in question (ie, 160)
% OPTIONAL name/value pairs
%     -iqFilePath: Full path to the IQ txt file containing the left/top (corner) information.
%         Format should be char. 
%     -convertToNm: default is 0, set to 1 to save output file in nm
%         rather than camera pixel units
%     -SplitPoint: The left/right split point of the DUV on the camera,
%         default is set to 256 for ixon on storm scope
%     -dontsave: deafult is 0 (off). In this case, the corrected loc
%         file will be saved. if set to 0, the corrected loc file will not
%         be saved.
%     -loc: if you defined filename as 'null' then this is the loc variable
%     -headerstring: if you defined filename as 'null' then this is the
%         headerstring variable
%     -savename: if you have filename as 'null' then you can optionally
%         define a name to save output file as - it will automatically append
%         _DUV.csv, so for example use 'savename','Bsn_PSD95'. Otherwise it will
%         save as 'null_DUV.csv'. If filename is not null it will save as
%         filename_DUV.csv
%     -figflag: default is 1 (on). Set to 0 to suppress figure showing DUV
%         corrected image
%     -loadinplace: default 0, set to 1 to use the folder of the absolute
%     path of the file you are DUV correcting as search directory for tform
%     file
% Output:
%     -returns headerstring and corrected loc in [px], this is for use in other programs.
%     -optional: csv file with shifted  side localizations, plus the added last column
%         will now be "Channel". Saved as csv file with all original x/y positions
%         adjusted, but other columns will be the same. If convertToNm is
%         set to 1, then output will be in nm, otherwise it will be in px
%      -Note that if output is in pixels and you load into ThunderSTORM,
%      it will convert to nm depending n what the camera settings are.
% Original code ADL 3/2020
% Updated to be bidirectional 1/2022
% Added figflag input 5/6/22 ADL

function [headerstring, loc] = ADLcorrectDUV_poly_bidirectional(filename,pixelsize,varargin)

% Parse the inputs
p = inputParser;
addRequired(p,'filename',@ischar);
addRequired(p,'pixelsize',@isfloat);
addParameter(p,'iqFilePath',@ischar);
addParameter(p,'convertToNm',0);
addParameter(p,'SplitPoint',256);
addParameter(p,'dontsave',0);
addParameter(p,'loc',[]);
addParameter(p,'headerstring','');
addParameter(p,'savename','',@ischar);
addParameter(p,'figflag',1,@isnumeric);
addParameter(p,'loadinplace',0,@isnumeric);

parse(p,filename,pixelsize,varargin{:});

filename = p.Results.filename;
pixelsize = p.Results.pixelsize;
iqFilePath = p.Results.iqFilePath; % only needed if correcting for corner which we generally are not.
convertToNm = p.Results.convertToNm;
splitpoint = p.Results.SplitPoint;
dontsave = p.Results.dontsave;
savename = p.Results.savename;
loadinplace = p.Results.loadinplace;


% Load tform and loc file, converting to pixels. Find the x and y columns.
%[fov_left, fov_top] = getCorner(iqFilePath);
if loadinplace == 0
    files = dir('tform*.mat');
    t = load(files.name);
elseif loadinplace == 1
    [folder,~,~] = fileparts(filename);
    files = dir(fullfile(folder,'tform*.mat'));
    t = load(fullfile(files.folder,files.name));
end

if contains(files.name,'l2r')
    direction = 'l2r';
elseif contains(files.name,'r2l')
    direction = 'r2l';
else
    S = 'You have some naming issue with your tform, it should be tforml2r or tformr2l.mat\n';
    fprintf(S)
end

if strcmp(filename,'null') == 1
    loc = p.Results.loc;
    headerstring = p.Results.headerstring;
    
    if ~ischar(headerstring) || ~isfloat(loc) || isempty(loc) || isempty(headerstring)
        warning(['You screwed up the loc/headerstring input loading\n'...
            'Either you are missing one or they are not formatted right. Try again.'])
        return
    end
    
    S = 'Using loc and headerstring that were provided in function input.\n';
    fprintf(S)
elseif strcmp(filename,'null') == 0
    [headerstring, loc] = Omniloader(filename,'convertToPixels',pixelsize,'verbose',0);
end

col = getColumns(headerstring);

S = 'Transforming locs...\n';
fprintf(S)

% shift the locs by FOV then assign channels to left and right
% loc(:,col.x) = loc(:,col.x) + fov_left; %uncomment these two lines if using corner
% loc(:,col.y) = loc(:,col.y) + fov_top;
left = loc(:,col.x) < splitpoint-1;
right = loc(:,col.x) > splitpoint;
if ~strcmp(headerstring, ['x [px],y [px],sigmax [px],amplitude [photons],offset [photons],'...
            'sigmay [px],ellipticity,frame,photconv,intensity [photons],'...
            'bkgintensity [photons],bkgstd [photons],error [nm],R^2'])
    loc(:,end+1) = nan;
    loc(left, end) = 1; %The left side is always channel 1
    loc(right, end) = 2; %The right side is always channel  2
    chcol = size(loc,2);
    headerstring = [headerstring, ',Channel'];
else
    loc(left,col.amp) = 1;
    loc(right,col.amp) = 2;
    chcol = col.amp;
    headers = strsplit(headerstring,',');
    headers{col.amp} = 'Channel';
    headerstring = strjoin(headers,',');
end
    

% apply the tform to spatially transform the locs, then replace the
% original values
switch direction
    case 'r2l'
        [xc, yc] = transformPointsInverse(t.tform, loc(right,col.x), loc(right,col.y)); 
        loc(right,col.x) = xc;
        loc(right,col.y) = yc;
        fprintf('Moved right side onto left side\n')
    case 'l2r'
        [xc, yc] = transformPointsInverse(t.tform, loc(left,col.x), loc(left,col.y)); 
        loc(left,col.x) = xc;
        loc(left,col.y) = yc;
        fprintf('Moved left side onto right side\n')
end
% subtracting off the minimum x/y position if desired
% lowx = min([min(xc) min(loc(left,col.x))]); 
% lowy = min([min(yc) min(loc(left,col.y))]); 
% loc(left,col.x) = (loc(left,col.x) - lowx); 
% loc(left,col.y) = (loc(left,col.y) - lowy); 
% loc(right,col.x) = (xc - lowx); 
% loc(right,col.y) = (yc - lowy);


%    
%Diagnostic figure
if p.Results.figflag
    figure
    hold on; axis equal;
    
    plot(loc(loc(:,chcol)==1,col.x),loc(loc(:,chcol)==1,col.y), 'Color', 'blue',... %plots x and y for 647 with characteristics below in blue
        'Markersize', 1,...
        'LineStyle', 'none',...
        'MarkerFaceColor', 'blue',...
        'Marker', 'o');
    plot(loc(loc(:,chcol)==2,col.x),loc(loc(:,chcol)==2,col.y), 'Color', 'red',... %plots x and y for 561 with characteristics below in red
        'Markersize', 1,...
        'LineStyle', 'none',...
        'MarkerFaceColor', 'red',...
        'Marker', 'o');
    title('DUV corrected without corner!')
end

S = 'Transformation done and locs are now DUV corrected.\n';
fprintf(S)

% Save file as CSV using TSwrite with appropriate header   

if isempty(savename)  && strcmp(filename,'null')
    [~,savename,~] = fileparts(filename);
end
savename = [savename '_DUV.csv'];

if dontsave == 0
    S = 'Saving...\n';
    fprintf(S)
    if convertToNm == 1 % if this option was set in input, save out CSV in nm.
        TSwrite(loc,savename,headerstring,'convertToNm',pixelsize,'verbose',0);
        S = 'Locs saved with nm units, DUV correction complete.\n';
        fprintf(S)
    elseif convertToNm == 0 % otherwise save it as pixels.
        TSwrite(loc,savename,headerstring,'verbose',0);
        S = 'Locs saved in pixel dimensions, DUV correction complete.\n';
        fprintf(S)
    end
end
    
    
end
