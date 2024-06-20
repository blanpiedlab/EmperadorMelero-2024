%% [headerstring,loc,col] = Omniloader(filename,varargin)
% This function will read in one of 3 localization filetypes: .csv, .txt, and .hdf5
% And extract a header string for it, defining all the columns.
% Optionally, it will also convert any column with the header including [nm]
% into camera pixels, or vice versa anything with [px] into [nm].
% REQUIRED INPUTS:
%     -filename: input as char ie, 'locfile.csv'
%     -notes: This function can handle:
%         -.csv files from ThunderSTORM (MLE, Phasor, or Picasso-localized)
%             -3d or 2d localization, with or without merging
%         -.txt files from Matlab
%             -'all_localizations_found.txt' will load and concatenate all of these in the
%                 active director
%             -'all_localization_found0.txt' will load just that file
%             -'localized_3d_calibrated.txt' will load an already processed 3d file
%             -'localized_2d_calibrated.txt' will load an already processed 2d file
%         -.hdf files from Picasso
%             -3d or 2d localization, with or without meging
% OPTIONAL INPUTS:
%     -convert to pixels: as a name/value pair, input 'convertToPixels', pixelsize
%         for ex, 'convertToPixels',160. This will convert xyz and sigma
%         columns to pixel space (ie, divide by pixel size). It will NOT
%         convert error columns to pixels.
%     -convert to nm: as a name/value pair, input 'convertToNm', pixelsize
%         for ex, 'convertToNm',160. This will convert xyz, sigma, and
%         lpx/y from picasso and convert to nm (ie, multiply by pixel
%         size).
%     -verbose: default 1, set to 0 to suppress columns text in cmd
% OUTPUTS:
%     -headerstring: This is a string (formatted for csv output to ThunderSTORM)
%     of the names of all the columns. If you use TSwrite, this can go straight
%     in as the "header" input, assuming you don't add any more columns
%     -loc: This is the loc table as an (n locs) x (m columns) array
%     -col: struct "col" containing the column numbers of each field, it is
%     the output of the "getColumns" function
%
% Created 3/26/20 by Aaron Levy during Coronacation 2K20
    

function [headerstring, loc, col] = Omniloader(filename,varargin)

p = inputParser;

addParameter(p,'convertToPixels',0);
addParameter(p,'convertToNm',0);
addParameter(p,'verbose',1);
parse(p,varargin{:});
convertToPixels = p.Results.convertToPixels;
convertToNm = p.Results.convertToNm;
verbose = p.Results.verbose;

if ~exist(filename,'file')
    [~,name,~] = fileparts(filename);
    if ~exist(fullfile([name '0.txt']),'file')
        error('Filename not found, please check spelling or directory and try again');
    end
end

% Load CSVs (ThunderSTORM, CSV output from Picasso)
if contains(filename,'.csv')
    
    loc = csvread(filename,1,0);
    fid = fopen(filename);
    headerline = textscan(fid, '%s' ,1,'Delimiter','\n');
    headerchar = headerline{1}{1};
    headerstring = replace(headerchar,'"','');
    fclose(fid);
    col = getColumns(headerstring);
    
    if contains(headerstring,'detections') == 1
        
        if isempty(col.z) == 0 && isempty(nonzeros(loc(:,col.z))) == 0
            
            S = 'Loaded a 3D, merged, ThunderSTORM-formatted loc table\n';
            fprintf(S)
            
        elseif (isempty(col.z) == 0 && isempty(nonzeros(loc(:,col.z))) == 1) || isempty(col.z) == 1
            
            S = 'Loaded a 2D, merged, ThunderSTORM-formatted loc table\n';
            fprintf(S)
            
        end
        
    elseif contains(headerstring,'detections') == 0
        
        if isempty(col.z) == 0 && isempty(nonzeros(loc(:,col.z))) == 0
            
            S = 'Loaded a 3D, un-merged, ThunderSTORM-formatted loc table\n';
            fprintf(S)
            
        elseif (isempty(col.z) == 0 && isempty(nonzeros(loc(:,col.z))) == 1) || isempty(col.z) == 1
            
            S = 'Loaded a 2D, un-merged, ThunderSTORM-formatted loc table\n';
            fprintf(S)
            
        end
        
    end
    
end


% Load matlab txt files; either all_localizations_found, or
% localized_3d/2d_calibrated
if contains(filename,'.txt')
    
    if contains(filename,'all_localizations_found')
        
        [~,name,~] = fileparts(filename);
        files = dir([name '*.txt']); % find the localization files
        
        if length(files) == 1
            
            loc = importdata(filename);
            
        elseif length(files) > 1
            
            loc = nan(20000000,14); % initialize the loc variable
            pos = 0;
            nfrm = 0;
            for i = 1:length(files) % loop for the number of files
                filename = fullfile([name num2str(i-1) '.txt']);
                fid = fopen(filename);
                temploc = textscan(fid,repmat('%f', 1, 14),'CollectOutput',1);
                temploc = temploc{:};
                fclose(fid);
                temploc(:,8) = temploc(:,8) + nfrm; % shift the frame numbers in loc1 to be sequential with the previous file
                maxloc = pos + size(temploc,1);
                pos = pos + 1;
                loc(pos:maxloc,:) = temploc;
                pos = maxloc;
                nfrm = max(temploc(:,8)); % shift nfrm to be the new max frame number
            end
            loc(isnan(loc(:,1)),:) = [];
            
        end
        
        headerstring = ['x [px],y [px],sigmax [px],amplitude [photons],offset [photons],'...
            'sigmay [px],ellipticity,frame,photconv,intensity [photons],'...
            'bkgintensity [photons],bkgstd [photons],error [nm],R^2'];
        S = 'Loaded and merged %d Matlab-localized .txt files\n';
        fprintf(S,length(files))
        
       
    end
    
    if contains(filename,'localized_3d_calibrated')
        
        filename = fullfile(filename);
        fid = fopen(filename);
        loc = textscan(fid,repmat('%f', 1, 10),'CollectOutput',1);
        loc = loc{:};
        fclose(fid);
        headerstring = ['x [px],y [px],channel,z [nm],frame,'...
            'intensity [photons],R^2,error [nm],sigmax [px],sigmay [px]'];
        S = 'Loaded a Matlab-localized localized_3d_calibrated file\n';
        fprintf(S)
        
    end
    
    if contains(filename,'localized_2d_calibrated')
        
        filename = fullfile(filename);
        fid = fopen(filename);
        loc = textscan(fid,repmat('%f', 1, 10),'CollectOutput',1);
        loc = loc{:};
        fclose(fid);
        headerstring = ['x [px],y [px],channel,frame,intensity [photons],'...
            'R^2,error [nm],sigmax [px],sigmay [px],ellipticity'];
        S = 'Loaded a Matlab-localized localized_2d_calibrated file\n';
        fprintf(S)
        
    end
    
end

% Load Picasso Hdf5 files
if contains(filename,'.hdf5')
    
    info = h5info(filename);
    infoname = info.Name;
    infotype = info.Datasets.Name;
    s = h5read(filename,[infoname infotype]);
    loc = structfun(@double,s,'UniformOutput',false);
    loc = struct2cell(loc);
    loc = permute(loc, [2 1]);
    loc = cell2mat(loc);
    
    % Parse the headerstring and add [px] or [nm] to appropriate columns
    headerstring = join(fieldnames(s),',');
    headerstring = headerstring{:};
    headers = strsplit(headerstring,',');
    columns = find(ismember(headers,{'x','y','sx','sy','lpx','lpy'}));
    for i = 1:length(columns)
        thiscolumn = columns(i);
        headers{thiscolumn} = strcat(headers{thiscolumn},' [px]');
    end
    if sum(strcmp(headers,'z')) == 1
        zcol = find(strcmp(headers,'z'));
        headers{zcol} = strcat(headers{zcol},' [nm]');
    end
    headerstring = strjoin(headers,',');
    
    col = getColumns(headerstring);
       
    if strcmp(infotype,'clusters') == 1
       
        S = 'Loaded a Picasso hdf5 dbclusters table\n';
        fprintf(S)
        
    elseif strcmp(infotype,'locs') == 1
        
        if contains(headerstring,'len') == 1
            mergetf = 'a merged';
        else
            mergetf = 'an unmerged';
        end
        
        if isempty(col.z) == 1
            ztf = '2D';
        else
            ztf = '3D';
        end
        
        if sum(strcmp(fieldnames(col),'groups')) == 1
            clusteredtf = 'with DBSCAN cluster info';
        else
            clusteredtf = 'with no cluster info';
        end
        
        S = 'Loaded %s, %s Picasso hdf5 loc table, %s.\n';
        fprintf(S,mergetf,ztf,clusteredtf)
        
    end
    
end

% Optionally convert xyz and sigmax/y to pixels
if convertToPixels ~= 0
    
    col = getColumns(headerstring);
    cols = [col.x col.y col.z col.sigmax col.sigmay];
    headers = strsplit(headerstring,',');
    for i = 1:length(cols)
        thiscol = cols(i);
        if ~contains(headers{thiscol},'[px]') %will prevent columns already in px from being double converted
            loc(:,thiscol) = loc(:,thiscol)./convertToPixels;
            headers(thiscol) = replace(headers(thiscol),'[nm]','[px]');
        end
    end
    
    headerstring = strjoin(headers,',');

    S = ['\nYou have converted xyz and sigma that were in nm to camera pixels '...
        'with pixel size %d nm.\n ***NOTE*** that error is still in nm, and if file is Phasor, sigma was not changed.\n'];
    fprintf(S,convertToPixels)
    
end

% Optionally convert xyz and lpx/y error to nm. Matlab and tstorm output
% error in nm already.
if convertToNm ~= 0
    
    col = getColumns(headerstring);
    cols = [col.x col.y col.z col.sigmax col.sigmay col.error col.lpx col.lpy];
    headers = strsplit(headerstring,',');
    for i = 1:length(cols)
        thiscol = cols(i);
        if ~contains(headers{thiscol},'[nm]') % prevent double conversion to nm
            loc(:,thiscol) = loc(:,thiscol).*convertToNm;
            headers(thiscol) = replace(headers(thiscol),'[px]','[nm]');
        end
    end
    
    headerstring = strjoin(headers,',');
    
    S = ['\nYou have converted xyz and sigma that were in camera pixels '...
        'to nm, with pixel size %d nm.\n***NOTE***Error should also be in nm, and if file is Phasor, sigma was not changed.\n'];
    fprintf(S,convertToNm)
    
end

% Display the column names
if verbose
    cnames = strsplit(headerstring,',');
    nums = strsplit(num2str(1:length(cnames)),' ');
    shortstring1 = '\t';
    shortstring2 = ': %s\n';
    string = strjoin(strcat(shortstring1,nums,shortstring2));
    fprintf('Your loaded columns are:\n');
    fprintf(string,cnames{:});
end

col = getColumns(headerstring);

end